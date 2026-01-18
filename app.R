# Gene Correlation App
# Based on Fredrik Wermeling's correlation script
# Standalone Shiny application

# ============================================================================
# SETUP - Install packages if needed
# ============================================================================

required_packages <- c("shiny", "data.table", "igraph", "ggplot2", "ggraph",
                       "ggrepel", "DT", "shinyjs", "httr", "jsonlite", "visNetwork")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Increase max upload size to 500MB
options(shiny.maxRequestSize = 500*1024^2)

library(shiny)
library(shinyjs)
library(data.table)
library(igraph)
library(ggplot2)
library(ggraph)
library(ggrepel)
library(DT)
library(httr)
library(jsonlite)
library(visNetwork)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Fetch genes from GO term using QuickGO API
fetch_go_genes <- function(go_id) {
  # Clean GO ID format
  go_id <- toupper(trimws(go_id))
  if (!grepl("^GO:", go_id)) {
    go_id <- paste0("GO:", go_id)
  }
  
  tryCatch({
    # Use QuickGO API
    url <- paste0(
      "https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?",
      "goId=", go_id,
      "&taxonId=9606",
      "&geneProductType=protein",
      "&downloadLimit=2000"
    )
    
    response <- GET(url, add_headers(Accept = "text/tsv"))
    
    if (status_code(response) == 200) {
      content_text <- content(response, "text", encoding = "UTF-8")
      lines <- strsplit(content_text, "\n")[[1]]
      
      if (length(lines) > 1) {
        # Parse TSV - gene symbols are typically in column 3
        genes <- unique(sapply(lines[-1], function(line) {
          parts <- strsplit(line, "\t")[[1]]
          if (length(parts) >= 3) parts[3] else NA
        }))
        genes <- genes[!is.na(genes) & genes != ""]
        return(list(success = TRUE, genes = genes, count = length(genes)))
      }
    }
    
    return(list(success = FALSE, error = "No genes found for this GO term"))
  }, error = function(e) {
    return(list(success = FALSE, error = paste("API error:", e$message)))
  })
}

# Run correlation analysis (adapted from your script)
run_correlation_analysis <- function(reference_data, gene_list, analysis_mode, correlation_cutoff) {
  
  # Clean gene list
  gene_list <- toupper(trimws(gene_list))
  gene_list <- gene_list[gene_list != ""]
  
  # Match genes
matched_columns <- colnames(reference_data)[colnames(reference_data) %in% gene_list]
  not_found_genes <- setdiff(gene_list, matched_columns)
  
  if (length(matched_columns) == 0) {
    return(list(
      success = FALSE,
      error = "No matching genes found in the reference data.",
      not_found = not_found_genes
    ))
  }
  
  # Compute correlations
  if (analysis_mode == "analysis") {
    filtered_data <- reference_data[, ..matched_columns]
    cor_matrix <- cor(filtered_data, use = "pairwise.complete.obs")
    
    genes_in_matrix <- colnames(cor_matrix)
    output_table <- data.table(
      Gene1 = rep(genes_in_matrix, times = length(genes_in_matrix)),
      Gene2 = rep(genes_in_matrix, each = length(genes_in_matrix)),
      Correlation = as.vector(cor_matrix)
    )
  } else {
    # Design mode
    filtered_data <- reference_data[, ..matched_columns]
    cor_matrix <- cor(reference_data, filtered_data, use = "pairwise.complete.obs")
    
    all_genes <- rownames(cor_matrix)
    matched_genes <- colnames(cor_matrix)
    output_table <- data.table(
      Gene1 = rep(matched_genes, each = length(all_genes)),
      Gene2 = rep(all_genes, times = length(matched_genes)),
      Correlation = as.vector(cor_matrix)
    )
  }
  
  # Filter by cutoff
  filtered_correlations <- output_table[
    abs(Correlation) > correlation_cutoff & (Gene1 != Gene2)
  ]
  filtered_correlations[, Correlation := round(Correlation, 3)]
  
  if (nrow(filtered_correlations) == 0) {
    return(list(
      success = FALSE,
      error = paste("No correlations found above cutoff of", correlation_cutoff),
      not_found = not_found_genes,
      matched = matched_columns
    ))
  }
  
  # Build graph and find clusters
  graph <- graph_from_data_frame(
    filtered_correlations[, .(Gene1, Gene2, Correlation)], 
    directed = FALSE
  )
  
  clusters_info <- components(graph)
  cluster_membership <- clusters_info$membership
  filtered_correlations[, Cluster := cluster_membership[Gene1]]
  
  # Add star markers for design mode
  if (analysis_mode == "design") {
    filtered_correlations[, InGeneList_Gene1 := ifelse(Gene1 %in% matched_columns, "*", "")]
    filtered_correlations[, InGeneList_Gene2 := ifelse(Gene2 %in% matched_columns, "*", "")]
  }
  
  # Build clusters table
  found_genes <- V(graph)$name
  found_clusters <- cluster_membership[found_genes]
  
  if (analysis_mode == "design") {
    cluster_dt <- data.table(
      Gene = found_genes,
      Cluster = found_clusters,
      In_GeneList = ifelse(found_genes %in% matched_columns, "*", "")
    )
  } else {
    cluster_dt <- data.table(
      Gene = found_genes,
      Cluster = found_clusters
    )
  }
  
  # Add effect size stats
  cluster_dt[, Mean_Effect := sapply(Gene, function(g) {
    if (g %in% colnames(reference_data)) {
      round(mean(reference_data[[g]], na.rm = TRUE), 2)
    } else NA
  })]
  
  cluster_dt[, SD_Effect := sapply(Gene, function(g) {
    if (g %in% colnames(reference_data)) {
      round(sd(reference_data[[g]], na.rm = TRUE), 2)
    } else NA
  })]
  
  setorder(cluster_dt, Cluster)
  setorder(filtered_correlations, Cluster)
  
  return(list(
    success = TRUE,
    correlations = filtered_correlations,
    clusters = cluster_dt,
    graph = graph,
    matched = matched_columns,
    not_found = not_found_genes,
    mode = analysis_mode,
    cutoff = correlation_cutoff
  ))
}

# Find gene synonyms using MyGene.info API
find_gene_synonyms <- function(gene_symbol, available_genes) {
  gene_symbol <- toupper(trimws(gene_symbol))

  tryCatch({
    # Query MyGene.info for the gene
    url <- paste0("https://mygene.info/v3/query?q=", gene_symbol,
                  "&species=human&fields=symbol,alias,other_names")

    response <- GET(url, add_headers(Accept = "application/json"))

    if (status_code(response) == 200) {
      result <- content(response, "parsed")

      if (result$total > 0 && length(result$hits) > 0) {
        # Collect all possible names
        all_names <- c()

        for (hit in result$hits) {
          if (!is.null(hit$symbol)) {
            all_names <- c(all_names, toupper(hit$symbol))
          }
          if (!is.null(hit$alias)) {
            if (is.list(hit$alias)) {
              all_names <- c(all_names, toupper(unlist(hit$alias)))
            } else {
              all_names <- c(all_names, toupper(hit$alias))
            }
          }
          if (!is.null(hit$other_names)) {
            if (is.list(hit$other_names)) {
              all_names <- c(all_names, toupper(unlist(hit$other_names)))
            } else {
              all_names <- c(all_names, toupper(hit$other_names))
            }
          }
        }

        all_names <- unique(all_names)

        # Find which synonyms are in the available genes
        matches <- intersect(all_names, available_genes)

        if (length(matches) > 0) {
          return(list(
            found = TRUE,
            original = gene_symbol,
            matches = matches
          ))
        }
      }
    }

    return(list(found = FALSE, original = gene_symbol, matches = character(0)))

  }, error = function(e) {
    return(list(found = FALSE, original = gene_symbol, matches = character(0)))
  })
}

# Run enrichment analysis using Enrichr API
run_enrichment_analysis <- function(gene_list, databases = c("GO_Biological_Process_2023",
                                                              "GO_Molecular_Function_2023",
                                                              "KEGG_2021_Human",
                                                              "Reactome_2022")) {
  gene_list <- unique(toupper(trimws(gene_list)))
  gene_list <- gene_list[gene_list != ""]

  if (length(gene_list) < 2) {
    return(list(success = FALSE, error = "Need at least 2 genes"))
  }

  tryCatch({
    # Step 1: Add gene list to Enrichr
    add_url <- "https://maayanlab.cloud/Enrichr/addList"
    gene_str <- paste(gene_list, collapse = "\n")

    add_response <- POST(
      add_url,
      body = list(list = gene_str, description = "Gene Correlation Analysis"),
      encode = "multipart"
    )

    if (status_code(add_response) != 200) {
      return(list(success = FALSE, error = "Failed to submit genes to Enrichr"))
    }

    add_result <- content(add_response, "parsed")
    user_list_id <- add_result$userListId

    # Step 2: Get enrichment results for each database
    all_results <- list()

    for (db in databases) {
      enrich_url <- paste0("https://maayanlab.cloud/Enrichr/enrich?userListId=",
                           user_list_id, "&backgroundType=", db)

      Sys.sleep(0.3)  # Rate limiting

      enrich_response <- GET(enrich_url)

      if (status_code(enrich_response) == 200) {
        enrich_data <- content(enrich_response, "parsed")

        if (length(enrich_data[[db]]) > 0) {
          # Parse results
          for (term in enrich_data[[db]]) {
            all_results[[length(all_results) + 1]] <- data.table(
              Database = db,
              Term = term[[2]],
              P_value = term[[3]],
              Adj_P_value = term[[7]],
              Odds_Ratio = term[[4]],
              Combined_Score = term[[5]],
              Genes = paste(term[[6]], collapse = ";"),
              Gene_Count = length(term[[6]])
            )
          }
        }
      }
    }

    if (length(all_results) == 0) {
      return(list(success = FALSE, error = "No enrichment results found"))
    }

    results_dt <- rbindlist(all_results)
    results_dt <- results_dt[order(Adj_P_value)]

    # Clean up term names
    results_dt[, Term := gsub(" \\(GO:\\d+\\)$", "", Term)]
    results_dt[, Term := gsub("_Homo sapiens.*$", "", Term)]

    return(list(success = TRUE, results = results_dt))

  }, error = function(e) {
    return(list(success = FALSE, error = paste("Enrichr API error:", e$message)))
  })
}

# ============================================================================
# UI
# ============================================================================

ui <- fluidPage(
  useShinyjs(),
  
  # Custom CSS for styling - Green Listed theme
  tags$head(
    tags$style(HTML("
      @import url('https://fonts.googleapis.com/css2?family=Source+Sans+Pro:wght@300;400;600;700&family=Roboto+Mono:wght@400;500&display=swap');

      :root {
        --green-900: #14532d;
        --green-800: #166534;
        --green-700: #15803d;
        --green-600: #16a34a;
        --green-500: #22c55e;
        --green-400: #4ade80;
        --green-300: #86efac;
        --green-200: #bbf7d0;
        --green-100: #dcfce7;
        --green-50: #f0fdf4;
        --gray-600: #4b5563;
        --gray-500: #6b7280;
        --gray-200: #e5e7eb;
        --gray-100: #f3f4f6;
      }

      body {
        font-family: 'Source Sans Pro', -apple-system, BlinkMacSystemFont, sans-serif;
        background-color: var(--green-50);
        min-height: 100vh;
        color: #1f2937;
      }

      .main-header {
        background: linear-gradient(135deg, var(--green-700) 0%, var(--green-600) 100%);
        color: white;
        padding: 25px 30px;
        margin: -15px -15px 25px -15px;
        border-bottom: 3px solid var(--green-400);
      }

      .main-header h1 {
        margin: 0;
        font-weight: 700;
        font-size: 1.8em;
        letter-spacing: -0.3px;
      }

      .main-header p {
        margin: 8px 0 0 0;
        opacity: 0.9;
        font-size: 1.05em;
        font-weight: 300;
      }

      .card {
        background: white;
        border-radius: 8px;
        padding: 20px;
        margin-bottom: 16px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.08);
        border: 1px solid var(--gray-200);
      }

      .card-title {
        font-weight: 600;
        font-size: 1em;
        color: var(--green-800);
        margin-bottom: 16px;
        padding-bottom: 10px;
        border-bottom: 2px solid var(--green-100);
      }

      .btn-primary {
        background: var(--green-600);
        border: none;
        padding: 10px 24px;
        font-weight: 600;
        border-radius: 6px;
        transition: all 0.2s ease;
        color: white;
      }

      .btn-primary:hover, .btn-primary:focus {
        background: var(--green-700);
        box-shadow: 0 2px 8px rgba(13, 148, 136, 0.3);
        color: white;
      }

      .btn-success {
        background: var(--green-500);
        border: none;
        border-radius: 6px;
        font-weight: 600;
      }

      .btn-success:hover, .btn-success:focus {
        background: var(--green-600);
      }

      .btn-default, .btn-secondary {
        background: var(--gray-100);
        border: 1px solid var(--gray-200);
        color: var(--gray-600);
        border-radius: 6px;
      }

      .form-control, .selectize-input {
        border-radius: 6px;
        border: 1px solid var(--gray-200);
        padding: 8px 12px;
        font-size: 0.95em;
      }

      .form-control:focus {
        border-color: var(--green-500);
        box-shadow: 0 0 0 3px rgba(20, 184, 166, 0.15);
        outline: none;
      }

      .nav-tabs {
        border-bottom: 2px solid var(--gray-200);
      }

      .nav-tabs .nav-link {
        border: none;
        color: var(--gray-500);
        font-weight: 500;
        padding: 10px 16px;
        border-radius: 6px 6px 0 0;
        font-size: 0.9em;
      }

      .nav-tabs .nav-link:hover {
        color: var(--green-600);
        background: var(--green-50);
      }

      .nav-tabs .nav-link.active {
        background: var(--green-600);
        color: white;
      }

      .status-box {
        padding: 12px 15px;
        border-radius: 6px;
        margin: 12px 0;
        font-family: 'Roboto Mono', monospace;
        font-size: 0.85em;
      }

      .status-success {
        background: var(--green-50);
        color: var(--green-800);
        border-left: 3px solid var(--green-500);
      }

      .status-warning {
        background: #fffbeb;
        color: #92400e;
        border-left: 3px solid #f59e0b;
      }

      .status-error {
        background: #fef2f2;
        color: #991b1b;
        border-left: 3px solid #ef4444;
      }

      .gene-count {
        display: inline-block;
        background: var(--green-600);
        color: white;
        padding: 3px 10px;
        border-radius: 12px;
        font-size: 0.8em;
        font-weight: 600;
        margin-left: 8px;
      }

      #gene_textarea {
        font-family: 'Roboto Mono', monospace;
        min-height: 180px;
        font-size: 0.9em;
      }

      .footer {
        text-align: center;
        padding: 20px;
        color: var(--gray-500);
        font-size: 0.85em;
        margin-top: 20px;
        border-top: 1px solid var(--gray-200);
      }

      .footer a {
        color: var(--green-600);
        text-decoration: none;
        font-weight: 500;
      }

      .footer a:hover {
        color: var(--green-700);
        text-decoration: underline;
      }

      /* Radio buttons and checkboxes */
      .radio label, .checkbox label {
        font-weight: 400;
      }

      /* Slider styling */
      .irs--shiny .irs-bar {
        background: var(--green-500);
        border-top: 1px solid var(--green-600);
        border-bottom: 1px solid var(--green-600);
      }

      .irs--shiny .irs-handle {
        background: var(--green-600);
        border: 2px solid white;
        box-shadow: 0 1px 3px rgba(0,0,0,0.2);
      }

      .irs--shiny .irs-from, .irs--shiny .irs-to, .irs--shiny .irs-single {
        background: var(--green-600);
      }

      /* DataTables styling */
      .dataTables_wrapper .dataTables_filter input {
        border: 1px solid var(--gray-200);
        border-radius: 6px;
        padding: 6px 10px;
      }

      .dataTables_wrapper .dataTables_length select {
        border: 1px solid var(--gray-200);
        border-radius: 6px;
      }

      table.dataTable thead th {
        background: var(--green-50);
        color: var(--green-800);
        font-weight: 600;
        border-bottom: 2px solid var(--green-200) !important;
      }

      /* Download buttons */
      .btn[id^='download'] {
        background: var(--green-600);
        color: white;
        border: none;
        border-radius: 6px;
      }

      .btn[id^='download']:hover {
        background: var(--green-700);
      }

      /* Help text */
      .help-block {
        color: var(--gray-500);
        font-size: 0.85em;
      }

      /* visNetwork styling */
      .visNetwork {
        border: 1px solid var(--gray-200);
        border-radius: 6px;
      }
    "))
  ),

  # Header
  div(class = "main-header",
      h1("Gene Correlation Explorer"),
      p("Analyze gene correlations from DepMap CRISPR screen data")
  ),
  
  fluidRow(
    # Left panel - Inputs
    column(4,
      # Reference data card
      div(class = "card",
        div(class = "card-title", "1. Load Reference Data"),

        tabsetPanel(id = "ref_input_method",
          tabPanel("Upload",
            br(),
            fileInput("reference_file",
                      "Upload DepMap CRISPRGeneEffect.csv",
                      accept = ".csv"),
            helpText("Max file size: 500MB")
          ),

          tabPanel("Local Path",
            br(),
            textInput("reference_path",
                      "Enter full path to file:",
                      value = "~/Downloads/CRISPRGeneEffect.csv"),
            actionButton("load_local", "Load File", class = "btn-success btn-sm")
          )
        ),

        helpText("Download from: ",
                 tags$a(href = "https://depmap.org/portal/data_page/?tab=allData",
                        "DepMap Portal", target = "_blank")),
        uiOutput("reference_status")
      ),
      
      # Gene input card
      div(class = "card",
        div(class = "card-title", "2. Input Genes"),
        
        tabsetPanel(id = "input_method",
          tabPanel("Paste/Type",
            br(),
            textAreaInput("gene_textarea",
                          "Enter gene symbols (one per line):",
                          rows = 8,
                          placeholder = "TP53\nBRCA1\nMYC\n..."),
            fluidRow(
              column(6, actionButton("clear_genes", "Clear", class = "btn-sm")),
              column(6, actionButton("load_test_genes", "Load Test Genes", class = "btn-sm btn-info"))
            )
          ),

          tabPanel("Upload File",
            br(),
            fileInput("gene_file", "Upload .txt file",
                      accept = c(".txt", ".csv")),
            helpText("One gene symbol per line")
          ),

          tabPanel("GO Term",
            br(),
            textInput("go_term", "Enter GO term ID:",
                      placeholder = "GO:0006955"),
            actionButton("fetch_go", "Fetch Genes", class = "btn-success btn-sm"),
            br(), br(),
            uiOutput("go_status"),
            helpText("Example: GO:0006955 (immune response)")
          )
        ),

        br(),
        uiOutput("gene_count_display"),
        uiOutput("gene_validation_display")
      ),
      
      # Parameters card
      div(class = "card",
        div(class = "card-title", "3. Set Parameters"),
        
        radioButtons("analysis_mode", "Analysis Mode:",
          choices = list(
            "Analysis (within gene list)" = "analysis",
            "Design (find correlated genes)" = "design"
          ),
          selected = "analysis"
        ),
        
        sliderInput("correlation_cutoff", "Correlation Cutoff:",
                    min = 0.2, max = 0.8, value = 0.5, step = 0.05),
        
        helpText("Higher cutoff = fewer, stronger correlations"),
        
        br(),
        actionButton("run_analysis", "Run Analysis", 
                     class = "btn-primary btn-lg", 
                     style = "width: 100%;"),
        
        br(), br(),
        uiOutput("analysis_status")
      )
    ),
    
    # Right panel - Results
    column(8,
      div(class = "card",
        div(class = "card-title", "Results"),
        
        tabsetPanel(id = "results_tabs",
          tabPanel("Network",
            br(),
            # Controls row
            fluidRow(
              column(3,
                sliderInput("net_font_size", "Font Size:",
                            min = 10, max = 40, value = 16, step = 2),
                checkboxInput("net_show_gene_effect",
                              "Show gene effect in labels",
                              value = FALSE),
                conditionalPanel(
                  condition = "input.net_show_gene_effect == true",
                  checkboxInput("net_show_sd", "Include +/- SD", value = FALSE)
                )
              ),
              column(3,
                sliderInput("net_node_size", "Node Size:",
                            min = 10, max = 60, value = 25, step = 5),
                sliderInput("net_edge_width", "Edge Width:",
                            min = 1, max = 10, value = 3, step = 1)
              ),
              column(3,
                checkboxInput("net_physics", "Enable physics", value = FALSE),
                actionButton("net_stabilize", "Auto-arrange", class = "btn-sm"),
                actionButton("net_fit", "Fit to View", class = "btn-sm")
              ),
              column(3,
                h5("Export Figure"),
                actionButton("export_png", "PNG", class = "btn-success btn-sm"),
                actionButton("export_svg", "SVG", class = "btn-success btn-sm"),
                actionButton("export_pdf", "PDF", class = "btn-success btn-sm"),
                br(),
                helpText("Uncheck physics, arrange nodes, then export")
              ),
              column(3,
                h5("Node Visibility"),
                actionButton("show_all_nodes", "Show All Hidden Nodes", class = "btn-outline-secondary btn-sm"),
                br(),
                checkboxInput("show_legend_nodes", "Show Legend in Network", value = TRUE),
                helpText("Click on a node to hide it")
              )
            ),
            hr(),
            # Hidden nodes info
            fluidRow(
              column(12,
                uiOutput("hidden_nodes_info")
              )
            ),
            # Network
            div(style = "background: white; border: 1px solid #ddd; border-radius: 8px;",
              visNetworkOutput("network_plot", height = "550px")
            ),
            # JavaScript for exports with white background (legend is in the network)
            tags$script(HTML("
              // Helper function to add white background
              function canvasWithWhiteBg(originalCanvas) {
                var newCanvas = document.createElement('canvas');
                newCanvas.width = originalCanvas.width;
                newCanvas.height = originalCanvas.height;
                var ctx = newCanvas.getContext('2d');
                ctx.fillStyle = '#FFFFFF';
                ctx.fillRect(0, 0, newCanvas.width, newCanvas.height);
                ctx.drawImage(originalCanvas, 0, 0);
                return newCanvas;
              }

              Shiny.addCustomMessageHandler('exportNetworkPNG', function(message) {
                var network = HTMLWidgets.find('#network_plot');
                if (network && network.network) {
                  var originalCanvas = network.network.canvas.frame.canvas;
                  var canvas = canvasWithWhiteBg(originalCanvas);
                  var link = document.createElement('a');
                  link.download = 'correlation_network_' + new Date().toISOString().slice(0,10) + '.png';
                  link.href = canvas.toDataURL('image/png');
                  link.click();
                }
              });

              Shiny.addCustomMessageHandler('exportNetworkSVG', function(message) {
                var network = HTMLWidgets.find('#network_plot');
                if (network && network.network) {
                  var originalCanvas = network.network.canvas.frame.canvas;
                  var canvas = canvasWithWhiteBg(originalCanvas);
                  var svgData = '<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"' + canvas.width + '\" height=\"' + canvas.height + '\">' +
                    '<rect width=\"100%\" height=\"100%\" fill=\"white\"/>' +
                    '<image href=\"' + canvas.toDataURL('image/png') + '\" width=\"' + canvas.width + '\" height=\"' + canvas.height + '\"/>' +
                    '</svg>';
                  var link = document.createElement('a');
                  link.download = 'correlation_network_' + new Date().toISOString().slice(0,10) + '.svg';
                  link.href = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(svgData);
                  link.click();
                }
              });
            "))
          ),

          tabPanel("Correlations",
            br(),
            downloadButton("download_correlations", "Download CSV"),
            br(), br(),
            DTOutput("correlations_table")
          ),

          tabPanel("Clusters",
            br(),
            downloadButton("download_clusters", "Download CSV"),
            br(), br(),
            DTOutput("clusters_table")
          ),

          tabPanel("Summary",
            br(),
            verbatimTextOutput("summary_text")
          ),

          tabPanel("Enrichment",
            br(),
            fluidRow(
              column(4,
                h5("Gene Set Enrichment Analysis"),
                helpText("Analyze your input genes against pathway databases"),
                checkboxGroupInput("enrichment_dbs", "Databases:",
                  choices = c(
                    "GO Biological Process" = "GO_Biological_Process_2023",
                    "GO Molecular Function" = "GO_Molecular_Function_2023",
                    "GO Cellular Component" = "GO_Cellular_Component_2023",
                    "KEGG Pathways" = "KEGG_2021_Human",
                    "Reactome" = "Reactome_2022",
                    "WikiPathways" = "WikiPathway_2023_Human"
                  ),
                  selected = c("GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022")
                ),
                actionButton("run_enrichment", "Run Enrichment Analysis",
                             class = "btn-primary"),
                br(), br(),
                uiOutput("enrichment_status")
              ),
              column(8,
                conditionalPanel(
                  condition = "output.enrichment_ready",
                  h5("Enrichment Results"),
                  helpText("Click a row to load pathway genes for correlation analysis"),
                  downloadButton("download_enrichment", "Download Results", class = "btn-sm"),
                  br(), br(),
                  DTOutput("enrichment_table")
                )
              )
            )
          )
        )
      )
    )
  ),

  # Footer
  div(class = "footer",
    p("Gene Correlation Explorer | Part of the ",
      tags$a(href = "https://greenlisted.cmm.se", "Green Listed"),
      " tool family"),
    p(tags$a(href = "https://wermelinglab.com", "Wermeling Lab"),
      " | Karolinska Institutet")
  )
)

# ============================================================================
# SERVER
# ============================================================================

server <- function(input, output, session) {
  
  # Reactive values to store data
  rv <- reactiveValues(
    reference_data = NULL,
    gene_list = character(0),
    results = NULL,
    enrichment_results = NULL
  )
  
  # -------------------------------------------------------------------------
  # Reference data loading
  # -------------------------------------------------------------------------
  
  # Helper function to load and process reference data
  load_reference_data <- function(filepath) {
    withProgress(message = "Loading reference data...", {
      tryCatch({
        ref_data <- fread(filepath, showProgress = FALSE)

        # Remove first column if it's sample IDs
        if (is.character(ref_data[[1]])) {
          ref_data <- ref_data[, -1, with = FALSE]
        }

        # Clean column names
        original_colnames <- colnames(ref_data)
        colnames(ref_data) <- toupper(trimws(gsub(" \\(.*\\)$", "", original_colnames)))

        rv$reference_data <- ref_data

        showNotification(
          paste("Loaded", ncol(ref_data), "genes from", nrow(ref_data), "cell lines"),
          type = "message"
        )
      }, error = function(e) {
        showNotification(paste("Error loading file:", e$message), type = "error")
      })
    })
  }

  observeEvent(input$reference_file, {
    req(input$reference_file)
    load_reference_data(input$reference_file$datapath)
  })

  # Load from local path
  observeEvent(input$load_local, {
    req(input$reference_path)

    filepath <- path.expand(trimws(input$reference_path))

    if (!file.exists(filepath)) {
      showNotification("File not found. Please check the path.", type = "error")
      return()
    }

    load_reference_data(filepath)
  })
  
  output$reference_status <- renderUI({
    if (is.null(rv$reference_data)) {
      div(class = "status-box status-warning",
          "⚠️ No reference data loaded")
    } else {
      div(class = "status-box status-success",
          paste("✓", ncol(rv$reference_data), "genes,", 
                nrow(rv$reference_data), "cell lines"))
    }
  })
  
  # -------------------------------------------------------------------------
  # Gene input handling
  # -------------------------------------------------------------------------
  
  # Update gene list from textarea
  observe({
    genes <- unlist(strsplit(input$gene_textarea, "\\s+"))
    genes <- genes[genes != ""]
    rv$gene_list <- genes
  })
  
  # Upload gene file
  observeEvent(input$gene_file, {
    req(input$gene_file)
    genes <- scan(input$gene_file$datapath, what = "character", quiet = TRUE)
    genes <- genes[genes != ""]
    updateTextAreaInput(session, "gene_textarea", value = paste(genes, collapse = "\n"))
  })
  
  # Clear genes
  observeEvent(input$clear_genes, {
    updateTextAreaInput(session, "gene_textarea", value = "")
  })
  
  # Fetch GO term genes
  observeEvent(input$fetch_go, {
    req(input$go_term)
    
    withProgress(message = "Fetching genes from GO term...", {
      result <- fetch_go_genes(input$go_term)
      
      if (result$success) {
        current_genes <- rv$gene_list
        new_genes <- unique(c(current_genes, result$genes))
        updateTextAreaInput(session, "gene_textarea", 
                            value = paste(new_genes, collapse = "\n"))
        
        output$go_status <- renderUI({
          div(class = "status-box status-success",
              paste("✓ Added", result$count, "genes from", input$go_term))
        })
      } else {
        output$go_status <- renderUI({
          div(class = "status-box status-error",
              paste("✗", result$error))
        })
      }
    })
  })
  
  # Load test genes
  observeEvent(input$load_test_genes, {
    test_genes <- c("TP53", "BRCA1", "BRCA2", "MYC", "KRAS", "EGFR", "PTEN",
                    "RB1", "APC", "VHL", "CDKN2A", "NOTCH1", "PIK3CA", "BRAF",
                    "ATM", "ERBB2", "CDK4", "MDM2", "NRAS", "ARID1A")
    updateTextAreaInput(session, "gene_textarea", value = paste(test_genes, collapse = "\n"))
    showNotification("Loaded 20 test genes (cancer-related)", type = "message")
  })

  # Gene count display
  output$gene_count_display <- renderUI({
    n <- length(rv$gene_list)
    if (n > 0) {
      div(
        strong("Genes entered:"),
        span(class = "gene-count", n)
      )
    }
  })

  # Gene validation display - check against reference data
  output$gene_validation_display <- renderUI({
    req(length(rv$gene_list) > 0)
    req(rv$reference_data)

    genes <- toupper(trimws(rv$gene_list))
    available_genes <- colnames(rv$reference_data)

    found <- genes[genes %in% available_genes]
    not_found <- genes[!genes %in% available_genes]

    if (length(not_found) == 0) {
      div(class = "status-box status-success",
          style = "font-size: 0.85em; padding: 10px;",
          paste0("✓ All ", length(found), " genes found in reference data"))
    } else {
      div(
        div(class = "status-box status-warning",
            style = "font-size: 0.85em; padding: 10px;",
            HTML(paste0(
              "<b>", length(found), " found</b>, <b>", length(not_found), " not found</b><br>",
              "<span style='color: #92400e;'>Not found: ", paste(head(not_found, 10), collapse = ", "),
              if(length(not_found) > 10) paste0(" (+", length(not_found) - 10, " more)") else "",
              "</span>"
            ))
        ),
        actionButton("find_synonyms", "Find Synonyms", class = "btn-sm btn-warning",
                     style = "margin-top: 5px;")
      )
    }
  })

  # Find synonyms for genes not found
  observeEvent(input$find_synonyms, {
    req(rv$reference_data)
    req(length(rv$gene_list) > 0)

    genes <- toupper(trimws(rv$gene_list))
    available_genes <- colnames(rv$reference_data)
    not_found <- genes[!genes %in% available_genes]

    if (length(not_found) == 0) {
      showNotification("All genes already found!", type = "message")
      return()
    }

    withProgress(message = "Looking up gene synonyms...", value = 0, {
      replacements <- list()
      n_total <- length(not_found)

      for (i in seq_along(not_found)) {
        gene <- not_found[i]
        incProgress(1/n_total, detail = gene)

        result <- find_gene_synonyms(gene, available_genes)

        if (result$found && length(result$matches) > 0) {
          # Use the first match as replacement
          replacements[[gene]] <- result$matches[1]
        }

        Sys.sleep(0.1)  # Rate limiting
      }
    })

    if (length(replacements) > 0) {
      # Update gene list with replacements
      updated_genes <- genes
      for (original in names(replacements)) {
        updated_genes[updated_genes == original] <- replacements[[original]]
      }

      updateTextAreaInput(session, "gene_textarea", value = paste(updated_genes, collapse = "\n"))

      msg <- paste0("Replaced ", length(replacements), " gene(s) with synonyms:\n",
                    paste(sapply(names(replacements), function(x)
                      paste0(x, " → ", replacements[[x]])), collapse = ", "))
      showNotification(msg, type = "message", duration = 8)
    } else {
      showNotification("No synonyms found for the missing genes", type = "warning")
    }
  })

  # -------------------------------------------------------------------------
  # Run analysis
  # -------------------------------------------------------------------------
  
  observeEvent(input$run_analysis, {
    # Validation
    if (is.null(rv$reference_data)) {
      showNotification("Please load reference data first", type = "error")
      return()
    }
    
    if (length(rv$gene_list) < 2) {
      showNotification("Please enter at least 2 genes", type = "error")
      return()
    }
    
    withProgress(message = "Running correlation analysis...", {
      rv$results <- run_correlation_analysis(
        reference_data = rv$reference_data,
        gene_list = rv$gene_list,
        analysis_mode = input$analysis_mode,
        correlation_cutoff = input$correlation_cutoff
      )
    })
    
    if (rv$results$success) {
      showNotification(
        paste("Found", nrow(rv$results$correlations), "correlations in",
              nrow(rv$results$clusters), "genes"),
        type = "message"
      )
    } else {
      showNotification(rv$results$error, type = "warning")
    }
  })
  
  output$analysis_status <- renderUI({
    if (is.null(rv$results)) return(NULL)
    
    if (rv$results$success) {
      div(class = "status-box status-success",
          paste("✓ Analysis complete\n",
                nrow(rv$results$correlations), "correlations found\n",
                nrow(rv$results$clusters), "genes in network\n",
                length(rv$results$not_found), "genes not found"))
    } else {
      div(class = "status-box status-error",
          paste("✗", rv$results$error))
    }
  })
  
  # -------------------------------------------------------------------------
  # Results display
  # -------------------------------------------------------------------------
  
  # Store network data for updates
  network_data <- reactiveValues(nodes = NULL, edges = NULL, hidden_nodes = character(0),
                                  all_nodes = NULL, all_edges = NULL)

  # Network visualization - only re-render when results change or legend toggle
  output$network_plot <- renderVisNetwork({
    req(rv$results$success)

    # React to legend toggle
    show_legend <- input$show_legend_nodes

    # Isolate display settings so they don't trigger re-render
    font_size <- isolate(input$net_font_size)
    node_size <- isolate(input$net_node_size)
    edge_width <- isolate(input$net_edge_width)
    show_effect <- isolate(input$net_show_gene_effect)
    show_sd <- isolate(input$net_show_sd)

    graph <- rv$results$graph
    correlations <- rv$results$correlations
    matched <- rv$results$matched
    mode <- rv$results$mode
    clusters_dt <- rv$results$clusters
    cutoff <- rv$results$cutoff

    # Create node labels (optionally with gene effect and SD)
    node_names <- V(graph)$name
    if (show_effect) {
      labels <- sapply(node_names, function(gene) {
        effect_row <- clusters_dt[Gene == gene]
        if (nrow(effect_row) > 0 && !is.na(effect_row$Mean_Effect[1])) {
          if (show_sd && !is.na(effect_row$SD_Effect[1])) {
            paste0(gene, " (", effect_row$Mean_Effect[1], " +/- ", effect_row$SD_Effect[1], ")")
          } else {
            paste0(gene, " (", effect_row$Mean_Effect[1], ")")
          }
        } else {
          gene
        }
      })
    } else {
      labels <- node_names
    }

    # Create nodes with all columns needed for both gene and legend nodes
    nodes <- data.frame(
      id = node_names,
      label = labels,
      size = node_size,
      font.size = font_size,
      group = "gene",
      shape = "dot",
      x = NA,
      y = NA,
      physics = TRUE,
      stringsAsFactors = FALSE
    )

    # Add tooltips with gene effect
    nodes$title <- sapply(node_names, function(gene) {
      effect_row <- clusters_dt[Gene == gene]
      if (nrow(effect_row) > 0 && !is.na(effect_row$Mean_Effect[1])) {
        paste0("<b>", gene, "</b><br>Mean Effect: ", effect_row$Mean_Effect[1],
               "<br>SD: ", effect_row$SD_Effect[1], "<br><i>Click to hide</i>")
      } else {
        paste0("<b>", gene, "</b><br><i>Click to hide</i>")
      }
    })

    # Color nodes based on whether they're in original gene list (green theme)
    if (mode == "design") {
      nodes$color <- ifelse(nodes$id %in% matched, "#16a34a", "#86efac")
      nodes$title <- paste0(nodes$title, "<br><i>",
                            ifelse(nodes$id %in% matched, "Input gene", "Correlated gene"), "</i>")
    } else {
      nodes$color <- "#16a34a"
    }

    # Store base correlation for edge width updates
    edges <- data.frame(
      id = paste0(correlations$Gene1, "-", correlations$Gene2),
      from = correlations$Gene1,
      to = correlations$Gene2,
      correlation = correlations$Correlation,
      width = abs(correlations$Correlation) * edge_width * 3,
      color = ifelse(correlations$Correlation > 0, "#3182ce", "#e53e3e"),
      title = paste("r =", round(correlations$Correlation, 3)),
      stringsAsFactors = FALSE
    )

    # Remove duplicate edges (since graph is undirected)
    edges <- edges[!duplicated(t(apply(edges[,c("from", "to")], 1, sort))),]

    # Add legend as movable nodes if enabled
    if (show_legend) {
      # Legend nodes - positioned at bottom right (same columns as gene nodes)
      legend_nodes <- data.frame(
        id = c("legend_title", "legend_pos_left", "legend_pos_right",
               "legend_neg_left", "legend_neg_right",
               "legend_thick_title",
               "legend_thin_left", "legend_thin_right",
               "legend_med_left", "legend_med_right",
               "legend_thick_left", "legend_thick_right",
               "legend_cutoff"),
        label = c(paste0("Correlation (cutoff: ", round(cutoff, 2), ")"),
                  "", "Positive",
                  "", "Negative",
                  "Edge thickness:",
                  "", as.character(round(cutoff, 1)),
                  "", "0.7",
                  "", "1.0",
                  paste0("* = Input gene")),
        size = c(1, 3, 1, 3, 1, 1, 3, 1, 3, 1, 3, 1, 1),
        font.size = c(14, 1, 11, 1, 11, 12, 1, 10, 1, 10, 1, 10, 10),
        group = rep("legend", 13),
        shape = c("text", "dot", "text", "dot", "text", "text",
                  "dot", "text", "dot", "text", "dot", "text", "text"),
        x = c(350, 350, 380, 350, 380, 350, 350, 380, 350, 380, 350, 380, 350),
        y = c(250, 275, 275, 300, 300, 330, 355, 355, 375, 375, 395, 395, 420),
        physics = rep(FALSE, 13),
        title = rep(NA, 13),
        color = c("#ffffff", "#3182ce", "#ffffff",
                  "#e53e3e", "#ffffff",
                  "#ffffff",
                  "#666666", "#ffffff",
                  "#666666", "#ffffff",
                  "#666666", "#ffffff",
                  "#ffffff"),
        stringsAsFactors = FALSE
      )

      # Legend edges (the lines showing positive/negative)
      legend_edges <- data.frame(
        id = c("legend_pos_edge", "legend_neg_edge",
               "legend_thin_edge", "legend_med_edge", "legend_thick_edge"),
        from = c("legend_pos_left", "legend_neg_left",
                 "legend_thin_left", "legend_med_left", "legend_thick_left"),
        to = c("legend_pos_right", "legend_neg_right",
               "legend_thin_right", "legend_med_right", "legend_thick_right"),
        correlation = c(NA, NA, NA, NA, NA),
        width = c(3, 3, 1, 4, 8),
        color = c("#3182ce", "#e53e3e", "#666666", "#666666", "#666666"),
        title = c("Positive correlation", "Negative correlation",
                  paste0("r = ", round(cutoff, 1)), "r = 0.7", "r = 1.0"),
        stringsAsFactors = FALSE
      )

      # Ensure columns match before rbind
      nodes <- rbind(nodes[, c("id", "label", "size", "font.size", "group", "shape", "x", "y", "physics", "title", "color")],
                     legend_nodes[, c("id", "label", "size", "font.size", "group", "shape", "x", "y", "physics", "title", "color")])
      edges <- rbind(edges, legend_edges)
    }

    # Store all nodes/edges for restoration
    network_data$all_nodes <- nodes
    network_data$all_edges <- edges
    network_data$nodes <- nodes
    network_data$edges <- edges
    network_data$clusters_dt <- clusters_dt
    network_data$hidden_nodes <- character(0)

    visNetwork(nodes, edges) %>%
      visNodes(borderWidth = 2,
               borderWidthSelected = 4,
               font = list(face = "arial", color = "#333333"),
               color = list(border = "#ffffff")) %>%
      visEdges(smooth = FALSE) %>%
      visGroups(groupname = "legend", physics = FALSE) %>%
      visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
                 nodesIdSelection = FALSE) %>%
      visPhysics(enabled = FALSE,
                 solver = "forceAtlas2Based",
                 forceAtlas2Based = list(gravitationalConstant = -50),
                 stabilization = list(iterations = 200)) %>%
      visInteraction(navigationButtons = TRUE,
                     dragNodes = TRUE,
                     dragView = TRUE,
                     zoomView = TRUE) %>%
      visEvents(click = "function(params) {
        if (params.nodes.length > 0) {
          var nodeId = params.nodes[0];
          if (!nodeId.startsWith('legend_')) {
            Shiny.setInputValue('clicked_node', nodeId, {priority: 'event'});
          }
        }
      }")
  })

  # Update font size without re-rendering
  observeEvent(input$net_font_size, {
    req(network_data$nodes)
    nodes <- network_data$nodes
    nodes$font.size <- input$net_font_size
    visNetworkProxy("network_plot") %>%
      visUpdateNodes(nodes)
  }, ignoreInit = TRUE)

  # Update node size without re-rendering
  observeEvent(input$net_node_size, {
    req(network_data$nodes)
    nodes <- network_data$nodes
    nodes$size <- input$net_node_size
    network_data$nodes <- nodes
    visNetworkProxy("network_plot") %>%
      visUpdateNodes(nodes)
  }, ignoreInit = TRUE)

  # Update edge width without re-rendering
  observeEvent(input$net_edge_width, {
    req(network_data$edges)
    edges <- network_data$edges
    edges$width <- abs(edges$correlation) * input$net_edge_width * 3
    network_data$edges <- edges
    visNetworkProxy("network_plot") %>%
      visUpdateEdges(edges)
  }, ignoreInit = TRUE)

  # Update labels with gene effect without re-rendering
  # Update labels when gene effect toggle changes
  observeEvent(input$net_show_gene_effect, {
    req(network_data$nodes)
    req(rv$results$success)

    nodes <- network_data$nodes
    clusters_dt <- network_data$clusters_dt
    show_sd <- input$net_show_sd

    # Only update gene nodes (not legend nodes)
    gene_nodes <- nodes$group == "gene"

    if (input$net_show_gene_effect) {
      nodes$label[gene_nodes] <- sapply(nodes$id[gene_nodes], function(gene) {
        effect_row <- clusters_dt[Gene == gene]
        if (nrow(effect_row) > 0 && !is.na(effect_row$Mean_Effect[1])) {
          if (show_sd && !is.na(effect_row$SD_Effect[1])) {
            paste0(gene, " (", effect_row$Mean_Effect[1], " +/- ", effect_row$SD_Effect[1], ")")
          } else {
            paste0(gene, " (", effect_row$Mean_Effect[1], ")")
          }
        } else {
          gene
        }
      })
    } else {
      nodes$label[gene_nodes] <- nodes$id[gene_nodes]
    }

    network_data$nodes <- nodes
    visNetworkProxy("network_plot") %>%
      visUpdateNodes(nodes)
  }, ignoreInit = TRUE)

  # Update labels when SD toggle changes
  observeEvent(input$net_show_sd, {
    req(network_data$nodes)
    req(rv$results$success)
    req(input$net_show_gene_effect)  # Only update if gene effect is shown

    nodes <- network_data$nodes
    clusters_dt <- network_data$clusters_dt
    show_sd <- input$net_show_sd

    # Only update gene nodes (not legend nodes)
    gene_nodes <- nodes$group == "gene"

    nodes$label[gene_nodes] <- sapply(nodes$id[gene_nodes], function(gene) {
      effect_row <- clusters_dt[Gene == gene]
      if (nrow(effect_row) > 0 && !is.na(effect_row$Mean_Effect[1])) {
        if (show_sd && !is.na(effect_row$SD_Effect[1])) {
          paste0(gene, " (", effect_row$Mean_Effect[1], " +/- ", effect_row$SD_Effect[1], ")")
        } else {
          paste0(gene, " (", effect_row$Mean_Effect[1], ")")
        }
      } else {
        gene
      }
    })

    network_data$nodes <- nodes
    visNetworkProxy("network_plot") %>%
      visUpdateNodes(nodes)
  }, ignoreInit = TRUE)

  # Export handlers - legend is now part of the network
  observeEvent(input$export_png, {
    session$sendCustomMessage("exportNetworkPNG", list())
  })

  observeEvent(input$export_svg, {
    session$sendCustomMessage("exportNetworkSVG", list())
  })

  observeEvent(input$export_pdf, {
    session$sendCustomMessage("exportNetworkPNG", list())
    showNotification("For PDF: Open the PNG in Preview and export as PDF", type = "message")
  })

  # Stabilize network
  observeEvent(input$net_stabilize, {
    visNetworkProxy("network_plot") %>%
      visStabilize()
  })

  # Fit to view
  observeEvent(input$net_fit, {
    visNetworkProxy("network_plot") %>%
      visFit()
  })

  # Toggle physics
  observeEvent(input$net_physics, {
    visNetworkProxy("network_plot") %>%
      visPhysics(enabled = input$net_physics)
  })

  # Hide node when clicked
  observeEvent(input$clicked_node, {
    req(input$clicked_node)
    node_id <- input$clicked_node

    # Add to hidden nodes
    network_data$hidden_nodes <- unique(c(network_data$hidden_nodes, node_id))

    # Remove the node from current display
    current_nodes <- network_data$nodes
    current_edges <- network_data$edges

    # Filter out the hidden node
    current_nodes <- current_nodes[!current_nodes$id %in% network_data$hidden_nodes, ]

    # Filter out edges connected to hidden nodes
    current_edges <- current_edges[
      !current_edges$from %in% network_data$hidden_nodes &
      !current_edges$to %in% network_data$hidden_nodes, ]

    network_data$nodes <- current_nodes
    network_data$edges <- current_edges

    # Update the network
    visNetworkProxy("network_plot") %>%
      visRemoveNodes(id = node_id)

    showNotification(paste("Hidden:", node_id), type = "message", duration = 2)
  })

  # Show all hidden nodes
  observeEvent(input$show_all_nodes, {
    req(length(network_data$hidden_nodes) > 0)

    # Restore all nodes and edges
    network_data$nodes <- network_data$all_nodes
    network_data$edges <- network_data$all_edges

    hidden_count <- length(network_data$hidden_nodes)
    network_data$hidden_nodes <- character(0)

    # Re-render the network with all nodes
    visNetworkProxy("network_plot") %>%
      visUpdateNodes(network_data$all_nodes) %>%
      visUpdateEdges(network_data$all_edges)

    showNotification(paste("Restored", hidden_count, "hidden nodes"), type = "message")
  })

  # Display hidden nodes info
  output$hidden_nodes_info <- renderUI({
    hidden <- network_data$hidden_nodes
    if (length(hidden) > 0) {
      div(style = "background: #fef3c7; padding: 8px 12px; border-radius: 6px; margin-bottom: 10px;",
        strong(paste(length(hidden), "hidden node(s): ")),
        span(paste(hidden, collapse = ", ")),
        actionLink("clear_one_hidden", " [Show all]", style = "margin-left: 10px;")
      )
    }
  })

  # Correlations table
  output$correlations_table <- renderDT({
    req(rv$results$success)
    datatable(rv$results$correlations,
              options = list(pageLength = 15, scrollX = TRUE),
              filter = "top")
  })
  
  # Clusters table
  output$clusters_table <- renderDT({
    req(rv$results$success)
    datatable(rv$results$clusters,
              options = list(pageLength = 15, scrollX = TRUE),
              filter = "top")
  })
  
  # Summary text
  output$summary_text <- renderText({
    req(rv$results)

    if (!rv$results$success) {
      return(paste("Analysis failed:", rv$results$error))
    }

    paste(
      "=== Analysis Summary ===",
      paste("Date:", Sys.time()),
      paste("Mode:", ifelse(rv$results$mode == "analysis",
                            "Analysis (within gene list)",
                            "Design (find correlated genes)")),
      paste("Correlation cutoff:", rv$results$cutoff),
      "",
      paste("Genes matched:", length(rv$results$matched)),
      paste("Genes not found:", length(rv$results$not_found)),
      "",
      paste("Correlations found:", nrow(rv$results$correlations)),
      paste("Genes in network:", nrow(rv$results$clusters)),
      paste("Number of clusters:", max(rv$results$clusters$Cluster)),
      "",
      "=== Genes Not Found ===",
      if (length(rv$results$not_found) > 0) {
        paste(rv$results$not_found, collapse = "\n")
      } else {
        "All genes were found"
      },
      sep = "\n"
    )
  })

  # -------------------------------------------------------------------------
  # Downloads
  # -------------------------------------------------------------------------
  
  output$download_correlations <- downloadHandler(
    filename = function() {
      paste0("correlations_", format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      req(rv$results$success)
      fwrite(rv$results$correlations, file)
    }
  )
  
  output$download_clusters <- downloadHandler(
    filename = function() {
      paste0("clusters_", format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      req(rv$results$success)
      fwrite(rv$results$clusters, file)
    }
  )

  # -------------------------------------------------------------------------
  # Enrichment Analysis
  # -------------------------------------------------------------------------

  # Run enrichment analysis
  observeEvent(input$run_enrichment, {
    if (length(rv$gene_list) < 2) {
      showNotification("Please enter at least 2 genes first", type = "error")
      return()
    }

    if (length(input$enrichment_dbs) == 0) {
      showNotification("Please select at least one database", type = "error")
      return()
    }

    withProgress(message = "Running enrichment analysis...", {
      rv$enrichment_results <- run_enrichment_analysis(
        rv$gene_list,
        databases = input$enrichment_dbs
      )
    })

    if (rv$enrichment_results$success) {
      showNotification(
        paste("Found", nrow(rv$enrichment_results$results), "enriched terms"),
        type = "message"
      )
    } else {
      showNotification(rv$enrichment_results$error, type = "warning")
    }
  })

  # Enrichment status
  output$enrichment_status <- renderUI({
    if (is.null(rv$enrichment_results)) return(NULL)

    if (rv$enrichment_results$success) {
      div(class = "status-box status-success",
          paste("Found", nrow(rv$enrichment_results$results), "enriched terms"))
    } else {
      div(class = "status-box status-error",
          paste("Error:", rv$enrichment_results$error))
    }
  })

  # Output flag for conditional panel
  output$enrichment_ready <- reactive({
    !is.null(rv$enrichment_results) && rv$enrichment_results$success
  })
  outputOptions(output, "enrichment_ready", suspendWhenHidden = FALSE)

  # Enrichment results table
  output$enrichment_table <- renderDT({
    req(rv$enrichment_results$success)

    dt <- rv$enrichment_results$results[, .(
      Database,
      Term,
      `P-value` = signif(P_value, 3),
      `Adj P-value` = signif(Adj_P_value, 3),
      `Gene Count` = Gene_Count,
      Genes
    )]

    datatable(
      dt,
      selection = "single",
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        order = list(list(3, "asc"))  # Order by Adj P-value
      ),
      filter = "top",
      rownames = FALSE
    ) %>%
      formatStyle("Adj P-value",
                  backgroundColor = styleInterval(
                    c(0.01, 0.05),
                    c("#dcfce7", "#fef9c3", "white")
                  ))
  })

  # Handle row click to load pathway genes
  observeEvent(input$enrichment_table_rows_selected, {
    req(input$enrichment_table_rows_selected)
    req(rv$enrichment_results$success)

    selected_row <- input$enrichment_table_rows_selected
    selected_term <- rv$enrichment_results$results[selected_row]

    # Get genes from the selected pathway
    pathway_genes <- unlist(strsplit(selected_term$Genes, ";"))
    pathway_genes <- trimws(pathway_genes)

    # Update the gene textarea with pathway genes
    updateTextAreaInput(session, "gene_textarea",
                        value = paste(pathway_genes, collapse = "\n"))

    showNotification(
      paste("Loaded", length(pathway_genes), "genes from:",
            substr(selected_term$Term, 1, 50)),
      type = "message"
    )

    # Switch to the gene input tab
    updateTabsetPanel(session, "input_method", selected = "Paste/Type")
  })

  # Download enrichment results
  output$download_enrichment <- downloadHandler(
    filename = function() {
      paste0("enrichment_", format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      req(rv$enrichment_results$success)
      fwrite(rv$enrichment_results$results, file)
    }
  )
}

# ============================================================================
# RUN APP
# ============================================================================

shinyApp(ui = ui, server = server)
