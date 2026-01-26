# Gene Correlation Explorer - v81
# Standalone Shiny application for analyzing DepMap CRISPR screen correlations
# 
# Author: Fredrik Wermeling / Wermeling Lab
# URL: https://wermelinglab.com | https://greenlisted.cmm.se
# Repository: https://github.com/fredrikwermeling/correlation-app
#
# Optimizations in v81:
# - Vectorized slope/N calculations for better performance
# - Cleaned up duplicate/dead code blocks
# - Added version constant
# - Consolidated helper functions

APP_VERSION <- "81"

# ============================================================================
# SETUP - Install packages if needed
# ============================================================================

required_packages <- c("shiny", "data.table", "igraph", "ggplot2", "ggraph",
                       "ggrepel", "DT", "shinyjs", "httr", "jsonlite", "visNetwork", "dplyr", "htmltools", "base64enc", "cowplot", "tidyr")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Increase max upload size to 500MB
options(shiny.maxRequestSize = 500*1024^2)

# Force use of period as decimal separator for consistent display
options(OutDec = ".")

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
library(dplyr)
library(htmltools)
library(base64enc)
library(cowplot)
library(tidyr)

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================


sanitize_limits <- function(lim) {
  if (is.null(lim) || length(lim) != 2) return(NULL)
  lim <- suppressWarnings(as.numeric(lim))
  if (any(is.na(lim)) || any(!is.finite(lim))) return(NULL)
  lim
}

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

    # In design mode, correlate the selected genes against ALL numeric gene-effect columns.
    # This avoids including non-numeric metadata columns (e.g., DepMap_ID) which can break downstream steps.
    reference_gene_cols <- names(reference_data)[vapply(reference_data, is.numeric, logical(1))]
    if (length(reference_gene_cols) == 0) {
      return(list(
        success = FALSE,
        error = "No numeric gene-effect columns found in the reference data.",
        not_found = not_found_genes
      ))
    }

    cor_matrix <- cor(reference_data[, ..reference_gene_cols], filtered_data, use = "pairwise.complete.obs")

    all_genes <- rownames(cor_matrix)
    matched_genes <- colnames(cor_matrix)
    output_table <- data.table(
      Gene1 = rep(matched_genes, each = length(all_genes)),
      Gene2 = rep(all_genes, times = length(matched_genes)),
      Correlation = as.vector(cor_matrix)
    )
  }
  
  # Filter by cutoff and remove duplicates (keep only Gene1 < Gene2)
  filtered_correlations <- output_table[
    abs(Correlation) >= correlation_cutoff & (Gene1 != Gene2)
  ]
  
  # Remove duplicate pairs (a-b and b-a are the same)
  filtered_correlations <- filtered_correlations[Gene1 < Gene2]
  
  # Check if we have any correlations after filtering
  if (nrow(filtered_correlations) == 0) {
    return(list(
      success = FALSE,
      error = paste("No correlations found above cutoff of", correlation_cutoff),
      not_found = not_found_genes,
      matched = matched_columns
    ))
  }
  
  filtered_correlations[, Correlation := round(Correlation, 3)]
  
  # OPTIMIZED: Calculate slope and N using vectorized approach
  gene1_vec <- filtered_correlations$Gene1
  gene2_vec <- filtered_correlations$Gene2
  n_pairs <- length(gene1_vec)
  
  slopes <- numeric(n_pairs)
  ns <- integer(n_pairs)
  
  for (i in seq_len(n_pairs)) {
    g1 <- gene1_vec[i]
    g2 <- gene2_vec[i]
    
    x <- reference_data[[g1]]
    y <- reference_data[[g2]]
    
    if (is.null(x) || is.null(y)) {
      slopes[i] <- NA_real_
      ns[i] <- NA_integer_
      next
    }
    
    valid <- complete.cases(x, y)
    n_valid <- sum(valid)
    ns[i] <- n_valid
    
    if (n_valid < 3) {
      slopes[i] <- NA_real_
    } else {
      slopes[i] <- round(coef(lm(y[valid] ~ x[valid]))[2], 3)
    }
  }
  
  filtered_correlations[, Slope := slopes]
  filtered_correlations[, N := ns]

  # Replace missing slopes/N with 0 to avoid downstream filtering issues
  filtered_correlations[is.na(Slope), Slope := 0]
  filtered_correlations[is.na(N), N := 0L]

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
    filtered_correlations[, InGeneList_Gene1 := fifelse(Gene1 %in% matched_columns, "*", "")]
    filtered_correlations[, InGeneList_Gene2 := fifelse(Gene2 %in% matched_columns, "*", "")]
  }
  
  # Build clusters table
  found_genes <- V(graph)$name
  found_clusters <- cluster_membership[found_genes]
  
  if (analysis_mode == "design") {
    cluster_dt <- data.table(
      Gene = found_genes,
      Cluster = found_clusters,
      In_GeneList = fifelse(found_genes %in% matched_columns, "*", "")
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

# Find gene synonyms using ortholog lookup (if available) then MyGene.info API
find_gene_synonyms <- function(gene_symbol, available_genes, ortholog_lookup = NULL) {
  gene_symbol <- toupper(trimws(gene_symbol))

  # STEP 1: Check ortholog lookup table first (fastest and most reliable)
  if (!is.null(ortholog_lookup) && nrow(ortholog_lookup) > 0) {
    # Look for this gene as a mouse symbol
    orthologs <- ortholog_lookup[mouse_symbol == gene_symbol, human_symbol]
    if (length(orthologs) > 0) {
      matches <- intersect(orthologs, available_genes)
      if (length(matches) > 0) {
        return(list(
          found = TRUE,
          original = gene_symbol,
          matches = matches,
          source = "ortholog"
        ))
      }
    }
  }
  
  # STEP 2: Simple case conversion (Mouse Abc1 -> Human ABC1)
  # Many mouse genes are just titlecase versions of human genes
  if (gene_symbol %in% available_genes) {
    return(list(
      found = TRUE,
      original = gene_symbol,
      matches = gene_symbol,
      source = "direct"
    ))
  }

  # STEP 3: Fall back to MyGene.info API
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
            matches = matches,
            source = "mygene"
          ))
        }
      }
    }

    return(list(found = FALSE, original = gene_symbol, matches = character(0), source = NA))

  }, error = function(e) {
    return(list(found = FALSE, original = gene_symbol, matches = character(0), source = NA))
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
    tags$script(src = "https://cdn.jsdelivr.net/npm/html2canvas@1.4.1/dist/html2canvas.min.js"),
    # Force period decimal separator in number inputs
    tags$script(HTML("
      $(document).on('shiny:sessioninitialized', function() {
        // Override toLocaleString for consistent decimal formatting
        Number.prototype.toLocaleString = function() {
          return this.toString();
        };
        
        // Force period decimals in number inputs
        $(document).on('input change', 'input[type=number]', function() {
          var val = $(this).val();
          if (val && val.includes(',')) {
            $(this).val(val.replace(',', '.'));
          }
        });
      });
    ")),
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
      p(paste0("Analyze gene correlations from DepMap CRISPR screen data (v", APP_VERSION, ")"))
  ),
  
  fluidRow(
    # Left panel - Inputs
    column(4,
      # Reference data card
      div(class = "card",
        div(class = "card-title", "1. Load Reference Data"),

        # File upload - no tabs needed
        br(),
        fileInput("reference_file",
                  "1. Upload CRISPRGeneEffect.csv",
                  accept = ".csv"),
        helpText("Required. Max file size: 500MB"),
        br(),
        fileInput("metadata_file",
                  "2. Upload Model.csv (optional, for searchable cell names)",
                  accept = ".csv"),
        helpText("Without this file, you can only search by ACH IDs (e.g., ACH-000001).",
                "With this file, you can search by cell line names (e.g., HELA, A549)."),
        br(),
        fileInput("hotspot_file",
                  "3. Upload OmicsSomaticMutationsMatrixHotspot.csv (optional, hotspot mutations)",
                  accept = ".csv"),
        helpText("Enables hotspot mutation overlays (0/1/2) in Inspect view and stratified correlations."),
        br(),
        fileInput("ortholog_file",
                  "4. Upload ortholog_lookup_best.csv (optional, mouse-human mapping)",
                  accept = ".csv"),
        helpText("Improves mouse→human gene mapping. Generate using create_ortholog_lookup.R script."),
        br(),
        helpText("Download files from: ",
                 tags$a(href = "https://depmap.org/portal/data_page/?tab=allData",
                        "DepMap Portal", target = "_blank"),
                 br(),
                 "Look for 'CRISPRGeneEffect.csv' (required), 'Model.csv' (optional), and 'OmicsSomaticMutationsMatrixHotspot.csv' (optional)"),
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

          tabPanel("Upload with Stats",
            br(),
            fileInput("gene_stats_file", "Upload CSV/TSV with gene statistics",
                      accept = c(".csv", ".tsv", ".txt")),
            helpText("File can contain: (1) only gene names, OR (2) gene names plus optional LFC and/or statistics columns"),
            downloadButton("download_example_stats", "Download Example File", class = "btn-sm btn-info"),
            br(), br(),
            conditionalPanel(
              condition = "output.stats_file_uploaded",
              selectInput("stats_gene_col", "Gene column:", choices = NULL),
              selectInput("stats_lfc_col", "Log fold change (LFC) column (optional):", choices = NULL),
              selectInput("stats_stat_col", "Statistics column (optional, FDR/p-value):", choices = NULL),
              actionButton("load_stats", "Load Gene Stats", class = "btn-success btn-sm"),
              br(), br(),
              uiOutput("stats_status")
            )
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
        
        numericInput("min_cell_lines", "Minimum Cell Lines per Correlation:",
                    value = 10, min = 3, max = 50, step = 1),
        
        sliderInput("min_slope", "Minimum Absolute Slope:",
                    min = 0, max = 2, value = 0.1, step = 0.05),
        
        helpText("Higher correlation cutoff = fewer, stronger correlations. Minimum N filters unreliable correlations. Minimum slope filters flat relationships."),
        
        # Oncotree filter section
        conditionalPanel(
          condition = "output.metadata_loaded == true",
          hr(),
          div(style = "background: #f0fdf4; padding: 10px; border-radius: 6px; margin-bottom: 10px;",
            strong(icon("filter"), " Cell Line Filter (optional)"),
            helpText("Limit analysis to cell lines from specific tissue types/cancer lineages. Note: The actual N in scatter plots may be lower than the filter count because some cell lines have missing data (NA) for specific genes.")
          ),
          selectInput("oncotree_lineage_filter", 
                     "Filter by Lineage:",
                     choices = c("All lineages" = ""),
                     selected = "",
                     multiple = TRUE),
          conditionalPanel(
            condition = "input.oncotree_lineage_filter != null && input.oncotree_lineage_filter.length > 0",
            selectInput("oncotree_disease_filter",
                       "Filter by Disease (optional):",
                       choices = c("All diseases" = ""),
                       selected = "",
                       multiple = TRUE)
          ),
          uiOutput("filter_status")
        ),
        
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
            # Controls - reorganized for cleaner layout
            div(style = "padding: 8px 0; margin-bottom: 0; font-size: 12px;",
              # Row 1: Font Size, Node Size, Edge Width (sliders)
              fluidRow(
                column(4, style = "padding-right: 10px;",
                  sliderInput("net_font_size", "Font Size:", min = 10, max = 40, value = 16, step = 2)
                ),
                column(4, style = "padding: 0 10px;",
                  sliderInput("net_node_size", "Node Size:", min = 10, max = 60, value = 25, step = 5)
                ),
                column(4, style = "padding-left: 10px;",
                  sliderInput("net_edge_width", "Edge Width:", min = 1, max = 10, value = 3, step = 1)
                )
              ),
              # Row 2: Gene Statistics (conditional) + Gene Effect options + View/Export controls
              fluidRow(
                column(5, style = "padding-right: 10px;",
                  conditionalPanel(
                    condition = "output.gene_stats_loaded",
                    div(style = "background: #f8fafc; padding: 8px; border-radius: 6px; border: 1px solid #e2e8f0;",
                      fluidRow(
                        column(6,
                          tags$b("Gene Statistics:", style = "font-size: 11px; display: block; margin-bottom: 4px;"),
                          radioButtons("net_stats_display", NULL, inline = TRUE,
                                      choices = list("None" = "none", "LFC" = "lfc", "FDR" = "fdr"),
                                      selected = "none")
                        ),
                        column(6,
                          checkboxInput("net_color_by_stats", "Color by stats", value = FALSE),
                          conditionalPanel(
                            condition = "input.net_color_by_stats == true",
                            radioButtons("net_color_stat", NULL, inline = TRUE,
                                        choices = list("Signed LFC" = "signed_lfc", "|LFC|" = "abs_lfc", "FDR" = "fdr"),
                                        selected = "signed_lfc"),
                            radioButtons("net_color_scale_scope", "Scale:", inline = TRUE,
                             choices = c("All genes" = "all", "Network only" = "network"),
                             selected = "all")
                          )
                        )
                      ),
                      conditionalPanel(
                        condition = "input.net_color_by_stats == true",
                        uiOutput("net_color_legend")
                      )
                    )
                  )
                ),
                column(3, style = "padding: 0 10px;",
                  div(style = "background: #f8fafc; padding: 8px; border-radius: 6px; border: 1px solid #e2e8f0;",
                    checkboxInput("net_show_gene_effect", "Show gene effect (GE)", value = FALSE),
                    conditionalPanel(
                      condition = "input.net_show_gene_effect == true",
                      checkboxInput("net_show_sd", "Include ± SD", value = FALSE)
                    )
                  )
                ),
                # Right column: View controls + Export (stacked)
                column(4, style = "padding-left: 10px;",
                  div(style = "background: #f8fafc; padding: 10px; border-radius: 6px; border: 1px solid #e2e8f0;",
                    # View Controls section
                    tags$b("View Controls:", style = "font-size: 11px; display: block; margin-bottom: 6px; color: #374151;"),
                    div(style = "display: flex; gap: 6px; margin-bottom: 10px;",
                      actionButton("net_fit", "Fit to View", class = "btn-sm btn-outline-secondary", style = "flex: 1; padding: 4px 8px;"),
                      actionButton("show_all_nodes", "Show Hidden", class = "btn-sm btn-outline-secondary", style = "flex: 1; padding: 4px 8px;")
                    ),
                    tags$small("Double-click node to hide", style = "color: #6b7280; font-size: 9px; display: block; margin-bottom: 8px;"),
                    
                    # Export Image section
                    tags$b("Export Image + Legend:", style = "font-size: 11px; display: block; margin-bottom: 6px; color: #374151;"),
                    div(style = "display: flex; gap: 6px;",
                      actionButton("export_network_png", "PNG", class = "btn-success btn-sm", style = "flex: 1; padding: 4px 8px;"),
                      actionButton("export_network_svg", "SVG", class = "btn-success btn-sm", style = "flex: 1; padding: 4px 8px;")
                    )
                  )
                )
              ),
              # Row 3: Download All (separate, for data export)
              fluidRow(style = "margin-top: 10px;",
                column(4, offset = 8,
                  downloadButton("download_all", "Download All Data", class = "btn-outline-success btn-sm", style = "width: 100%; padding: 5px 10px;"),
                  tags$small("ZIP: correlation tables, clusters, summary, legend", style = "color: #6b7280; font-size: 9px; display: block; margin-top: 2px;")
                )
              )
            ),
            # Hidden nodes info
            fluidRow(
              column(12,
                uiOutput("hidden_nodes_info")
              )
            ),
            # Network + legend (legend always shown below)
            div(id = "network_export_container",
              div(style = "background: white; border: 1px solid #ddd; border-radius: 8px;",
                visNetworkOutput("network_plot", height = "600px")
              ),
              # Static legend below network (always visible)
              br(),
              uiOutput("network_legend_static")
            ),
            # JavaScript for exports with white background (includes legend)
            tags$script(HTML("
              // Hide visNetwork navigation buttons for clean export
              function hideNavButtons() {
                var navBtns = document.querySelectorAll('.vis-navigation, .vis-button');
                navBtns.forEach(function(btn) { btn.style.display = 'none'; });
              }
              function showNavButtons() {
                var navBtns = document.querySelectorAll('.vis-navigation, .vis-button');
                navBtns.forEach(function(btn) { btn.style.display = ''; });
              }

              Shiny.addCustomMessageHandler('exportNetworkPNG', function(message) {
                hideNavButtons();
                setTimeout(function() {
                  var container = document.getElementById('network_export_container');
                  if (container) {
                    html2canvas(container, {
                      backgroundColor: '#FFFFFF',
                      scale: 2
                    }).then(function(canvas) {
                      var link = document.createElement('a');
                      link.download = 'correlation_network_' + new Date().toISOString().slice(0,10) + '.png';
                      link.href = canvas.toDataURL('image/png');
                      try { Shiny.setInputValue('last_network_png', link.href, {priority: 'event'}); } catch(e) {}
                      link.click();
                      showNavButtons();
                    }).catch(function(err) {
                      console.error('Export failed:', err);
                      showNavButtons();
                    });
                  } else {
                    showNavButtons();
                  }
                }, 100);
              });

              Shiny.addCustomMessageHandler('exportNetworkSVG', function(message) {
                hideNavButtons();
                setTimeout(function() {
                  var container = document.getElementById('network_export_container');
                  if (container) {
                    html2canvas(container, {
                      backgroundColor: '#FFFFFF',
                      scale: 2
                    }).then(function(canvas) {
                      var svgData = '<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"' + canvas.width + '\" height=\"' + canvas.height + '\">' +
                        '<rect width=\"100%\" height=\"100%\" fill=\"white\"/>' +
                        '<image href=\"' + canvas.toDataURL('image/png') + '\" width=\"' + canvas.width + '\" height=\"' + canvas.height + '\"/>' +
                        '</svg>';
                      var link = document.createElement('a');
                      link.download = 'correlation_network_' + new Date().toISOString().slice(0,10) + '.svg';
                      link.href = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(svgData);
                      link.click();
                      showNavButtons();
                    }).catch(function(err) {
                      console.error('Export failed:', err);
                      showNavButtons();
                    });
                  } else {
                    showNavButtons();
                  }
                }, 100);
              });
            "))
          ),

          tabPanel("Correlations",
            br(),
            div(style = "background: #dcfce7; padding: 10px 15px; border-radius: 6px; margin-bottom: 15px; border-left: 4px solid #16a34a;",
              icon("info-circle"),
              strong(" Tip: "),
              "Click the ", tags$strong("Inspect"), " button to view a detailed scatter plot showing all cell line data points. You can adjust axis ranges, highlight specific cell lines, and download as PNG, SVG, or Excel."
            ),
            downloadButton("download_correlations", "Download CSV"),
            br(), br(),
            DTOutput("correlations_table")
          ),

          tabPanel("Clusters",
            br(),
            helpText("Mean_Effect and SD_Effect are from DepMap CRISPR data. LFC and FDR are from your uploaded dataset (if provided)."),
            downloadButton("download_clusters", "Download CSV"),
            br(), br(),
            DTOutput("clusters_table")
          ),

          tabPanel("Summary",
            br(),
            downloadButton("download_summary", "Download Summary"),            br(), br(),
            verbatimTextOutput("summary_text")
          )

        )
      )
    )
  ),

  # Footer
  div(class = "footer",
    p("Gene Correlation Explorer v", APP_VERSION, " | Part of the ",
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

  # ---- Helper: map statistics to node colors + legend ----
  map_abs_lfc_color <- function(x, max_x) {
    if (is.na(x)) return("#d1d5db")  # light gray for missing
    if (is.na(max_x) || max_x <= 0) max_x <- 1
    t <- min(max(abs(x) / max_x, 0), 1)
    # Standard sequential gradient: light gray -> orange -> red
    pal <- grDevices::colorRamp(c("#f5f5f5", "#fdae61", "#d7191c"))
    rgbv <- pal(t) / 255
    grDevices::rgb(rgbv[1], rgbv[2], rgbv[3])
  }

  

  
map_signed_lfc_color <- function(x, min_x, max_x) {
  # Blue -> white -> red, scaled by the observed min/max LFC (not symmetric by abs()).
  if (is.na(x)) return("#d1d5db")  # light gray for missing

  if (!is.finite(min_x) || is.na(min_x)) min_x <- x
  if (!is.finite(max_x) || is.na(max_x)) max_x <- x
  if (!is.finite(min_x)) min_x <- -1
  if (!is.finite(max_x)) max_x <- 1

  # Handle degenerate case
  if (max_x == min_x) {
    return("#f7f7f7")
  }

  # Piecewise scaling so 0 maps to white even when range is asymmetric.
  if (min_x < 0 && max_x > 0) {
    if (x <= 0) {
      u <- 0.5 * (x - min_x) / (0 - min_x)   # [0..0.5]
    } else {
      u <- 0.5 + 0.5 * (x - 0) / (max_x - 0) # [0.5..1]
    }
    u <- min(max(u, 0), 1)
    pal <- grDevices::colorRamp(c("#2166ac", "#f7f7f7", "#b2182b"))
    rgbv <- pal(u) / 255
    return(grDevices::rgb(rgbv[1], rgbv[2], rgbv[3]))
  }

  # If all positive: white -> red (lowest -> highest)
  if (min_x >= 0 && max_x > 0) {
    u <- (x - min_x) / (max_x - min_x)
    u <- min(max(u, 0), 1)
    pal <- grDevices::colorRamp(c("#f7f7f7", "#b2182b"))
    rgbv <- pal(u) / 255
    return(grDevices::rgb(rgbv[1], rgbv[2], rgbv[3]))
  }

  # If all negative: blue -> white (lowest (most negative) -> highest (closest to 0))
  u <- (x - min_x) / (max_x - min_x)
  u <- min(max(u, 0), 1)
  pal <- grDevices::colorRamp(c("#2166ac", "#f7f7f7"))
  rgbv <- pal(u) / 255
  grDevices::rgb(rgbv[1], rgbv[2], rgbv[3])
}
map_fdr_color <- function(fdr, min_fdr) {
    if (is.na(fdr)) return("#d1d5db")
    if (!is.finite(fdr) || fdr <= 0) fdr <- .Machine$double.xmin
    if (is.na(min_fdr) || !is.finite(min_fdr) || min_fdr <= 0) min_fdr <- .Machine$double.xmin
    # Use -log10(FDR) for intensity, scaled to the most significant value in the stats table
    logv <- -log10(fdr)
    maxlog <- -log10(min_fdr)
    if (!is.finite(maxlog) || maxlog <= 0) maxlog <- 1
    t <- min(max(logv / maxlog, 0), 1)
    # Standard sequential gradient: light gray -> light blue -> dark blue
    pal <- grDevices::colorRamp(c("#f5f5f5", "#9ecae1", "#08519c"))
    rgbv <- pal(t) / 255
    grDevices::rgb(rgbv[1], rgbv[2], rgbv[3])
  }

  legend_bar_html <- function(title, left_lab, right_lab, colors) {
    # colors: vector of 2 or 3 hex colors used in CSS linear-gradient
    grad <- paste(colors, collapse = ", ")
    shiny::HTML(paste0(
      "<div style='padding:4px 0;'>",
      "<div style='font-size:11px; font-weight:600; margin-bottom:4px;'>", title, "</div>",
      "<div style='height:10px; width:180px; border-radius:6px; background: linear-gradient(90deg, ", grad, "); border:1px solid #e5e7eb;'></div>",
      "<div style='display:flex; justify-content:space-between; width:180px; font-size:10px; color:#374151; margin-top:2px;'>",
      "<span>", left_lab, "</span><span>", right_lab, "</span>",
      "</div>",
      "<div style='font-size:10px; color:#6b7280; margin-top:2px;'>Missing: <span style='display:inline-block; width:8px; height:8px; background:#d1d5db; border:1px solid #e5e7eb; vertical-align:middle; border-radius:2px;'></span></div>",
      "</div>"
    ))
  }

  
  # Reactive values to store data
  rv <- reactiveValues(
    synonym_map = NULL,
    reverse_syn_map = NULL,

    reference_data = NULL,
    cell_line_names = NULL,  # Store cell line names (extracted or from metadata)
    cell_line_ids = NULL,  # Store original cell line IDs (ACH-xxx)
    cell_line_name_lookup = NULL,  # Lookup table from metadata file (ACH -> StrippedCellLineName)
    cell_line_search_lookup = NULL,  # Lookup for searching (includes both CellLineName and StrippedCellLineName)
    full_metadata = NULL,  # Store full Model.csv data for Excel export
    gene_list = character(0),
    gene_stats = NULL,  # Gene statistics (gene, LFC, FDR)
    gene_stats_raw = NULL,  # Raw uploaded file for column mapping
    results = NULL,
    enrichment_results = NULL,
    current_plot_data = NULL,  # Store plot data with cell line names
    current_scatter_genes = NULL,  # Store gene names for download
    current_cor_value = NULL,  # Store correlation value
    scatter_xlim = NULL,  # Store x-axis limits
    scatter_ylim = NULL,  # Store y-axis limits
    scatter_width = 8,  # Default width
    scatter_height = 8,  # Default height (square)
    scatter_cell_search = "",  # Cell line search terms
    hotspot_file_path = NULL,
    hotspot_genes = NULL,
    current_hotspot_gene = NULL,
    current_hotspot_mode = "none",
    scatter_plot_trigger = 0,  # Trigger for plot updates
    clicked_cells = character(0),  # Cell lines identified by clicking
    mutation_compare_data = NULL,  # Data for mutation comparison table download
    ortholog_lookup = NULL,  # Pre-computed ortholog lookup table (mouse -> human)
    ortholog_source_map = NULL  # Track which genes came from ortholog mapping
  )
  
  # -------------------------------------------------------------------------
  # Reference data loading
  # -------------------------------------------------------------------------
  
  # Helper function to load and process reference data
  load_reference_data <- function(filepath) {
    withProgress(message = "Loading reference data...", {
      tryCatch({
        ref_data <- fread(filepath, showProgress = FALSE)

        # Store cell line IDs from first column
        if (is.character(ref_data[[1]])) {
          rv$cell_line_ids <- ref_data[[1]]  # Store ACH IDs
          
          # Use metadata lookup if available, otherwise use ACH IDs
          if (!is.null(rv$cell_line_name_lookup)) {
            rv$cell_line_names <- sapply(rv$cell_line_ids, function(id) {
              result <- rv$cell_line_name_lookup[id]
              if (is.na(result)) id else result
            }, USE.NAMES = FALSE)
          } else {
            rv$cell_line_names <- rv$cell_line_ids  # Just use ACH IDs
          }
          
          ref_data <- ref_data[, -1, with = FALSE]
        } else {
          rv$cell_line_ids <- NULL
          rv$cell_line_names <- NULL
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

  # Cache last exported network figure (PNG data URL) so Download all can include it
  observeEvent(input$last_network_png, {
    req(input$last_network_png)
    url <- input$last_network_png
    if (!is.character(url) || !nzchar(url)) return()
    if (!grepl("^data:image/png;base64,", url)) return()

    b64 <- sub("^data:image/png;base64,", "", url)
    raw <- tryCatch(base64enc::base64decode(b64), error = function(e) NULL)
    if (is.null(raw)) return()

    fn <- file.path(tempdir(), paste0("network_export_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"))
    writeBin(raw, fn)
    rv$last_network_png_file <- fn
  }, ignoreInit = TRUE)



  # Load metadata file for cell line names (Model.csv from DepMap)
  load_metadata <- function(filepath) {
    withProgress(message = "Loading cell line metadata...", {
      tryCatch({
        metadata <- fread(filepath, showProgress = FALSE)
        
        # Store full metadata for Excel export
        # Column names from DepMap Model.csv
        if (ncol(metadata) >= 4) {
          rv$full_metadata <- as.data.frame(metadata)
          
          id_col <- 1  # Column A = ModelID (ACH-...)
          name_col <- 4  # Column D = StrippedCellLineName (HELA, A549, etc.)
          cellline_col <- 3  # Column C = CellLineName (HE-LA, A549, etc.)
          
          # Create main lookup table (ACH -> StrippedCellLineName)
          lookup <- setNames(as.character(metadata[[name_col]]), as.character(metadata[[id_col]]))
          rv$cell_line_name_lookup <- lookup
          
          # Create search lookup that maps BOTH CellLineName and StrippedCellLineName to ACH ID
          # This allows users to search by either name format
          search_lookup <- list()
          for (i in seq_len(nrow(metadata))) {
            ach_id <- as.character(metadata[[id_col]][i])
            stripped_name <- toupper(as.character(metadata[[name_col]][i]))
            cell_name <- toupper(as.character(metadata[[cellline_col]][i]))
            
            # Map both names to ACH ID (for reverse lookup during search)
            if (!is.na(stripped_name) && stripped_name != "") {
              search_lookup[[stripped_name]] <- ach_id
            }
            if (!is.na(cell_name) && cell_name != "" && cell_name != stripped_name) {
              search_lookup[[cell_name]] <- ach_id
            }
          }
          rv$cell_line_search_lookup <- search_lookup
          
          # Update cell line names if we have reference data loaded
          # Use single brackets [] which return NA for missing keys (not error)
          if (!is.null(rv$cell_line_ids)) {
            rv$cell_line_names <- sapply(rv$cell_line_ids, function(id) {
              result <- lookup[id]
              if (is.na(result)) id else result
            }, USE.NAMES = FALSE)
          }
          
          showNotification(
            paste("✓ Loaded metadata for", length(lookup), "cell lines"),
            type = "message",
            duration = 3
          )
        } else {
          showNotification("Model.csv format not recognized", type = "error")
        }
      }, error = function(e) {
        showNotification(paste("Error loading Model.csv:", e$message), type = "error")
      })
    })
  }
  
  observeEvent(input$metadata_file, {
    req(input$metadata_file)
    load_metadata(input$metadata_file$datapath)
  })

  # Hotspot mutation matrix upload (OmicsSomaticMutationsMatrixHotspot.csv)
  observeEvent(input$hotspot_file, {
    req(input$hotspot_file)
    rv$hotspot_file_path <- input$hotspot_file$datapath

    # Read header only to get available columns (keep memory low)
    hdr <- tryCatch({
      data.table::fread(rv$hotspot_file_path, nrows = 0, showProgress = FALSE)
    }, error = function(e) NULL)

    if (is.null(hdr)) {
      rv$hotspot_genes <- NULL
      rv$hotspot_genes_clean <- NULL
      rv$hotspot_gene_counts <- NULL
      rv$hotspot_id_col <- NULL
      showNotification("Failed to read hotspot file header. Please confirm it is a valid CSV.", type = "error")
      return()
    }

    cn <- names(hdr)

    # Detect an ID column robustly (prefer DepMap_ID, then any depmap*id, then ModelID, then any column containing ACH- IDs)
    id_col <- NULL

    # 1) direct matches
    if ("DepMap_ID" %in% cn) id_col <- "DepMap_ID"
    if (is.null(id_col)) {
      cand <- cn[grepl("depmap", cn, ignore.case = TRUE) & grepl("id", cn, ignore.case = TRUE)]
      if (length(cand) >= 1) id_col <- cand[1]
    }
    if (is.null(id_col) && "ModelID" %in% cn) id_col <- "ModelID"
    if (is.null(id_col) && "model_id" %in% cn) id_col <- "model_id"

    # 2) content-based detection: scan a small sample and pick the first column that looks like ACH- IDs
    samp <- tryCatch({
      data.table::fread(input$hotspot_file$datapath, nrows = 200, showProgress = FALSE)
    }, error = function(e) NULL)

    if (is.null(id_col) && !is.null(samp) && ncol(samp) >= 1) {
      for (cc in names(samp)) {
        v <- samp[[cc]]
        if (is.character(v) && any(grepl("^ACH-", v))) {
          id_col <- cc
          break
        }
      }
    }

    if (is.null(id_col) || !(id_col %in% cn)) {
      rv$hotspot_genes <- NULL
      rv$hotspot_genes_clean <- NULL
      rv$hotspot_gene_counts <- NULL
      rv$hotspot_id_col <- NULL
      showNotification("Hotspot file is missing a recognizable DepMap ID column (expected DepMap_ID, ModelID, or a column containing ACH- IDs). Please use the DepMap OmicsSomaticMutationsMatrixHotspot.csv file.", type = "error")
      return()
    }

    rv$hotspot_id_col <- id_col
    
    # Filter out non-gene columns (common metadata columns) - be very strict
    non_gene_cols <- c(id_col, "V1", "SequencingID", "ModelID", "model_id", "CellLineName", 
                       "StrippedCellLineName", "OncotreeLineage", "OncotreePrimaryDisease",
                       "OncotreeSubtype", "Age", "Sex", "PrimaryOrMetastasis",
                       "ModelConditionID", "IsDefaultEntryForModel", "IsDefaultEntryForMC",
                       "ProfileID", "PatientID", "CatalogNumber", "CCLEName", "COSMICID",
                       "PublicComments", "WTSIIssue", "RRID", "EngineeredModel")
    gene_cols <- setdiff(cn, non_gene_cols)
    
    # Further filter: only columns that look like gene names
    # Must start with a letter, and either:
    # - Be all caps with optional numbers (e.g., TP53, BRCA1, HLA-A)
    # - Have format "GENE (number)" like "TP53 (7157)"
    # Exclude columns that look like metadata (contain words like "ID", "Entry", "Model", "Default", "Condition")
    gene_cols <- gene_cols[grepl("^[A-Z][A-Z0-9-]*($|\\s*\\([0-9]+\\)$)", gene_cols, ignore.case = FALSE)]
    gene_cols <- gene_cols[!grepl("(ID|Entry|Model|Default|Condition|Profile|Patient|Catalog|Comment|Issue)", gene_cols, ignore.case = TRUE)]
    
    if (length(gene_cols) == 0) {
      rv$hotspot_genes <- NULL
      rv$hotspot_genes_clean <- NULL
      rv$hotspot_gene_counts <- NULL
      showNotification("No gene columns found in hotspot file.", type = "error")
      return()
    }
    
    # Read a larger sample to count mutations per gene (for filtering genes with ≥10 mutations)
    # Important: Only count mutations in cell lines that are present in the gene effect file
    withProgress(message = "Scanning hotspot mutations...", {
      full_data <- tryCatch({
        data.table::fread(rv$hotspot_file_path, select = c(id_col, gene_cols), showProgress = FALSE)
      }, error = function(e) NULL)
    })
    
    if (!is.null(full_data)) {
      # Filter to only cell lines present in gene effect file (if available)
      if (!is.null(rv$cell_line_ids) && length(rv$cell_line_ids) > 0) {
        full_data <- full_data[full_data[[id_col]] %in% rv$cell_line_ids, ]
      }
      
      # Count non-zero mutations per gene (only in overlapping cell lines)
      mutation_counts <- sapply(gene_cols, function(g) {
        vals <- full_data[[g]]
        sum(vals > 0, na.rm = TRUE)
      })
      
      # Filter genes with ≥10 mutations
      genes_with_mutations <- gene_cols[mutation_counts >= 10]
      counts_filtered <- mutation_counts[mutation_counts >= 10]
      
      # Remove HLA genes (these are HLA typing, not cancer driver mutations)
      clean_names <- gsub("\\s*\\([^)]+\\)$", "", genes_with_mutations)
      is_hla <- grepl("^HLA-", clean_names, ignore.case = TRUE)
      genes_with_mutations <- genes_with_mutations[!is_hla]
      counts_filtered <- counts_filtered[!is_hla]
      clean_names <- clean_names[!is_hla]
      
      # Sort by mutation count (descending)
      sort_order <- order(-counts_filtered)
      genes_with_mutations <- genes_with_mutations[sort_order]
      counts_filtered <- counts_filtered[sort_order]
      clean_names <- clean_names[sort_order]
      
      rv$hotspot_genes <- genes_with_mutations
      rv$hotspot_gene_counts <- counts_filtered
      rv$hotspot_genes_clean <- clean_names
      
      n_filtered <- length(gene_cols) - length(genes_with_mutations)
      showNotification(paste0("Hotspot file loaded: ", length(genes_with_mutations), " genes shown (", n_filtered, " filtered). Search for others."), type = "message")
    } else {
      # Fallback: just use all gene columns without filtering
      rv$hotspot_genes <- gene_cols
      rv$hotspot_genes_clean <- gsub("\\s*\\([^)]+\\)$", "", gene_cols)
      rv$hotspot_gene_counts <- NULL
      showNotification(paste0("Hotspot file loaded (header only): ", length(gene_cols), " gene columns available."), type = "message")
    }
  })
  
  # Handle ortholog lookup file upload
  observeEvent(input$ortholog_file, {
    req(input$ortholog_file)
    
    tryCatch({
      # Read the ortholog lookup file
      ortholog_data <- data.table::fread(input$ortholog_file$datapath, showProgress = FALSE)
      
      # Check for required columns
      required_cols <- c("mouse_symbol", "human_symbol")
      if (!all(required_cols %in% names(ortholog_data))) {
        showNotification("Ortholog file must have 'mouse_symbol' and 'human_symbol' columns", type = "error")
        return()
      }
      
      # Store the lookup table (uppercase for matching)
      ortholog_data[, mouse_symbol := toupper(mouse_symbol)]
      ortholog_data[, human_symbol := toupper(human_symbol)]
      
      # If there's a confidence column, sort by it
      if ("confidence" %in% names(ortholog_data)) {
        ortholog_data <- ortholog_data[order(-confidence)]
      }
      
      rv$ortholog_lookup <- ortholog_data
      
      n_mappings <- nrow(ortholog_data)
      n_unique_mouse <- length(unique(ortholog_data$mouse_symbol))
      
      showNotification(
        paste0("Ortholog lookup loaded: ", n_unique_mouse, " mouse genes mapped to human orthologs"),
        type = "message",
        duration = 5
      )
      
    }, error = function(e) {
      showNotification(paste("Error loading ortholog file:", e$message), type = "error")
    })
  })
  
  # Display reference status
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
  
  # Output flag for conditional panel
  output$metadata_loaded <- reactive({
    !is.null(rv$full_metadata)
  })
  outputOptions(output, "metadata_loaded", suspendWhenHidden = FALSE)
  # Ensure results tables update even when their tab is not active
  # Update Oncotree lineage filter choices when metadata is loaded
  observe({
    req(rv$full_metadata)
    
    # Get unique lineages from OncotreeLineage column (column 6)
    if ("OncotreeLineage" %in% colnames(rv$full_metadata)) {
      lineages <- sort(unique(na.omit(rv$full_metadata$OncotreeLineage)))
      lineages <- lineages[lineages != ""]
      
      updateSelectInput(session, "oncotree_lineage_filter",
                       choices = c("All lineages" = "", lineages),
                       selected = "")
    }
  })
  
  # Update disease filter based on selected lineage
  observeEvent(input$oncotree_lineage_filter, {
    req(rv$full_metadata)
    
    if (length(input$oncotree_lineage_filter) > 0 && any(input$oncotree_lineage_filter != "")) {
      # Filter to selected lineages
      selected_lineages <- input$oncotree_lineage_filter[input$oncotree_lineage_filter != ""]
      filtered_meta <- rv$full_metadata[rv$full_metadata$OncotreeLineage %in% selected_lineages, ]
      
      if ("OncotreePrimaryDisease" %in% colnames(filtered_meta)) {
        diseases <- sort(unique(na.omit(filtered_meta$OncotreePrimaryDisease)))
        diseases <- diseases[diseases != ""]
        
        updateSelectInput(session, "oncotree_disease_filter",
                         choices = c("All diseases" = "", diseases),
                         selected = "")
      }
    }
  }, ignoreNULL = FALSE)
  
  # Show filter status
  output$filter_status <- renderUI({
    if (is.null(rv$full_metadata) || is.null(rv$cell_line_ids)) return(NULL)
    
    # Calculate how many cell lines match the filter
    n_total <- length(rv$cell_line_ids)
    
    if (length(input$oncotree_lineage_filter) == 0 || all(input$oncotree_lineage_filter == "")) {
      return(div(style = "font-size: 0.85em; color: #666;",
                 paste("Using all", n_total, "cell lines")))
    }
    
    # Get ACH IDs that match the filter
    selected_lineages <- input$oncotree_lineage_filter[input$oncotree_lineage_filter != ""]
    filtered_meta <- rv$full_metadata[rv$full_metadata$OncotreeLineage %in% selected_lineages, ]
    
    if (length(input$oncotree_disease_filter) > 0 && any(input$oncotree_disease_filter != "")) {
      selected_diseases <- input$oncotree_disease_filter[input$oncotree_disease_filter != ""]
      filtered_meta <- filtered_meta[filtered_meta$OncotreePrimaryDisease %in% selected_diseases, ]
    }
    
    matching_ids <- intersect(rv$cell_line_ids, filtered_meta$ModelID)
    n_filtered <- length(matching_ids)
    
    div(class = "status-box status-success", style = "font-size: 0.85em; padding: 8px;",
        paste("✓ Filter active:", n_filtered, "of", n_total, "cell lines selected"))
  })
  
  # Helper function to get filtered cell line indices
  get_filtered_indices <- reactive({
    if (is.null(rv$cell_line_ids)) return(NULL)
    
    # If no metadata or no filter selected, return all indices
    if (is.null(rv$full_metadata) || 
        length(input$oncotree_lineage_filter) == 0 || 
        all(input$oncotree_lineage_filter == "")) {
      return(seq_along(rv$cell_line_ids))
    }
    
    # Filter by lineage
    selected_lineages <- input$oncotree_lineage_filter[input$oncotree_lineage_filter != ""]
    filtered_meta <- rv$full_metadata[rv$full_metadata$OncotreeLineage %in% selected_lineages, ]
    
    # Filter by disease if selected
    if (length(input$oncotree_disease_filter) > 0 && any(input$oncotree_disease_filter != "")) {
      selected_diseases <- input$oncotree_disease_filter[input$oncotree_disease_filter != ""]
      filtered_meta <- filtered_meta[filtered_meta$OncotreePrimaryDisease %in% selected_diseases, ]
    }
    
    # Return indices of matching cell lines
    matching_ids <- filtered_meta$ModelID
    which(rv$cell_line_ids %in% matching_ids)
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
  # Download example gene stats file
  output$download_example_stats <- downloadHandler(
    filename = function() { "example_gene_stats.csv" },
    content = function(file) {
      example_data <- data.frame(
        gene = c("TP53", "MDM2", "BRCA1", "MYC", "KRAS", "EGFR", "PTEN", "AKT1", "PIK3CA", "BRAF",
                "ATM", "CHEK2", "RB1", "CDKN2A", "SMAD4", "ARID1A", "NFE2L2", "KEAP1", "STK11", "CREBBP"),
        LFC = c(1.8, -2.3, 0.3, -1.5, 2.1, -0.8, 1.2, -0.4, 0.9, -1.9,
               0.6, -1.1, 1.4, -0.5, 0.7, -1.7, 1.0, -0.9, 1.6, -1.3),
        FDR = c(0.00001, 0.000001, 0.15, 0.0001, 0.000005, 0.02, 0.001, 0.08, 0.005, 0.00002,
               0.04, 0.003, 0.0002, 0.06, 0.01, 0.0003, 0.002, 0.007, 0.0001, 0.0004)
      )
      write.csv(example_data, file, row.names = FALSE)
    }
  )
  
  # Handle gene stats file upload
  observeEvent(input$gene_stats_file, {
    req(input$gene_stats_file)
    
    tryCatch({
      # Auto-detect separator (fread handles this well)
      stats_data <- fread(input$gene_stats_file$datapath, sep = "auto")
      rv$gene_stats_raw <- stats_data
      
      col_choices <- colnames(stats_data)
      
      # Try to auto-select columns based on common names
      gene_col <- col_choices[1]  # Default to first column
      if (any(grepl("^gene$", col_choices, ignore.case = TRUE))) {
        gene_col <- col_choices[grepl("^gene$", col_choices, ignore.case = TRUE)][1]
      }
      
      lfc_col <- if(length(col_choices) >= 2) col_choices[2] else col_choices[1]
      if (any(grepl("lfc|log.*fold", col_choices, ignore.case = TRUE))) {
        lfc_col <- col_choices[grepl("lfc|log.*fold", col_choices, ignore.case = TRUE)][1]
      }
      
      fdr_col <- if(length(col_choices) >= 3) col_choices[3] else col_choices[1]
      if (any(grepl("fdr|statistic|pval|p.value", col_choices, ignore.case = TRUE))) {
        fdr_col <- col_choices[grepl("fdr|statistic|pval|p.value", col_choices, ignore.case = TRUE)][1]
      }
      
      updateSelectInput(session, "stats_gene_col", choices = col_choices, selected = gene_col)
      updateSelectInput(session, "stats_lfc_col", choices = col_choices, selected = lfc_col)
      updateSelectInput(session, "stats_stat_col", choices = col_choices, selected = fdr_col)
      
      showNotification("File uploaded. Map columns and click 'Load Gene Stats'", type = "message")
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Load gene stats after column mapping
  observeEvent(input$load_stats, {
    req(rv$gene_stats_raw, input$stats_gene_col)
    
    stats_data <- rv$gene_stats_raw
    
    # Check if LFC and FDR columns are selected (not empty/NA)
    # Check if LFC and/or statistics columns are selected (not empty/NA)
    has_lfc <- !is.null(input$stats_lfc_col) && input$stats_lfc_col != "" && input$stats_lfc_col %in% colnames(stats_data)
    has_fdr <- !is.null(input$stats_stat_col) && input$stats_stat_col != "" && input$stats_stat_col %in% colnames(stats_data)

    if (has_lfc || has_fdr) {
      # Stats mode (LFC and/or FDR)
      gene_stats <- data.table(
        gene = toupper(trimws(as.character(stats_data[[input$stats_gene_col]]))),
        LFC = if (has_lfc) as.numeric(stats_data[[input$stats_lfc_col]]) else NA_real_,
        FDR = if (has_fdr) as.numeric(stats_data[[input$stats_stat_col]]) else NA_real_
      )
      rv$gene_stats <- gene_stats[!is.na(gene) & gene != ""]
    } else {
      # Genes only mode - no stats
      # Genes only mode - no stats
      genes <- toupper(trimws(as.character(stats_data[[input$stats_gene_col]])))
      genes <- genes[!is.na(genes) & genes != ""]
      rv$gene_stats <- NULL  # No stats
      updateTextAreaInput(session, "gene_textarea", value = paste(genes, collapse = "\n"))
      showNotification(paste("Loaded", length(genes), "genes (no statistics)"), type = "message")
      return()
    }
    
    updateTextAreaInput(session, "gene_textarea", value = paste(rv$gene_stats$gene, collapse = "\n"))
    
    output$stats_status <- renderUI({
      div(class = "status-box status-success",
          paste("✓ Loaded stats for", nrow(rv$gene_stats), "genes"))
    })
    
    showNotification(paste("Loaded stats for", nrow(rv$gene_stats), "genes"), type = "message")
  })
  
  # Output flags
  output$stats_file_uploaded <- reactive({ !is.null(rv$gene_stats_raw) })
  outputOptions(output, "stats_file_uploaded", suspendWhenHidden = FALSE)
  
  output$gene_stats_loaded <- reactive({ !is.null(rv$gene_stats) })
  outputOptions(output, "gene_stats_loaded", suspendWhenHidden = FALSE)
  
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
        actionButton("find_synonyms", "Find Synonyms/Orthologs", class = "btn-sm btn-warning",
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

    # Check if ortholog lookup is available
    has_ortholog_lookup <- !is.null(rv$ortholog_lookup) && nrow(rv$ortholog_lookup) > 0
    if (has_ortholog_lookup) {
      showNotification(paste0("Looking up ", length(not_found), " genes using ortholog table + MyGene.info..."), 
                       type = "message", duration = 3)
    }

    withProgress(message = "Looking up gene synonyms/orthologs...", value = 0, {
      replacements <- list()
      ortholog_sources <- list()  # Track which source found each gene
      n_total <- length(not_found)

      for (i in seq_along(not_found)) {
        gene <- not_found[i]
        incProgress(1/n_total, detail = gene)

        # Pass ortholog lookup to the function
        result <- find_gene_synonyms(gene, available_genes, rv$ortholog_lookup)

        if (result$found && length(result$matches) > 0) {
          # Use the first match as replacement
          replacements[[gene]] <- result$matches[1]
          # Track the source
          if (!is.null(result$source)) {
            ortholog_sources[[gene]] <- result$source
          }
        }

        # Only rate limit for API calls (when no ortholog lookup or ortholog lookup failed)
        if (is.null(rv$ortholog_lookup) || !result$found || result$source == "mygene") {
          Sys.sleep(0.1)  # Rate limiting for MyGene API
        }
      }
    })

    if (length(replacements) > 0) {

      # Store synonym mappings for display (replacement -> original) and reporting
      rv$synonym_map <- replacements
      rv$reverse_syn_map <- setNames(names(replacements), unname(replacements))
      rv$ortholog_source_map <- ortholog_sources
      
      # Update gene list with replacements
      updated_genes <- genes
      for (original in names(replacements)) {
        updated_genes[updated_genes == original] <- replacements[[original]]
      }

      updateTextAreaInput(session, "gene_textarea", value = paste(updated_genes, collapse = "\n"))

      # Create detailed message showing source of each mapping
      msg_parts <- sapply(names(replacements), function(x) {
        source_label <- if (!is.null(ortholog_sources[[x]])) paste0(" [", ortholog_sources[[x]], "]") else ""
        paste0(x, " → ", replacements[[x]], source_label)
      })
      
      # Count by source
      n_ortholog <- sum(sapply(ortholog_sources, function(s) !is.null(s) && s == "ortholog"))
      n_mygene <- sum(sapply(ortholog_sources, function(s) !is.null(s) && s == "mygene"))
      
      source_summary <- paste0(
        "Replaced ", length(replacements), " gene(s): ",
        if (n_ortholog > 0) paste0(n_ortholog, " via ortholog table") else "",
        if (n_ortholog > 0 && n_mygene > 0) ", " else "",
        if (n_mygene > 0) paste0(n_mygene, " via MyGene.info") else ""
      )
      
      msg <- paste0(source_summary, "\n", paste(msg_parts, collapse = ", "))
      showNotification(msg, type = "message", duration = 10)
    } else {
      showNotification("No synonyms/orthologs found for the missing genes", type = "warning")
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
    
    if (length(rv$gene_list) < 1) {
      showNotification("Please enter at least 1 gene", type = "error")
      return()
    }
    
    # In analysis mode, need at least 2 genes
    if (input$analysis_mode == "analysis" && length(rv$gene_list) < 2) {
      showNotification("Analysis mode requires at least 2 genes (or use Design mode for single gene)", type = "error")
      return()
    }
    
    # Get filtered indices
    filtered_idx <- get_filtered_indices()
    
    # Check if we have enough cell lines after filtering
    if (length(filtered_idx) < 10) {
      showNotification("Too few cell lines match the filter (need at least 10)", type = "error")
      return()
    }
    
    # Filter reference data if needed
    if (length(filtered_idx) < nrow(rv$reference_data)) {
      filtered_data <- rv$reference_data[filtered_idx, ]
      filter_msg <- paste0(" (filtered to ", length(filtered_idx), " cell lines)")
    } else {
      filtered_data <- rv$reference_data
      filter_msg <- ""
    }
    
    withProgress(message = "Running correlation analysis...", {
      rv$results <- run_correlation_analysis(
        reference_data = filtered_data,
        gene_list = rv$gene_list,
        analysis_mode = input$analysis_mode,
        correlation_cutoff = input$correlation_cutoff
      )
      # Store filter info for later use
      rv$results$filtered_indices <- filtered_idx
      rv$results$filter_msg <- filter_msg
    })
    
    if (rv$results$success) {
      corr_count <- nrow(rv$results$correlations)
      showNotification(
        paste("Found", corr_count, "correlations in", nrow(rv$results$clusters), "genes", filter_msg),
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
  
  # Render static legend below network - updates with edge width
  output$network_legend_static <- renderUI({
  req(rv$results$success)

  cutoff <- rv$results$cutoff
  mode <- rv$results$mode
  edge_width <- if (!is.null(input$net_edge_width)) input$net_edge_width else 3

  # Optional node color legend (included in exported PNG)
  node_color_block <- NULL
  if (isTRUE(input$net_color_by_stats) && !is.null(rv$gene_stats)) {
    if (input$net_color_stat == "signed_lfc") {
      max_abs <- suppressWarnings(max(abs(scale_stats_dt()$LFC), na.rm = TRUE))
      if (!is.finite(max_abs) || is.na(max_abs) || max_abs <= 0) max_abs <- 1
      node_color_block <- legend_bar_html(
        title = "Node color: Signed LFC",
        left_lab = paste0("-", round(max_abs, 2)),
        right_lab = paste0("+", round(max_abs, 2)),
        colors = c("#2166ac", "#f7f7f7", "#b2182b")
      )
} else if (input$net_color_stat == "abs_lfc") {
      max_abs <- suppressWarnings(max(abs(scale_stats_dt()$LFC), na.rm = TRUE))
      if (!is.finite(max_abs) || max_abs <= 0) max_abs <- 1
      node_color_block <- legend_bar_html(
        title = "Node color: Absolute LFC",
        left_lab = "0",
        right_lab = round(max_abs, 2),
        colors = c("#f5f5f5", "#fdae61", "#d7191c")
      )
    } else if (input$net_color_stat == "fdr") {
      min_fdr <- suppressWarnings(min(scale_stats_dt()$FDR[scale_stats_dt()$FDR > 0], na.rm = TRUE))
      if (!is.finite(min_fdr) || min_fdr <= 0) min_fdr <- .Machine$double.xmin
      node_color_block <- legend_bar_html(
        title = "Node color: Smallest FDR",
        left_lab = format(min_fdr, scientific = TRUE, digits = 2),
        right_lab = "1",
        colors = c("#1d4ed8", "#93c5fd", "#f5f5f5")
      )
    }
  }

  div(
    style = "background: #f9fafb; border: 1px solid #e5e7eb; border-radius: 8px; padding: 12px 15px;",
    fluidRow(
      # Column 1: Correlation type
      column(
        2,
        div(
          style = "font-size: 11px;",
          strong("Correlation:"),
          div(
            style = "margin-top: 6px;",
            div(
              style = "margin-bottom: 4px; display: flex; align-items: center;",
              span(style = "display: inline-block; width: 25px; height: 3px; background: #3182ce; margin-right: 6px;"),
              span("Positive")
            ),
            div(
              style = "display: flex; align-items: center;",
              span(style = "display: inline-block; width: 25px; height: 3px; background: #e53e3e; margin-right: 6px;"),
              span("Negative")
            )
          )
        )
      ),
      # Column 2: Edge Thickness
      column(
        3,
        div(
          style = "font-size: 11px;",
          strong(paste0("Edge Thickness (cutoff: ", round(cutoff, 2), "):")),
          div(
            style = "margin-top: 6px; display: flex; align-items: center; gap: 10px;",
            div(
              style = "display: flex; align-items: center;",
              span(style = "display: inline-block; width: 20px; height: 1px; background: #666; margin-right: 4px;"),
              span(paste0("r=", round(cutoff, 2)))
            ),
            div(
              style = "display: flex; align-items: center;",
              span(style = paste0("display: inline-block; width: 20px; height: ", round(edge_width * 1.3), "px; background: #666; margin-right: 4px;")),
              span("r=0.7")
            ),
            div(
              style = "display: flex; align-items: center;",
              span(style = paste0("display: inline-block; width: 20px; height: ", round(edge_width * 2.7), "px; background: #666; margin-right: 4px;")),
              span("r=1.0")
            )
          )
        )
      ),
      # Column 3: Node Type (only for design mode)
      if (mode == "design") {
        column(
          2,
          div(
            style = "font-size: 11px;",
            strong("Node Type:"),
            div(
              style = "margin-top: 6px;",
              div(
                style = "margin-bottom: 4px; display: flex; align-items: center;",
                span(style = "display: inline-block; width: 10px; height: 10px; background: #16a34a; border-radius: 50%; margin-right: 6px;"),
                span("Input")
              ),
              div(
                style = "display: flex; align-items: center;",
                span(style = "display: inline-block; width: 10px; height: 10px; background: #86efac; border-radius: 50%; margin-right: 6px;"),
                span("Correlated")
              )
            )
          )
        )
      } else {
        NULL
      },
      # Column 4: Node color legend (if stats coloring is enabled)
      if (!is.null(node_color_block)) {
        column(
          if (mode == "design") 3 else 5,
          node_color_block
        )
      } else {
        NULL
      },
      # Column 5: Synonym indicator (if synonyms were used)
      if (!is.null(rv$synonym_map) && length(rv$synonym_map) > 0) {
        column(
          2,
          div(
            style = "font-size: 11px;",
            strong("Note:"),
            div(
              style = "margin-top: 6px; display: flex; align-items: center;",
              span("*", style = "font-weight: bold; margin-right: 4px;"),
              span("synonym/ortholog")
            )
          )
        )
      } else {
        NULL
      }
    )
  )
})
  
  # Store network data for updates
  network_data <- reactiveValues(nodes = NULL, edges = NULL, hidden_nodes = character(0),
                                  all_nodes = NULL, all_edges = NULL)

  # Network visualization - only re-render when results change
  output$network_plot <- renderVisNetwork({
    req(rv$results$success)

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
    
    # First, create base labels with asterisks for synonym-replaced genes
    base_labels <- sapply(node_names, function(gene) {
      if (!is.null(rv$reverse_syn_map) && gene %in% names(rv$reverse_syn_map)) {
        paste0(gene, "*")  # Show TP53* for genes that were mapped from synonyms/orthologs
      } else {
        gene
      }
    })
    
    if (show_effect) {
      labels <- sapply(seq_along(node_names), function(i) {
        gene <- node_names[i]
        display_gene <- base_labels[i]
        effect_row <- clusters_dt[Gene == gene]
        if (nrow(effect_row) > 0 && !is.na(effect_row$Mean_Effect[1])) {
          if (show_sd && !is.na(effect_row$SD_Effect[1])) {
            paste0(display_gene, "\n(GE:", effect_row$Mean_Effect[1], " ± ", effect_row$SD_Effect[1], ")")
          } else {
            paste0(display_gene, "\n(GE:", effect_row$Mean_Effect[1], ")")
          }
        } else {
          display_gene
        }
      })
    } else {
      labels <- base_labels
    }
    
    # Add LFC/FDR to labels if gene stats are loaded
    stats_display <- isolate(input$net_stats_display)
    
    if (!is.null(rv$gene_stats) && !is.null(stats_display) && stats_display != "none") {
      new_labels <- character(length(node_names))
      for (i in seq_along(node_names)) {
        gene <- node_names[i]
        base_label <- labels[i]
        gene_up <- toupper(trimws(gene)); stats_row <- rv$gene_stats[rv$gene_stats$gene == gene_up, ]
        
        if (nrow(stats_row) > 0) {
          if (stats_display == "lfc" && !is.na(stats_row$LFC[1])) {
            new_labels[i] <- paste0(base_label, "\n[LFC:", round(stats_row$LFC[1], 1), "]")
          } else if (stats_display == "fdr" && !is.na(stats_row$FDR[1])) {
            new_labels[i] <- paste0(base_label, "\n[FDR:", formatC(stats_row$FDR[1], format = "e", digits = 1), "]")
          } else {
            new_labels[i] <- base_label
          }
        } else {
          new_labels[i] <- base_label
        }
      }
      labels <- new_labels
    }

    # Use igraph layout for initial positioning - spread nodes nicely
    set.seed(42)  # For reproducible layout
    n_nodes <- length(node_names)
    
    # Choose layout based on network size
    if (n_nodes <= 10) {
      layout_coords <- layout_with_fr(graph, niter = 500)
    } else if (n_nodes <= 30) {
      layout_coords <- layout_with_fr(graph, niter = 300)
    } else {
      layout_coords <- layout_with_kk(graph)
    }
    
    # Scale layout to reasonable coordinates
    layout_coords <- scale(layout_coords) * 200

    # Create nodes with pre-computed positions
    nodes <- data.frame(
      id = node_names,
      label = labels,
      size = node_size,
      font.size = font_size,
      font.multi = "html",
      group = "gene",
      shape = "dot",
      x = layout_coords[, 1],
      y = layout_coords[, 2],
      physics = FALSE,  # Disable physics since we have good initial layout
      stringsAsFactors = FALSE
    )

    # Add tooltips with gene effect
    nodes$title <- sapply(node_names, function(gene) {
      effect_row <- clusters_dt[Gene == gene]
      if (nrow(effect_row) > 0 && !is.na(effect_row$Mean_Effect[1])) {
        paste0("<b>", gene, "</b><br>Mean Effect: ", effect_row$Mean_Effect[1],
               "<br>SD: ", effect_row$SD_Effect[1], "<br><i>Double-click to hide</i>")
      } else {
        paste0("<b>", gene, "</b><br><i>Double-click to hide</i>")
      }
    })

    # Color nodes based on gene stats or original gene list
    color_by_stats <- isolate(input$net_color_by_stats)
    color_stat_type <- isolate(input$net_color_stat)
    
    if (!is.null(rv$gene_stats) && !is.null(color_by_stats) && color_by_stats) {
      # Color by statistics
      color_vec <- sapply(node_names, function(gene) {
        gene_up <- toupper(trimws(gene))
        stats_row <- rv$gene_stats[rv$gene_stats$gene == gene_up, ]

        need_col <- if (!is.null(color_stat_type) && color_stat_type %in% c("abs_lfc", "signed_lfc")) "LFC" else "FDR"
        if (nrow(stats_row) == 0 || is.na(stats_row[[need_col]][1])) {
          return("#cccccc")  # Gray for missing data
        }

        if (!is.null(color_stat_type) && color_stat_type == "signed_lfc") {
          lfc_val <- stats_row$LFC[1]
          min_lfc <- suppressWarnings(min(rv$gene_stats$LFC, na.rm = TRUE))
          max_lfc <- suppressWarnings(max(rv$gene_stats$LFC, na.rm = TRUE))
          map_signed_lfc_color(lfc_val, min_lfc, max_lfc)
        } else if (!is.null(color_stat_type) && color_stat_type == "abs_lfc") {
          lfc_val <- stats_row$LFC[1]
          max_abs <- suppressWarnings(max(abs(rv$gene_stats$LFC), na.rm = TRUE))
          map_abs_lfc_color(lfc_val, max_abs)
        } else {
          fdr_val <- stats_row$FDR[1]
          min_fdr <- suppressWarnings(min(rv$gene_stats$FDR[rv$gene_stats$FDR > 0], na.rm = TRUE))
          map_fdr_color(fdr_val, min_fdr)
        }
      })
      nodes$color <- lapply(unname(color_vec), function(clr) list(background = clr, border = "#ffffff"))
    } else if (mode == "design") {
      # Original coloring by input gene
      nodes$color <- lapply(ifelse(nodes$id %in% matched, "#16a34a", "#86efac"), function(clr) list(background = clr, border = "#ffffff"))
      nodes$title <- paste0(nodes$title, "<br><i>",
                            ifelse(nodes$id %in% matched, "Input gene", "Correlated gene"), "</i>")
    } else {
      nodes$color <- lapply(rep("#16a34a", nrow(nodes)), function(clr) list(background = clr, border = "#ffffff"))
    }

    # Create edges with width proportional to correlation
    # Scale edge width more dramatically for visibility
    min_width <- 1
    max_width <- edge_width * 4
    
    edges <- data.frame(
      id = paste0(correlations$Gene1, "-", correlations$Gene2),
      from = correlations$Gene1,
      to = correlations$Gene2,
      correlation = correlations$Correlation,
      width = min_width + (abs(correlations$Correlation) - cutoff) / (1 - cutoff) * (max_width - min_width),
      color = ifelse(correlations$Correlation > 0, "#3182ce", "#e53e3e"),
      title = paste("r =", round(correlations$Correlation, 3)),
      smooth = FALSE,
      stringsAsFactors = FALSE
    )

    # Remove duplicate edges (since graph is undirected)
    edges <- edges[!duplicated(t(apply(edges[,c("from", "to")], 1, sort))),]

    # Store all nodes/edges for restoration
    network_data$all_nodes <- nodes
    network_data$all_edges <- edges
    network_data$nodes <- nodes
    network_data$edges <- edges
    network_data$clusters_dt <- clusters_dt
    network_data$hidden_nodes <- character(0)
    network_data$cutoff <- cutoff
    network_data$mode <- mode

    # Build the network
    vis <- visNetwork(nodes, edges) %>%
      visNodes(borderWidth = 2,
               borderWidthSelected = 4,
               font = list(face = "arial", color = "#333333", multi = "html")) %>%
      visEdges(smooth = FALSE) %>%
      visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
                 nodesIdSelection = FALSE) %>%
      visPhysics(enabled = FALSE) %>%
      visInteraction(navigationButtons = TRUE,
                     dragNodes = TRUE,
                     dragView = TRUE,
                     zoomView = TRUE) %>%
      visEvents(doubleClick = "function(params) {
        if (params.nodes.length > 0) {
          var nodeId = params.nodes[0];
          Shiny.setInputValue('clicked_node', nodeId, {priority: 'event'});
        }
      }")
    
    vis
  })

  # Update font size without re-rendering
  observeEvent(input$net_font_size, {
    req(network_data$nodes)
    nodes <- data.frame(network_data$nodes)
    nodes$font.size <- input$net_font_size
    # Preserve current node size
    if (!is.null(input$net_node_size)) {
      nodes$size <- input$net_node_size
    }
    visNetworkProxy("network_plot") %>%
      visUpdateNodes(nodes)
  }, ignoreInit = TRUE)

  # Update node size without re-rendering
  observeEvent(input$net_node_size, {
    req(network_data$nodes)
    nodes <- data.frame(network_data$nodes)
    nodes$size <- input$net_node_size
    # Preserve current font size
    if (!is.null(input$net_font_size)) {
      nodes$font.size <- input$net_font_size
    }
    visNetworkProxy("network_plot") %>%
      visUpdateNodes(nodes)
  }, ignoreInit = TRUE)

  # Update edge width without re-rendering
  observeEvent(input$net_edge_width, {
    req(network_data$edges)
    req(network_data$cutoff)
    edges <- data.frame(network_data$edges)
    cutoff <- network_data$cutoff
    min_width <- 1
    max_width <- input$net_edge_width * 4
    edges$width <- min_width + (abs(edges$correlation) - cutoff) / (1 - cutoff) * (max_width - min_width)
    visNetworkProxy("network_plot") %>%
      visUpdateEdges(edges)
  }, ignoreInit = TRUE)

  # ---- Helper: build node labels (gene + optional gene effect + optional LFC/FDR) ----
build_node_label <- function(gene, clusters_dt, gene_stats_dt) {
  base_label <- gene

  # If this gene was a synonym replacement, display human gene with asterisk
  # gene is the human gene (e.g., TP53), reverse_syn_map maps human -> original input (e.g., Trp53)
  display_gene <- gene
  if (!is.null(rv$reverse_syn_map) && gene %in% names(rv$reverse_syn_map)) {
    display_gene <- paste0(gene, "*")  # Show TP53* not Trp53*
  }
  base_label <- display_gene


  # Gene effect label
  if (isTRUE(input$net_show_gene_effect)) {
    effect_row <- clusters_dt[Gene == gene]
    if (nrow(effect_row) > 0 && !is.na(effect_row$Mean_Effect[1])) {
      if (isTRUE(input$net_show_sd) && !is.na(effect_row$SD_Effect[1])) {
        base_label <- paste0(display_gene, "\n(GE:", effect_row$Mean_Effect[1], " ± ", effect_row$SD_Effect[1], ")")
      } else {
        base_label <- paste0(display_gene, "\n(GE:", effect_row$Mean_Effect[1], ")")
      }
    }
  }

  # Gene statistics (uploaded)
  if (!is.null(gene_stats_dt) && !is.null(input$net_stats_display) && input$net_stats_display != "none") {
    gene_up <- toupper(trimws(gene))
    stats_row <- gene_stats_dt[gene == gene_up]  # keyed DT supports this; fallback works too
    if (nrow(stats_row) == 0) stats_row <- gene_stats_dt[gene_stats_dt$gene == gene_up, ]
    if (nrow(stats_row) > 0) {
      if (input$net_stats_display == "lfc" && !is.na(stats_row$LFC[1])) {
        base_label <- paste0(base_label, "\n[LFC:", round(stats_row$LFC[1], 2), "]")
      } else if (input$net_stats_display == "fdr" && !is.na(stats_row$FDR[1])) {
        base_label <- paste0(base_label, "\n[FDR:", formatC(stats_row$FDR[1], format = "e", digits = 1), "]")
      }
    }
  }

  base_label
}

# Update labels in one place so options act independently (no overwrite-by-click-order)
# Update labels in one place so options act independently (no overwrite-by-click-order)
  observeEvent(
    list(input$net_show_gene_effect, input$net_show_sd, input$net_stats_display),
    {
      req(network_data$nodes)
      nodes <- network_data$nodes
      clusters_dt <- network_data$clusters_dt
      gene_stats_dt <- rv$gene_stats

      # Rebuild labels
      nodes$label <- vapply(nodes$id, build_node_label, character(1),
                            clusters_dt = clusters_dt, gene_stats_dt = gene_stats_dt)

      # IMPORTANT: include existing colors so label updates do not reset node colors
      vis_nodes <- data.frame(
        id = nodes$id,
        label = nodes$label,
        color = I(nodes$color),
        stringsAsFactors = FALSE
      )

      visNetworkProxy("network_plot") %>% visUpdateNodes(vis_nodes)

      # Keep stored nodes in sync (so future updates have the correct color/label)
      network_data$nodes$label <- nodes$label
      if (!is.null(network_data$all_nodes)) {
        idx <- match(network_data$all_nodes$id, nodes$id)
        network_data$all_nodes$label <- nodes$label[idx]
      }
    },
    ignoreInit = TRUE
  )
  


  
  

  
  # Helper: determine which genes to use for color scale ranges
  scale_genes <- reactive({
    req(rv$gene_stats)
    scope <- input$net_color_scale_scope
    if (is.null(scope) || scope == "") scope <- "all"

    all_genes <- unique(toupper(trimws(rv$gene_stats$gene)))

    if (scope == "network") {
      if (!is.null(network_data$edges) && nrow(network_data$edges) > 0) {
        net_genes <- unique(toupper(trimws(c(network_data$edges$from, network_data$edges$to))))
        return(intersect(all_genes, net_genes))
      } else {
        return(character(0))
      }
    }

    all_genes
  })

  # Helper: fetch stats subset for scale computations
  scale_stats_dt <- reactive({
    req(rv$gene_stats)
    genes <- scale_genes()
    if (length(genes) == 0) return(rv$gene_stats[0, ])
    rv$gene_stats[rv$gene_stats$gene %in% genes, ]
  })

# Legend for node coloring by statistics
  output$net_color_legend <- renderUI({
    req(input$net_color_by_stats)
    if (!isTRUE(input$net_color_by_stats)) return(NULL)
    req(rv$gene_stats)

    if (input$net_color_stat == "signed_lfc") {
      min_lfc <- suppressWarnings(min(scale_stats_dt()$LFC, na.rm = TRUE))
max_lfc <- suppressWarnings(max(scale_stats_dt()$LFC, na.rm = TRUE))
if (!is.finite(min_lfc) || is.na(min_lfc)) min_lfc <- -1
if (!is.finite(max_lfc) || is.na(max_lfc)) max_lfc <- 1
legend_bar_html(
  title = "Node color: Signed LFC",
  left_lab = round(min_lfc, 2),
  right_lab = round(max_lfc, 2),
  colors = c("#2166ac", "#f7f7f7", "#b2182b")
)
    } else if (input$net_color_stat == "abs_lfc") {
      max_lfc <- suppressWarnings(max(abs(scale_stats_dt()$LFC), na.rm = TRUE))
      if (!is.finite(max_lfc) || max_lfc <= 0) max_lfc <- 1
      legend_bar_html(
        title = "Node color: Absolute LFC",
        left_lab = "0",
        right_lab = paste0(round(max_lfc, 2)),
        colors = c("#f5f5f5", "#fdae61", "#d7191c")
      )
    } else {
      min_fdr <- suppressWarnings(min(scale_stats_dt()$FDR[rv$gene_stats$FDR > 0], na.rm = TRUE))
      if (!is.finite(min_fdr) || is.na(min_fdr) || min_fdr <= 0) min_fdr <- .Machine$double.xmin
      legend_bar_html(
        title = "Node color: Smallest FDR",
        left_lab = format(min_fdr, scientific = TRUE, digits = 2),
        right_lab = "1",
        # Low FDR should be strongest (dark) on the left
        colors = c("#08519c", "#9ecae1", "#f5f5f5")
      )
    }
  })
# Update colors when color_by_stats checkbox changes
  observeEvent(list(input$net_color_by_stats, input$net_color_stat, input$net_color_scale_scope), {
    req(network_data$nodes)
    req(rv$results$success)
    req(rv$gene_stats)
    
    nodes <- data.frame(network_data$nodes)
    color_by_stats <- input$net_color_by_stats
    color_stat_type <- input$net_color_stat
    mode <- rv$results$mode
    matched <- rv$results$matched
    
    if (color_by_stats) {
      # Color by statistics
      color_vec <- sapply(nodes$id, function(gene) {
        gene_up <- toupper(trimws(gene)); stats_row <- rv$gene_stats[rv$gene_stats$gene == gene_up, ]
        
        if (nrow(stats_row) == 0 || is.na(stats_row[[if(color_stat_type %in% c("abs_lfc","signed_lfc")) "LFC" else "FDR"]][1])) {
          return("#cccccc")  # Gray for missing data
        }
        
        if (color_stat_type == "signed_lfc") {
        lfc_val <- stats_row$LFC[1]
        max_abs <- suppressWarnings(max(abs(scale_stats_dt()$LFC), na.rm = TRUE))
        if (!is.finite(max_abs) || is.na(max_abs) || max_abs <= 0) max_abs <- 1
        min_lfc <- -max_abs
        max_lfc <- max_abs
        
        map_signed_lfc_color(lfc_val, min_lfc, max_lfc)
      } else if (color_stat_type == "abs_lfc") {
        lfc_val <- stats_row$LFC[1]
        max_lfc <- suppressWarnings(max(abs(scale_stats_dt()$LFC), na.rm = TRUE))
        map_abs_lfc_color(lfc_val, max_lfc)
      } else {
        fdr_val <- stats_row$FDR[1]
        min_fdr <- suppressWarnings(min(scale_stats_dt()$FDR[scale_stats_dt()$FDR > 0], na.rm = TRUE))
        map_fdr_color(fdr_val, min_fdr)
      }
      })
      nodes$color <- lapply(unname(color_vec), function(clr) list(background = clr, border = "#ffffff"))
    } else {
      # Revert to original colors
      if (mode == "design") {
        nodes$color <- lapply(ifelse(nodes$id %in% matched, "#16a34a", "#86efac"), function(clr) list(background = clr, border = "#ffffff"))
      } else {
        nodes$color <- lapply(rep("#16a34a", nrow(nodes)), function(clr) list(background = clr, border = "#ffffff"))
      }
    }
    
    # Persist colors so other UI updates (e.g., label toggles) do not reset styling
    network_data$nodes$color <- nodes$color
    if (!is.null(network_data$all_nodes)) {
      idx <- match(network_data$all_nodes$id, nodes$id)
      network_data$all_nodes$color <- nodes$color[idx]
    }

    vis_nodes <- data.frame(id = nodes$id, color = I(nodes$color), stringsAsFactors = FALSE)
    visNetworkProxy("network_plot") %>%
      visUpdateNodes(vis_nodes)
  }, ignoreInit = TRUE)
  
  # Export handlers - network uses JavaScript to capture canvas
  observeEvent(input$export_network_png, {
    session$sendCustomMessage("exportNetworkPNG", list())
  })
  
  observeEvent(input$export_network_svg, {
    session$sendCustomMessage("exportNetworkSVG", list())
  })
  
  # Export legend as PNG (keep for backward compatibility)
    draw_combined_legend <- function(cutoff, edge_width) {
    par(mar = c(0.8, 0.8, 0.8, 0.8))
    plot.new()
    
    # Column anchors (spread to avoid overlap)
    x1 <- 0.06; x2 <- 0.36; x3 <- 0.72
    
    # --- Correlation sign (column 1) ---
    text(x1, 0.86, "Correlation", font = 2, cex = 1.1, adj = 0)
    
    segments(x1, 0.72, x1 + 0.08, 0.72, col = "#3182ce", lwd = 3)
    text(x1 + 0.10, 0.72, "Positive", adj = 0, cex = 1)
    
    segments(x1, 0.62, x1 + 0.08, 0.62, col = "#e53e3e", lwd = 3)
    text(x1 + 0.10, 0.62, "Negative", adj = 0, cex = 1)
    
    # --- Edge thickness (column 2) ---
    text(x2, 0.86, paste0("Edge thickness (cutoff r=", round(cutoff, 2), ")"), font = 2, cex = 1.1, adj = 0)
    
    segments(x2, 0.72, x2 + 0.08, 0.72, col = "gray30", lwd = edge_width * 0.8)
    text(x2 + 0.10, 0.72, paste0("r = ", round(cutoff, 2)), adj = 0, cex = 1)
    
    segments(x2, 0.62, x2 + 0.08, 0.62, col = "gray30", lwd = edge_width * 1.7)
    text(x2 + 0.10, 0.62, paste0("r = ", round((1 + cutoff) / 2, 2)), adj = 0, cex = 1)
    
    segments(x2, 0.52, x2 + 0.08, 0.52, col = "gray30", lwd = edge_width * 2.7)
    text(x2 + 0.10, 0.52, "r = 1.0", adj = 0, cex = 1)
    
    # --- Design mode indicators (under column 1) ---
    if (mode == "design") {
      text(x1, 0.42, "Node type", font = 2, cex = 1.2, adj = 0)
    
      points(x1 + 0.03, 0.30, pch = 19, col = "#16a34a", cex = 2)
      text(x1 + 0.07, 0.30, "Input gene", adj = 0, cex = 1)
    
      points(x1 + 0.03, 0.20, pch = 19, col = "#86efac", cex = 2)
      text(x1 + 0.07, 0.20, "Correlated gene", adj = 0, cex = 1)
    }
    
    # --- Optional node color bar (column 3) ---
    if (isTRUE(input$net_color_by_stats) && !is.null(rv$gene_stats)) {
      text(x3, 0.86, "Node color", font = 2, cex = 1.1, adj = 0)
    
      if (input$net_color_stat == "signed_lfc") {
      max_abs <- suppressWarnings(max(abs(scale_stats_dt()$LFC), na.rm = TRUE))
      if (!is.finite(max_abs) || is.na(max_abs) || max_abs <= 0) max_abs <- 1
      cols <- colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(100)
      left_lab <- paste0("-", round(max_abs, 2))
      right_lab <- paste0("+", round(max_abs, 2))
      subtitle <- "Signed LFC"
    } else if (input$net_color_stat == "abs_lfc") {
        max_lfc <- suppressWarnings(max(abs(scale_stats_dt()$LFC), na.rm = TRUE))
        if (!is.finite(max_lfc) || max_lfc <= 0) max_lfc <- 1
        cols <- colorRampPalette(c("#f5f5f5", "#fdae61", "#d7191c"))(100)
        left_lab <- "0"
        right_lab <- paste0(round(max_lfc, 2))
        subtitle <- "Absolute LFC"
      } else {
        min_fdr <- suppressWarnings(min(scale_stats_dt()$FDR[scale_stats_dt()$FDR > 0], na.rm = TRUE))
        if (!is.finite(min_fdr) || is.na(min_fdr) || min_fdr <= 0) min_fdr <- .Machine$double.xmin
        # Low FDR should be strongest (dark) on the left
        cols <- colorRampPalette(c("#1d4ed8", "#93c5fd", "#f5f5f5"))(100)
        left_lab <- format(min_fdr, scientific = TRUE, digits = 2)
        right_lab <- "1"
        subtitle <- "Smallest FDR"
      }
    
      # Subtitle
      text(x3, 0.74, subtitle, adj = 0, cex = 1)
    
      # Gradient bar
      x0 <- x3
      x1b <- x3 + 0.22
      y0 <- 0.61
      y1b <- 0.67
      xs <- seq(x0, x1b, length.out = 101)
      for (i in 1:100) {
        rect(xs[i], y0, xs[i + 1], y1b, col = cols[i], border = NA)
      }
      rect(x0, y0, x1b, y1b, border = "gray60")
    
      # Labels
      text(x0, 0.56, left_lab, adj = c(0, 0.5), cex = 0.95)
      text(x1b, 0.56, right_lab, adj = c(1, 0.5), cex = 0.95)
    }
    
    
  }

output$export_legend <- downloadHandler(
    filename = function() {
      paste0("correlation_legend_", format(Sys.Date(), "%Y%m%d"), ".png")
    },
    content = function(file) {req(rv$results$success)

cutoff <- rv$results$cutoff
edge_width <- if (!is.null(input$net_edge_width)) input$net_edge_width else 3

png(file, width = 1300, height = 420, res = 150, bg = "white")
par(mar = c(0.8, 0.8, 0.8, 0.8))
draw_combined_legend(cutoff = cutoff, edge_width = edge_width)
dev.off()
    }
  )


output$export_legend_svg <- downloadHandler(
    filename = function() {
      paste0("correlation_legend_", format(Sys.Date(), "%Y%m%d"), ".svg")
    },
    content = function(file) {req(rv$results$success)

cutoff <- rv$results$cutoff
edge_width <- if (!is.null(input$net_edge_width)) input$net_edge_width else 3

svg(file, width = 1300/150, height = 420/150)
par(mar = c(0.8, 0.8, 0.8, 0.8))
draw_combined_legend(cutoff = cutoff, edge_width = edge_width)
dev.off()
    }
  )

  # Fit to view
  observeEvent(input$net_fit, {
    visNetworkProxy("network_plot") %>%
      visFit()
  })

  # Hide node when DOUBLE-CLICKED
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
  
  # Observer for [Show all] link in the yellow banner
  observeEvent(input$clear_one_hidden, {
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

  # Correlations table with Inspect button
  
  # Correlations table with Inspect button
  output$correlations_table <- renderDT({
    req(rv$results)

    if (is.null(rv$results$success) || !isTRUE(rv$results$success)) {
      return(DT::datatable(
        data.frame(Message = "Run analysis to populate correlations."),
        options = list(dom = 't'),
        rownames = FALSE
      ))
    }

    if (is.null(rv$results$correlations) || nrow(rv$results$correlations) == 0) {
      return(DT::datatable(
        data.frame(Message = "No correlations found (table is empty)."),
        options = list(dom = 't'),
        rownames = FALSE
      ))
    }

    # Always work on a plain data.frame to avoid any data.table edge cases in DT
    display_data <- as.data.frame(rv$results$correlations, stringsAsFactors = FALSE)

    # Add Inspect button
    display_data$Inspect <- sprintf(
      paste0(
        '<button class="btn btn-sm btn-outline-primary" ',
        'onclick="Shiny.setInputValue(\'inspect_row\', %d, {priority: \'event\'})">',
        'Inspect</button>'
      ),
      seq_len(nrow(display_data))
    )

    # Add By tissue button (always present; disabled if Model.csv not loaded)
    disabled_attr <- if (is.null(rv$full_metadata)) 'disabled="disabled" title="Upload Model.csv to enable By tissue"' else ''
    display_data$Oncotree <- sprintf(
      paste0(
        '<button class="btn btn-sm btn-outline-success" %s ',
        'onclick="Shiny.setInputValue(\'oncotree_row\', %d, {priority: \'event\'})">',
        'By tissue</button>'
      ),
      disabled_attr,
      seq_len(nrow(display_data))
    )

    # Select and order columns for display
    base_cols <- c("Gene1", "Gene2", "Correlation", "Slope", "N", "Cluster")
    available_base <- base_cols[base_cols %in% names(display_data)]
    display_data <- display_data[, c(available_base, "Inspect", "Oncotree"), drop = FALSE]

    DT::datatable(
      display_data,
      escape = FALSE,
      rownames = FALSE,
      options = list(
        pageLength = 10,
        autoWidth = TRUE,
        scrollX = TRUE,
        deferRender = TRUE,
        order = list(list(2, "desc"))
      ),
      filter = "top"
    )
  })
# Handle clicking on Inspect button
  observeEvent(input$inspect_row, {
    req(input$inspect_row)
    req(rv$reference_data)
    req(rv$results$success)
    
    selected_row <- input$inspect_row
    correlation_data <- rv$results$correlations[selected_row]
    
    gene1 <- correlation_data$Gene1
    gene2 <- correlation_data$Gene2
    cor_value <- correlation_data$Correlation
    
    # Use stored cell line names if available
    cell_line_names <- if (!is.null(rv$cell_line_names)) {
      rv$cell_line_names
    } else {
      paste0("Cell_", seq_len(nrow(rv$reference_data)))
    }
    
    # Also get cell line IDs if available (for Excel export)
    cell_line_ids <- if (!is.null(rv$cell_line_ids)) {
      rv$cell_line_ids
    } else {
      NULL
    }
    
    # Get original CellLineName (with hyphens etc.) for searching, if metadata available
    cell_line_original_names <- NULL
    cell_line_lineages <- NULL
    if (!is.null(rv$full_metadata) && !is.null(cell_line_ids)) {
      # Create lookup from full metadata (ACH ID -> CellLineName from column 3)
      orig_lookup <- setNames(
        as.character(rv$full_metadata[[3]]),  # CellLineName
        as.character(rv$full_metadata[[1]])   # ModelID
      )
      cell_line_original_names <- sapply(cell_line_ids, function(id) {
        result <- orig_lookup[id]
        if (is.na(result)) "" else result
      }, USE.NAMES = FALSE)
      
      # Also get lineage for highlighting filtered cells
      if ("OncotreeLineage" %in% colnames(rv$full_metadata)) {
        lineage_lookup <- setNames(
          as.character(rv$full_metadata$OncotreeLineage),
          as.character(rv$full_metadata[[1]])
        )
        cell_line_lineages <- sapply(cell_line_ids, function(id) {
          result <- lineage_lookup[id]
          if (is.na(result)) "" else result
        }, USE.NAMES = FALSE)
      }
    }
    
    # Extract data for these two genes with cell line names (use ALL data, not filtered)
    plot_data <- data.frame(
      cell_line_id = cell_line_ids,
      cell_line = cell_line_names,
      cell_line_original = cell_line_original_names,
      lineage = cell_line_lineages,
      OncotreeLineage = cell_line_lineages,  # Also store as OncotreeLineage for filter
      gene1_effect = rv$reference_data[[gene1]],
      gene2_effect = rv$reference_data[[gene2]],
      stringsAsFactors = FALSE
    )
    
    # Remove NA values (this filters all columns together)
    plot_data <- plot_data[complete.cases(plot_data[, c("gene1_effect", "gene2_effect")]), ]
    
    # Calculate correlation for ALL data
    cor_all <- round(cor(plot_data$gene1_effect, plot_data$gene2_effect, use = "complete.obs"), 3)
    
    # Check if a filter is active and calculate filtered correlation
    filter_active <- FALSE
    filter_label <- ""
    cor_filtered <- NULL
    filtered_lineages <- NULL
    
    if (!is.null(rv$full_metadata) && 
        length(input$oncotree_lineage_filter) > 0 && 
        any(input$oncotree_lineage_filter != "")) {
      
      filter_active <- TRUE
      filtered_lineages <- input$oncotree_lineage_filter[input$oncotree_lineage_filter != ""]
      filter_label <- paste(filtered_lineages, collapse = ", ")
      
      # Mark which cells are in the filter
      plot_data$in_filter <- plot_data$lineage %in% filtered_lineages
      
      # Calculate correlation for filtered subset
      filtered_data <- plot_data[plot_data$in_filter, ]
      if (nrow(filtered_data) >= 3) {
        cor_filtered <- round(cor(filtered_data$gene1_effect, filtered_data$gene2_effect, use = "complete.obs"), 3)
      }
    } else {
      plot_data$in_filter <- FALSE
    }
    
    # Store data for reactive plotting
    rv$current_plot_data <- plot_data
    rv$current_scatter_genes <- c(gene1, gene2)
    rv$current_cor_value <- cor_value  # This is from the analysis (may be filtered)
    rv$current_cor_all <- cor_all
    rv$current_cor_filtered <- cor_filtered
    rv$current_filter_active <- filter_active
    rv$current_filter_label <- filter_label
    rv$current_filtered_lineages <- filtered_lineages
    
    # Initialize plot settings with auto ranges
    xlim <- range(plot_data$gene1_effect, na.rm = TRUE)
    ylim <- range(plot_data$gene2_effect, na.rm = TRUE)
    rv$scatter_xlim <- xlim
    rv$scatter_ylim <- ylim
    rv$scatter_width <- 8
    rv$scatter_height <- 8  # Square by default
    rv$scatter_cell_search <- ""
    rv$scatter_plot_trigger <- 0  # Trigger for plot updates
    
    # Build subtitle with correlation info
    # Calculate N for filtered data
    n_all <- nrow(plot_data)
    n_filtered <- if (filter_active) sum(plot_data$in_filter) else n_all
    
    # Calculate slope for all data
    if (n_all >= 3) {
      lm_all <- lm(gene2_effect ~ gene1_effect, data = plot_data)
      slope_all <- round(coef(lm_all)[2], 3)
    } else {
      slope_all <- NA
    }
    
    # Calculate slope for filtered data if applicable
    if (filter_active && n_filtered >= 3) {
      filtered_plot_data <- plot_data[plot_data$in_filter, ]
      lm_filtered <- lm(gene2_effect ~ gene1_effect, data = filtered_plot_data)
      slope_filtered <- round(coef(lm_filtered)[2], 3)
    } else {
      slope_filtered <- NULL
    }
    
    # Build title with N counts, r, slope, and warning for small samples
    if (filter_active && !is.null(cor_filtered)) {
      modal_title <- paste0("Correlation: ", gene1, " vs ", gene2, 
                           " | All: r=", cor_all, ", slope=", slope_all, " (n=", n_all, ")",
                           " | ", filter_label, ": r=", cor_filtered, ", slope=", slope_filtered, " (n=", n_filtered, ")",
                           if(n_filtered < 10) " ⚠️ Small sample!" else "")
    } else {
      modal_title <- paste0("Correlation: ", gene1, " vs ", gene2, 
                           " | r=", cor_all, ", slope=", slope_all, " (n=", n_all, ")",
                           if(n_all < 10) " ⚠️ Small sample!" else "")
    }
    
    # Show modal with plot and controls - using custom CSS for extra width
    showModal(modalDialog(
      title = modal_title,
      size = "xl",
      easyClose = TRUE,
      footer = tagList(
        actionButton("scatter_by_tissue", "By tissue", class = "btn btn-info"),
        downloadButton("download_scatter_png", "Download PNG"),
        downloadButton("download_scatter_svg", "Download SVG"),
        downloadButton("download_scatter_csv", "Download CSV"),
        modalButton("Close")
      ),
      tags$style(HTML("
        .modal-xl { 
          max-width: 1400px !important; 
          width: 95% !important;
        }
        .control-box {
          background: #f9fafb;
          padding: 10px 12px;
          border-radius: 8px;
          border: 1px solid #e5e7eb;
          margin-bottom: 10px;
        }
        .control-box-compact {
          background: #f9fafb;
          padding: 8px 10px;
          border-radius: 6px;
          border: 1px solid #e5e7eb;
          margin-bottom: 8px;
        }
        .control-label {
          font-size: 12px;
          font-weight: 600;
          color: #374151;
          margin-bottom: 6px;
          display: block;
        }
        .control-input {
          width: 100% !important;
          font-size: 13px !important;
          padding: 6px 8px !important;
          height: 32px !important;
          text-align: center !important;
          border: 1.5px solid #d1d5db !important;
          border-radius: 6px !important;
          background: white !important;
        }
        .control-input:focus {
          border-color: #16a34a !important;
          outline: none !important;
          box-shadow: 0 0 0 3px rgba(22, 163, 74, 0.1) !important;
        }
        .control-input::-webkit-inner-spin-button,
        .control-input::-webkit-outer-spin-button {
          -webkit-appearance: none;
          margin: 0;
        }
        .control-input {
          -moz-appearance: textfield;
        }
        .axis-label {
          font-size: 11px;
          color: #6b7280;
          margin-bottom: 2px;
        }
      ")),
      fluidRow(
        column(4,
          div(style = "padding-right: 15px;",
            # Axis Ranges - X and Y on same row
            div(class = "control-box-compact",
              div(class = "control-label", "Axis Ranges"),
              fluidRow(
                column(6,
                  tags$div(class = "axis-label", "X-axis"),
                  fluidRow(
                    column(6, style = "padding-right: 2px;",
                      tags$input(
                        id = "scatter_xmin",
                        type = "text",
                        class = "control-input",
                        value = sprintf("%.1f", round(xlim[1], 1)),
                        pattern = "-?[0-9]*\\.?[0-9]*",
                        title = "X min",
                        onchange = "Shiny.setInputValue('scatter_xmin', parseFloat(this.value.replace(',', '.')))",
                        oninput = "this.value = this.value.replace(',', '.')"
                      )
                    ),
                    column(6, style = "padding-left: 2px;",
                      tags$input(
                        id = "scatter_xmax",
                        type = "text",
                        class = "control-input",
                        value = sprintf("%.1f", round(xlim[2], 1)),
                        pattern = "-?[0-9]*\\.?[0-9]*",
                        title = "X max",
                        onchange = "Shiny.setInputValue('scatter_xmax', parseFloat(this.value.replace(',', '.')))",
                        oninput = "this.value = this.value.replace(',', '.')"
                      )
                    )
                  )
                ),
                column(6,
                  tags$div(class = "axis-label", "Y-axis"),
                  fluidRow(
                    column(6, style = "padding-right: 2px;",
                      tags$input(
                        id = "scatter_ymin",
                        type = "text",
                        class = "control-input",
                        value = sprintf("%.1f", round(ylim[1], 1)),
                        pattern = "-?[0-9]*\\.?[0-9]*",
                        title = "Y min",
                        onchange = "Shiny.setInputValue('scatter_ymin', parseFloat(this.value.replace(',', '.')))",
                        oninput = "this.value = this.value.replace(',', '.')"
                      )
                    ),
                    column(6, style = "padding-left: 2px;",
                      tags$input(
                        id = "scatter_ymax",
                        type = "text",
                        class = "control-input",
                        value = sprintf("%.1f", round(ylim[2], 1)),
                        pattern = "-?[0-9]*\\.?[0-9]*",
                        title = "Y max",
                        onchange = "Shiny.setInputValue('scatter_ymax', parseFloat(this.value.replace(',', '.')))",
                        oninput = "this.value = this.value.replace(',', '.')"
                      )
                    )
                  )
                )
              ),
              actionButton("scatter_reset_axes", "Reset", class = "btn-sm btn-outline-secondary", style = "width: 100%; margin-top: 6px; padding: 4px;")
            ),

            # Filter by Cancer Type (new)
            if (!is.null(rv$full_metadata) && "OncotreeLineage" %in% colnames(rv$full_metadata)) {
              lineage_choices <- c("All cancer types" = "", sort(unique(rv$full_metadata$OncotreeLineage[rv$full_metadata$OncotreeLineage != ""])))
              div(class = "control-box-compact",
                div(class = "control-label", "Filter by Cancer Type"),
                selectInput("scatter_cancer_filter", NULL,
                            choices = lineage_choices,
                            selected = "",
                            width = "100%")
              )
            },

            # Cell Line Highlighting - merged text search and click identify
            div(class = "control-box-compact",
              div(class = "control-label", "Highlight Cell Lines"),
              fluidRow(
                column(8,
                  textAreaInput("scatter_cell_search", NULL,
                                value = "", height = "60px",
                                placeholder = "Type names or click plot...")
                ),
                column(4,
                  numericInput("highlight_label_size", "Font", value = 3, min = 1, max = 6, step = 0.5, width = "70px")
                )
              ),
              div(style = "margin-top: 4px;",
                actionButton("clear_highlights", "Clear All", class = "btn-xs btn-outline-secondary", style = "padding: 2px 8px; font-size: 10px;")
              ),
              tags$small(style = "color: #6b7280; font-size: 10px; margin-top: 4px; display: block;",
                if (is.null(rv$cell_line_name_lookup)) "Enter ACH IDs or click points on plot" else "Name, ACH ID, or click points")
            ),

            # Hotspot mutation overlay - with color by mutation mode
            div(class = "control-box",
              div(class = "control-label", "Hotspot Mutation Overlay"),
              selectizeInput(
                "hotspot_gene",
                NULL,
                choices = NULL,
                selected = "",
                options = list(
                  placeholder = "Type gene (e.g., TP53)...",
                  maxOptions = 50
                )
              ),
              tags$small(style = "color: #6b7280; font-size: 10px; margin-top: -10px; display: block;", 
                "Common mutations shown; type to search all genes"),
              selectInput("hotspot_mode", NULL,
                          choices = list(
                            "None" = "none",
                            "Color by hotspot mutation (0/1/2)" = "color_by_mut",
                            "3-panel (0 / 1 / 2 hotspot mutations)" = "three_panel",
                            "Highlight 1+2 hotspot mutations" = "12",
                            "Highlight 1 hotspot mutation" = "1",
                            "Highlight 2 hotspot mutations" = "2",
                            "Compare 0 vs 1+2 (table)" = "compare_table"
                          ),
                          selected = "color_by_mut"),
              conditionalPanel(
                condition = "input.hotspot_mode == 'three_panel'",
                radioButtons("pie_scale_mode", "Pie chart scaling:", 
                            choices = list("By count (proportional)" = "count", 
                                         "Equal size (100%)" = "percent"),
                            selected = "count", inline = TRUE),
                selectInput("pie_highlight_cancer", "Highlight cancer type:",
                           choices = c("None" = ""),
                           selected = "",
                           width = "100%")
              ),
              tags$small(style = "color: #6b7280; font-size: 10px;", 
                "Upload OmicsSomaticMutationsMatrixHotspot.csv to enable")
            )
          )
        ),

        column(8,
          # Conditional: show table for compare_table mode, otherwise show plot
          conditionalPanel(
            condition = "input.hotspot_mode == 'compare_table'",
            div(style = "padding: 10px;",
              h4("Mutation Effect on Correlation by Cancer Type"),
              p(style = "color: #6b7280; font-size: 12px;", 
                "Comparing correlation between WT (0 mutations) vs Mutant (1+2 mutations) cells, stratified by cancer type."),
              DTOutput("mutation_compare_table"),
              br(),
              downloadButton("download_mutation_compare", "Download CSV", class = "btn-sm btn-outline-success")
            )
          ),
          conditionalPanel(
            condition = "input.hotspot_mode != 'compare_table'",
            div(
              verbatimTextOutput("clicked_cell_info", placeholder = TRUE),
              plotOutput("scatter_plot_modal", height = "620px", click = "scatter_click")
            )
          )
        )
      )
    ))
    
    # Initialize clicked cells storage
    if (is.null(rv$clicked_cells)) rv$clicked_cells <- character(0)

    # Populate hotspot gene choices (server-side) once the modal exists
    if (!is.null(rv$hotspot_genes) && length(rv$hotspot_genes) > 0) {
      # Use clean gene names for display, but keep original column names as values
      if (!is.null(rv$hotspot_genes_clean) && !is.null(rv$hotspot_gene_counts)) {
        # Create named vector: display "TP53 (847 mut)" to indicate mutation count in file
        display_names <- paste0(rv$hotspot_genes_clean, " (", rv$hotspot_gene_counts, " mut)")
        choices <- setNames(rv$hotspot_genes, display_names)
      } else if (!is.null(rv$hotspot_genes_clean)) {
        choices <- setNames(rv$hotspot_genes, rv$hotspot_genes_clean)
      } else {
        choices <- rv$hotspot_genes
      }
      updateSelectizeInput(session, "hotspot_gene", choices = choices, selected = "", server = TRUE)
    } else {
      updateSelectizeInput(session, "hotspot_gene", choices = character(0), selected = "", server = TRUE)
    }
    
    # Populate cancer type highlight choices for pie charts
    if (!is.null(rv$current_plot_data) && "OncotreeLineage" %in% names(rv$current_plot_data)) {
      cancer_types <- sort(unique(rv$current_plot_data$OncotreeLineage[!is.na(rv$current_plot_data$OncotreeLineage) & rv$current_plot_data$OncotreeLineage != ""]))
      updateSelectInput(session, "pie_highlight_cancer", 
                       choices = c("None" = "", cancer_types),
                       selected = "")
    }
  })
  
  # Handle click on scatter plot to identify cell lines (always enabled now)
  observeEvent(input$scatter_click, {
    req(rv$current_plot_data)
    
    # Don't handle clicks in 3-panel mode (faceted plots don't work with clicks)
    if (isTRUE(input$hotspot_mode == "three_panel")) return()
    
    click <- input$scatter_click
    plot_data <- rv$current_plot_data
    
    # Find nearest point to click
    distances <- sqrt((plot_data$gene1_effect - click$x)^2 + (plot_data$gene2_effect - click$y)^2)
    nearest_idx <- which.min(distances)
    
    # Only identify if click is reasonably close (within 0.2 units)
    if (distances[nearest_idx] < 0.2) {
      cell_name <- plot_data$cell_line[nearest_idx]
      cell_id <- plot_data$cell_line_id[nearest_idx]
      cancer_type <- if ("OncotreeLineage" %in% names(plot_data)) plot_data$OncotreeLineage[nearest_idx] else ""
      
      # Create label with cancer type
      cell_label <- if (!is.null(cancer_type) && !is.na(cancer_type) && cancer_type != "") {
        paste0(cell_name, " (", cancer_type, ")")
      } else {
        cell_name
      }
      
      # Toggle - if already in list, remove it; otherwise add it
      # Store both the simple name (for matching) and the full label (for display)
      if (cell_name %in% names(rv$clicked_cells)) {
        rv$clicked_cells <- rv$clicked_cells[names(rv$clicked_cells) != cell_name]
      } else {
        new_entry <- setNames(cell_label, cell_name)
        rv$clicked_cells <- c(rv$clicked_cells, new_entry)
      }
      
      # Update the highlight in plot data
      rv$scatter_plot_trigger <- Sys.time()
    }
  })
  
  # Clear all highlights button (merged clear function)
  observeEvent(input$clear_highlights, {
    rv$clicked_cells <- character(0)
    updateTextAreaInput(session, "scatter_cell_search", value = "")
    rv$scatter_plot_trigger <- Sys.time()
  })
  
  # Output clicked cell info
  output$clicked_cell_info <- renderText({
    # Check if in 3-panel mode
    is_three_panel <- isTRUE(input$hotspot_mode == "three_panel")
    
    # Count text search matches
    search_terms <- if (!is.null(input$scatter_cell_search) && input$scatter_cell_search != "") {
      trimws(strsplit(input$scatter_cell_search, "\n")[[1]])
    } else {
      character(0)
    }
    search_terms <- search_terms[search_terms != ""]
    
    if (length(rv$clicked_cells) == 0 && length(search_terms) == 0) {
      if (is_three_panel) {
        "3-panel mode: Type cell line names in the box on the left to highlight them."
      } else {
        "Click points to highlight, or type names in the search box. Click again to remove."
      }
    } else {
      parts <- c()
      if (length(rv$clicked_cells) > 0) {
        parts <- c(parts, paste("Clicked:", paste(rv$clicked_cells, collapse = ", ")))
      }
      if (length(search_terms) > 0) {
        parts <- c(parts, paste("Search:", paste(search_terms, collapse = ", ")))
      }
      paste(parts, collapse = " | ")
    }
  })
  
  # Handle clicking on Oncotree breakdown button
  observeEvent(input$oncotree_row, {
    req(input$oncotree_row)
    req(rv$reference_data)
    req(rv$results$success)
    req(rv$full_metadata)
    
    selected_row <- input$oncotree_row
    correlation_data <- rv$results$correlations[selected_row]
    
    gene1 <- correlation_data$Gene1
    gene2 <- correlation_data$Gene2
    
    # Get cell line IDs
    cell_line_ids <- rv$cell_line_ids
    if (is.null(cell_line_ids)) {
      showNotification("Cell line IDs not available", type = "error")
      return()
    }
    
    # Create data frame with gene effects and metadata
    breakdown_data <- data.frame(
      ModelID = cell_line_ids,
      gene1_effect = rv$reference_data[[gene1]],
      gene2_effect = rv$reference_data[[gene2]],
      stringsAsFactors = FALSE
    )
    
    # Merge with metadata
    meta_cols <- c("ModelID", "OncotreeLineage", "OncotreePrimaryDisease")
    meta_subset <- rv$full_metadata[, intersect(meta_cols, colnames(rv$full_metadata)), drop = FALSE]
    breakdown_data <- merge(breakdown_data, meta_subset, by = "ModelID", all.x = TRUE)
    
    # Remove rows with missing gene effects
    breakdown_data <- breakdown_data[complete.cases(breakdown_data[, c("gene1_effect", "gene2_effect")]), ]
    
    # Calculate statistics by lineage
    lineage_stats <- breakdown_data %>%
      dplyr::group_by(OncotreeLineage) %>%
      dplyr::summarise(
        n = dplyr::n(),
        correlation = if(dplyr::n() >= 3) round(cor(gene1_effect, gene2_effect, use = "complete.obs"), 3) else NA,
        mean_gene1 = round(mean(gene1_effect, na.rm = TRUE), 3),
        sd_gene1 = round(sd(gene1_effect, na.rm = TRUE), 3),
        mean_gene2 = round(mean(gene2_effect, na.rm = TRUE), 3),
        sd_gene2 = round(sd(gene2_effect, na.rm = TRUE), 3),
        .groups = "drop"
      ) %>%
      dplyr::filter(!is.na(OncotreeLineage) & OncotreeLineage != "") %>%
      dplyr::arrange(dplyr::desc(n))
    
    # Store for plotting
    rv$oncotree_breakdown_data <- breakdown_data
    rv$oncotree_lineage_stats <- lineage_stats
    rv$oncotree_genes <- c(gene1, gene2)
    
    # Show modal
    showModal(modalDialog(
      title = paste("Tissue Breakdown:", gene1, "vs", gene2),
      size = "xl",
      easyClose = TRUE,
      footer = tagList(
        downloadButton("download_oncotree_breakdown", "Download CSV"),
        modalButton("Close")
      ),
      
      # Custom CSS for wider modal
      tags$style(HTML("
        .modal-xl { max-width: 1400px !important; width: 95% !important; }
        .oncotree-table .dataTables_wrapper { font-size: 12px; }
      ")),
      
      fluidRow(
        column(5,
          h4("Correlation by Tissue Type"),
          downloadButton("download_oncotree_bar_png", "PNG", class = "btn-sm btn-success"),
          downloadButton("download_oncotree_bar_svg", "SVG", class = "btn-sm btn-success"),
          br(), br(),
          plotOutput("oncotree_bar_plot", height = "550px")
        ),
        column(7,
          h4("Statistics by Lineage"),
          div(class = "oncotree-table", style = "max-height: 550px; overflow-y: auto; overflow-x: auto;",
            DTOutput("oncotree_stats_table")
          )
        )
      ),
      
      hr(),
      
      fluidRow(
        column(12,
          h4("Scatter Plot by Selected Lineage"),
          fluidRow(
            column(4,
              selectInput("oncotree_scatter_lineage", "Select lineage to view:",
                         choices = c("All lineages", as.character(lineage_stats$OncotreeLineage)),
                         selected = "All lineages",
                         width = "100%")
            ),
            column(4,
              br(),
              downloadButton("download_oncotree_scatter_png", "PNG", class = "btn-sm btn-success"),
              downloadButton("download_oncotree_scatter_svg", "SVG", class = "btn-sm btn-success")
            )
          ),
          plotOutput("oncotree_scatter_plot", height = "400px")
        )
      )
    ))
  })
  
  # Render Oncotree bar plot
  output$oncotree_bar_plot <- renderPlot({
    req(rv$oncotree_lineage_stats)
    
    stats <- rv$oncotree_lineage_stats
    stats <- stats[!is.na(stats$correlation), ]
    
    if (nrow(stats) == 0) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, label = "No data available") +
               theme_void())
    }
    
    # Order by correlation
    stats$OncotreeLineage <- factor(stats$OncotreeLineage, 
                                     levels = stats$OncotreeLineage[order(stats$correlation)])
    
    # Create label with n
    stats$label <- paste0("n=", stats$n)
    
    ggplot(stats, aes(x = OncotreeLineage, y = correlation, fill = correlation)) +
      geom_col(width = 0.7) +
      geom_text(aes(label = label, 
                    y = ifelse(correlation >= 0, correlation + 0.03, correlation - 0.03)),
                hjust = ifelse(stats$correlation >= 0, 0, 1), 
                size = 3) +
      scale_fill_gradient2(low = "#e53e3e", mid = "#f5f5f5", high = "#16a34a", 
                           midpoint = 0, limits = c(-1, 1)) +
      coord_flip(clip = "off") +
      labs(
        title = paste(rv$oncotree_genes[1], "vs", rv$oncotree_genes[2]),
        x = NULL,
        y = "Correlation",
        fill = "r"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        legend.position = "right",
        plot.title = element_text(face = "bold", size = 13),
        plot.margin = margin(10, 40, 10, 10),  # Extra right margin for labels
        axis.text.y = element_text(size = 9)
      ) +
      scale_y_continuous(limits = c(min(stats$correlation, -0.1) - 0.1, 
                                    max(stats$correlation, 0.1) + 0.15),
                        expand = c(0, 0))
  })
  
  # Render Oncotree stats table
  output$oncotree_stats_table <- renderDT({
    req(rv$oncotree_lineage_stats)
    
    stats <- rv$oncotree_lineage_stats
    gene1 <- rv$oncotree_genes[1]
    gene2 <- rv$oncotree_genes[2]
    
    # Format for display - clear column names
    display_stats <- data.frame(
      Lineage = stats$OncotreeLineage,
      N = stats$n,
      Correlation = stats$correlation,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    
    # Add gene effect columns with clear names
    display_stats[[paste0(gene1, " Effect (mean)")]] <- stats$mean_gene1
    display_stats[[paste0(gene1, " Effect (SD)")]] <- stats$sd_gene1
    display_stats[[paste0(gene2, " Effect (mean)")]] <- stats$mean_gene2
    display_stats[[paste0(gene2, " Effect (SD)")]] <- stats$sd_gene2
    
    datatable(display_stats,
              options = list(
                pageLength = 25, 
                dom = 't',
                scrollX = TRUE,
                autoWidth = FALSE,
                order = list(list(2, 'desc')),  # Sort by Correlation (column index 2) descending
                columnDefs = list(
                  list(width = '100px', targets = 0),   # Lineage
                  list(width = '35px', targets = 1),    # N
                  list(width = '70px', targets = 2),    # Correlation
                  list(width = '70px', targets = 3:6)   # gene columns
                )
              ),
              rownames = FALSE,
              class = 'compact stripe') %>%
      formatRound(columns = 'Correlation', digits = 2) %>%
      formatRound(columns = c(4:7), digits = 2) %>%
      formatStyle('Correlation',
                 backgroundColor = styleInterval(
                   c(-0.3, 0, 0.3),
                   c("#fee2e2", "#fef3c7", "#dcfce7", "#bbf7d0")
                 ))
  })
  
  # Render Oncotree scatter plot
  output$oncotree_scatter_plot <- renderPlot({
    req(rv$oncotree_breakdown_data)
    req(input$oncotree_scatter_lineage)
    
    data <- rv$oncotree_breakdown_data
    gene1 <- rv$oncotree_genes[1]
    gene2 <- rv$oncotree_genes[2]
    
    # Filter by selected lineage
    if (input$oncotree_scatter_lineage != "All lineages") {
      data <- data[data$OncotreeLineage == input$oncotree_scatter_lineage, ]
      subtitle_text <- paste0(input$oncotree_scatter_lineage, " (n = ", nrow(data), ")")
    } else {
      subtitle_text <- paste0("All lineages (n = ", nrow(data), ")")
    }
    
    if (nrow(data) < 2) {
      return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, label = "Not enough data points") +
               theme_void())
    }
    
    # Calculate correlation for this subset
    cor_val <- round(cor(data$gene1_effect, data$gene2_effect, use = "complete.obs"), 3)
    
    ggplot(data, aes(x = gene1_effect, y = gene2_effect)) +
      geom_point(alpha = 0.6, color = "#16a34a", size = 2) +
      geom_smooth(method = "lm", se = TRUE, color = "#059669", fill = "#86efac") +
      labs(
        title = paste(gene1, "vs", gene2),
        subtitle = paste0(subtitle_text, ", r = ", cor_val),
        x = paste(gene1, "Effect"),
        y = paste(gene2, "Effect")
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold"),
        aspect.ratio = 1
      )
  })
  
  # Download Oncotree breakdown
  output$download_oncotree_breakdown <- downloadHandler(
    filename = function() {
      paste0(rv$oncotree_genes[1], "_vs_", rv$oncotree_genes[2], "_tissue_breakdown_", 
             format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      req(rv$oncotree_lineage_stats)
      fwrite(rv$oncotree_lineage_stats, file)
    }
  )
  
  # Download oncotree bar plot as PNG
  output$download_oncotree_bar_png <- downloadHandler(
    filename = function() {
      paste0(rv$oncotree_genes[1], "_vs_", rv$oncotree_genes[2], "_by_tissue_", 
             format(Sys.Date(), "%Y%m%d"), ".png")
    },
    content = function(file) {
      req(rv$oncotree_lineage_stats)
      
      stats <- rv$oncotree_lineage_stats
      stats <- stats[!is.na(stats$correlation), ]
      stats$OncotreeLineage <- factor(stats$OncotreeLineage, 
                                       levels = stats$OncotreeLineage[order(stats$correlation)])
      stats$label <- paste0("n=", stats$n)
      
      p <- ggplot(stats, aes(x = OncotreeLineage, y = correlation, fill = correlation)) +
        geom_col(width = 0.7) +
        geom_text(aes(label = label, 
                      y = ifelse(correlation >= 0, correlation + 0.03, correlation - 0.03)),
                  hjust = ifelse(stats$correlation >= 0, 0, 1), 
                  size = 3) +
        scale_fill_gradient2(low = "#e53e3e", mid = "#f5f5f5", high = "#16a34a", 
                             midpoint = 0, limits = c(-1, 1)) +
        coord_flip(clip = "off") +
        labs(
          title = paste(rv$oncotree_genes[1], "vs", rv$oncotree_genes[2], "- By Tissue"),
          x = NULL,
          y = "Correlation",
          fill = "r"
        ) +
        theme_minimal(base_size = 11) +
        theme(
          legend.position = "right",
          plot.title = element_text(face = "bold", size = 13),
          plot.margin = margin(10, 40, 10, 10),
          axis.text.y = element_text(size = 9)
        ) +
        scale_y_continuous(limits = c(min(stats$correlation, -0.1) - 0.1, 
                                      max(stats$correlation, 0.1) + 0.15),
                          expand = c(0, 0))
      
      ggsave(file, plot = p, width = 8, height = max(6, nrow(stats) * 0.3), dpi = 300)
    }
  )
  
  # Download oncotree bar plot as SVG
  output$download_oncotree_bar_svg <- downloadHandler(
    filename = function() {
      paste0(rv$oncotree_genes[1], "_vs_", rv$oncotree_genes[2], "_by_tissue_", 
             format(Sys.Date(), "%Y%m%d"), ".svg")
    },
    content = function(file) {
      req(rv$oncotree_lineage_stats)
      
      stats <- rv$oncotree_lineage_stats
      stats <- stats[!is.na(stats$correlation), ]
      stats$OncotreeLineage <- factor(stats$OncotreeLineage, 
                                       levels = stats$OncotreeLineage[order(stats$correlation)])
      stats$label <- paste0("n=", stats$n)
      
      p <- ggplot(stats, aes(x = OncotreeLineage, y = correlation, fill = correlation)) +
        geom_col(width = 0.7) +
        geom_text(aes(label = label, 
                      y = ifelse(correlation >= 0, correlation + 0.03, correlation - 0.03)),
                  hjust = ifelse(stats$correlation >= 0, 0, 1), 
                  size = 3) +
        scale_fill_gradient2(low = "#e53e3e", mid = "#f5f5f5", high = "#16a34a", 
                             midpoint = 0, limits = c(-1, 1)) +
        coord_flip(clip = "off") +
        labs(
          title = paste(rv$oncotree_genes[1], "vs", rv$oncotree_genes[2], "- By Tissue"),
          x = NULL,
          y = "Correlation",
          fill = "r"
        ) +
        theme_minimal(base_size = 11) +
        theme(
          legend.position = "right",
          plot.title = element_text(face = "bold", size = 13),
          plot.margin = margin(10, 40, 10, 10),
          axis.text.y = element_text(size = 9)
        ) +
        scale_y_continuous(limits = c(min(stats$correlation, -0.1) - 0.1, 
                                      max(stats$correlation, 0.1) + 0.15),
                          expand = c(0, 0))
      
      ggsave(file, plot = p, width = 8, height = max(6, nrow(stats) * 0.3), device = "svg")
    }
  )
  
  # Download oncotree scatter plot as PNG
  output$download_oncotree_scatter_png <- downloadHandler(
    filename = function() {
      lineage <- if (!is.null(input$oncotree_scatter_lineage) && input$oncotree_scatter_lineage != "All lineages") {
        paste0("_", gsub("[^A-Za-z0-9]", "_", input$oncotree_scatter_lineage))
      } else {
        ""
      }
      paste0(rv$oncotree_genes[1], "_vs_", rv$oncotree_genes[2], lineage, "_scatter_", 
             format(Sys.Date(), "%Y%m%d"), ".png")
    },
    content = function(file) {
      req(rv$oncotree_breakdown_data)
      
      data <- rv$oncotree_breakdown_data
      
      # Filter by selected lineage
      if (!is.null(input$oncotree_scatter_lineage) && input$oncotree_scatter_lineage != "All lineages") {
        data <- data[data$OncotreeLineage == input$oncotree_scatter_lineage, ]
      }
      
      if (nrow(data) >= 3) {
        cor_val <- cor(data$gene1_effect, data$gene2_effect, use = "complete.obs")
        title_text <- paste0(rv$oncotree_genes[1], " vs ", rv$oncotree_genes[2],
                             ifelse(input$oncotree_scatter_lineage != "All lineages", 
                                    paste0(" (", input$oncotree_scatter_lineage, ")"), ""),
                             " | r = ", round(cor_val, 3), " | n = ", nrow(data))
      } else {
        title_text <- paste0(rv$oncotree_genes[1], " vs ", rv$oncotree_genes[2])
      }
      
      p <- ggplot(data, aes(x = gene1_effect, y = gene2_effect)) +
        geom_point(color = "#16a34a", alpha = 0.6, size = 2) +
        geom_smooth(method = "lm", se = TRUE, color = "#3182ce", linewidth = 1, alpha = 0.2) +
        labs(
          title = title_text,
          x = paste(rv$oncotree_genes[1], "Effect"),
          y = paste(rv$oncotree_genes[2], "Effect")
        ) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, color = "gray70")
        )
      
      ggsave(file, plot = p, width = 8, height = 6, dpi = 300)
    }
  )
  
  # Download oncotree scatter plot as SVG
  output$download_oncotree_scatter_svg <- downloadHandler(
    filename = function() {
      lineage <- if (!is.null(input$oncotree_scatter_lineage) && input$oncotree_scatter_lineage != "All lineages") {
        paste0("_", gsub("[^A-Za-z0-9]", "_", input$oncotree_scatter_lineage))
      } else {
        ""
      }
      paste0(rv$oncotree_genes[1], "_vs_", rv$oncotree_genes[2], lineage, "_scatter_", 
             format(Sys.Date(), "%Y%m%d"), ".svg")
    },
    content = function(file) {
      req(rv$oncotree_breakdown_data)
      
      data <- rv$oncotree_breakdown_data
      
      # Filter by selected lineage
      if (!is.null(input$oncotree_scatter_lineage) && input$oncotree_scatter_lineage != "All lineages") {
        data <- data[data$OncotreeLineage == input$oncotree_scatter_lineage, ]
      }
      
      if (nrow(data) >= 3) {
        cor_val <- cor(data$gene1_effect, data$gene2_effect, use = "complete.obs")
        title_text <- paste0(rv$oncotree_genes[1], " vs ", rv$oncotree_genes[2],
                             ifelse(input$oncotree_scatter_lineage != "All lineages", 
                                    paste0(" (", input$oncotree_scatter_lineage, ")"), ""),
                             " | r = ", round(cor_val, 3), " | n = ", nrow(data))
      } else {
        title_text <- paste0(rv$oncotree_genes[1], " vs ", rv$oncotree_genes[2])
      }
      
      p <- ggplot(data, aes(x = gene1_effect, y = gene2_effect)) +
        geom_point(color = "#16a34a", alpha = 0.6, size = 2) +
        geom_smooth(method = "lm", se = TRUE, color = "#3182ce", linewidth = 1, alpha = 0.2) +
        labs(
          title = title_text,
          x = paste(rv$oncotree_genes[1], "Effect"),
          y = paste(rv$oncotree_genes[2], "Effect")
        ) +
        theme_minimal(base_size = 12) +
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(fill = NA, color = "gray70")
        )
      
      ggsave(file, plot = p, width = 8, height = 6, device = "svg")
    }
  )
  
  # Update plot when search box changes or axes change
  observeEvent(c(input$scatter_cell_search, input$scatter_xmin, input$scatter_xmax, 
                 input$scatter_ymin, input$scatter_ymax, input$scatter_cancer_filter,
                 input$hotspot_gene, input$hotspot_mode,
                 input$pie_scale_mode, input$highlight_label_size, input$pie_highlight_cancer), {
    rv$scatter_plot_trigger <- rv$scatter_plot_trigger + 1
  }, ignoreInit = TRUE)
  
  # Reset axes button
  observeEvent(input$scatter_reset_axes, {
    req(rv$current_plot_data)
    xlim <- range(rv$current_plot_data$gene1_effect, na.rm = TRUE)
    ylim <- range(rv$current_plot_data$gene2_effect, na.rm = TRUE)
    
    updateNumericInput(session, "scatter_xmin", value = round(xlim[1], 1))
    updateNumericInput(session, "scatter_xmax", value = round(xlim[2], 1))
    updateNumericInput(session, "scatter_ymin", value = round(ylim[1], 1))
    updateNumericInput(session, "scatter_ymax", value = round(ylim[2], 1))
    
    rv$scatter_plot_trigger <- rv$scatter_plot_trigger + 1
  })
  
  # Render scatter plot in modal with all controls
  output$scatter_plot_modal <- renderPlot({
    req(rv$current_plot_data)
    req(rv$current_scatter_genes)
    
    # Trigger dependency
    rv$scatter_plot_trigger
    
    plot_data <- rv$current_plot_data

    # Apply cancer type filter from scatter_cancer_filter if set
    if (!is.null(input$scatter_cancer_filter) && nzchar(input$scatter_cancer_filter) && 
        "OncotreeLineage" %in% colnames(plot_data)) {
      plot_data <- plot_data[plot_data$OncotreeLineage == input$scatter_cancer_filter, ]
    }

    # Optional hotspot mutation overlay (loaded from OmicsSomaticMutationsMatrixHotspot.csv)
    plot_data$hotspot_status <- NA_integer_
    plot_data$hotspot_highlight <- FALSE
    hotspot_gene_input <- if (!is.null(input$hotspot_gene)) trimws(input$hotspot_gene) else ""
    hotspot_mode <- if (!is.null(input$hotspot_mode)) input$hotspot_mode else "none"
    
    # The dropdown value IS the column name (e.g., "TP53 (7157)"), display shows "TP53 (n=1409)"
    hotspot_col <- hotspot_gene_input
    
    # Extract just the gene symbol for display purposes (remove parenthetical part)
    hotspot_gene <- gsub("\\s*\\([^)]+\\)$", "", hotspot_gene_input)

    if (!is.null(rv$hotspot_file_path) && nzchar(hotspot_col) && hotspot_mode != "none") {

      if (is.null(rv$hotspot_id_col) || !nzchar(rv$hotspot_id_col)) {
        showNotification("Hotspot file loaded without a detectable ID column. Please re-upload the hotspot file.", type = "error")
        hotspot_mode <- "none"
      } else if (is.null(rv$hotspot_genes) || !(hotspot_col %in% rv$hotspot_genes)) {
        showNotification(paste0("Hotspot gene '", hotspot_gene, "' not found in hotspot file."), type = "warning")
        hotspot_mode <- "none"
      }
    }

    if (!is.null(rv$hotspot_file_path) && nzchar(hotspot_col) && hotspot_mode != "none") {
      # Read only ID + the requested gene column to keep memory low

hs_dt <- tryCatch({
        data.table::fread(rv$hotspot_file_path,
                          select = c(rv$hotspot_id_col, hotspot_col),
                          showProgress = FALSE)
      }, error = function(e) NULL)

      if (!is.null(hs_dt) && !is.null(rv$hotspot_id_col) && rv$hotspot_id_col %in% names(hs_dt) && hotspot_col %in% names(hs_dt)) {
        data.table::setnames(hs_dt, c(rv$hotspot_id_col, hotspot_col), c("cell_line_id", "hotspot_status"))
        
        # Deduplicate hotspot data - keep first occurrence of each cell line
        hs_dt <- hs_dt[!duplicated(hs_dt$cell_line_id), ]
        
        # Use match-based lookup instead of join to prevent row duplication
        match_idx <- match(plot_data$cell_line_id, hs_dt$cell_line_id)
        plot_data$hotspot_status <- hs_dt$hotspot_status[match_idx]
        plot_data$hotspot_status <- suppressWarnings(as.integer(plot_data$hotspot_status))

        if (hotspot_mode == "0") {
          plot_data$hotspot_highlight <- !is.na(plot_data$hotspot_status) & plot_data$hotspot_status == 0
        } else if (hotspot_mode == "1") {
          plot_data$hotspot_highlight <- !is.na(plot_data$hotspot_status) & plot_data$hotspot_status == 1
        } else if (hotspot_mode == "2") {
          plot_data$hotspot_highlight <- !is.na(plot_data$hotspot_status) & plot_data$hotspot_status == 2
        } else if (hotspot_mode %in% c("12", "compare", "color_by_mut", "three_panel", "compare_table")) {
          plot_data$hotspot_highlight <- !is.na(plot_data$hotspot_status) & plot_data$hotspot_status %in% c(1, 2)
        }
      }
    }
    gene1 <- rv$current_scatter_genes[1]
    gene2 <- rv$current_scatter_genes[2]
    
    # Get correlation values
    cor_all <- rv$current_cor_all
    cor_filtered <- rv$current_cor_filtered
    filter_active <- rv$current_filter_active
    filter_label <- rv$current_filter_label
    
    # Get axis limits from inputs
    xlim <- c(input$scatter_xmin, input$scatter_xmax)
    ylim <- c(input$scatter_ymin, input$scatter_ymax)
    
    # Get cell line search terms
    search_terms <- if (!is.null(input$scatter_cell_search) && input$scatter_cell_search != "") {
      toupper(trimws(strsplit(input$scatter_cell_search, "\n")[[1]]))
    } else {
      character(0)
    }
    
    # Mark highlighted cell lines - search cell_line names, ACH IDs, AND original names
    plot_data$highlight <- FALSE
    plot_data$search_label <- NA_character_  # Store label for search-matched cells
    if (length(search_terms) > 0) {
      search_terms <- search_terms[search_terms != ""]  # Remove empty strings
      for (term in search_terms) {
        # Search in displayed cell line names (StrippedCellLineName)
        matches_name <- grepl(term, toupper(plot_data$cell_line), fixed = FALSE)
        # Search in ACH IDs
        matches_id <- if (!is.null(plot_data$cell_line_id)) {
          grepl(term, toupper(plot_data$cell_line_id), fixed = FALSE)
        } else {
          FALSE
        }
        # Search in original CellLineName (with hyphens etc.)
        matches_original <- if (!is.null(plot_data$cell_line_original)) {
          grepl(term, toupper(plot_data$cell_line_original), fixed = FALSE)
        } else {
          FALSE
        }
        plot_data$highlight <- plot_data$highlight | matches_name | matches_id | matches_original
      }
      
      # Create labels with cancer type for search-matched cells
      for (i in which(plot_data$highlight)) {
        cell_name <- plot_data$cell_line[i]
        cancer_type <- if ("OncotreeLineage" %in% names(plot_data)) plot_data$OncotreeLineage[i] else ""
        if (!is.null(cancer_type) && !is.na(cancer_type) && cancer_type != "") {
          plot_data$search_label[i] <- paste0(cell_name, " (", cancer_type, ")")
        } else {
          plot_data$search_label[i] <- cell_name
        }
      }
    }
    
    # Initialize click_label column (needed even if no cells clicked, for search highlights)
    plot_data$click_label <- NA_character_
    
    # Also highlight clicked cells (use names which are the cell_line values)
    if (length(rv$clicked_cells) > 0) {
      clicked_names <- names(rv$clicked_cells)
      plot_data$highlight <- plot_data$highlight | (plot_data$cell_line %in% clicked_names)
      # Store the full labels for display
      for (i in seq_along(clicked_names)) {
        plot_data$click_label[plot_data$cell_line == clicked_names[i]] <- rv$clicked_cells[i]
      }
    }
    
    # Build subtitle with correlation info
    if (filter_active && !is.null(cor_filtered)) {
      n_filtered <- sum(plot_data$in_filter, na.rm = TRUE)
        # Slope values for any displayed r
    slope_all <- if (nrow(plot_data) >= 3) round(coef(lm(gene2_effect ~ gene1_effect, data = plot_data))[2], 3) else NA
    
    subtitle_text <- paste0("All cells: r=", cor_all, ", slope=", slope_all, " (n=", nrow(plot_data), ")",
                             " | ", filter_label, ": r=", cor_filtered, " (n=", n_filtered, ", shown in orange)")
    } else {
      slope_all <- if (nrow(plot_data) >= 3) round(coef(lm(gene2_effect ~ gene1_effect, data = plot_data))[2], 3) else NA
    subtitle_text <- paste0("Correlation: ", cor_all, ", slope=", slope_all, " (n=", nrow(plot_data), " cell lines)")
    }
    
    if (sum(plot_data$highlight) > 0) {
      subtitle_text <- paste0(subtitle_text, " | ", sum(plot_data$highlight), " search matches")
    }

    # Append hotspot mutation info (if active)
    if (hotspot_mode == "three_panel" && any(!is.na(plot_data$hotspot_status))) {
      n_wt <- sum(is.na(plot_data$hotspot_status) | plot_data$hotspot_status == 0, na.rm = TRUE)
      n_mut1 <- sum(!is.na(plot_data$hotspot_status) & plot_data$hotspot_status == 1, na.rm = TRUE)
      n_mut2 <- sum(!is.na(plot_data$hotspot_status) & plot_data$hotspot_status == 2, na.rm = TRUE)
      subtitle_text <- paste0(subtitle_text, "\n", hotspot_gene, ": WT (n=", n_wt, 
                              "), 1 mut (n=", n_mut1, "), 2 mut (n=", n_mut2, ")")
    } else if (hotspot_mode == "color_by_mut" && any(!is.na(plot_data$hotspot_status))) {
      n_total <- nrow(plot_data)
      n_wt <- sum(is.na(plot_data$hotspot_status) | plot_data$hotspot_status == 0, na.rm = TRUE)
      n_mut1 <- sum(!is.na(plot_data$hotspot_status) & plot_data$hotspot_status == 1, na.rm = TRUE)
      n_mut2 <- sum(!is.na(plot_data$hotspot_status) & plot_data$hotspot_status == 2, na.rm = TRUE)
      
      # Calculate percentages
      pct_wt <- round(100 * n_wt / n_total, 1)
      pct_mut1 <- round(100 * n_mut1 / n_total, 1)
      pct_mut2 <- round(100 * n_mut2 / n_total, 1)
      
      # Calculate r and slope for each group
      wt_data <- plot_data[is.na(plot_data$hotspot_status) | plot_data$hotspot_status == 0, ]
      mut1_data <- plot_data[!is.na(plot_data$hotspot_status) & plot_data$hotspot_status == 1, ]
      mut2_data <- plot_data[!is.na(plot_data$hotspot_status) & plot_data$hotspot_status == 2, ]
      
      r_wt <- if (nrow(wt_data) >= 3) round(cor(wt_data$gene1_effect, wt_data$gene2_effect, use = "complete.obs"), 3) else NA
      r_mut1 <- if (nrow(mut1_data) >= 3) round(cor(mut1_data$gene1_effect, mut1_data$gene2_effect, use = "complete.obs"), 3) else NA
      r_mut2 <- if (nrow(mut2_data) >= 3) round(cor(mut2_data$gene1_effect, mut2_data$gene2_effect, use = "complete.obs"), 3) else NA
      
      slope_wt <- if (nrow(wt_data) >= 3) round(coef(lm(gene2_effect ~ gene1_effect, data = wt_data))[2], 3) else NA
      slope_mut1 <- if (nrow(mut1_data) >= 3) round(coef(lm(gene2_effect ~ gene1_effect, data = mut1_data))[2], 3) else NA
      slope_mut2 <- if (nrow(mut2_data) >= 3) round(coef(lm(gene2_effect ~ gene1_effect, data = mut2_data))[2], 3) else NA
      
      # Build detailed subtitle with line breaks for each group
      subtitle_text <- paste0(subtitle_text, "\n", hotspot_gene, " mutations:")
      subtitle_text <- paste0(subtitle_text, "\n  WT (gray): n=", n_wt, " (", pct_wt, "%), r=", r_wt, ", slope=", slope_wt)
      subtitle_text <- paste0(subtitle_text, "\n  1 mut (blue): n=", n_mut1, " (", pct_mut1, "%), r=", r_mut1, ", slope=", slope_mut1)
      subtitle_text <- paste0(subtitle_text, "\n  2 mut (red): n=", n_mut2, " (", pct_mut2, "%), r=", r_mut2, ", slope=", slope_mut2)
    } else if (sum(plot_data$hotspot_highlight, na.rm = TRUE) >= 3) {
      hs_df <- plot_data[plot_data$hotspot_highlight, ]
      hs_r <- suppressWarnings(round(cor(hs_df$gene1_effect, hs_df$gene2_effect, use = "complete.obs"), 3))
      hs_slope <- suppressWarnings(round(coef(lm(gene2_effect ~ gene1_effect, data = hs_df))[2], 3))
      subtitle_text <- paste0(subtitle_text, "\nHotspot (", hotspot_gene, "): r=", hs_r, ", slope=", hs_slope, " (n=", nrow(hs_df), ")")
    }

    if (hotspot_mode == "compare") {
      g0 <- plot_data[!is.na(plot_data$hotspot_status) & plot_data$hotspot_status == 0, ]
      g12 <- plot_data[!is.na(plot_data$hotspot_status) & plot_data$hotspot_status %in% c(1,2), ]
      if (nrow(g0) >= 3 && nrow(g12) >= 3) {
        r0 <- suppressWarnings(cor(g0$gene1_effect, g0$gene2_effect, use = "complete.obs"))
        r12 <- suppressWarnings(cor(g12$gene1_effect, g12$gene2_effect, use = "complete.obs"))
        s0 <- suppressWarnings(coef(lm(gene2_effect ~ gene1_effect, data = g0))[2])
        s12 <- suppressWarnings(coef(lm(gene2_effect ~ gene1_effect, data = g12))[2])

        # Fisher r-to-z difference test
        z <- (atanh(r0) - atanh(r12)) / sqrt(1/(nrow(g0)-3) + 1/(nrow(g12)-3))
        p_diff <- 2 * pnorm(-abs(z))

        subtitle_text <- paste0(subtitle_text,
                                "\nCompare 0 vs 1+2: r0=", round(r0,3), ", slope0=", round(s0,3), " (n=", nrow(g0), 
                                "); r1+2=", round(r12,3), ", slope1+2=", round(s12,3), " (n=", nrow(g12),
                                "); p(diff)=", format(p_diff, scientific = TRUE, digits = 2))
      }
    }
    
    # Determine point colors based on filter status
    if (hotspot_mode == "three_panel" && any(!is.na(plot_data$hotspot_status))) {
      # 3-panel mode: create faceted plot by mutation status with pie charts below
      plot_data$mut_panel <- factor(
        ifelse(is.na(plot_data$hotspot_status), "WT (0 mut)", 
               ifelse(plot_data$hotspot_status == 0, "WT (0 mut)",
                      ifelse(plot_data$hotspot_status == 1, "1 mutation", "2 mutations"))),
        levels = c("WT (0 mut)", "1 mutation", "2 mutations")
      )
      
      # Calculate stats for each panel (including medians)
      panel_stats <- plot_data %>%
        dplyr::group_by(mut_panel) %>%
        dplyr::summarise(
          n = dplyr::n(),
          r = round(cor(gene1_effect, gene2_effect, use = "complete.obs"), 3),
          slope = round(coef(lm(gene2_effect ~ gene1_effect))[2], 3),
          median_x = round(median(gene1_effect, na.rm = TRUE), 2),
          median_y = round(median(gene2_effect, na.rm = TRUE), 2),
          .groups = "drop"
        )
      
      # Check if a cancer type is highlighted (for both scatter and pie)
      highlight_cancer <- if (!is.null(input$pie_highlight_cancer) && nzchar(input$pie_highlight_cancer)) {
        input$pie_highlight_cancer
      } else {
        NULL
      }
      
      # If a cancer type is highlighted, calculate stats for that cancer type too
      if (!is.null(highlight_cancer) && "OncotreeLineage" %in% names(plot_data)) {
        cancer_data <- plot_data[!is.na(plot_data$OncotreeLineage) & plot_data$OncotreeLineage == highlight_cancer, ]
        
        if (nrow(cancer_data) > 0) {
          cancer_stats <- cancer_data %>%
            dplyr::group_by(mut_panel) %>%
            dplyr::summarise(
              n_cancer = dplyr::n(),
              r_cancer = if (dplyr::n() >= 3) round(cor(gene1_effect, gene2_effect, use = "complete.obs"), 3) else NA,
              slope_cancer = if (dplyr::n() >= 3) round(coef(lm(gene2_effect ~ gene1_effect))[2], 3) else NA,
              .groups = "drop"
            )
          
          # Merge with panel_stats
          panel_stats <- panel_stats %>%
            dplyr::left_join(cancer_stats, by = "mut_panel")
          
          # Fill NAs with 0 for n_cancer
          panel_stats$n_cancer[is.na(panel_stats$n_cancer)] <- 0
          
          # Create labels showing both All and highlighted cancer stats
          # Regression line is always for ALL cells
          panel_labels <- setNames(
            paste0(panel_stats$mut_panel, 
                   "\nAll: n=", panel_stats$n, ", r=", panel_stats$r, ", slope=", panel_stats$slope,
                   "\n", highlight_cancer, ": n=", panel_stats$n_cancer, 
                   ", r=", ifelse(is.na(panel_stats$r_cancer), "NA", panel_stats$r_cancer),
                   ", slope=", ifelse(is.na(panel_stats$slope_cancer), "NA", panel_stats$slope_cancer),
                   "\n(Regression line: All cells)"),
            panel_stats$mut_panel
          )
        } else {
          # No cells for this cancer type
          panel_labels <- setNames(
            paste0(panel_stats$mut_panel, "\n(n=", panel_stats$n, ", r=", panel_stats$r, ", slope=", panel_stats$slope, 
                   ")\nMedian ", gene1, "=", panel_stats$median_x, ", ", gene2, "=", panel_stats$median_y),
            panel_stats$mut_panel
          )
        }
      } else {
        # Create labels for facets (include slope and medians) - no cancer highlighting
        panel_labels <- setNames(
          paste0(panel_stats$mut_panel, "\n(n=", panel_stats$n, ", r=", panel_stats$r, ", slope=", panel_stats$slope, 
                 ")\nMedian ", gene1, "=", panel_stats$median_x, ", ", gene2, "=", panel_stats$median_y),
          panel_stats$mut_panel
        )
      }
      
      # Include highlighted cells in the plot data for 3-panel
      plot_data_3p <- plot_data
      
      # Add cancer type highlighting flag
      plot_data_3p$cancer_highlight <- FALSE
      if (!is.null(highlight_cancer) && "OncotreeLineage" %in% names(plot_data_3p)) {
        plot_data_3p$cancer_highlight <- !is.na(plot_data_3p$OncotreeLineage) & 
                                          plot_data_3p$OncotreeLineage == highlight_cancer
      }
      
      # We need to compute the cancer colors consistently with the pie chart
      # First, get the filtered cancer types (>= 2% representation) like the pie chart does
      cancer_counts_for_colors <- NULL
      scatter_cancer_colors <- NULL
      highlight_cancer_color <- "#f97316"  # Default fallback
      
      if ("OncotreeLineage" %in% names(plot_data) && !all(is.na(plot_data$OncotreeLineage))) {
        cancer_counts_for_colors <- plot_data %>%
          dplyr::filter(!is.na(OncotreeLineage) & OncotreeLineage != "") %>%
          dplyr::group_by(mut_panel, OncotreeLineage) %>%
          dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
          dplyr::group_by(mut_panel) %>%
          dplyr::mutate(total = sum(n), pct = n / total * 100) %>%
          dplyr::ungroup()
        
        # Keep highlighted cancer even if < 2%
        if (!is.null(highlight_cancer)) {
          cancer_counts_for_colors <- cancer_counts_for_colors %>%
            dplyr::filter(pct >= 2 | OncotreeLineage == highlight_cancer)
        } else {
          cancer_counts_for_colors <- cancer_counts_for_colors %>%
            dplyr::filter(pct >= 2)
        }
        
        # Use the same sorted list as pie chart
        all_cancer_types_sorted <- sort(unique(cancer_counts_for_colors$OncotreeLineage))
        n_cancer_colors <- length(all_cancer_types_sorted)
        
        color_palette_3p <- c("#e11d48", "#f97316", "#eab308", "#22c55e", "#14b8a6", 
                              "#0ea5e9", "#6366f1", "#a855f7", "#ec4899", "#78716c",
                              "#dc2626", "#84cc16", "#06b6d4", "#8b5cf6", "#f43f5e",
                              "#64748b", "#059669", "#7c3aed", "#db2777", "#ca8a04",
                              "#991b1b", "#166534", "#1e40af", "#6b21a8", "#be185d",
                              "#854d0e", "#0f766e", "#4338ca", "#9d174d", "#365314")
        if (n_cancer_colors > length(color_palette_3p)) {
          extra_cols <- colorRampPalette(c("#e11d48", "#0ea5e9", "#22c55e", "#a855f7"))(n_cancer_colors - length(color_palette_3p))
          color_palette_3p <- c(color_palette_3p, extra_cols)
        }
        scatter_cancer_colors <- setNames(color_palette_3p[seq_len(n_cancer_colors)], all_cancer_types_sorted)
        
        # Get the highlight color for the selected cancer type
        if (!is.null(highlight_cancer) && highlight_cancer %in% names(scatter_cancer_colors)) {
          highlight_cancer_color <- scatter_cancer_colors[highlight_cancer]
        }
      }
      
      # Main scatter plot with 3 panels
      if (!is.null(highlight_cancer) && sum(plot_data_3p$cancer_highlight) > 0) {
        # When a cancer type is highlighted: gray out non-highlighted cells, color highlighted with cancer color
        non_highlighted <- plot_data_3p[!plot_data_3p$cancer_highlight, ]
        highlighted <- plot_data_3p[plot_data_3p$cancer_highlight, ]
        
        p_scatter <- ggplot(plot_data_3p, aes(x = gene1_effect, y = gene2_effect)) +
          geom_smooth(method = "lm", se = TRUE, color = "#059669", fill = "#86efac", linewidth = 0.8) +
          # Gray points for non-highlighted cells
          geom_point(data = non_highlighted, color = "#d1d5db", alpha = 0.4, size = 1.5) +
          # Colored points for highlighted cancer type (use the cancer type's color from palette)
          geom_point(data = highlighted, color = highlight_cancer_color, alpha = 0.8, size = 2.5) +
          facet_wrap(~mut_panel, labeller = labeller(mut_panel = panel_labels), nrow = 1) +
          coord_cartesian(xlim = sanitize_limits(xlim), ylim = sanitize_limits(ylim)) +
          labs(
            title = paste0(gene1, " vs ", gene2, " - ", hotspot_gene, " hotspot mutation stratification"),
            subtitle = paste0("Highlighting: ", highlight_cancer),
            x = paste0(gene1, " CRISPR Effect"),
            y = paste0(gene2, " CRISPR Effect")
          ) +
          theme_minimal(base_size = 11) +
          theme(
            legend.position = "none",
            plot.title = element_text(face = "bold", size = 13),
            plot.subtitle = element_text(color = highlight_cancer_color, size = 10),
            strip.text = element_text(face = "bold", size = 9),
            panel.spacing = unit(0.5, "lines"),
            panel.border = element_rect(color = "#374151", fill = NA, linewidth = 1),
            aspect.ratio = 1
          )
      } else {
        # No cancer highlight - use mutation status colors
        p_scatter <- ggplot(plot_data_3p, aes(x = gene1_effect, y = gene2_effect)) +
          geom_smooth(method = "lm", se = TRUE, color = "#059669", fill = "#86efac", linewidth = 0.8) +
          geom_point(aes(color = mut_panel), alpha = 0.6, size = 1.8) +
          scale_color_manual(values = c("WT (0 mut)" = "#9ca3af", "1 mutation" = "#3b82f6", "2 mutations" = "#dc2626")) +
          facet_wrap(~mut_panel, labeller = labeller(mut_panel = panel_labels), nrow = 1) +
          coord_cartesian(xlim = sanitize_limits(xlim), ylim = sanitize_limits(ylim)) +
          labs(
            title = paste0(gene1, " vs ", gene2, " - ", hotspot_gene, " hotspot mutation stratification"),
            x = paste0(gene1, " CRISPR Effect"),
            y = paste0(gene2, " CRISPR Effect")
          ) +
          theme_minimal(base_size = 11) +
          theme(
            legend.position = "none",
            plot.title = element_text(face = "bold", size = 13),
            strip.text = element_text(face = "bold", size = 9),
            panel.spacing = unit(0.5, "lines"),
            panel.border = element_rect(color = "#374151", fill = NA, linewidth = 1),
            aspect.ratio = 1
          )
      }
      
      # Add highlighted points if any (from text search or click)
      if (sum(plot_data_3p$highlight) > 0) {
        # Use highlight_label_size for both text search and clicked points
        highlight_font_size <- if (!is.null(input$highlight_label_size)) input$highlight_label_size else 3
        highlight_3p <- plot_data_3p[plot_data_3p$highlight, ]
        highlight_3p$display_label <- ifelse(
          !is.na(highlight_3p$click_label),
          highlight_3p$click_label,
          ifelse(!is.na(highlight_3p$search_label),
                 highlight_3p$search_label,
                 highlight_3p$cell_line)
        )
        
        p_scatter <- p_scatter + 
          geom_point(data = highlight_3p, color = "#dc2626", size = 3, alpha = 0.9) +
          ggrepel::geom_text_repel(
            data = highlight_3p,
            aes(label = display_label),
            size = highlight_font_size,
            color = "#dc2626",
            max.overlaps = 20,
            seed = 42
          )
      }
      
      # Create pie charts for cancer type distribution in each panel
      if ("OncotreeLineage" %in% names(plot_data) && !all(is.na(plot_data$OncotreeLineage))) {
        
        # highlight_cancer is already defined above for scatter plot highlighting
        
        # Get ALL cancer counts first (before filtering by percentage)
        cancer_counts_all <- plot_data %>%
          dplyr::filter(!is.na(OncotreeLineage) & OncotreeLineage != "") %>%
          dplyr::group_by(mut_panel, OncotreeLineage) %>%
          dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
          dplyr::group_by(mut_panel) %>%
          dplyr::mutate(
            total = sum(n),
            pct = n / total * 100
          ) %>%
          dplyr::ungroup()
        
        # Filter to >= 2% representation, BUT always keep the highlighted cancer type
        if (!is.null(highlight_cancer)) {
          cancer_counts <- cancer_counts_all %>%
            dplyr::filter(pct >= 2 | OncotreeLineage == highlight_cancer)
        } else {
          cancer_counts <- cancer_counts_all %>%
            dplyr::filter(pct >= 2)
        }
        
        # Calculate total N for each panel (for center label)
        panel_totals <- plot_data %>%
          dplyr::filter(!is.na(OncotreeLineage) & OncotreeLineage != "") %>%
          dplyr::group_by(mut_panel) %>%
          dplyr::summarise(total_n = dplyr::n(), .groups = "drop")
        
        # If a cancer type is highlighted, update center labels to show that cancer's count/percent
        if (!is.null(highlight_cancer)) {
          # Get the highlighted cancer's stats per panel from the unfiltered data
          highlight_stats <- cancer_counts_all %>%
            dplyr::filter(OncotreeLineage == highlight_cancer)
          
          # Merge with panel_totals to get stats for each panel
          panel_totals <- panel_totals %>%
            dplyr::left_join(highlight_stats %>% dplyr::select(mut_panel, n, pct), by = "mut_panel")
          
          # Check pie scaling mode for label format
          pie_mode <- if (!is.null(input$pie_scale_mode)) input$pie_scale_mode else "count"
          
          # Create center label based on mode
          if (pie_mode == "percent") {
            panel_totals$center_label <- ifelse(
              is.na(panel_totals$n),
              "0%",
              paste0(round(panel_totals$pct, 1), "%")
            )
          } else {
            panel_totals$center_label <- ifelse(
              is.na(panel_totals$n),
              "n=0",
              paste0("n=", panel_totals$n)
            )
          }
        } else {
          # No cancer highlighted - show total n in center
          panel_totals$center_label <- paste0("n=", panel_totals$total_n)
        }
        
        # Get ALL unique cancer types across all panels for consistent legend
        all_cancer_types <- sort(unique(cancer_counts$OncotreeLineage))
        n_colors <- length(all_cancer_types)
        
        # Create consistent color palette (expanded to 30 colors)
        color_palette <- c("#e11d48", "#f97316", "#eab308", "#22c55e", "#14b8a6", 
                          "#0ea5e9", "#6366f1", "#a855f7", "#ec4899", "#78716c",
                          "#dc2626", "#84cc16", "#06b6d4", "#8b5cf6", "#f43f5e",
                          "#64748b", "#059669", "#7c3aed", "#db2777", "#ca8a04",
                          "#991b1b", "#166534", "#1e40af", "#6b21a8", "#be185d",
                          "#854d0e", "#0f766e", "#4338ca", "#9d174d", "#365314")
        # If we need more colors, generate them
        if (n_colors > length(color_palette)) {
          extra_colors <- colorRampPalette(c("#e11d48", "#0ea5e9", "#22c55e", "#a855f7"))(n_colors - length(color_palette))
          color_palette <- c(color_palette, extra_colors)
        }
        cancer_colors <- setNames(color_palette[1:n_colors], all_cancer_types)
        
        # Calculate total N per cancer type for legend
        cancer_totals <- cancer_counts %>%
          dplyr::group_by(OncotreeLineage) %>%
          dplyr::summarise(total = sum(n), .groups = "drop") %>%
          dplyr::arrange(desc(total))
        
        # Check pie scaling mode (count = proportional, percent = equal size)
        pie_mode <- if (!is.null(input$pie_scale_mode)) input$pie_scale_mode else "count"
        
        # For equal-size mode, convert counts to percentages within each panel
        if (pie_mode == "percent") {
          cancer_counts <- cancer_counts %>%
            dplyr::group_by(mut_panel) %>%
            dplyr::mutate(
              n_display = n / sum(n) * 100  # Convert to percentage
            ) %>%
            dplyr::ungroup()
        } else {
          cancer_counts$n_display <- cancer_counts$n
        }
        
        # Create legend labels with N (always show counts in legend)
        legend_labels <- setNames(
          paste0(cancer_totals$OncotreeLineage, " (n=", cancer_totals$total, ")"),
          cancer_totals$OncotreeLineage
        )
        
        # Modify colors if highlighting a specific cancer type
        if (!is.null(highlight_cancer) && highlight_cancer %in% all_cancer_types) {
          # Store the original color for the highlighted cancer type
          highlight_color <- cancer_colors[highlight_cancer]
          # Gray out all colors
          cancer_colors[] <- "#e5e7eb"
          # Restore the highlighted cancer type's original color
          cancer_colors[highlight_cancer] <- highlight_color
        }
        
        p_pie <- ggplot(cancer_counts, aes(x = "", y = n_display, fill = OncotreeLineage)) +
          geom_bar(stat = "identity", width = 1) +
          coord_polar("y", start = 0) +
          facet_wrap(~mut_panel, nrow = 1) +
          scale_fill_manual(values = cancer_colors, labels = legend_labels, name = NULL) +
          geom_text(data = panel_totals, aes(x = 0, y = 0, label = center_label), 
                    inherit.aes = FALSE, size = 3.5, fontface = "bold") +
          theme_void() +
          theme(
            strip.text = element_blank(),
            legend.position = "bottom",
            legend.text = element_text(size = 7),
            legend.key.size = unit(0.3, "cm"),
            legend.spacing.x = unit(0.1, "cm"),
            legend.spacing.y = unit(0.05, "cm"),
            legend.margin = margin(t = 0, b = 0),
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA)
          ) +
          guides(fill = guide_legend(nrow = 3, byrow = TRUE))
        
        # Combine scatter and pie plots - reduce scatter height for more square panels
        p <- cowplot::plot_grid(p_scatter, p_pie, ncol = 1, rel_heights = c(1.8, 1.2))
      } else {
        p <- p_scatter
      }
      
      # Return early since we've already built the complete plot
      return(p)
        
    } else if (hotspot_mode == "color_by_mut" && any(!is.na(plot_data$hotspot_status))) {
      # Color by mutation mode: WT=gray/small, 1 mut=blue/large, 2 mut=red/largest
      # Separate data into layers
      wt_data <- plot_data[is.na(plot_data$hotspot_status) | plot_data$hotspot_status == 0, ]
      mut1_data <- plot_data[!is.na(plot_data$hotspot_status) & plot_data$hotspot_status == 1, ]
      mut2_data <- plot_data[!is.na(plot_data$hotspot_status) & plot_data$hotspot_status == 2, ]
      
      # Remove highlighted from these to handle separately
      wt_data_no_hl <- wt_data[!wt_data$highlight, ]
      mut1_data_no_hl <- mut1_data[!mut1_data$highlight, ]
      mut2_data_no_hl <- mut2_data[!mut2_data$highlight, ]
      
      # Create legend labels with N and %
      n_wt <- nrow(wt_data)
      n_mut1 <- nrow(mut1_data)
      n_mut2 <- nrow(mut2_data)
      n_total <- nrow(plot_data)
      pct_wt <- round(100 * n_wt / n_total, 1)
      pct_mut1 <- round(100 * n_mut1 / n_total, 1)
      pct_mut2 <- round(100 * n_mut2 / n_total, 1)
      
      # Add a mutation_label column for legend (only if data exists)
      if (nrow(wt_data_no_hl) > 0) wt_data_no_hl$mutation_label <- paste0("WT (n=", n_wt, ", ", pct_wt, "%)")
      if (nrow(mut1_data_no_hl) > 0) mut1_data_no_hl$mutation_label <- paste0("1 mut (n=", n_mut1, ", ", pct_mut1, "%)")
      if (nrow(mut2_data_no_hl) > 0) mut2_data_no_hl$mutation_label <- paste0("2 mut (n=", n_mut2, ", ", pct_mut2, "%)")
      
      # Combine for legend - only include non-empty dataframes
      dfs_to_bind <- list()
      if (nrow(wt_data_no_hl) > 0) dfs_to_bind <- c(dfs_to_bind, list(wt_data_no_hl))
      if (nrow(mut1_data_no_hl) > 0) dfs_to_bind <- c(dfs_to_bind, list(mut1_data_no_hl))
      if (nrow(mut2_data_no_hl) > 0) dfs_to_bind <- c(dfs_to_bind, list(mut2_data_no_hl))
      
      if (length(dfs_to_bind) > 0) {
        legend_data <- do.call(rbind, dfs_to_bind)
        legend_data$mutation_label <- factor(legend_data$mutation_label, 
          levels = c(paste0("WT (n=", n_wt, ", ", pct_wt, "%)"), 
                     paste0("1 mut (n=", n_mut1, ", ", pct_mut1, "%)"), 
                     paste0("2 mut (n=", n_mut2, ", ", pct_mut2, "%)")))
      } else {
        # No data to plot - create empty legend_data
        legend_data <- plot_data[0, ]
        legend_data$mutation_label <- factor(character(0))
      }
      p <- ggplot(plot_data, aes(x = gene1_effect, y = gene2_effect)) +
        geom_smooth(method = "lm", se = TRUE, color = "#059669", fill = "#86efac") +
        # Plot with legend
        geom_point(data = legend_data, aes(color = mutation_label, size = mutation_label, alpha = mutation_label)) +
        scale_color_manual(
          name = hotspot_gene,
          values = c("#9ca3af", "#3b82f6", "#dc2626")
        ) +
        scale_size_manual(
          name = hotspot_gene,
          values = c(2, 3, 3.5)
        ) +
        scale_alpha_manual(
          name = hotspot_gene,
          values = c(0.3, 0.5, 0.6)
        ) +
        guides(
          color = guide_legend(title = hotspot_gene, override.aes = list(size = 4)),
          size = "none",
          alpha = "none"
        ) +
        coord_cartesian(xlim = sanitize_limits(xlim), ylim = sanitize_limits(ylim)) +
        theme(
          legend.position = "right",
          legend.title = element_text(face = "bold", size = 11),
          legend.text = element_text(size = 10)
        )
        
    } else if (filter_active && "in_filter" %in% colnames(plot_data)) {
      # Non-highlighted points: gray for non-filtered, orange for filtered
      non_highlight_data <- plot_data[!plot_data$highlight, ]
      
      p <- ggplot(plot_data, aes(x = gene1_effect, y = gene2_effect)) +
        geom_smooth(method = "lm", se = TRUE, color = "#059669", fill = "#86efac") +
        # Non-filtered cells (gray)
        geom_point(data = non_highlight_data[!non_highlight_data$in_filter, ], 
                   alpha = 0.3, size = 2, color = "#9ca3af") +
        # Filtered cells (orange) - the selected cancer type
        geom_point(data = non_highlight_data[non_highlight_data$in_filter, ], 
                   alpha = 0.7, size = 2.5, color = "#f97316") +
        coord_cartesian(xlim = sanitize_limits(xlim), ylim = sanitize_limits(ylim))
    } else {
      # No filter - all points green
      p <- ggplot(plot_data, aes(x = gene1_effect, y = gene2_effect)) +
        geom_smooth(method = "lm", se = TRUE, color = "#059669", fill = "#86efac") +
        geom_point(data = plot_data[!plot_data$highlight, ], 
                   alpha = 0.5, size = 2, color = "#16a34a") +
        coord_cartesian(xlim = sanitize_limits(xlim), ylim = sanitize_limits(ylim))
    }
    
    p <- p + labs(
        title = paste0(gene1, " vs ", gene2),
        subtitle = subtitle_text,
        x = paste0(gene1, " CRISPR Effect"),
        y = paste0(gene2, " CRISPR Effect")
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "gray30", size = 11),
        aspect.ratio = if (hotspot_mode == "three_panel") NULL else 1  # Allow flexible ratio for facets
      )
    
    
    # Add hotspot overlay points (blue) below search highlights
    # Skip if using color_by_mut or three_panel mode (already rendered with color coding)
    if (!hotspot_mode %in% c("color_by_mut", "three_panel") && sum(plot_data$hotspot_highlight, na.rm = TRUE) > 0) {
      p <- p +
        geom_point(
          data = plot_data[plot_data$hotspot_highlight & !plot_data$highlight, ],
          alpha = 0.85, size = 3.2, color = "#2563eb"
        )
    }

# Add highlighted points on top (search results and clicked cells in red)
    if (sum(plot_data$highlight) > 0) {
      highlight_font_size <- if (!is.null(input$highlight_label_size)) input$highlight_label_size else 3
      
      # For clicked cells use click_label, for search-matched use search_label, otherwise just cell_line
      highlight_data <- plot_data[plot_data$highlight, ]
      highlight_data$display_label <- ifelse(
        !is.na(highlight_data$click_label),
        highlight_data$click_label,
        ifelse(!is.na(highlight_data$search_label),
               highlight_data$search_label,
               highlight_data$cell_line)
      )
      
      p <- p + 
        geom_point(data = highlight_data, 
                   alpha = 0.9, size = 4, color = "#dc2626") +
        ggrepel::geom_text_repel(
          data = highlight_data,
          aes(label = display_label),
          size = highlight_font_size, 
          color = "#dc2626",
          max.overlaps = 30,
          box.padding = 0.3,
          point.padding = 0.2,
          seed = 42
        )
    }
    
    p
  })
  
  # Mutation comparison table (compare WT vs Mutant by cancer type)
  output$mutation_compare_table <- renderDT({
    req(rv$current_plot_data)
    req(rv$current_scatter_genes)
    req(input$hotspot_gene)
    req(rv$hotspot_file_path)
    
    plot_data <- rv$current_plot_data
    gene1 <- rv$current_scatter_genes[1]
    gene2 <- rv$current_scatter_genes[2]
    hotspot_col <- trimws(input$hotspot_gene)
    hotspot_gene <- gsub("\\s*\\([^)]+\\)$", "", hotspot_col)
    
    # Load hotspot data
    hs_dt <- tryCatch({
      data.table::fread(rv$hotspot_file_path,
                        select = c(rv$hotspot_id_col, hotspot_col),
                        showProgress = FALSE)
    }, error = function(e) NULL)
    
    if (is.null(hs_dt)) {
      return(DT::datatable(data.frame(Message = "Could not load hotspot data."), options = list(dom = 't')))
    }
    
    data.table::setnames(hs_dt, c(rv$hotspot_id_col, hotspot_col), c("cell_line_id", "hotspot_status"))
    hs_dt <- hs_dt[!duplicated(hs_dt$cell_line_id), ]
    
    # Merge hotspot status
    match_idx <- match(plot_data$cell_line_id, hs_dt$cell_line_id)
    plot_data$hotspot_status <- suppressWarnings(as.integer(hs_dt$hotspot_status[match_idx]))
    
    # Need OncotreeLineage for stratification
    if (!"OncotreeLineage" %in% names(plot_data) || all(is.na(plot_data$OncotreeLineage))) {
      return(DT::datatable(data.frame(Message = "Upload Model.csv to enable cancer type stratification."), options = list(dom = 't')))
    }
    
    # Calculate stats by cancer type
    results <- plot_data %>%
      dplyr::filter(!is.na(hotspot_status) & !is.na(OncotreeLineage) & OncotreeLineage != "") %>%
      dplyr::mutate(mut_group = ifelse(hotspot_status == 0, "WT", "Mutant")) %>%
      dplyr::group_by(OncotreeLineage, mut_group) %>%
      dplyr::summarise(
        n = dplyr::n(),
        r = cor(gene1_effect, gene2_effect, use = "complete.obs"),
        slope = if(dplyr::n() >= 3) coef(lm(gene2_effect ~ gene1_effect))[2] else NA_real_,
        slope_se = if(dplyr::n() >= 3) summary(lm(gene2_effect ~ gene1_effect))$coefficients[2, 2] else NA_real_,
        .groups = "drop"
      ) %>%
      tidyr::pivot_wider(
        names_from = mut_group,
        values_from = c(n, r, slope, slope_se),
        names_sep = "_"
      ) %>%
      dplyr::filter(!is.na(n_WT) & !is.na(n_Mutant) & n_WT >= 3 & n_Mutant >= 3) %>%
      dplyr::mutate(
        r_diff = r_Mutant - r_WT,
        slope_diff = slope_Mutant - slope_WT,
        # Fisher z-test for correlation difference
        z_r = (atanh(r_Mutant) - atanh(r_WT)) / sqrt(1/(n_Mutant-3) + 1/(n_WT-3)),
        p_r = 2 * pnorm(-abs(z_r)),
        # t-test for slope difference (using pooled SE)
        se_diff = sqrt(slope_se_WT^2 + slope_se_Mutant^2),
        t_slope = slope_diff / se_diff,
        df_slope = n_WT + n_Mutant - 4,
        p_slope = 2 * pt(-abs(t_slope), df_slope)
      ) %>%
      dplyr::arrange(p_r) %>%
      dplyr::select(
        `Cancer Type` = OncotreeLineage,
        `N (WT)` = n_WT,
        `r (WT)` = r_WT,
        `slope (WT)` = slope_WT,
        `N (Mut)` = n_Mutant,
        `r (Mut)` = r_Mutant,
        `slope (Mut)` = slope_Mutant,
        `Δr` = r_diff,
        `p(Δr)` = p_r,
        `Δslope` = slope_diff,
        `p(Δslope)` = p_slope
      ) %>%
      dplyr::mutate(
        `r (WT)` = round(`r (WT)`, 3),
        `slope (WT)` = round(`slope (WT)`, 3),
        `r (Mut)` = round(`r (Mut)`, 3),
        `slope (Mut)` = round(`slope (Mut)`, 3),
        `Δr` = round(`Δr`, 3),
        `Δslope` = round(`Δslope`, 3)
        # Keep p-values as numeric for proper sorting
      )
    
    # Store for download (format p-values for CSV)
    results_for_download <- results %>%
      dplyr::mutate(
        `p(Δr)` = format(`p(Δr)`, scientific = TRUE, digits = 2),
        `p(Δslope)` = format(`p(Δslope)`, scientific = TRUE, digits = 2)
      )
    rv$mutation_compare_data <- results_for_download
    
    DT::datatable(
      results,
      rownames = FALSE,
      options = list(
        pageLength = 15,
        order = list(list(8, "asc")),  # Sort by p(Δr) ascending (column index 8, most significant first)
        dom = 'frtip',
        scrollX = TRUE,
        columnDefs = list(
          list(targets = c(8, 10), render = DT::JS(
            "function(data, type, row) {",
            "  if (type === 'display' && data !== null) {",
            "    return parseFloat(data).toExponential(1);",
            "  }",
            "  return data;",
            "}"
          ))
        )
      )
    ) %>%
      DT::formatStyle('Δr', 
        background = DT::styleInterval(c(-0.2, 0, 0.2), c('#fee2e2', '#fef3c7', '#fef3c7', '#dcfce7')),
        fontWeight = 'bold'
      ) %>%
      DT::formatStyle('Δslope', 
        background = DT::styleInterval(c(-0.3, 0, 0.3), c('#fee2e2', '#fef3c7', '#fef3c7', '#dcfce7')),
        fontWeight = 'bold'
      )
  })
  
  # Download mutation comparison table
  output$download_mutation_compare <- downloadHandler(
    filename = function() {
      gene1 <- rv$current_scatter_genes[1]
      gene2 <- rv$current_scatter_genes[2]
      hotspot_gene <- gsub("\\s*\\([^)]+\\)$", "", trimws(input$hotspot_gene))
      paste0(gene1, "_vs_", gene2, "_", hotspot_gene, "_mutation_comparison.csv")
    },
    content = function(file) {
      write.csv(rv$mutation_compare_data, file, row.names = FALSE)
    }
  )
  
  # Store current plot for download - mirrors the renderPlot logic
  current_scatter_plot <- reactive({
    req(rv$current_plot_data)
    req(rv$current_scatter_genes)
    
    plot_data <- rv$current_plot_data
    gene1 <- rv$current_scatter_genes[1]
    gene2 <- rv$current_scatter_genes[2]
    cor_all <- rv$current_cor_all
    
    xlim <- c(input$scatter_xmin %||% rv$scatter_xlim[1], 
              input$scatter_xmax %||% rv$scatter_xlim[2])
    ylim <- c(input$scatter_ymin %||% rv$scatter_ylim[1], 
              input$scatter_ymax %||% rv$scatter_ylim[2])
    
    # Apply cancer type filter
    if (!is.null(input$scatter_cancer_filter) && nzchar(input$scatter_cancer_filter) && 
        "OncotreeLineage" %in% colnames(plot_data)) {
      plot_data <- plot_data[plot_data$OncotreeLineage == input$scatter_cancer_filter, ]
    }
    
    # Handle hotspot mutation overlay
    plot_data$hotspot_status <- NA_integer_
    plot_data$hotspot_highlight <- FALSE
    hotspot_gene_input <- if (!is.null(input$hotspot_gene)) trimws(input$hotspot_gene) else ""
    hotspot_mode <- if (!is.null(input$hotspot_mode)) input$hotspot_mode else "none"
    hotspot_col <- hotspot_gene_input
    hotspot_gene <- gsub("\\s*\\([^)]+\\)$", "", hotspot_gene_input)
    
    if (!is.null(rv$hotspot_file_path) && nzchar(hotspot_col) && hotspot_mode != "none" &&
        !is.null(rv$hotspot_id_col) && nzchar(rv$hotspot_id_col)) {
      hs_dt <- tryCatch({
        data.table::fread(rv$hotspot_file_path,
                          select = c(rv$hotspot_id_col, hotspot_col),
                          showProgress = FALSE)
      }, error = function(e) NULL)
      
      if (!is.null(hs_dt) && hotspot_col %in% names(hs_dt)) {
        data.table::setnames(hs_dt, c(rv$hotspot_id_col, hotspot_col), c("cell_line_id", "hotspot_status"))
        hs_dt <- hs_dt[!duplicated(hs_dt$cell_line_id), ]
        match_idx <- match(plot_data$cell_line_id, hs_dt$cell_line_id)
        plot_data$hotspot_status <- hs_dt$hotspot_status[match_idx]
        plot_data$hotspot_status <- suppressWarnings(as.integer(plot_data$hotspot_status))
        
        if (hotspot_mode %in% c("12", "compare", "color_by_mut", "three_panel")) {
          plot_data$hotspot_highlight <- !is.na(plot_data$hotspot_status) & plot_data$hotspot_status %in% c(1, 2)
        }
      }
    }
    
    # Handle search/highlight
    search_terms <- if (!is.null(input$scatter_cell_search) && input$scatter_cell_search != "") {
      toupper(trimws(strsplit(input$scatter_cell_search, "\n")[[1]]))
    } else {
      character(0)
    }
    
    plot_data$highlight <- FALSE
    if (length(search_terms) > 0) {
      search_terms <- search_terms[search_terms != ""]
      for (term in search_terms) {
        matches_name <- grepl(term, toupper(plot_data$cell_line), fixed = FALSE)
        matches_id <- if (!is.null(plot_data$cell_line_id)) {
          grepl(term, toupper(plot_data$cell_line_id), fixed = FALSE)
        } else {
          FALSE
        }
        plot_data$highlight <- plot_data$highlight | matches_name | matches_id
      }
    }
    
    # Include clicked cells
    if (length(rv$clicked_cells) > 0) {
      clicked_names <- names(rv$clicked_cells)
      plot_data$highlight <- plot_data$highlight | (plot_data$cell_line %in% clicked_names)
    }
    
    # Calculate slope
    slope_all <- if (nrow(plot_data) >= 3) round(coef(lm(gene2_effect ~ gene1_effect, data = plot_data))[2], 3) else NA
    subtitle_text <- paste0("r=", cor_all, ", slope=", slope_all, " (n=", nrow(plot_data), ")")

    # --- FIX: force numeric for ggplot / geom_smooth ---
    plot_data$gene1_effect <- suppressWarnings(as.numeric(plot_data$gene1_effect))
    plot_data$gene2_effect <- suppressWarnings(as.numeric(plot_data$gene2_effect))
    plot_data <- plot_data[is.finite(plot_data$gene1_effect) & is.finite(plot_data$gene2_effect), ]

    
    # Get highlight label size
    highlight_size <- if (!is.null(input$highlight_label_size)) input$highlight_label_size else 3
    
    # Build plot based on hotspot mode
    if (hotspot_mode == "three_panel" && any(!is.na(plot_data$hotspot_status))) {
      # 3-panel mode with pie charts
      plot_data$mut_panel <- factor(
        ifelse(is.na(plot_data$hotspot_status), "WT (0 mut)", 
               ifelse(plot_data$hotspot_status == 0, "WT (0 mut)",
                      ifelse(plot_data$hotspot_status == 1, "1 mutation", "2 mutations"))),
        levels = c("WT (0 mut)", "1 mutation", "2 mutations")
      )
      
      panel_stats <- plot_data %>%
        dplyr::group_by(mut_panel) %>%
        dplyr::summarise(
          n = dplyr::n(),
          r = round(cor(gene1_effect, gene2_effect, use = "complete.obs"), 3),
          slope = round(coef(lm(gene2_effect ~ gene1_effect))[2], 3),
          .groups = "drop"
        )
      
      panel_labels <- setNames(
        paste0(panel_stats$mut_panel, "\n(n=", panel_stats$n, ", r=", panel_stats$r, ", slope=", panel_stats$slope, ")"),
        panel_stats$mut_panel
      )
      
      # Check for cancer type highlighting (same as display)
      highlight_cancer <- if (!is.null(input$pie_highlight_cancer) && nzchar(input$pie_highlight_cancer)) {
        input$pie_highlight_cancer
      } else {
        NULL
      }
      
      # Add cancer highlighting flag
      plot_data$cancer_highlight <- FALSE
      if (!is.null(highlight_cancer) && "OncotreeLineage" %in% names(plot_data)) {
        plot_data$cancer_highlight <- !is.na(plot_data$OncotreeLineage) & 
                                       plot_data$OncotreeLineage == highlight_cancer
      }
      
      # Compute cancer colors consistently with pie chart
      scatter_cancer_colors <- NULL
      highlight_cancer_color <- "#f97316"
      
      if ("OncotreeLineage" %in% names(plot_data) && !all(is.na(plot_data$OncotreeLineage))) {
        cancer_counts_for_colors <- plot_data %>%
          dplyr::filter(!is.na(OncotreeLineage) & OncotreeLineage != "") %>%
          dplyr::group_by(mut_panel, OncotreeLineage) %>%
          dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
          dplyr::group_by(mut_panel) %>%
          dplyr::mutate(total = sum(n), pct = n / total * 100) %>%
          dplyr::ungroup()
        
        if (!is.null(highlight_cancer)) {
          cancer_counts_for_colors <- cancer_counts_for_colors %>%
            dplyr::filter(pct >= 2 | OncotreeLineage == highlight_cancer)
        } else {
          cancer_counts_for_colors <- cancer_counts_for_colors %>%
            dplyr::filter(pct >= 2)
        }
        
        all_cancer_types_sorted <- sort(unique(cancer_counts_for_colors$OncotreeLineage))
        n_cancer_colors <- length(all_cancer_types_sorted)
        
        color_palette_exp <- c("#e11d48", "#f97316", "#eab308", "#22c55e", "#14b8a6", 
                               "#0ea5e9", "#6366f1", "#a855f7", "#ec4899", "#78716c",
                               "#dc2626", "#84cc16", "#06b6d4", "#8b5cf6", "#f43f5e",
                               "#64748b", "#059669", "#7c3aed", "#db2777", "#ca8a04",
                               "#991b1b", "#166534", "#1e40af", "#6b21a8", "#be185d",
                               "#854d0e", "#0f766e", "#4338ca", "#9d174d", "#365314")
        if (n_cancer_colors > length(color_palette_exp)) {
          extra_cols <- colorRampPalette(c("#e11d48", "#0ea5e9", "#22c55e", "#a855f7"))(n_cancer_colors - length(color_palette_exp))
          color_palette_exp <- c(color_palette_exp, extra_cols)
        }
        scatter_cancer_colors <- setNames(color_palette_exp[seq_len(n_cancer_colors)], all_cancer_types_sorted)
        
        if (!is.null(highlight_cancer) && highlight_cancer %in% names(scatter_cancer_colors)) {
          highlight_cancer_color <- scatter_cancer_colors[highlight_cancer]
        }
      }
      
      # Build scatter plot with or without cancer highlighting
      if (!is.null(highlight_cancer) && sum(plot_data$cancer_highlight) > 0) {
        non_highlighted <- plot_data[!plot_data$cancer_highlight, ]
        highlighted <- plot_data[plot_data$cancer_highlight, ]
        
        p_scatter <- ggplot(plot_data, aes(x = gene1_effect, y = gene2_effect)) +
          geom_smooth(method = "lm", se = TRUE, color = "#059669", fill = "#86efac", linewidth = 0.8) +
          geom_point(data = non_highlighted, color = "#d1d5db", alpha = 0.4, size = 1.5) +
          geom_point(data = highlighted, color = highlight_cancer_color, alpha = 0.8, size = 2.5) +
          facet_wrap(~mut_panel, labeller = labeller(mut_panel = panel_labels), nrow = 1) +
          coord_cartesian(xlim = sanitize_limits(xlim), ylim = sanitize_limits(ylim)) +
          labs(
            title = paste0(gene1, " vs ", gene2, " - ", hotspot_gene, " hotspot mutation stratification"),
            subtitle = paste0("Highlighting: ", highlight_cancer),
            x = paste0(gene1, " CRISPR Effect"),
            y = paste0(gene2, " CRISPR Effect")
          ) +
          theme_minimal(base_size = 11) +
          theme(legend.position = "none", plot.title = element_text(face = "bold", size = 13),
                plot.subtitle = element_text(color = highlight_cancer_color, size = 10),
                strip.text = element_text(face = "bold", size = 9), 
                panel.border = element_rect(color = "#374151", fill = NA),
                aspect.ratio = 1)
      } else {
        p_scatter <- ggplot(plot_data, aes(x = gene1_effect, y = gene2_effect)) +
          geom_smooth(method = "lm", se = TRUE, color = "#059669", fill = "#86efac", linewidth = 0.8) +
          geom_point(aes(color = mut_panel), alpha = 0.6, size = 1.8) +
          scale_color_manual(values = c("WT (0 mut)" = "#9ca3af", "1 mutation" = "#3b82f6", "2 mutations" = "#dc2626")) +
          facet_wrap(~mut_panel, labeller = labeller(mut_panel = panel_labels), nrow = 1) +
          coord_cartesian(xlim = sanitize_limits(xlim), ylim = sanitize_limits(ylim)) +
          labs(
            title = paste0(gene1, " vs ", gene2, " - ", hotspot_gene, " hotspot mutation stratification"),
            x = paste0(gene1, " CRISPR Effect"),
            y = paste0(gene2, " CRISPR Effect")
          ) +
          theme_minimal(base_size = 11) +
          theme(legend.position = "none", plot.title = element_text(face = "bold", size = 13),
                strip.text = element_text(face = "bold", size = 9), 
                panel.border = element_rect(color = "#374151", fill = NA),
                aspect.ratio = 1)
      }
      
      # Add highlighted points
      if (sum(plot_data$highlight) > 0) {
        highlight_data <- plot_data[plot_data$highlight, ]
        p_scatter <- p_scatter + geom_point(data = highlight_data, color = "#dc2626", size = 3, alpha = 0.9) +
          ggrepel::geom_text_repel(data = highlight_data, aes(label = cell_line), size = highlight_size, color = "#dc2626")
      }
      
      # Create pie charts for cancer type distribution
      if ("OncotreeLineage" %in% names(plot_data) && !all(is.na(plot_data$OncotreeLineage))) {
        # Get all cancer counts first (before filtering)
        cancer_counts_all <- plot_data %>%
          dplyr::filter(!is.na(OncotreeLineage) & OncotreeLineage != "") %>%
          dplyr::group_by(mut_panel, OncotreeLineage) %>%
          dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
          dplyr::group_by(mut_panel) %>%
          dplyr::mutate(total = sum(n), pct = n / total * 100) %>%
          dplyr::ungroup()
        
        # Filter to >= 2% but keep highlighted cancer
        if (!is.null(highlight_cancer)) {
          cancer_counts <- cancer_counts_all %>%
            dplyr::filter(pct >= 2 | OncotreeLineage == highlight_cancer)
        } else {
          cancer_counts <- cancer_counts_all %>%
            dplyr::filter(pct >= 2)
        }
        
        panel_totals <- plot_data %>%
          dplyr::filter(!is.na(OncotreeLineage) & OncotreeLineage != "") %>%
          dplyr::group_by(mut_panel) %>%
          dplyr::summarise(total_n = dplyr::n(), .groups = "drop")
        
        # If highlighting, add stats for center label
        if (!is.null(highlight_cancer)) {
          highlight_stats <- cancer_counts_all %>%
            dplyr::filter(OncotreeLineage == highlight_cancer)
          panel_totals <- panel_totals %>%
            dplyr::left_join(highlight_stats %>% dplyr::select(mut_panel, n, pct), by = "mut_panel")
          
          pie_mode <- if (!is.null(input$pie_scale_mode)) input$pie_scale_mode else "count"
          if (pie_mode == "percent") {
            panel_totals$center_label <- ifelse(
              is.na(panel_totals$n), "0%", paste0(round(panel_totals$pct, 1), "%"))
          } else {
            panel_totals$center_label <- ifelse(
              is.na(panel_totals$n), "n=0", paste0("n=", panel_totals$n))
          }
        } else {
          panel_totals$center_label <- paste0("n=", panel_totals$total_n)
        }
        
        all_cancer_types <- sort(unique(cancer_counts$OncotreeLineage))
        n_colors <- length(all_cancer_types)
        
        color_palette <- c("#e11d48", "#f97316", "#eab308", "#22c55e", "#14b8a6", 
                          "#0ea5e9", "#6366f1", "#a855f7", "#ec4899", "#78716c",
                          "#dc2626", "#84cc16", "#06b6d4", "#8b5cf6", "#f43f5e",
                          "#64748b", "#059669", "#7c3aed", "#db2777", "#ca8a04",
                          "#991b1b", "#166534", "#1e40af", "#6b21a8", "#be185d",
                          "#854d0e", "#0f766e", "#4338ca", "#9d174d", "#365314")
        if (n_colors > length(color_palette)) {
          extra_colors <- colorRampPalette(c("#e11d48", "#0ea5e9", "#22c55e", "#a855f7"))(n_colors - length(color_palette))
          color_palette <- c(color_palette, extra_colors)
        }
        cancer_colors <- setNames(color_palette[1:n_colors], all_cancer_types)
        
        cancer_totals <- cancer_counts %>%
          dplyr::group_by(OncotreeLineage) %>%
          dplyr::summarise(total = sum(n), .groups = "drop") %>%
          dplyr::arrange(desc(total))
        
        pie_mode <- if (!is.null(input$pie_scale_mode)) input$pie_scale_mode else "count"
        if (pie_mode == "percent") {
          cancer_counts <- cancer_counts %>%
            dplyr::group_by(mut_panel) %>%
            dplyr::mutate(n_display = n / sum(n) * 100) %>%
            dplyr::ungroup()
        } else {
          cancer_counts$n_display <- cancer_counts$n
        }
        
        legend_labels <- setNames(
          paste0(cancer_totals$OncotreeLineage, " (n=", cancer_totals$total, ")"),
          cancer_totals$OncotreeLineage
        )
        
        # Gray out colors if highlighting
        if (!is.null(highlight_cancer) && highlight_cancer %in% all_cancer_types) {
          highlight_color_pie <- cancer_colors[highlight_cancer]
          cancer_colors[] <- "#e5e7eb"
          cancer_colors[highlight_cancer] <- highlight_color_pie
        }
        
        p_pie <- ggplot(cancer_counts, aes(x = "", y = n_display, fill = OncotreeLineage)) +
          geom_bar(stat = "identity", width = 1) +
          coord_polar("y", start = 0) +
          facet_wrap(~mut_panel, nrow = 1) +
          scale_fill_manual(values = cancer_colors, labels = legend_labels, name = NULL) +
          geom_text(data = panel_totals, aes(x = 0, y = 0, label = center_label), 
                    inherit.aes = FALSE, size = 3.5, fontface = "bold") +
          theme_void() +
          theme(strip.text = element_blank(), legend.position = "bottom",
                legend.text = element_text(size = 7), legend.key.size = unit(0.3, "cm"),
                plot.background = element_rect(fill = "white", color = NA),
                panel.background = element_rect(fill = "white", color = NA)) +
          guides(fill = guide_legend(nrow = 3, byrow = TRUE))
        
        p <- cowplot::plot_grid(p_scatter, p_pie, ncol = 1, rel_heights = c(1.8, 1.2))
      } else {
        p <- p_scatter
      }
      
      return(p)
      
    } else if (hotspot_mode == "color_by_mut" && any(!is.na(plot_data$hotspot_status))) {
      # Color by mutation mode
      wt_data <- plot_data[is.na(plot_data$hotspot_status) | plot_data$hotspot_status == 0, ]
      mut1_data <- plot_data[!is.na(plot_data$hotspot_status) & plot_data$hotspot_status == 1, ]
      mut2_data <- plot_data[!is.na(plot_data$hotspot_status) & plot_data$hotspot_status == 2, ]
      
      n_wt <- nrow(wt_data)
      n_mut1 <- nrow(mut1_data)
      n_mut2 <- nrow(mut2_data)
      n_total <- nrow(plot_data)
      
      # Calculate percentages
      pct_wt <- round(100 * n_wt / n_total, 1)
      pct_mut1 <- round(100 * n_mut1 / n_total, 1)
      pct_mut2 <- round(100 * n_mut2 / n_total, 1)
      
      # Calculate r and slope for each group
      r_wt <- if (n_wt >= 3) round(cor(wt_data$gene1_effect, wt_data$gene2_effect, use = "complete.obs"), 3) else NA
      r_mut1 <- if (n_mut1 >= 3) round(cor(mut1_data$gene1_effect, mut1_data$gene2_effect, use = "complete.obs"), 3) else NA
      r_mut2 <- if (n_mut2 >= 3) round(cor(mut2_data$gene1_effect, mut2_data$gene2_effect, use = "complete.obs"), 3) else NA
      
      slope_wt <- if (n_wt >= 3) round(coef(lm(gene2_effect ~ gene1_effect, data = wt_data))[2], 3) else NA
      slope_mut1 <- if (n_mut1 >= 3) round(coef(lm(gene2_effect ~ gene1_effect, data = mut1_data))[2], 3) else NA
      slope_mut2 <- if (n_mut2 >= 3) round(coef(lm(gene2_effect ~ gene1_effect, data = mut2_data))[2], 3) else NA
      
      # Build detailed subtitle with line breaks
      subtitle_text <- paste0(subtitle_text, "\n", hotspot_gene, " mutations:")
      subtitle_text <- paste0(subtitle_text, "\n  WT (gray): n=", n_wt, " (", pct_wt, "%), r=", r_wt, ", slope=", slope_wt)
      subtitle_text <- paste0(subtitle_text, "\n  1 mut (blue): n=", n_mut1, " (", pct_mut1, "%), r=", r_mut1, ", slope=", slope_mut1)
      subtitle_text <- paste0(subtitle_text, "\n  2 mut (red): n=", n_mut2, " (", pct_mut2, "%), r=", r_mut2, ", slope=", slope_mut2)
      
      wt_data_no_hl <- wt_data[!wt_data$highlight, ]
      mut1_data_no_hl <- mut1_data[!mut1_data$highlight, ]
      mut2_data_no_hl <- mut2_data[!mut2_data$highlight, ]
      
      if (nrow(wt_data_no_hl) > 0) wt_data_no_hl$mutation_label <- paste0("WT (n=", n_wt, ", ", pct_wt, "%)")
      if (nrow(mut1_data_no_hl) > 0) mut1_data_no_hl$mutation_label <- paste0("1 mut (n=", n_mut1, ", ", pct_mut1, "%)")
      if (nrow(mut2_data_no_hl) > 0) mut2_data_no_hl$mutation_label <- paste0("2 mut (n=", n_mut2, ", ", pct_mut2, "%)")
      
      dfs_to_bind <- list()
      if (nrow(wt_data_no_hl) > 0) dfs_to_bind <- c(dfs_to_bind, list(wt_data_no_hl))
      if (nrow(mut1_data_no_hl) > 0) dfs_to_bind <- c(dfs_to_bind, list(mut1_data_no_hl))
      if (nrow(mut2_data_no_hl) > 0) dfs_to_bind <- c(dfs_to_bind, list(mut2_data_no_hl))
      
      if (length(dfs_to_bind) > 0) {
        legend_data <- do.call(rbind, dfs_to_bind)
        legend_data$mutation_label <- factor(legend_data$mutation_label, 
          levels = c(paste0("WT (n=", n_wt, ", ", pct_wt, "%)"), 
                     paste0("1 mut (n=", n_mut1, ", ", pct_mut1, "%)"), 
                     paste0("2 mut (n=", n_mut2, ", ", pct_mut2, "%)")))
      } else {
        legend_data <- plot_data[0, ]
        legend_data$mutation_label <- factor(character(0))
      }
      
      p <- ggplot(plot_data, aes(x = gene1_effect, y = gene2_effect)) +
        geom_smooth(method = "lm", se = TRUE, color = "#059669", fill = "#86efac") +
        geom_point(data = legend_data, aes(color = mutation_label, size = mutation_label, alpha = mutation_label)) +
        scale_color_manual(name = hotspot_gene, values = c("#9ca3af", "#3b82f6", "#dc2626")) +
        scale_size_manual(name = hotspot_gene, values = c(2, 3, 3.5)) +
        scale_alpha_manual(name = hotspot_gene, values = c(0.3, 0.5, 0.6)) +
        guides(color = guide_legend(title = hotspot_gene, override.aes = list(size = 4)), size = "none", alpha = "none") +
        coord_cartesian(xlim = sanitize_limits(xlim), ylim = sanitize_limits(ylim)) +
        labs(title = paste0(gene1, " vs ", gene2), subtitle = subtitle_text,
             x = paste0(gene1, " CRISPR Effect"), y = paste0(gene2, " CRISPR Effect")) +
        theme_minimal(base_size = 14) +
        theme(plot.title = element_text(face = "bold", size = 16), plot.subtitle = element_text(color = "gray30"),
              legend.position = "right", aspect.ratio = 1)
              
    } else {
      # Default mode - all points green
      p <- ggplot(plot_data, aes(x = gene1_effect, y = gene2_effect)) +
        geom_smooth(method = "lm", se = TRUE, color = "#059669", fill = "#86efac") +
        geom_point(data = plot_data[!plot_data$highlight, ], alpha = 0.5, size = 2, color = "#16a34a") +
        coord_cartesian(xlim = sanitize_limits(xlim), ylim = sanitize_limits(ylim)) +
        labs(title = paste0(gene1, " vs ", gene2), subtitle = subtitle_text,
             x = paste0(gene1, " CRISPR Effect"), y = paste0(gene2, " CRISPR Effect")) +
        theme_minimal(base_size = 14) +
        theme(plot.title = element_text(face = "bold", size = 16), plot.subtitle = element_text(color = "gray30"), aspect.ratio = 1)
    }
    
    # Add highlighted points for non-3panel modes
    if (!hotspot_mode %in% c("three_panel") && sum(plot_data$highlight) > 0) {
      p <- p + 
        geom_point(data = plot_data[plot_data$highlight, ], alpha = 0.9, size = 4, color = "#dc2626") +
        ggrepel::geom_text_repel(data = plot_data[plot_data$highlight, ], aes(label = cell_line),
                                 size = highlight_size, color = "#dc2626", max.overlaps = 20)
    }
    
    p
  })
  
  # Download scatter plot as PNG
# ---- Inspect modal: correlations by tissue/lineage ----
observeEvent(input$scatter_by_tissue, {
  req(rv$current_plot_data)
  req(rv$current_scatter_genes)
  df <- rv$current_plot_data

  # Try common lineage column names
  lineage_col <- NULL
  for (cand in c("lineage", "OncotreeLineage", "oncotree_lineage", "Tissue", "tissue")) {
    if (cand %in% colnames(df)) { lineage_col <- cand; break }
  }
  if (is.null(lineage_col)) {
    showNotification("No tissue/lineage column available for this correlation.", type = "warning")
    return()
  }

  df$.__lineage <- df[[lineage_col]]
  df <- df[!is.na(df$.__lineage) & df$.__lineage != "", ]

  if (nrow(df) < 3) {
    showNotification("Too few samples to compute tissue-specific correlations.", type = "warning")
    return()
  }

  res_list <- lapply(split(df, df$.__lineage), function(d) {
    n <- nrow(d)
    if (n < 3) return(NULL)
    r <- suppressWarnings(cor(d$gene1_effect, d$gene2_effect, use = "complete.obs"))
    slope <- suppressWarnings(coef(lm(gene2_effect ~ gene1_effect, data = d))[2])
    data.frame(
      lineage = as.character(d$.__lineage[1]),
      n = n,
      r = r,
      slope = slope,
      stringsAsFactors = FALSE
    )
  })
  res <- do.call(rbind, res_list)
  if (is.null(res) || nrow(res) == 0) {
    showNotification("No tissue/lineage groups have n >= 3 for this gene pair.", type = "warning")
    return()
  }
  res$r <- round(res$r, 3)
  res$slope <- round(res$slope, 3)
  res <- res[order(-abs(res$r)), ]

  showModal(modalDialog(
    title = paste0("By tissue/lineage: ", rv$current_scatter_genes[1], " vs ", rv$current_scatter_genes[2]),
    size = "l",
    easyClose = TRUE,
    DT::DTOutput("by_tissue_table"),
    footer = tagList(
      downloadButton("download_by_tissue_csv", "Download CSV"),
      modalButton("Close")
    )
  ))

  output$by_tissue_table <- DT::renderDT(res, options = list(pageLength = 15))
  output$download_by_tissue_csv <- downloadHandler(
    filename = function() paste0(rv$current_scatter_genes[1], "_vs_", rv$current_scatter_genes[2], "_by_tissue.csv"),
    content = function(file) write.csv(res, file, row.names = FALSE)
  )
}, ignoreInit = TRUE)

output$download_scatter_png <- downloadHandler(
    filename = function() {
      paste0(rv$current_scatter_genes[1], "_vs_", rv$current_scatter_genes[2], "_", 
             format(Sys.Date(), "%Y%m%d"), ".png")
    },
    content = function(file) {
      # Use wider dimensions for 3-panel mode with pie charts
      hotspot_mode <- if (!is.null(input$hotspot_mode)) input$hotspot_mode else "none"
      if (hotspot_mode == "three_panel") {
        width <- 14
        height <- 10
      } else {
        width <- input$scatter_width_input %||% 8
        height <- input$scatter_height_input %||% 8
      }
      ggsave(file, plot = current_scatter_plot(), width = width, height = height, dpi = 300)
    }
  )
  
  # Download scatter plot as SVG
  output$download_scatter_svg <- downloadHandler(
    filename = function() {
      paste0(rv$current_scatter_genes[1], "_vs_", rv$current_scatter_genes[2], "_", 
             format(Sys.Date(), "%Y%m%d"), ".svg")
    },
    content = function(file) {
      # Use wider dimensions for 3-panel mode with pie charts
      hotspot_mode <- if (!is.null(input$hotspot_mode)) input$hotspot_mode else "none"
      if (hotspot_mode == "three_panel") {
        width <- 14
        height <- 10
      } else {
        width <- input$scatter_width_input %||% 8
        height <- input$scatter_height_input %||% 8
      }
      ggsave(file, plot = current_scatter_plot(), width = width, height = height)
    }
  )
  
  # Download scatter plot data as CSV (includes mutation status if hotspot selected)
  output$download_scatter_csv <- downloadHandler(
    filename = function() {
      paste0(rv$current_scatter_genes[1], "_vs_", rv$current_scatter_genes[2], "_data_", 
             format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      req(rv$current_plot_data)
      req(rv$current_scatter_genes)
      
      # Start with basic export data
      export_data <- data.frame(
        ModelID = rv$current_plot_data$cell_line_id,
        CellLineName = rv$current_plot_data$cell_line,
        stringsAsFactors = FALSE
      )
      
      # Add gene effect columns with proper names
      export_data[[rv$current_scatter_genes[1]]] <- rv$current_plot_data$gene1_effect
      export_data[[rv$current_scatter_genes[2]]] <- rv$current_plot_data$gene2_effect
      
      # Add hotspot mutation status - read from file if selected
      hotspot_gene_input <- if (!is.null(input$hotspot_gene)) trimws(input$hotspot_gene) else ""
      hotspot_col <- hotspot_gene_input  # The dropdown value IS the column name
      hotspot_gene <- gsub("\\s*\\([^)]+\\)$", "", hotspot_gene_input)  # Clean name for column header
      
      if (!is.null(rv$hotspot_file_path) && nzchar(hotspot_col) && !is.null(rv$hotspot_id_col)) {
        # Read hotspot data for the selected gene
        hs_dt <- tryCatch({
          data.table::fread(rv$hotspot_file_path,
                            select = c(rv$hotspot_id_col, hotspot_col),
                            showProgress = FALSE)
        }, error = function(e) NULL)
        
        if (!is.null(hs_dt) && hotspot_col %in% names(hs_dt)) {
          data.table::setnames(hs_dt, c(rv$hotspot_id_col, hotspot_col), c("cell_line_id", "hotspot_status"))
          hs_dt <- hs_dt[!duplicated(hs_dt$cell_line_id), ]
          
          # Match to export data
          match_idx <- match(export_data$ModelID, hs_dt$cell_line_id)
          col_name <- paste0(hotspot_gene, "_mutation")
          export_data[[col_name]] <- hs_dt$hotspot_status[match_idx]
        }
      }
      
      # Merge in full metadata if available
      if (!is.null(rv$full_metadata) && !is.null(rv$current_plot_data$cell_line_id)) {
        # Select useful columns from Model.csv (exclude very long or redundant ones)
        metadata_cols <- c("ModelID", "CellLineName", "StrippedCellLineName", 
                          "OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype",
                          "Age", "Sex", "PrimaryOrMetastasis", "SampleCollectionSite",
                          "SourceType", "SourceDetail")
        
        # Only keep columns that exist in the metadata
        available_cols <- intersect(metadata_cols, colnames(rv$full_metadata))
        
        if (length(available_cols) > 0) {
          metadata_subset <- rv$full_metadata[, available_cols, drop = FALSE]
          
          # Merge by ModelID
          export_data <- merge(
            export_data, 
            metadata_subset, 
            by = "ModelID", 
            all.x = TRUE,
            suffixes = c("", "_meta")
          )
          
          # Remove duplicate CellLineName column if created
          if ("CellLineName_meta" %in% colnames(export_data)) {
            export_data$CellLineName_meta <- NULL
          }
          
          # Reorder columns: ModelID, CellLineName, gene effects, mutation status, then metadata
          gene_cols <- rv$current_scatter_genes
          mut_col_pattern <- paste0(hotspot_gene, "_mutation")
          mut_col <- if (mut_col_pattern %in% colnames(export_data)) mut_col_pattern else NULL
          priority_cols <- c("ModelID", "CellLineName", gene_cols)
          if (!is.null(mut_col)) {
            priority_cols <- c(priority_cols, mut_col)
          }
          other_cols <- setdiff(colnames(export_data), priority_cols)
          export_data <- export_data[, c(priority_cols, other_cols), drop = FALSE]
        }
      }
      
      # Write clean CSV
      fwrite(export_data, file = file, row.names = FALSE)
    }
  )
  
  # Clusters table
  output$clusters_table <- renderDT({
    req(rv$results$success)
    
    clusters_display <- copy(rv$results$clusters)
    
    # Add LFC/FDR if gene stats are loaded
    if (!is.null(rv$gene_stats)) {
      gene_stats_copy <- copy(rv$gene_stats)
      setnames(gene_stats_copy, c("gene", "LFC", "FDR"), c("Gene", "LFC", "FDR"))
      clusters_display <- merge(clusters_display, gene_stats_copy, by = "Gene", all.x = TRUE)
      
      # Round values
      if ("LFC" %in% colnames(clusters_display)) {
        clusters_display[, LFC := round(LFC, 1)]
        clusters_display[, FDR := sapply(FDR, function(x) if(is.na(x)) NA_character_ else formatC(x, format = "e", digits = 1))]
      }
    }
    
    datatable(clusters_display,
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
      "=== Input Genes ===",
      paste("Total input genes:", length(rv$gene_list)),
      paste(rv$gene_list, collapse = ", "),
      "",
            "=== Synonyms Used ===",
      (if (!is.null(rv$synonym_map) && length(rv$synonym_map) > 0) paste(sapply(names(rv$synonym_map), function(k) paste0(k, " -> ", rv$synonym_map[[k]])), collapse = ", ") else "None"),
      "",
"=== Analysis Results ===",
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
  
  # Download summary
  output$download_summary <- downloadHandler(
    filename = function() {
      paste0("correlation_summary_", format(Sys.Date(), "%Y%m%d"), ".txt")
    },
    content = function(file) {
      req(rv$results$success)
      
      summary_text <- paste(
        "=== Analysis Summary ===",
        paste("Date:", Sys.time()),
        paste("Mode:", ifelse(rv$results$mode == "analysis",
                              "Analysis (within gene list)",
                              "Design (find correlated genes)")),
        paste("Correlation cutoff:", rv$results$cutoff),
        "",
        "=== Input Genes ===",
        paste("Total input genes:", length(rv$gene_list)),
        paste(rv$gene_list, collapse = ", "),
                "",
        "=== Synonyms Used ===",
        (if (!is.null(rv$synonym_map) && length(rv$synonym_map) > 0) paste(sapply(names(rv$synonym_map), function(k) paste0(k, " -> ", rv$synonym_map[[k]])), collapse = ", ") else "None"),
"",
        "=== Analysis Results ===",
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
      
      writeLines(summary_text, file)
    }
  )

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
  
  # Download all files as zip
  output$download_all <- downloadHandler(
    filename = function() {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      paste0("correlation_analysis_", timestamp, ".zip")
    },
    content = function(file) {
      req(rv$results$success)

      # Create temporary directory
      temp_dir <- tempdir()
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      uid <- paste0(sample(c(letters, LETTERS, 0:9), 6, replace = TRUE), collapse = "")
      export_folder <- paste0("correlation_export_", timestamp, "_", uid)
      export_dir <- file.path(temp_dir, export_folder)
      dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)
# Export legend files (standalone)
legend_png <- file.path(export_dir, paste0("legend_", timestamp, ".png"))
legend_svg <- file.path(export_dir, paste0("legend_", timestamp, ".svg"))
try({
  png(legend_png, width = 1300, height = 420, res = 150, bg = "white")
  par(mar = c(0.8, 0.8, 0.8, 0.8))
  draw_combined_legend(cutoff = rv$results$cutoff, edge_width = if (!is.null(input$net_edge_width)) input$net_edge_width else 3)
  dev.off()
}, silent = TRUE)
try({
  svg(legend_svg, width = 1300/150, height = 420/150)
  par(mar = c(0.8, 0.8, 0.8, 0.8))
  draw_combined_legend(cutoff = rv$results$cutoff, edge_width = if (!is.null(input$net_edge_width)) input$net_edge_width else 3)
  dev.off()
}, silent = TRUE)


      # Save correlations
      corr_file <- file.path(export_dir, paste0("correlations_", timestamp, ".csv"))
      fwrite(rv$results$correlations, corr_file)
      
      # Save clusters
      clusters_file <- file.path(export_dir, paste0("clusters_", timestamp, ".csv"))
      clusters_display <- copy(rv$results$clusters)
      if (!is.null(rv$gene_stats)) {
        gene_stats_copy <- copy(rv$gene_stats)
        setnames(gene_stats_copy, c("gene", "LFC", "FDR"), c("Gene", "LFC", "FDR"))
        clusters_display <- merge(clusters_display, gene_stats_copy, by = "Gene", all.x = TRUE)
      }
      fwrite(clusters_display, clusters_file)
      
      # Save summary
      summary_file <- file.path(export_dir, paste0("summary_", timestamp, ".txt"))
      summary_text <- paste(
        "=== Analysis Summary ===",
        paste("Date:", Sys.time()),
        paste("Mode:", ifelse(rv$results$mode == "analysis",
                              "Analysis (within gene list)",
                              "Design (find correlated genes)")),
        paste("Correlation cutoff:", rv$results$cutoff),
        "",
        "=== Input Genes ===",
        paste("Total input genes:", length(rv$gene_list)),
        paste(rv$gene_list, collapse = ", "),
        "",
        "=== Analysis Results ===",
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
      writeLines(summary_text, summary_file)
      
      # Create list of files to zip
      files_to_zip <- c(corr_file, clusters_file, summary_file)
      
      # Include last exported network figure (if available)
      if (!is.null(rv$last_network_png_file) && file.exists(rv$last_network_png_file)) {
        file.copy(rv$last_network_png_file, file.path(export_dir, paste0("network_figure_", timestamp, ".png")), overwrite = TRUE)
      }
      
      # Verify legend files were created and add to export
      files_in_dir <- list.files(export_dir, full.names = FALSE)

      # Zip export folder (preserve directory structure)
      old_wd <- getwd()
      setwd(temp_dir)
      on.exit(setwd(old_wd), add = TRUE)
      zip(file, export_folder, flags = "-r9X")
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
