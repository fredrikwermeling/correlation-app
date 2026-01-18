# Gene Correlation Explorer

A Shiny app for analyzing gene correlations using DepMap CRISPR screen data.

Based on Fredrik Wermeling's correlation script.

## Features

- **Multiple input methods**: Paste genes, upload file, or fetch from GO terms
- **Two analysis modes**: 
  - Analysis (correlations within your gene list)
  - Design (find new genes correlated with your list)
- **Interactive network visualization**: Zoom, pan, click nodes
- **Downloadable results**: CSV files for correlations and clusters

## Installation

### 1. Install R and RStudio

If not already installed:
- R: https://cran.r-project.org/
- RStudio: https://posit.co/download/rstudio-desktop/

### 2. Install required packages

Open RStudio and run:

```r
install.packages(c("shiny", "data.table", "igraph", "ggplot2", "ggraph", 
                   "DT", "shinyjs", "httr", "jsonlite", "visNetwork"))
```

### 3. Download DepMap reference data

1. Go to: https://depmap.org/portal/data_page/?tab=allData
2. Find "CRISPRGeneEffect.csv" (latest release, ~450MB)
3. Download and save somewhere accessible

## Running the App

### Option 1: Run from RStudio

1. Open `app.R` in RStudio
2. Click "Run App" button (top right of editor)

### Option 2: Run from R console

```r
shiny::runApp("path/to/correlation_app")
```

### Option 3: Run from terminal

```bash
Rscript -e "shiny::runApp('path/to/correlation_app')"
```

## Usage

1. **Load Reference Data**: Upload the DepMap CRISPRGeneEffect.csv file
2. **Input Genes**: 
   - Paste gene symbols (one per line)
   - Upload a .txt file
   - Enter a GO term ID to fetch genes
3. **Set Parameters**:
   - Choose Analysis or Design mode
   - Set correlation cutoff (0.3-0.5 recommended)
4. **Run Analysis**: Click the button and wait
5. **Explore Results**:
   - Interactive network visualization
   - Downloadable correlation and cluster tables

## GO Term Import

You can fetch genes associated with any GO term. Examples:
- `GO:0006955` - immune response
- `GO:0006915` - apoptotic process  
- `GO:0007049` - cell cycle
- `GO:0016310` - phosphorylation

## Tips

- Start with a correlation cutoff of 0.5 and adjust as needed
- Design mode generates larger outputs (genome-wide search)
- Not all gene symbols may match - check the Summary tab for "not found" genes
- Use the [Green Listed synonym tool](https://greenlisted.cmm.se/) to find alternative symbols

## Future Development

Planned features:
- KEGG pathway import
- MSigDB gene set import
- Cell line filtering
- Integration with Green Listed website

## Credits

- DepMap data: https://depmap.org/
- Wermeling Lab: https://wermelinglab.com/
- Green Listed: https://greenlisted.cmm.se/

## License

For research use. Please follow DepMap's terms and conditions when using their data.
