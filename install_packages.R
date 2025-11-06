#----------------------------------------------
# R Package Installation Script
#----------------------------------------------
#
# This script installs all R packages required
# to run the Conserved Epitope Shiny App.
# Run this from the command line:
# $ Rscript install_packages.R
#
# Or run it interactively in an R session.
#----------------------------------------------

# --- CRAN Packages ---
message("Installing CRAN packages...")
cran_packages <- c(
  "shiny", 
  "httr", 
  "dplyr", 
  "plotly", 
  "DT", 
  "shinyjs", 
  "shinycssloaders", 
  "writexl", 
  "data.table"
)

install.packages(cran_packages, repos = "https://cloud.r-project.org/")

# --- Bioconductor Packages ---
message("Installing Bioconductor packages...")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "msa", 
  "Biostrings"
)

BiocManager::install(bioc_packages)

message("\n[✓] All R packages installed successfully.")
