#----------------------------------------------
# R Package Installation Script (Windows-Friendly)
#----------------------------------------------
#
# This script installs all R packages required
# to run the Conserved Epitope Shiny App.
# It forces the installation of pre-compiled "binaries"
# to avoid requiring RTools on Windows.
#
# Run this from the command line:
# $ Rscript install_packages.R
#
# Or run it interactively in an R session.
#----------------------------------------------

# --- CRAN Packages ---
message("Installing CRAN packages (binaries only)...")
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

# Add type = "binary" to skip source builds
install.packages(
  cran_packages, 
  repos = "https://cloud.r-project.org/",
  type = "binary"
)

# --- Bioconductor Packages ---
message("Installing Bioconductor packages (binaries only)...")

# First, install BiocManager itself as a binary
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", type = "binary")
}

bioc_packages <- c(
  "msa", 
  "Biostrings"
)

# Add type = "binary" to skip source builds
BiocManager::install(
  bioc_packages,
  type = "binary",
  update = FALSE # Prevents it from trying to update to a source version
)

message("\n[✓] All R packages installed successfully.")
