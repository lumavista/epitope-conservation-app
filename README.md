# Conserved Human MHC-I Epitope Predictor

This Shiny application provides an alignment-consistent consensus prediction pipeline for conserved human MHC-I epitopes.

## Copyright and License

**Copyright © 2025 LumaVista Bio. All rights reserved.**

This software is proprietary. Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated data files (the "Software"), to use the Software for **non-commercial, academic, and research purposes only**.

This permission **does not** include the right to modify, merge, publish, distribute, sublicense, and/or sell copies of the Software. The Software is provided "as is," without warranty of any kind.

## Data Source and Citation

The data used in this application (specifically the `combined_db.rds` file) is compiled from the **Immune Epitope Database (IEDB)**, a free public resource. All predictions are powered by the IEDB's public API.
Per the IEDB's terms of use, any use of this tool for research that results in a publication should cite the IEDB's main reference paper. As of late 2024, the primary citation is:
> Vita R, et al. The Immune Epitope Database (IEDB) in 2024. Nucleic Acids Res. 2024 Jan 5;52(D1):D1263-D1269. (doi: 10.1093/nar/gkad1046)

## How to Run (Local Development)

1.  **Clone this repository:**
    ```bash
    git clone https://github.com/lumavista/epitope-conservation-app.git
cd epitope-conservation-app
    ```

2.  **Install Dependencies:**
    Open R/RStudio and run the installation script:
    ```r
    source("install_packages.R")
    ```

3.  **Run the App:**
    In RStudio, click the **"Run App"** button. Or, from the console:
    ```r
    shiny::runApp()
    ```

---

## How to Deploy on a Shiny Server (Linux)

These instructions are for a typical Linux server (e.g., Ubuntu/Debian).

### 1. Install System Dependencies

You must install R and the system libraries needed by the R packages.

```bash
# Update package lists
sudo apt-get update

# Install R and essential build tools
sudo apt-get install -y r-base r-base-dev build-essential

# Install system libraries required by R packages
sudo apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev
```

### 2. Install R Packages

As root (or a user with permission), navigate to the app directory and run the installation script:

```bash
cd epitope-conservation-app
sudo Rscript install_packages.R
```

### 3. Run the App

You can now run the app from the R console:

```r
shiny::runApp(port = 3838, host = "0.0.0.0")
```

(Or, if you have Shiny Server installed, copy the app directory to `/srv/shiny-server/`.)
