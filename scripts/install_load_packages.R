# Install remotes if not already installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Define a function to install and load packages
install_and_load <- function(packages) {
  for (pkg in packages) {
    # Check if the package is installed
    if (!requireNamespace(pkg, quietly = TRUE)) {
      # Install dive package from GitHub if it is not installed
      if (pkg == "dive") {
        remotes::install_github("AMahajna/dive")
      } else {
        # Handle Bioconductor packages
        if (pkg %in% c("Biostrings", "miaViz", "mia", "scater")) {
          BiocManager::install(pkg)
        } else {
          install.packages(pkg)
        }
      }
    }
    library(pkg, character.only = TRUE)
  }
}

# List of packages to install and load, including the GitHub package
packages <- c(
  "dive", # Include dive in the list
  "corrplot", "plotly", "ggcorrplot", "BiocManager", "Biostrings", 
  "GGally", "miaViz", "mia", "scater", "devtools", "tidyverse", 
  "readxl", "astsa", "readr", "dplyr", "ggplot2", "forecast", 
  "RColorBrewer", "gridExtra", "writexl", "cowplot"
)

# Install and load the packages
install_and_load(packages)