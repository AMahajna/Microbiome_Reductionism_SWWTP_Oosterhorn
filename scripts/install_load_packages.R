###############################################################################
# Define a function to install and load packages
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (pkg %in% c("Biostrings", "miaViz", "mia", "scater")) {
        BiocManager::install(pkg)
      } else {
        install.packages(pkg)
      }
    }
    library(pkg, character.only = TRUE)
  }
}

# List of packages to install and load
packages <- c(
  "corrplot", "plotly", "ggcorrplot", "BiocManager", "Biostrings", 
  "GGally", "miaViz", "mia", "scater", "devtools", "tidyverse", 
  "readxl", "astsa", "readr", "dplyr", "ggplot2", "forecast", 
  "RColorBrewer", "gridExtra", "writexl"
)

# Install and load the packages
install_and_load(packages)


