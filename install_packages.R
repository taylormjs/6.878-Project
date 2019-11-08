# Download instructions: https://bioconductor.org/packages/release/bioc/html/EpiDISH.html

# NOTE(milo): Had to run rstudio as sudo to be able to install without errors.
install.packages("BiocManager")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EpiDISH")
BiocManager::install("GEOquery")
BiocManager::install("Biobase")

# Make sure we can load in the libraries.
library("EpiDISH")
library("Biobase")
library("GEOquery")