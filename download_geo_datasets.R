library(Biobase)    # Boxplot for selected GEO samples.
library(GEOquery)   # For downloading GEO2R data.

#================= COPIED FROM GEO2R ======================
# See: https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE59999
# This will use a cached version if it is found! Still has to unzip though (takes about a minute).

# NOTE(milo): Only need to run this once!
download_and_save_data_to_analysis_folder = function() {
  # Load MARTINO 2018.
  gset.martino2018 <- getGEO("GSE114134", GSEMatrix =TRUE, getGPL=FALSE)
  if (length(gset.martino2018) > 1) idx <- grep("GPL23976", attr(gset.martino2018, "names")) else idx <- 1
  gset.martino2018 <- gset.martino2018[[idx]]
  save(gset.martino2018, "./analysis/gset.martino2018.Rda")
  
  # Load MARTINO 2015.
  gset.martino2015 <- getGEO("GSE59999", GSEMatrix =TRUE, getGPL=FALSE)
  if (length(gset.martino2015) > 1) idx <- grep("GPL13534", attr(gset.martino2015, "names")) else idx <- 1
  gset.martino2015 <- gset.martino2015[[idx]]
  save(gset.martino2015, file="./analysis/gset.martino2015.Rda")
}