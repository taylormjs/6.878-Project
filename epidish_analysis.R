library(Biobase)    # Boxplot for selected GEO samples.
library(GEOquery)   # For downloading GEO2R data.
library(EpiDISH)

source("epidish_helpers.R")
source("modified_celldmc.R")

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


#================= MARTINO 2018 ===================
# Load in the dataframes if possible to avoid the parsing time.
load("./analysis/gset.martino2018.Rda")
pheno.martino2018 = makeBinaryPhenotypesMartino2018(gset.martino2018)
beta.m.martino2018 = getBetaMatrixMartino2018(gset.martino2018)
cellfrac.m.martino2018 = getEpidishCellFrac(beta.m.martino2018)

# NOTE(milo): This should be entirely CD4T cells!
boxplot(cellfrac.m.martino2018)

# Need to drop everything except CD4T and CD8T otherwise this will fail!
# https://adv-r.hadley.nz/subsetting.html
cellfrac.m.martino2018 = cellfrac.m.martino2018[,c("CD4T", "CD8T")]
boxplot(cellfrac.m.martino2018)

martino2018.celldmc.o <- CellDMC(beta.m.martino2018, pheno.martino2018, cellfrac.m.martino2018)
# save(martino2018.celldmc.o, file="./analysis/martino2018_celldmc.o")

load("./analysis/martino2018_celldmc.o")
summarizeDMCTs(martino2018.celldmc.o)

# Write things to .csv so we can escape R.
dmct = martino2018.celldmc.o$dmct
coe = martino2018.celldmc.o$coe
write.csv(dmct, file="./analysis/martino2018_dmct.csv")
write.csv(coe, file="./analysis/martino2018_coe.csv")
#==================================================


#================= MARTINO 2015 ===================
# Load in the dataframes if possible to avoid the parsing time.
load("./analysis/gset.martino2015.Rda")
pheno.martino2015 = makeBinaryPhenotypesMartino2015(gset.martino2015)
beta.m.martino2015 = getBetaMatrixMartino2015(gset.martino2015)
cellfrac.m.martino2015 = getEpidishCellFrac(beta.m.martino2015)

# beta.m.martino2015 = beta.m.martino2015[1:100,]

# martino2015.celldmc.o <- CellDMC(beta.m.martino2015, pheno.martino2015, cellfrac.m.martino2015)
martino2015.celldmc.o <- ModifiedCellDMC(beta.m.martino2015, pheno.martino2015, cellfrac.m.martino2015,
                                         mc.cores=6)
save(martino2015.celldmc.o, file="./analysis/martino2015_celldmc.o")

load("./analysis/martino2015_celldmc.o")
summarizeDMCTs(martino2015.celldmc.o)

# Write things to .csv so we can escape R.
write.csv(martino2015.celldmc.o$dmct, file="./analysis/martino2015_dmct.csv")
write.csv(martino2015.celldmc.o$coe.change, file="./analysis/martino2015_coe_change.csv")
write.csv(martino2015.celldmc.o$coe.control, file="./analysis/martino2015_coe_control.csv")
#==================================================
