library(Biobase)    # Boxplot for selected GEO samples.
library(GEOquery)   # For downloading GEO2R data.
library(EpiDISH)
library(scales)
library(LICORS)

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

# NOTE(milo): This should be entirely CD4T cells! It looks like mostly CD4 cells are detected
# with some CD8T cells detected too. Could the cell sorting that Martino does have some error?
# Maybe there actually are some CD8T cells in his samples, which we can use to argue that its
# better to do deconvolution in software during analysis.
boxplot(cellfrac.m.martino2018)

# Need to drop everything except CD4T and CD8T otherwise this will fail!
# https://adv-r.hadley.nz/subsetting.html
cellfrac.m.martino2018 = cellfrac.m.martino2018[,c("CD4T", "CD8T")]
boxplot(cellfrac.m.martino2018)

# martino2018.celldmc.o <- CellDMC(beta.m.martino2018, pheno.martino2018, cellfrac.m.martino2018,
#                                  mc.cores=6)
martino2018.celldmc.o <- ModifiedCellDMC(beta.m.martino2018, pheno.martino2018, cellfrac.m.martino2018,
                                 mc.cores=6)
save(martino2018.celldmc.o, file="./analysis/martino2018_celldmc.o")

load("./analysis/martino2018_celldmc.o")
summarizeDMCTs(martino2018.celldmc.o)

# Write things to .csv so we can escape R.
write.csv(martino2018.celldmc.o$dmct, file="./analysis/martino2018_dmct.csv")
write.csv(martino2018.celldmc.o$coe.change, file="./analysis/martino2018_coe_change.csv")
write.csv(martino2018.celldmc.o$coe.control, file="./analysis/martino2018_coe_control.csv")
#==================================================


#================= MARTINO 2015 ===================
# Load in the dataframes if possible to avoid the parsing time.
load("./analysis/gset.martino2015.Rda")
pheno.martino2015 = makeBinaryPhenotypesMartino2015(gset.martino2015)
beta.m.martino2015 = getBetaMatrixMartino2015(gset.martino2015)
cellfrac.m.martino2015 = getEpidishCellFrac(beta.m.martino2015)
boxplot(cellfrac.m.martino2015)

# NOTE(milo): We seem to run into linear regression problems when all of the cell types
# are included. If we take this subset of cell types, the ModifiedCellDMC seems to work.
cellfrac.m.martino2015 = cellfrac.m.martino2015[,c("B", "NK", "CD4T", "CD8T", "Mono", "Eosino")]
# cellfrac.m.martino2015 = cellfrac.m.martino2015 + 1e-2
# cellfrac.m.martino2015 = normalize(cellfrac.m.martino2015, byrow = TRUE, tol = 1e-06)
boxplot(cellfrac.m.martino2015)

# NOTE(milo): Confirmed that my ModifiedCellDMC gives the same outputs as original CellDMC.
martino2015.celldmc.o <- CellDMC(beta.m.martino2015, pheno.martino2015, cellfrac.m.martino2015,
                                 mc.cores=6)
# martino2015.celldmc.o <- ModifiedCellDMC(beta.m.martino2015, pheno.martino2015, cellfrac.m.martino2015,
                                         # mc.cores=6)
# save(martino2015.celldmc.o, file="./analysis/martino2015_celldmc.o")

load("./analysis/martino2015_celldmc.o")
summarizeDMCTs(martino2015.celldmc.o)

m2015.dmct = martino2015.celldmc.o$dmct
m2015.signif = m2015.dmct[m2015.dmct[,1] == 1,]

# Write things to .csv so we can escape R.
write.csv(martino2015.celldmc.o$dmct, file="./analysis/martino2015_dmct.csv")
write.csv(martino2015.celldmc.o$coe.change, file="./analysis/martino2015_coe_change.csv")
write.csv(martino2015.celldmc.o$coe.control, file="./analysis/martino2015_coe_control.csv")
write.csv(cellfrac.m.martino2015, file="./analysis/martino2015_cellfrac.csv")

martino2015.phenotypes = gset.martino2015@phenoData@data
write.csv(martino2015.phenotypes, file="./analysis/martino2015_phenotypes.csv")

write.csv(beta.m.martino2015, file="./analysis/martino2015_beta.csv")
#==================================================
