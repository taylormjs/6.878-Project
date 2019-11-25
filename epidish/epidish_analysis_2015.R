library(Biobase)    # Boxplot for selected GEO samples.
library(GEOquery)   # For downloading GEO2R data.
library(EpiDISH)
library(scales)
library(LICORS)

source("epidish/epidish_helpers.R")
source("epidish/modified_celldmc.R")


#============== SENSITIZED VS. ALLERGIC =================
# NOTE(milo): This didn't work!
# Load in the dataframes if possible to avoid the parsing time.
load("./analysis/gset.martino2015.Rda")

pheno.martino2015 = makeBinaryPhenotypesSensitizedVsAllergicMartino2015(gset.martino2015)
beta.m.martino2015 = getBetaMatrixSensitizedVsAllergicMartino2015(gset.martino2015)
cellfrac.m.martino2015 = getEpidishCellFrac(beta.m.martino2015)
boxplot(cellfrac.m.martino2015)

# NOTE(milo): We seem to run into linear regression problems when all of the cell types
# are included. If we take this subset of cell types, the ModifiedCellDMC seems to work.
# cellfrac.m.martino2015 = cellfrac.m.martino2015[,c("B", "NK", "CD4T", "CD8T", "Mono", "Eosino")]
boxplot(cellfrac.m.martino2015)

# NOTE(milo): Confirmed that my ModifiedCellDMC gives the same outputs as original CellDMC.
martino2015.celldmc.o <- ModifiedCellDMC(beta.m.martino2015, pheno.martino2015, cellfrac.m.martino2015,
                                         mc.cores=6)
# save(martino2015.celldmc.o, file="./analysis/martino2015/alternate/martino2015_celldmc.o")

# load("./analysis/martino2015_celldmc.o")
summarizeDMCTs(martino2015.celldmc.o)

m2015.dmct = martino2015.celldmc.o$dmct
m2015.signif = m2015.dmct[m2015.dmct[,1] == 1,]

# Write things to .csv so we can escape R.
write.csv(martino2015.celldmc.o$dmct, file="./analysis/martino2015/alternate/martino2015_dmct.csv")
write.csv(martino2015.celldmc.o$coe.change, file="./analysis/martino2015/alternate/martino2015_coe_change.csv")
write.csv(martino2015.celldmc.o$coe.control, file="./analysis/martino2015/alternate/martino2015_coe_control.csv")
write.csv(cellfrac.m.martino2015, file="./analysis/martino2015/alternate/martino2015_cellfrac.csv")

martino2015.phenotypes = gset.martino2015@phenoData@data
write.csv(martino2015.phenotypes, file="./analysis/martino2015/alternate/martino2015_phenotypes.csv")

write.csv(beta.m.martino2015, file="./analysis/martino2015/alternate/martino2015_beta.csv")
#==================================================


#================= NONALLERGIC VS. ALLERGIC ===================
# Load in the dataframes if possible to avoid the parsing time.
load("./analysis/gset.martino2015.Rda")
pheno.martino2015 = makeBinaryPhenotypesNonallergicVsAllergicMartino2015(gset.martino2015)
beta.m.martino2015 = getBetaMatrixNonallergicVsAllergicMartino2015(gset.martino2015)
cellfrac.m.martino2015 = getEpidishCellFrac(beta.m.martino2015)
boxplot(cellfrac.m.martino2015)

# NOTE(milo): We seem to run into linear regression problems when all of the cell types
# are included. If we take this subset of cell types, the ModifiedCellDMC seems to work.
cellfrac.m.martino2015 = cellfrac.m.martino2015[,c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro")]
# cellfrac.m.martino2015 = normalize(cellfrac.m.martino2015, byrow = TRUE, tol = 1e-06)
boxplot(cellfrac.m.martino2015)

# NOTE(milo): Confirmed that my ModifiedCellDMC gives the same outputs as original CellDMC.
# martino2015.celldmc.o <- CellDMC(beta.m.martino2015, pheno.martino2015, cellfrac.m.martino2015,
                                 # mc.cores=6)
martino2015.celldmc.o <- ModifiedCellDMC(beta.m.martino2015, pheno.martino2015, cellfrac.m.martino2015,
                                         mc.cores=6)
# save(martino2015.celldmc.o, file="./analysis/martino2015/nonallergic_vs_allergic/martino2015_celldmc.o")

load("./analysis/martino2015_celldmc.o")
summarizeDMCTs(martino2015.celldmc.o)

m2015.dmct = martino2015.celldmc.o$dmct
m2015.signif = m2015.dmct[m2015.dmct[,1] == 1,]

# Write things to .csv so we can escape R.
write.csv(martino2015.celldmc.o$dmct, file="./analysis/martino2015/nonallergic_vs_allergic/martino2015_dmct.csv")
write.csv(martino2015.celldmc.o$coe.change, file="./analysis/martino2015/nonallergic_vs_allergic/martino2015_coe_change.csv")
write.csv(martino2015.celldmc.o$coe.control, file="./analysis/martino2015/nonallergic_vs_allergic/martino2015_coe_control.csv")
write.csv(cellfrac.m.martino2015, file="./analysis/martino2015/nonallergic_vs_allergic/martino2015_cellfrac.csv")

martino2015.phenotypes = gset.martino2015@phenoData@data
write.csv(martino2015.phenotypes, file="./analysis/martino2015/nonallergic_vs_allergic/martino2015_phenotypes.csv")

write.csv(beta.m.martino2015, file="./analysis/martino2015/nonallergic_vs_allergic/martino2015_beta.csv")
#==================================================

