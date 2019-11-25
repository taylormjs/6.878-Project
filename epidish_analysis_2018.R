library(Biobase)    # Boxplot for selected GEO samples.
library(GEOquery)   # For downloading GEO2R data.
library(EpiDISH)
library(scales)
library(LICORS)

source("epidish/epidish_helpers.R")
source("epidish/modified_celldmc.R")

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
