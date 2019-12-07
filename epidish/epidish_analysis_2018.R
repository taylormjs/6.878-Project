library(Biobase)    # Boxplot for selected GEO samples.
library(GEOquery)   # For downloading GEO2R data.
library(EpiDISH)
library(scales)
library(LICORS)

source("epidish/epidish_helpers.R")
source("epidish/modified_celldmc.R")
source("epidish/bulk_celldmc.R")

# Load in the dataframes if possible to avoid the parsing time.
load("./analysis/gset.martino2018.Rda")
load("./analysis/Mvalues2018.Rda")
pheno.martino2018 = makeBinaryPhenotypesMartino2018(gset.martino2018)
beta.m.martino2018 = getBetaMatrixMartino2018(gset.martino2018)

# IMPORTANT - choose the training set here!!!
martino2018.phenotypes = gset.martino2018@phenoData@data
martino2018.phenotypes = martino2018.phenotypes[martino2018.phenotypes$`allergy status:ch1` != "resolved",]
patients = rownames(martino2018.phenotypes)

set.seed(123)
train_idx <- sample(seq_len(length(patients)), size=113)
training_set <- patients[train_idx]
test_set <- patients[-train_idx]

lapply(training_set, write, "./analysis/martino2018/training_set.txt", append=TRUE, ncolumns=1)
lapply(test_set, write, "./analysis/martino2018/test_set.txt", append=TRUE, ncolumns=1)

pheno_train = martino2018.phenotypes[training_set,]
pheno_train = pheno_train$`allergy status:ch1`
pheno_train[pheno_train == "allergic"] = 1
pheno_train[pheno_train == "control"] = 0
pheno_train = as.integer(pheno_train)

# write.csv(martino2018.phenotypes, file="./analysis/martino2018/phenotypes.csv")
# write.csv(Mvalues2018.filt, file="./analysis/martino2018/Mvalues.csv")

# ======================== CELL SPECIFIC ============================
# NOTE(milo): This should be entirely CD4T cells! It looks like mostly CD4 cells are detected
# with some CD8T cells detected too. Could the cell sorting that Martino does have some error?
# Maybe there actually are some CD8T cells in his samples, which we can use to argue that its
# better to do deconvolution in software during analysis.
cellfrac.m.martino2018 = getEpidishCellFrac(beta.m.martino2018)
boxplot(cellfrac.m.martino2018)

# Need to drop everything except CD4T and CD8T otherwise this will fail!
# https://adv-r.hadley.nz/subsetting.html
cellfrac.m.martino2018 = cellfrac.m.martino2018[,c("CD4T", "CD8T")]
cellfrac.m.martino2018 = cellfrac.m.martino2018[training_set,]
boxplot(cellfrac.m.martino2018)

Mvalues2018.filt = getMvalues2018(Mvalues2018, gset.martino2018)
Mvalues_train = Mvalues2018.filt[,training_set]

martino2018.celldmc.o <- ModifiedCellDMC(Mvalues_train, pheno_train, cellfrac.m.martino2018, mc.cores=6)

folder = "./analysis/martino2018/Mvalues_control_vs_allergic/"
save(martino2018.celldmc.o, file=str_c(folder, "celldmc.o"))

load(str_c(folder, "celldmc.o"))
summarizeDMCTs(martino2018.celldmc.o)

# Write things to .csv so we can escape R.
write.csv(martino2018.celldmc.o$dmct, file=str_c(folder, "dmct.csv"))
write.csv(martino2018.celldmc.o$coe.change, file=str_c(folder, "coe_change.csv"))
write.csv(martino2018.celldmc.o$coe.control, file=str_c(folder, "coe_control.csv"))
write.csv(cellfrac.m.martino2018, file=str_c(folder, "cellfrac.csv"))


#================= MARTINO 2018 (M-VALUES) BULK DATA ===================
# Load in the dataframes if possible to avoid the parsing time.
load("./analysis/gset.martino2018.Rda")
load("./analysis/Mvalues2018.Rda")
pheno.martino2018 = makeBinaryPhenotypesMartino2018(gset.martino2018)
beta.m.martino2018 = getBetaMatrixMartino2018(gset.martino2018)

Mvalues2018.filt = getMvalues2018(Mvalues2018, gset.martino2018)
martino2018.bulkdmc.o <- BulkCellDMC(Mvalues2018.filt, pheno.martino2018, mc.cores=6)

folder = "./analysis/martino2018/Mvalues_control_vs_allergic_bulk/"
save(martino2018.bulkdmc.o, file=str_c(folder, "bulkdmc.o"))

# load(str_c(folder, "bulkdmc.o"))

dmct = martino2018.bulkdmc.o$dmct
signif = dmct[dmct[,1] != 0,]
dim(signif)

# Write things to .csv so we can escape R.
write.csv(martino2018.bulkdmc.o$dmct, file=str_c(folder, "dmct.csv"))
write.csv(martino2018.bulkdmc.o$coe.change, file=str_c(folder, "coe_change.csv"))
write.csv(martino2018.bulkdmc.o$coe.control, file=str_c(folder, "coe_control.csv"))

martino2018.phenotypes = gset.martino2018@phenoData@data
write.csv(martino2018.phenotypes, file=str_c(folder, "phenotypes.csv"))
write.csv(Mvalues2018.filt, file=str_c(folder, "Mvalues.csv"))


