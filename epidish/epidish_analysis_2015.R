library(Biobase)    # Boxplot for selected GEO samples.
library(GEOquery)   # For downloading GEO2R data.
library(EpiDISH)
library(scales)
library(LICORS)

source("epidish/epidish_helpers.R")
source("epidish/modified_celldmc.R")
source("epidish/bulk_celldmc.R")

# Load in the dataframes if possible to avoid the parsing time.
load("./analysis/gset.martino2015.Rda")
load("./analysis/Mvalues2015.Rda")

beta.m.martino2015 = getBetaMatrixNonallergicVsAllergicMartino2015(gset.martino2015)

# IMPORTANT - choose the training set here!!!
martino2015.phenotypes = gset.martino2015@phenoData@data
martino2015.phenotypes = martino2015.phenotypes[martino2015.phenotypes$`challenge outcome:ch1` != "sensitized",]
patients = rownames(martino2015.phenotypes)

set.seed(123)
train_idx <- sample(seq_len(length(patients)), size=32)
training_set <- patients[train_idx]
test_set <- patients[-train_idx]

lapply(training_set, write, "./analysis/martino2015/training_set.txt", append=TRUE, ncolumns=1)
lapply(test_set, write, "./analysis/martino2015/test_set.txt", append=TRUE, ncolumns=1)

# Mvalues2015.filt = Mvalues2015.filt[1:100,]
Mvalues2015.filt = getMvalues2015NonallergicVsAllergic(Mvalues2015, gset.martino2015)
# write.csv(Mvalues2015.filt, file="./analysis/martino2015/Mvalues.csv")
# write.csv(martino2015.phenotypes, file="./analysis/martino2015/phenotypes.csv")

# TRAIN ON ONLY THE TRAINING SUBSET
Mvalues_train = Mvalues2015.filt[, training_set]

pheno_train = martino2015.phenotypes[training_set,]
pheno_train = pheno_train$`challenge outcome:ch1`
pheno_train[pheno_train == "allergic"] = 1
pheno_train[pheno_train == "nonallergic"] = 0
pheno_train = as.integer(pheno_train)

#========================== BULK ANALYSIS ========================
martino2015.bulkdmc.o <- BulkCellDMC(Mvalues_train, pheno_train, mc.cores=6)

# NOTE(milo): MAKE SURE THIS FOLDER IS CORRECT!
folder = "./analysis/martino2015/Mvalues_nonallergic_vs_allergic_bulk/"
save(martino2015.bulkdmc.o, file=str_c(folder, "bulkdmc.o"))

# load(str_c(folder, "bulkdmc.o"))
dmct = martino2015.bulkdmc.o$dmct
coe.change = martino2015.bulkdmc.o$coe.change[[1]]
signif = coe.change[coe.change$adjP < 0.7,]

# Write things to .csv so we can escape R.
write.csv(martino2015.bulkdmc.o$dmct, file=str_c(folder, "dmct.csv"))
write.csv(martino2015.bulkdmc.o$coe.change, file=str_c(folder, "coe_change.csv"))
write.csv(martino2015.bulkdmc.o$coe.control, file=str_c(folder, "coe_control.csv"))
#==================================================


#========================= CELL SPECIFIC ===================
cellfrac.m.martino2015 = getEpidishCellFrac(beta.m.martino2015)
cellfrac.m.martino2015 = cellfrac.m.martino2015[training_set,]
boxplot(cellfrac.m.martino2015)

# NOTE(milo): We seem to run into linear regression problems when all of the cell types
# are included. If we take this subset of cell types, the ModifiedCellDMC seems to work.
# NOTE(milo): Martino 2015 uses PBMC (Peripheral Blood Mononuclear Cells) which doesn't
# include neutrophils and eosinophils.

# cellfrac.m.martino2015 = cellfrac.m.martino2015[,c("B", "NK", "CD4T", "CD8T", "Mono")]
# cellfrac.m.martino2015 = cellfrac.m.martino2015[,c("B", "NK", "CD4T", "CD8T", "Mono", "Eosino")]
cellfrac.m.martino2015 = cellfrac.m.martino2015[,c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro")]
boxplot(cellfrac.m.martino2015)

# NOTE(milo): Confirmed that my ModifiedCellDMC gives the same outputs as original CellDMC.
martino2015.celldmc.o <- ModifiedCellDMC(Mvalues_train, pheno_train, cellfrac.m.martino2015, mc.cores=6)

# NOTE(milo): MAKE SURE THIS FOLDER IS CORRECT!
# folder = "./analysis/martino2015/Mvalues_nonallergic_vs_allergic_all/"
# folder = "./analysis/martino2015/Mvalues_nonallergic_vs_allergic_only_pbmc/"
# folder = "./analysis/martino2015/Mvalues_nonallergic_vs_allergic_with_eosino/"
folder = "./analysis/martino2015/Mvalues_nonallergic_vs_allergic_with_neutro/"
save(martino2015.celldmc.o, file=str_c(folder, "celldmc.o"))

# load(str_c(folder, "celldmc.o"))
summarizeDMCTs(martino2015.celldmc.o)

# coe.change = martino2015.celldmc.o$coe.change
# head(coe.change$Eosino)

# Write things to .csv so we can escape R.
write.csv(martino2015.celldmc.o$dmct, file=str_c(folder, "dmct.csv"))
write.csv(martino2015.celldmc.o$coe.change, file=str_c(folder, "coe_change.csv"))
write.csv(martino2015.celldmc.o$coe.control, file=str_c(folder, "coe_control.csv"))
write.csv(cellfrac.m.martino2015, file=str_c(folder, "cellfrac.csv"))
#==================================================
