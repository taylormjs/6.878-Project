library(Biobase)    # Boxplot for selected GEO samples.
library(GEOquery)   # For downloading GEO2R data.

#================= COPIED FROM GEO2R ======================
# See: https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE59999
# This will use a cached version if it is found! Still has to unzip though (takes about a minute).
gset <- getGEO("GSE59999", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL13534", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

phenotypes = gset@phenoData@data
phenotypes$`gender:ch1`
phenotypes$`phenotype:ch1`
phenotypes$`challenge outcome:ch1`

beta.m = gset@assayData$exprs
colnames(beta.m)
rownames(beta.m)

#===================== Run EPIDISH =======================
# A reference-based function to infer the fractions of a priori known cell subtypes present in a sample
# representing a mixture of such cell-types. Inference proceeds via one of 3 methods
# (Robust PartialCorrelations-RPC, Cibersort-CBS, Constrained Projection-CP), as determined by the user

# Whole blood reference of 333 tsDHS-DMCs and 7 blood cell subtypes
data(centDHSbloodDMC.m)
epidish_out <- epidish(beta.m, centDHSbloodDMC.m, metod ='RPC')
cellfrac.m = epidish_out$estF
boxplot(cellfrac.m)

# Get allergic phenotypes: 0=nonallergic, 1=sensitized, 2=allergic.
# NOTE(milo): Just using the text phenotypes for now - celldmc seems OK with it?
pheno.v = phenotypes$`challenge outcome:ch1`

# Run celldmc to determine differentially methylated CpG locations.
celldmc.o <- CellDMC(beta.m, pheno.v, cellfrac.m)

head(celldmc.o$dmct)
head(celldmc.o$coe$CD4)

save(celldmc.o, file="martino2015_celldmc.o")
load("martino2015_celldmc.o")
