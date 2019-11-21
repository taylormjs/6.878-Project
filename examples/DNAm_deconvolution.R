library(Biobase)    # Boxplot for selected GEO samples.
library(GEOquery)   # For downloading GEO2R data.
library(EpiDISH)

#================= COPIED FROM GEO2R ======================
# See: https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE59999
# This will use a cached version if it is found! Still has to unzip though (takes about a minute).
gset <- getGEO("GSE59999", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL13534", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# View some of the relevant phenotype data.
phenotypes = gset@phenoData@data
phenotypes$`gender:ch1`
phenotypes$`phenotype:ch1`
phenotypes$`challenge outcome:ch1`

beta.m = gset@assayData$exprs
colnames(beta.m)
rownames(beta.m)

#===================== Run EPIDISH =======================
# Create binary phenotype labels for each patient.
# 0 = CONTROL = NONALLERGIC
# 1 = DISEASE = ALLERGIC
pheno.v = phenotypes$`challenge outcome:ch1`
pheno.v.binary = pheno.v
pheno.v.binary = pheno.v.binary[pheno.v.binary != "sensitized"] # Remove "sensitized" cases.
pheno.v.binary[pheno.v.binary == "allergic"] = 1
pheno.v.binary[pheno.v.binary == "nonallergic"] = 0
pheno.v.binary = as.integer(pheno.v.binary)
pheno.v.binary

# Run celldmc to determine differentially methylated CpG locations.
beta.m = beta.m[,pheno.v != "sensitized"]

# Whole blood reference of 333 tsDHS-DMCs and 7 blood cell subtypes
data(centDHSbloodDMC.m)
epidish_out <- epidish(beta.m, centDHSbloodDMC.m, method ='RPC')
cellfrac.m = epidish_out$estF
boxplot(cellfrac.m)

# Run celldmc.
# NOTE(milo): This will take a really long time to run! Load precomputed results if possible.
# celldmc.o <- CellDMC(beta.m, pheno.v.binary, cellfrac.m)

# save(celldmc.o, file="martino2015_celldmc.o")
load("../analysis/martino2015_celldmc.o")

# See: https://bioconductor.org/packages/release/bioc/manuals/EpiDISH/man/EpiDISH.pdf
# DMC = Differentially Methylated Cytosines
# Matrix gives wheter the input CpGs are DMCTs and DMCs. The first column tells whether
# a CpG is a DMC or not. If the CpG is called asDMC, the value will be 1, otherwise it
# is 0. The following columns give DMCTs for each cell-type.  If a CpG is a DMCT, the
# value will be 1 (hypermethylated for case compared to control) or -1  (hypomethylated
# for case compared to control).  Otherwise, the value is 0 (non-DMCT). The rows of this
# matrix are ordered as the same as that of the input beta.m.
dmct = celldmc.o$dmct
out.DMC = dmct[dmct[,1] == 1,]
message(sprintf("===== EPIDISH METHYLATION RESULTS ====="))
message(sprintf("EPIDISH found: %d DMCs", nrow(out.DMC)))
for (i in 2:ncol(out.DMC)) {
  cell_type_name = colnames(out.DMC)[i]
  num_dmc_celltype = nrow(dmct[dmct[,i] != 0,])
  message(sprintf("Cell type %s has %d DMCs", cell_type_name, num_dmc_celltype))
}

# This  list  contains  several  data frames,  which  correspond  to  each  cell-type
# infrac.m. Each data frame contains all CpGs in input beta.m.  All data frames contain
# estimated DNAm changes (Estimate), standard error (SE), estimated t statistics (t),
# raw P values (p), and multiple hypothesis corrected P values (adjP).

# NOTE(milo): It looks like with a p-value threshold of 0.05, we find exactly as many
# significant CpG locations as EPIDISH reported DMC locations, so they're probably
# using this threshold internally.
adj_p_thresh = 0.05

for (i in 1:length(out.coe)) {
  cell.name = names(out.coe)[[i]]
  cell.coe = out.coe[[i]] # NOTE(milo): Indexing by an integer requires [[]] double brackets!
  cell.signif_coe = cell_type_coe[cell.coe$adjP < adj_p_thresh,]
  message(sprintf("Cell type %s has %d CpG locations with adjusted p-value < %.2f",
                  cell.name, nrow(cell.signif_coe), adj_p_thresh))
}


testFunction = function(a) {
  return(2 * a)
}


