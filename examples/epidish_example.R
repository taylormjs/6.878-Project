library(EpiDISH)    # For DNAm analysis.

#================= EXAMPLE FROM EPIDISH DOCUMENTATION ================
# Guide: https://bioconductor.org/packages/release/bioc/vignettes/EpiDISH/inst/doc/EpiDISH.html
# Matrix of reference "centroids" that define typical methylation levels for each cell type.
# Each row is a CpG location (i.e cg16509569) and columns are "Epi", "Fib", "IC". I think
# these are the different cell type names?
data(centEpiFibIC.m)

# Beta values - each row is CpG location, each col is a patient.
data(DummyBeta.m)

# A reference-based function to infer the fractions of a priori known cell subtypes present in a sample
# representing a mixture of such cell-types. Inference proceeds via one of 3 methods
# (Robust PartialCorrelations-RPC, Cibersort-CBS, Constrained Projection-CP), as determined by the user
epidish_out <- epidish(DummyBeta.m, centEpiFibIC.m, method ='RPC')

# This is a matrix with cell fractions for each patient sample.
# Each patient is a row and each col is a cell type (rows should sum to 1.0).
cellfrac.m = epidish_out$estF

# Show the mean cell fractions across all patients, along with stdev.
boxplot(cellfrac.m)

pheno.v = integer(10) # Vector of phenotype values.
pheno.v[1:5] = 1 # Make the first 5 patients have phenotype "1".

celldmc.o <- CellDMC(DummyBeta.m, pheno.v, cellfrac.m)
head(celldmc.o$dmct)
head(celldmc.o$coe$Epi)
