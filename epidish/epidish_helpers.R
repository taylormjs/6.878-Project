# Create binary phenotype labels for each patient.
# 0 = CONTROL = NONALLERGIC
# 1 = DISEASE = ALLERGIC
makeBinaryPhenotypesMartino2015 = function(gset) {
  phenotypes = gset@phenoData@data
  pheno.v = phenotypes$`challenge outcome:ch1`
  pheno.v.binary = pheno.v
  pheno.v.binary = pheno.v.binary[pheno.v.binary != "sensitized"] # Remove "sensitized" cases.
  pheno.v.binary[pheno.v.binary == "allergic"] = 1
  pheno.v.binary[pheno.v.binary == "nonallergic"] = 0
  pheno.v.binary = as.integer(pheno.v.binary)
  
  return(pheno.v.binary)
}


# Create binary phenotype labels for each patient.
# 0 = CONTROL = NONALLERGIC
# 1 = DISEASE = ALLERGIC
makeBinaryPhenotypesMartino2018 = function(gset) {
  phenotypes = gset@phenoData@data
  phenotypes = gset.martino2018@phenoData@data
  pheno.v = phenotypes$`allergy status:ch1`
  pheno.v.binary = pheno.v
  pheno.v.binary = pheno.v.binary[pheno.v.binary != "resolved"] # Remove "resolved" cases.
  pheno.v.binary[pheno.v.binary == "allergic"] = 1
  pheno.v.binary[pheno.v.binary == "control"] = 0
  pheno.v.binary = as.integer(pheno.v.binary)
 
  return(pheno.v.binary) 
}


getBetaMatrixMartino2015 = function(gset) {
  phenotypes = gset@phenoData@data
  pheno.v = phenotypes$`challenge outcome:ch1`
  beta.m = gset@assayData$exprs
  beta.m = beta.m[,pheno.v != "sensitized"]

  return(beta.m)
}


getBetaMatrixMartino2018 = function(gset) {
  phenotypes = gset@phenoData@data
  phenotypes = gset.martino2018@phenoData@data
  pheno.v = phenotypes$`allergy status:ch1`
  beta.m = gset@assayData$exprs
  beta.m = beta.m[,pheno.v != "resolved"]

  return(beta.m)
}


# Get EPIDISH estimated cell fractions from a beta matrix.
getEpidishCellFrac = function(beta.m) {
  data(centDHSbloodDMC.m)
  epidish_out <- epidish(beta.m, centDHSbloodDMC.m, method ='RPC')
  cellfrac.m = epidish_out$estF

  return(cellfrac.m)
}

# See: https://bioconductor.org/packages/release/bioc/manuals/EpiDISH/man/EpiDISH.pdf
# DMC = Differentially Methylated Cytosines
# Matrix gives wheter the input CpGs are DMCTs and DMCs. The first column tells whether
# a CpG is a DMC or not. If the CpG is called asDMC, the value will be 1, otherwise it
# is 0. The following columns give DMCTs for each cell-type.  If a CpG is a DMCT, the
# value will be 1 (hypermethylated for case compared to control) or -1  (hypomethylated
# for case compared to control).  Otherwise, the value is 0 (non-DMCT). The rows of this
# matrix are ordered as the same as that of the input beta.m.
summarizeDMCTs = function(celldmc.o) {
  dmct = celldmc.o$dmct
  out.DMC = dmct[dmct[,1] == 1,]
  message(sprintf("===== EPIDISH METHYLATION RESULTS ====="))
  message(sprintf("EPIDISH found: %d DMCs", nrow(out.DMC)))
  for (i in 2:ncol(out.DMC)) {
    cell_type_name = colnames(out.DMC)[i]
    num_dmc_celltype = sum(dmct[,i] != 0)
    
    # NOTE(milo): This line below seems unreliable for the case where there is ONE row
    # or ZERO rows that meet the condition. R will turn a matrix with one row into a
    # different data type, and it will turn a matrix with just a header into something
    # weird too.
    # num_dmc_celltype = nrow(data.matrix(dmct[dmct[,i] != 0,]))
    message(sprintf("Cell type %s has %d DMCs", cell_type_name, num_dmc_celltype))
  }
}
