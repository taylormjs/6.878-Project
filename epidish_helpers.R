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
