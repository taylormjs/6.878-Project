# install.packages("matrixStats")
# install.packages("parallel")
# install.packages("stringr")
# install.packages("stats")
library("matrixStats")
library("parallel")
library("stringr")
library("stats")

# https://rdrr.io/github/sjczheng/EpiDISH/src/R/CellDMC.R
BulkCellDMC <- function(beta.m, pheno.v,
                        adjPMethod = "fdr", adjPThresh = 0.05, cov.mod = NULL, 
                        sort = FALSE, mc.cores = 1) {
  ### check input
  if (sum(is.na(pheno.v)) > 0) stop("No NA allowed in pheno.v!")
  if (ncol(beta.m) != length(pheno.v)) 
    stop("Number of columns of beta.m should equal to length of pheno.v!")
  
  ### check whether input is beta value matrix
  is.beta <- ((min(beta.m) >= 0) & (max(beta.m) <= 1))
  if (is.beta) {
    message("Detected a BETA matrix input")
  }
  ### guess factor input
  if (nlevels(factor(pheno.v)) == 2) {
    message("Binary phenotype detected. Predicted change will be 1 - 0.")
    pheno.v <- factor(pheno.v)
  }
  if (!is.factor(pheno.v) & !is.character(pheno.v)) 
    message("pheno.v is not factor or character. Treating as continuous variables.")

  # NOTE(milo): No interactions terms or frac.m terms now! Just pheno.v.
  print("Making design")
  design <- model.matrix(~ pheno.v)[, -1]
  if (!is.null(cov.mod)) design <- cbind(design, cov.mod[, -1])
  
  IntNames.v <- c("(Intercept)", "Pheno")

  print(sprintf("==> Fitting these %d covariates:", length(colnames(design))))
  print(colnames(design))

  AllNames.v = c("(Intercept)", colnames(design))

  allCoe.m <- do.call(rbind, mclapply(seq_len(nrow(beta.m)), function(i) {
    beta.v <- beta.m[i, ]
    Int.o <- lm(beta.v ~ ., data = data.frame(design))

    coef.matrix = summary(Int.o)$coefficients

    IntCoe.v <- unlist(apply(coef.matrix, 1, function(x) list(x)))
    names(IntCoe.v) <- NULL

    return(IntCoe.v)
  }, mc.preschedule = TRUE, mc.cores = mc.cores, mc.allow.recursive = TRUE))
  
  # Dimension of allCoe.m: num_cpg x (4 intercept + 4 pheno)
  print("==> Dimensions of allCoe.m:")
  print(dim(allCoe.m))
  
  # Get CONTROL coefficients for each cell type.
  print("==> GETTING CONTROL COEFFICIENTS ...")
  num_control_coef = 1 # Just the intercept.
  
  coe.control = lapply(seq_len(num_control_coef), function(j) {
    idx <- ((j - 1)*4 + 1):((j - 1)*4 + 4)
    tmp.m <- allCoe.m[, idx]
    tmp.m <- cbind(tmp.m, p.adjust(tmp.m[, 4], method = adjPMethod))
    if (is.beta) {
      tmp.m[which(tmp.m[,1] > 1),1] <- 1
      tmp.m[which(tmp.m[,1] < -1),1] <- -1
    }  ### if input is a beta values matrix, bound the estimated changes
    
    colnames(tmp.m) <- c("Estimate", "SE", "t", "p", "adjP")
    rownames(tmp.m) <- rownames(beta.m)
    return(data.frame(tmp.m))
  })
  
  # The control coefficients are put into a list of dataframes.
  # We want to concatenate all of them into a single matrix.
  coe.control.m = matrix(nrow=nrow(beta.m), ncol=0)
  for (i in 1:length(coe.control)) {
    var_name = AllNames.v[i]
    var_coe = coe.control[[i]]
    colnames(var_coe) = str_c(var_name, ".", colnames(var_coe))
    coe.control.m = cbind(coe.control.m, var_coe)
  }
  print("Concatenated coe.control ==> coe.control.m")
  print(dim(coe.control.m))
  print(colnames(coe.control.m))
  print("==> DONE")
  
  # Get the DISEASE coefficients (differential methylation) for each cell type.
  print("==> GETTING DISEASE COEFFICIENTS ...")
  offset_to_disease_coe = num_control_coef
  
  # NOTE(milo): Hacky, just replaced num cell types with 1.
  coe.ld <- lapply(seq_len(1) + offset_to_disease_coe, function(j) {
    idx <- ((j - 1)*4 + 1):((j - 1)*4 + 4)
    tmp.m <- allCoe.m[, idx]
    
    # Adjust the p-value of the 4th column.
    tmp.m <- cbind(tmp.m, p.adjust(tmp.m[, 4], method = adjPMethod))
    if (is.beta) {
      tmp.m[which(tmp.m[,1] > 1),1] <- 1
      tmp.m[which(tmp.m[,1] < -1),1] <- -1
    }  ### if input is a beta values matrix, bound the estimated changes
    
    colnames(tmp.m) <- c("Estimate", "SE", "t", "p", "adjP")
    rownames(tmp.m) <- rownames(beta.m)
    return(data.frame(tmp.m))
  })
  # names(coe.ld) <- colnames(frac.m)
  print("==> DONE")
  
  ### get dmct matrix
  print("==> MAKING DMCT MATRIX ...")
  dmct.m <- matrix(rep(0, nrow(beta.m)), ncol = 1)
  dmct.idx <- which(sapply(coe.ld, "[[", "adjP") < adjPThresh)
  dmct.m[dmct.idx] <- sign(sapply(coe.ld, "[[", "Estimate")[dmct.idx])
  dmc.v <- ifelse(rowAlls(dmct.m == 0), 0, 1)
  dmct.m <- cbind(dmc.v, dmct.m)
  colnames(dmct.m) <- c("DMC", "DUMMYCELL")
  rownames(dmct.m) <- rownames(beta.m)
  
  if(sort) coe.ld <- lapply(coe.ld, function(x) x[order(x$p),] )
  print("==> DONE")
  
  return(list(dmct = dmct.m, coe.change = coe.ld, coe.control = coe.control.m))
}