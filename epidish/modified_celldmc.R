# install.packages("matrixStats")
# install.packages("parallel")
# install.packages("stringr")
# install.packages("stats")
library("matrixStats")
library("parallel")
library("stringr")
library("stats")

# https://rdrr.io/github/sjczheng/EpiDISH/src/R/CellDMC.R
ModifiedCellDMC <- function(beta.m, pheno.v, frac.m, 
                    adjPMethod = "fdr", adjPThresh = 0.05, cov.mod = NULL, 
                    sort = FALSE, mc.cores = 1) {
  ### check input
  if (sum(is.na(pheno.v)) > 0) stop("No NA allowed in pheno.v!")
  if (ncol(beta.m) != length(pheno.v)) 
    stop("Number of columns of beta.m should equal to length of pheno.v!")
  if (ncol(beta.m) != nrow(frac.m)) 
    stop("Number of columns of beta.m should equal to number of rows of frac.m!")
  if (length(colnames(frac.m)) != ncol(frac.m)) 
    stop("Pls assign correct names of cell-type to frac.m")
  
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
  
  ### Fit model
  design <- model.matrix(~ frac.m + pheno.v:frac.m)[, -1]
  if (!is.null(cov.mod)) design <- cbind(design, cov.mod[, -1])
  IntNames.v <- str_c(colnames(frac.m), "Pheno")
  colnames(design)[(1 + ncol(frac.m)):(2*ncol(frac.m))] <- IntNames.v 
  
  # Design matrix has these columns:
  # [1] "frac.mB"      "frac.mNK"     "frac.mCD4T"   "frac.mCD8T"   "frac.mMono"   "frac.mNeutro"
  # [7] "frac.mEosino" "BPheno"       "NKPheno"      "CD4TPheno"    "CD8TPheno"    "MonoPheno"
  # The first 7 (just cell types) can be thought of as the "intercepts", or the contribution to the
  # beta value for each cell type in a control individual. The other 7 terms like BPheno are the
  # delta in beta values due to the disease or control phenotype. I think if we return the coefficients
  # for the frac.* terms, this is the beta values for the average control person.
  print(sprintf("==> Fitting these %d covariates:", length(colnames(design))))
  print(colnames(design))
  
  # We want to grab the INTERCEPT column, the frac.m.CELLTYPE columns
  # and the CELLTYPEPheno columns.
  AllNames.v = c("(Intercept)", colnames(design))

  # NOTE(milo): We do a linear regression for each CpG and then row-bind (rbind) it to the others.
  allCoe.m <- do.call(rbind, mclapply(seq_len(nrow(beta.m)), function(i) {
    beta.v <- beta.m[i, ]
    Int.o <- lm(beta.v ~ ., data = data.frame(design))
    
    # 4 columns: Estimate Std. Error    t value  Pr(>|t|)
    coef.matrix = matrix(nrow=length(AllNames.v), ncol=4)
    
    # This will mask out any colinear variables.
    valid_coeffs = summary(Int.o)$coefficients
    
    # Set the row and column names of our coeff matrix.
    colnames(coef.matrix) = colnames(valid_coeffs)
    rownames(coef.matrix) = AllNames.v
    
    # Fill in the non-NA valid coefficients.
    coef.matrix[rownames(valid_coeffs),] = valid_coeffs
    
    # Fill in the NA coeffs with zeros.
    coef.matrix[is.na(coef.matrix)] = 0
    
    ### get coe
    IntCoe.v <- unlist(apply(coef.matrix, 1, function(x) list(x)))
    names(IntCoe.v) <- NULL
    
    return(IntCoe.v)
  }, mc.preschedule = TRUE, mc.cores = mc.cores, mc.allow.recursive = TRUE))

  # Dimension of allCoe.m: (num_cpg x (1 + 2*num_cell_types)*4)
  # So for 7 cell types this would be 60 columns.
  print("==> Dimensions of allCoe.m:")
  print(dim(allCoe.m))
  
  # Get CONTROL coefficients for each cell type.
  # The first 4*(1+num_cell_types) things are the CONTROL coefficients.
  print("==> GETTING CONTROL COEFFICIENTS ...")
  num_cell_types = ncol(frac.m)
  num_control_coef = 1 + num_cell_types # Intercept + cell types.

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
  # NOTE(milo): Because we're keeping more coefficients in allCoe.m, we need to offset
  # by the number of cell types + 1 so that the old epidish stuff still works.
  print("==> GETTING DISEASE COEFFICIENTS ...")
  offset_to_disease_coe = num_control_coef

  coe.ld <- lapply(seq_len(ncol(frac.m)) + offset_to_disease_coe, function(j) {
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
  names(coe.ld) <- colnames(frac.m)
  print("==> DONE")
  
  ### get dmct matrix
  print("==> MAKING DMCT MATRIX ...")
  dmct.m <- matrix(rep(0, ncol(frac.m)*nrow(beta.m)), ncol = ncol(frac.m))
  dmct.idx <- which(sapply(coe.ld, "[[", "adjP") < adjPThresh)
  dmct.m[dmct.idx] <- sign(sapply(coe.ld, "[[", "Estimate")[dmct.idx])
  dmc.v <- ifelse(rowAlls(dmct.m == 0), 0, 1)
  dmct.m <- cbind(dmc.v, dmct.m)
  colnames(dmct.m) <- c("DMC", colnames(frac.m))
  rownames(dmct.m) <- rownames(beta.m)
  
  if(sort) coe.ld <- lapply(coe.ld, function(x) x[order(x$p),] )
  print("==> DONE")
  
  return(list(dmct = dmct.m, coe.change = coe.ld, coe.control = coe.control.m))
}