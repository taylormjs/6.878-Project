source("epidish/epidish_helpers.R")

# m-value = log_2 ( beta / (1 - beta) )
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-587
convertBetaToMValue = function(beta) {
  M = log2(beta / (1 - beta))
  return(M)
}

# Process Martino 2018.
load("./analysis/gset.martino2018.Rda")
beta.2018 = getBetaMatrixMartino2018(gset.martino2018)
Mvalues2018 = log2(beta.2018 / (1 - beta.2018))
save(Mvalues2018, file="./analysis/Mvalues2018.Rda")

# Process Martino 2015.
load("./analysis/gset.martino2015.Rda")
beta.2015 = getBetaMatrixNonallergicVsAllergicMartino2015(gset.martino2015)
Mvalues2015 = log2(beta.2015 / (1 - beta.2015))
save(Mvalues2015, file="./analysis/Mvalues2015.Rda")
