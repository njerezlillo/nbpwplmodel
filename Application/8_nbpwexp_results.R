library(stringr)
library(dplyr)
library(maxLik)
library(xtable)
library(survival)
library(survminer)
source("nbpwexp.R")

# Dataset -----------------------------------------------------------------

load("./Application/preliminar2.RData")

# Final estimation --------------------------------------------------------

target <- function(w) {
  arg <- rep(0, 16)
  arg[!zero_position] <- w
  -loglik_nbpwexp_reg_model(arg, df, p)
}

fit <- optim(fn = target, par = rep(0, 9), method = "CG", hessian = T)
round(2 * 9 + 2 * fit$value, 2) # AIC