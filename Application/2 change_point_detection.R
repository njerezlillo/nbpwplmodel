library(dplyr)
library(survival)
library(maxLik)
library(xtable)
library(survminer)
source("nbpwpl.R")

# Dataset -----------------------------------------------------------------

load("./Application/preliminar.RData")

t0 <- 1
n <- sum(df$status)

# k = 1 -------------------------------------------------------------------

profile <- function(y) profile_loglik_cens_pwpowerlaw(p = c(t0, y), df)

fit_p1 <-
  maxSANN(
    profile,
    start = c(3.5),
    constraints = list(
      ineqA = matrix(c(1, -1), ncol = 1, byrow = T),
      ineqB = c(-1.2, 8)
    )
  )

p_1 <- c(t0, fit_p1$estimate)

# k = 2 -------------------------------------------------------------------

profile <- function(y) profile_loglik_cens_pwpowerlaw(p = c(t0, y), df)

fit_p2 <-
  maxSANN(
    profile,
    start = c(3.4, 4.5),
    constraints = list(
      ineqA = matrix(c(1, 0, -1, 1, 0, -1), ncol = 2, byrow = T),
      ineqB = c(-1.2, -1, 8)
    )
  )

p_2 <- c(t0, fit_p2$estimate)

# k = 3 -------------------------------------------------------------------

profile <- function(y) profile_loglik_cens_pwpowerlaw(p = c(t0, y), df)

fit_p3 <-
  maxSANN(
    profile,
    start = c(2.7, 3.4, 4.5),
    constraints = list(
      ineqA = matrix(c(1, 0, 0, -1, 1, 0, 0, -1, 1, 0, 0, -1), ncol = 3, byrow = T),
      ineqB = c(-1.2, -0.5, -0.5, 8)
    )
  )

p_3 <- c(t0, fit_p3$estimate)

# Seleccion ---------------------------------------------------------------

alpha_1 <- mle_cens_pwpowerlaw(df, p_1)
alpha_2 <- mle_cens_pwpowerlaw(df, p_2)
alpha_3 <- mle_cens_pwpowerlaw(df, p_3)

# AIC ---------------------------------------------------------------------

AIC_1 <- 2 * (length(c(p_1, alpha_1)) - 1) - 2 * loglik_cens_pwpowerlaw(alpha_1, df, p_1)
AIC_2 <- 2 * (length(c(p_2, alpha_2)) - 1) - 2 * loglik_cens_pwpowerlaw(alpha_2, df, p_2)
AIC_3 <- 2 * (length(c(p_3, alpha_3)) - 1) - 2 * loglik_cens_pwpowerlaw(alpha_3, df, p_3)

# BIC ---------------------------------------------------------------------

BIC_1 <- log(n) * (length(c(p_1, alpha_1)) - 1) - 2 * loglik_cens_pwpowerlaw(alpha_1, df, p_1)
BIC_2 <- log(n) * (length(c(p_2, alpha_2)) - 1) - 2 * loglik_cens_pwpowerlaw(alpha_2, df, p_2)
BIC_3 <- log(n) * (length(c(p_3, alpha_3)) - 1) - 2 * loglik_cens_pwpowerlaw(alpha_3, df, p_3)

# Table -------------------------------------------------------------------

print(xtable(data.frame(
  k = 1:3,
  AIC = c(AIC_1, AIC_2, AIC_3),
  BIC = c(BIC_1, BIC_2, BIC_3)
)), include.rownames = F)

p <- p_2
save(df, p, q, file = "./Application/preliminar.RData")
save.image("./Application/run_script_2.RData")

