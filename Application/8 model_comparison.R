library(stringr)
library(dplyr)
library(maxLik)
library(xtable)
library(survival)
library(survminer)
source("nbpwexp.R")
source("nbcfwei.R")
source("nbpwpexp.R")

# Dataset -----------------------------------------------------------------

load("./Application/preliminar2.RData")
# load("./Application/run_script_8.RData")

# Cure fraction weibull model ---------------------------------------------

l_wei <- function(w) loglik_nbpwwei(w, df, q)

fit_wei <- 
  maxCG(
    l_wei,
    start = c(2, 2, 1, 1, 1, 1),
    constraints = list(
      ineqA = cbind(diag(c(1, 1)), matrix(0, ncol = 4, nrow = 2)),
      ineqB = c(-0.01, -0.01)
    )
  )
  
aic_wei <- 2 * 6 - 2 * fit_wei$maximum
bic_wei <- log(nrow(df)) * 6 - 2 * fit_wei$maximum

# PW Exponencial model ----------------------------------------------------

profile <- function(y) profile_loglik_cens_pwexp(p = c(1, y), df)

fit_p2 <-
  maxSANN(
    profile,
    start = c(3.4, 4.5),
    constraints = list(
      ineqA = matrix(c(1, 0, -1, 1, 0, -1), ncol = 2, byrow = T),
      ineqB = c(-1.2, -1, 8)
    )
  )

p_exp <- c(1, fit_p2$estimate)

l_pwexp <- function(w) {
  sum(apply(df, 1, function(z) loglik_nbpwexp(w, z, p_exp, q)))
}

fit_pwexp <- 
  maxCG(
    l_pwexp,
    start = c(1, 1, 1, 0, 0, 0, 0),
    constraints = list(
      ineqA = cbind(diag(c(1, 1, 1)), matrix(0, ncol = 4, nrow = 3)),
      ineqB = c(-0.01, -0.01, -0.01)
    )
  )

aic_pwexp <- 2 * 7 - 2 * fit_pwexp$maximum
bic_pwexp <- log(nrow(df)) * 7 - 2 * fit_pwexp$maximum

# PW Power-Exponencial model ----------------------------------------------

l_pwpexp <- function(w) {
  sum(apply(df, 1, function(z) loglik_nbpwpexp(w, z, p_exp, q)))
}

fit_pwpexp <- 
  maxCG(
    l_pwpexp,
    start = c(1, 1, 1, 1, 0, 0, 0, 0),
    constraints = list(
      ineqA = cbind(diag(c(1, 1, 1, 1)), matrix(0, ncol = 4, nrow = 4)),
      ineqB = c(-0.01, -0.01, -0.01, -0.01)
    )
  )

aic_pwpexp <- 2 * 8 - 2 * fit_pwpexp$maximum
bic_pwpexp <- log(nrow(df)) * 8 - 2 * fit_pwpexp$maximum

# save.image("./Application/run_script_8.RData")

# Figures -----------------------------------------------------------------

fit_km <-
  survfit(Surv(time, status) ~ disease + age + gender, data = df %>% mutate(time = time - 1))

###

theta_wei <- fit_wei$estimate

s_1 <- Vectorize(function(x) snbpwwei(x, q, plogis(theta_wei[3]), theta_wei[1:2]), "x")
s_2 <- Vectorize(function(x) snbpwwei(x, q, plogis(sum(theta_wei[c(3, 6)])), theta_wei[1:2]), "x")
s_3 <- Vectorize(function(x) snbpwwei(x, q, plogis(sum(theta_wei[c(3, 4, 6)])), theta_wei[1:2]), "x")
s_4 <- Vectorize(function(x) snbpwwei(x, q, plogis(sum(theta_wei[3:6])), theta_wei[1:2]), "x")

ggsurvplot(fit_km, conf.int = F,
           ggtheme = theme_bw(), size = 0.5,
           legend.title = "Disease", censor = F)$plot +
  labs(x = "Duration (years)", y = "Survival Rate") + 
  geom_function(fun = s_1) +
  geom_function(fun = s_2) +
  geom_function(fun = s_3) +
  geom_function(fun = s_4) +
  theme(legend.position = "top")

###

theta_pwexp <- fit_pwexp$estimate

s_1 <- Vectorize(function(x) snbpwexp(x + 1, q, plogis(theta_pwexp[4]), theta_pwexp[1:3], p_exp), "x")
s_2 <- Vectorize(function(x) snbpwexp(x + 1, q, plogis(sum(theta_pwexp[c(4, 7)])), theta_pwexp[1:3], p_exp), "x")
s_3 <- Vectorize(function(x) snbpwexp(x + 1, q, plogis(sum(theta_pwexp[c(4, 5, 7)])), theta_pwexp[1:3], p_exp), "x")
s_4 <- Vectorize(function(x) snbpwexp(x + 1, q, plogis(sum(theta_pwexp[4:7])), theta_pwexp[1:3], p_exp), "x")

ggsurvplot(fit_km, conf.int = F,
           ggtheme = theme_bw(), size = 0.5,
           legend.title = "Disease", censor = F)$plot +
  labs(x = "Duration (years)", y = "Survival Rate") + 
  geom_function(fun = s_1) +
  geom_function(fun = s_2) +
  geom_function(fun = s_3) +
  geom_function(fun = s_4) +
  theme(legend.position = "top")

###

theta_pwpexp <- fit_pwpexp$estimate

s_1 <- Vectorize(function(x) snbpwpexp(x + 1, q, plogis(theta_pwpexp[5]), theta_pwpexp[1:3], theta_pwpexp[4], p_exp), "x")
s_2 <- Vectorize(function(x) snbpwpexp(x + 1, q, plogis(sum(theta_pwpexp[c(5, 8)])), theta_pwpexp[1:3], theta_pwpexp[4], p_exp), "x")
s_3 <- Vectorize(function(x) snbpwpexp(x + 1, q, plogis(sum(theta_pwpexp[c(5, 6, 8)])), theta_pwpexp[1:3], theta_pwpexp[4], p_exp), "x")
s_4 <- Vectorize(function(x) snbpwpexp(x + 1, q, plogis(sum(theta_pwpexp[5:8])), theta_pwpexp[1:3], theta_pwpexp[4], p_exp), "x")

ggsurvplot(fit_km, conf.int = F,
           ggtheme = theme_bw(), size = 0.5,
           legend.title = "Disease", censor = F)$plot +
  labs(x = "Duration (years)", y = "Survival Rate") + 
  geom_function(fun = s_1) +
  geom_function(fun = s_2) +
  geom_function(fun = s_3) +
  geom_function(fun = s_4) +
  theme(legend.position = "top")
