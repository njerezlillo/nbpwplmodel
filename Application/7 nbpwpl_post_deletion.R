library(stringr)
library(dplyr)
library(maxLik)
library(survival)
library(survminer)
library(xtable)
source("nbpwpl.R")

# Dataset -----------------------------------------------------------------

load("./Application/preliminar3.RData")
df <- arrange(df, time)

groups <- vector(length = n)
groups[df$age == 0 & df$disease == 0] <-  "Group I"
groups[df$age == 0 & df$disease == 1] <-  "Group II"
groups[df$age == 1 & df$disease == 1] <-  "Group III"

# Influence ---------------------------------------------------------------

I_cw <- Influence(df, coef_final, p, H, "cw")
I_rs <- Influence(df, coef_final, p, H, "rs")

# cw: deleting influential observations -----------------------------------

## all
index_cw <- which(I_cw$Bi > mean(I_cw$Bi) + 2 * sd(I_cw$Bi))
df_cw <- df[-index_cw, ]

target_cw <- function(w) {
  arg <- rep(0, 16)
  arg[!zero_position] <- w
  -loglik_nbpwpl_reg_model(arg, df_cw, p)
}

fit_cw <- optim(fn = target_cw, par = rep(0, 9),
                method = "CG", hessian = T)

est <- fit_cw$par
se <- sqrt(diag(solve(fit_cw$hessian)))

results <- data.frame(est = round(est, 2),
                      se = round(se, 2),
                      pval = round(2 * pnorm(abs(est/se), lower.tail = F), 2))

print(xtable(t(results)), include.rownames = F)

## max
df_cw <- df[-which.max(I_cw$Bi), ]

target_cw <- function(w) {
  arg <- rep(0, 16)
  arg[!zero_position] <- w
  -loglik_nbpwpl_reg_model(arg, df_cw, p)
}

fit_cw <- optim(fn = target_cw, par = rep(0, 9),
                method = "CG", hessian = T)

est <- fit_cw$par
se <- sqrt(diag(solve(fit_cw$hessian)))

results <- data.frame(est = round(est, 2),
                      se = round(se, 2),
                      pval = round(2 * pnorm(abs(est/se), lower.tail = F), 2))

print(xtable(t(results)), include.rownames = F)

## group I
df_cw <- df[-na.omit(index_cw[groups[index_cw] == "Group I"]), ]

target_cw <- function(w) {
  arg <- rep(0, 16)
  arg[!zero_position] <- w
  -loglik_nbpwpl_reg_model(arg, df_cw, p)
}

fit_cw <- optim(fn = target_cw, par = rep(0, 9),
                method = "CG", hessian = T)

est <- fit_cw$par
se <- sqrt(diag(solve(fit_cw$hessian)))

results <- data.frame(est = round(est, 2),
                      se = round(se, 2),
                      pval = round(2 * pnorm(abs(est/se), lower.tail = F), 2))

print(xtable(t(results)), include.rownames = F)

## group II
df_cw <- df[-na.omit(index_cw[groups[index_cw] == "Group II"]), ]

target_cw <- function(w) {
  arg <- rep(0, 16)
  arg[!zero_position] <- w
  -loglik_nbpwpl_reg_model(arg, df_cw, p)
}

fit_cw <- optim(fn = target_cw, par = rep(0, 9),
                method = "CG", hessian = T)

est <- fit_cw$par
se <- sqrt(diag(solve(fit_cw$hessian)))

results <- data.frame(est = round(est, 2),
                      se = round(se, 2),
                      pval = round(2 * pnorm(abs(est/se), lower.tail = F), 2))

print(xtable(t(results)), include.rownames = F)

## group III
df_cw <- df[-na.omit(index_cw[groups[index_cw] == "Group III"]), ]

target_cw <- function(w) {
  arg <- rep(0, 16)
  arg[!zero_position] <- w
  -loglik_nbpwpl_reg_model(arg, df_cw, p)
}

fit_cw <- optim(fn = target_cw, par = rep(0, 9),
                method = "CG", hessian = T)

est <- fit_cw$par
se <- sqrt(diag(solve(fit_cw$hessian)))

results <- data.frame(est = round(est, 2),
                      se = round(se, 2),
                      pval = round(2 * pnorm(abs(est/se), lower.tail = F), 2))

print(xtable(t(results)), include.rownames = F)

# rs: deleting influential observations -----------------------------------

## all
index_rs <- which(I_rs$Bi > mean(I_rs$Bi) + 2 * sd(I_rs$Bi))
df_rs <- df[-index_rs, ]

target_rs <- function(w) {
  arg <- rep(0, 16)
  arg[!zero_position] <- w
  -loglik_nbpwpl_reg_model(arg, df_rs, p)
}

fit_rs <- optim(fn = target_rs, par = rep(0, 9),
                method = "CG", hessian = T)

est <- fit_rs$par
se <- sqrt(diag(solve(fit_rs$hessian)))

results <- data.frame(est = round(est, 2),
                      se = round(se, 2),
                      pval = round(2 * pnorm(abs(est/se), lower.tail = F), 2))

print(xtable(t(results)), include.rownames = F)

## max
df_rs <- df[-which.max(I_rs$Bi), ]

target_rs <- function(w) {
  arg <- rep(0, 16)
  arg[!zero_position] <- w
  -loglik_nbpwpl_reg_model(arg, df_rs, p)
}

fit_rs <- optim(fn = target_rs, par = rep(0, 9),
                method = "CG", hessian = T)

est <- fit_rs$par
se <- sqrt(diag(solve(fit_rs$hessian)))

results <- data.frame(est = round(est, 2),
                      se = round(se, 2),
                      pval = round(2 * pnorm(abs(est/se), lower.tail = F), 2))

print(xtable(t(results)), include.rownames = F)

## group I
df_rs <- df[-na.omit(index_rs[groups[index_rs] == "Group I"]), ]

target_rs <- function(w) {
  arg <- rep(0, 16)
  arg[!zero_position] <- w
  -loglik_nbpwpl_reg_model(arg, df_rs, p)
}

fit_rs <- optim(fn = target_rs, par = rep(0, 9),
                method = "CG", hessian = T)

est <- fit_rs$par
se <- sqrt(diag(solve(fit_rs$hessian)))

results <- data.frame(est = round(est, 2),
                      se = round(se, 2),
                      pval = round(2 * pnorm(abs(est/se), lower.tail = F), 2))

print(xtable(t(results)), include.rownames = F)

## group II
df_rs <- df[-na.omit(index_rs[groups[index_rs] == "Group II"]), ]

target_rs <- function(w) {
  arg <- rep(0, 16)
  arg[!zero_position] <- w
  -loglik_nbpwpl_reg_model(arg, df_rs, p)
}

fit_rs <- optim(fn = target_rs, par = rep(0, 9),
                method = "CG", hessian = T)

est <- fit_rs$par
se <- sqrt(diag(solve(fit_rs$hessian)))

results <- data.frame(est = round(est, 2),
                      se = round(se, 2),
                      pval = round(2 * pnorm(abs(est/se), lower.tail = F), 2))

print(xtable(t(results)), include.rownames = F)


# AAA ---------------------------------------------------------------------

yyy = na.omit(index_cw[groups[index_cw] == "Group II"])

fit_km_ori <- survfit(Surv(time, status) ~ age, data = df)
fit_km_del <- survfit(Surv(time, status) ~ age, data = df[-yyy,])

plot(fit_km_ori)
abline(v = p)

plot(fit_km_del)
abline(v = p)


