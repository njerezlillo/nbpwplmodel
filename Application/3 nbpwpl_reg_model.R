library(dplyr)
library(survival)
library(maxLik)
library(xtable)
library(survminer)
source("nbpwpl.R")

# Dataset -----------------------------------------------------------------

load("./Application/preliminar.RData")

# Estimation --------------------------------------------------------------

target <- function(w) -loglik_nbpwpl_reg_model(w, df, p)

fit <- optim(fn = target, par = rep(0, 16), method = "CG")
fit$par

theta_matrix(fit$par, df, p) # puntual estimation

Theta <- alpha_matrix(fit$par, df, p)

s_A <- Vectorize(function(x) snbpwpl(x + 1, q, Theta[1, 4], Theta[1, 1:3], p))
s_B <- Vectorize(function(x) snbpwpl(x + 1, q, Theta[2, 4], Theta[2, 1:3], p))
s_C <- Vectorize(function(x) snbpwpl(x + 1, q, Theta[3, 4], Theta[3, 1:3], p))
s_D <- Vectorize(function(x) snbpwpl(x + 1, q, Theta[4, 4], Theta[4, 1:3], p))
s_E <- Vectorize(function(x) snbpwpl(x + 1, q, Theta[5, 4], Theta[5, 1:3], p))
s_F <- Vectorize(function(x) snbpwpl(x + 1, q, Theta[6, 4], Theta[6, 1:3], p))
s_G <- Vectorize(function(x) snbpwpl(x + 1, q, Theta[7, 4], Theta[7, 1:3], p))
s_H <- Vectorize(function(x) snbpwpl(x + 1, q, Theta[8, 4], Theta[8, 1:3], p))

fit_km <- 
  survfit(Surv(time, status) ~ gender + disease + age,
          data = df[-c(39, 40, 50, 172),] %>% 
            mutate(time = time - 1))

plot(fit_km, xlab = "Duration (years)", ylab = "Survival Rate")
curve(s_A, 0, 16, col = "steelblue", add = T)
# curve(s_B, 0, 18, col = "red", add = T) #(NO EXISTE)
curve(s_C, 0, 16, col = "steelblue", add = T)
curve(s_D, 0, 18, col = "darkgreen", add = T)
# curve(s_E, 0, 16, col = "", add = T)  #(4 OBS)
curve(s_F, 0, 16, col = "red", add = T)
curve(s_G, 0, 16, col = "darkgreen", add = T)
curve(s_H, 0, 16, col = "red", add = T)

save.image("./Application/run_script_3.RData")
