library(stringr)
library(dplyr)
library(maxLik)
library(xtable)
library(survival)
library(survminer)
library(parallel)
library(doParallel)
source("nbpwpl.R")

# Dataset -----------------------------------------------------------------

load("./Application/preliminar.RData")

# Lasso -------------------------------------------------------------------

loglik_nbpwpl_reg_model_lasso <- function(w, lambda) {
  -loglik_nbpwpl_reg_model(w, df, p) +
    lambda * sum(abs(w[-c(1, 5, 9, 13)]))
}

FF <- function (z) {
  fit <- optim(par = rep(0, 16),
               fn = loglik_nbpwpl_reg_model_lasso,
               method = "CG",
               lambda = z)
  
  est <- round(fit$par, 3)
  est[which(!(abs(est) > 0.01))] <- 0
  n_par <- sum(est != 0)
  n_beta <- ncol(df[, -(1:2)])
  n_link <- length(p) 
  n <- nrow(df)
  
  # mean square error
  Theta <- alpha_matrix(est, df, p)
  
  fit_km <- 
    survfit(Surv(time, status) ~ age + gender + disease,
            data = df)
  
  xxx <- summary(survfit(Surv(time, status) ~ age + gender + disease,
            data = df), times=df$time)
  
  mse_index <- function (t) {
    DF <- fit_km$time[t]
    SURV <- fit_km$surv[t]
    group <- rep(names(fit_km$strata), fit_km$strata)[t]
    X <- c(1, as.numeric(str_extract_all(group, "\\d+")[[1]]))
    
    Alpha <- vector(length = n_beta)
    
    for (k in 0:(n_link - 2)) {
      Alpha[k + 1] <- 
        1.001 + exp(sum(est[which(X == 1) + n_beta * k]))
    }
    
    Alpha[n_link] <- 
      3.001 + exp(sum(est[which(X == 1) + n_beta * (n_link - 1)]))
    Alpha[n_link + 1] <- 
      sum(est[which(X == 1) + n_beta * n_link])
    
    (SURV - snbpwpl(DF, q, plogis(Alpha[n_link + 1]), Alpha[1:n_link], p))^2
  }
  
  MSE <- sum(vapply(1:length(fit_km$time), mse_index, numeric(1)))
  
  # output
  c(lambda = z, 
    MSE = MSE,
    AIC = 2 * (n_par + 2) - 2 * loglik_nbpwpl_reg_model(est, df, p), 
    BIC = log(n) * (n_par + 2) - 2 * loglik_nbpwpl_reg_model(est, df, p),
    npar = as.integer(n_par))
}

number_of_cores <- detectCores()
clusters <- parallel::makeCluster(number_of_cores)
doParallel::registerDoParallel(clusters)
registerDoParallel(cores = number_of_cores)

ptm <- proc.time() #start

estimates <- 
  foreach(zz = 0:20, .multicombine = TRUE) %dopar% FF(zz) 

proc.time() - ptm #final

estimates

lambdas <- data.frame(t(sapply(estimates, c)))

# Figure ------------------------------------------------------------------

# load("./Application/run_script_4.RData")

scaleFUN <- function(x) sprintf("%.2f", x)

fig1 <- ggplot(lambdas) + aes(x = lambda, y = BIC) +
  labs(x = expression(lambda), y = "Bayesian Information Criterion") +
  geom_segment(aes(x = 5, y = 2500, xend = 5, yend = 2580), 
               col = "red", lty = 2) +
  geom_line() + geom_point() + theme_bw()
  
fig2 <- ggplot(lambdas) + aes(x = lambda, y = npar) +
  labs(x = expression(lambda), y = "Number of parameters") +
  geom_segment(aes(x = 5, y = 0, xend = 5, yend = 16), 
               col = "red", lty = 2) +
  geom_line() + geom_point() + theme_bw() + 
  scale_y_continuous(breaks = seq(0, 16, by = 2),
                     labels = scaleFUN) 

plot_final <- gridExtra::grid.arrange(fig1, fig2, ncol = 2)
ggsave(plot = plot_final, "./Application/bic_lambda.pdf", height = 4, width = 12)

# Out ---------------------------------------------------------------------

fit_lasso <- optim(
  par = rep(0, 16),
  fn = loglik_nbpwpl_reg_model_lasso,
  method = "CG",
  lambda = 5
)

est <- round(fit_lasso$par, 3)
est[which(!(abs(est) > 0.01))] <- 0
zero_position <- (est == 0)

save(df, p, q, zero_position, file = "./Application/preliminar2.RData")
save.image("./Application/run_script_4.RData")
