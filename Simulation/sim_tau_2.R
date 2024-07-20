library(maxLik)
library(xtable)
library(numDeriv)
library(parallel)
library(doParallel)
source("nbpwpl.R")

Table_tau <- function(n) {
  set.seed(2024)
  Rep <- 300
  
  Beta_real <- c(-0.74,  0.00,  0.00,
                 0.67, 0.00, 0.00,
                 -1.81, 0.00, 0.00,
                 -5.71, -0.55, 4.81)
  p <- c(1, 2, 5)
  m <- 2
  
  PMLE <- matrix(ncol = Rep, nrow = length(p) - 1)
  
  for (i in 1:Rep) {
    fit_tau <- 1
    while (all(fit_tau == 1)) {
      x <- rnbpwplx(n, Beta_real, p = p, m = m)
      
      x_min <- quantile(sort(x$time), 0.05)
      x_max <- quantile(sort(x$time), 0.8)
      d <-  2.6
      A_2 <- matrix(c(1, -1, 0, 0, 1, -1), ncol = 2)
      d_2 <- c(-x_min, -d, x_max)
      
      profile <- function(y) profile_loglik_pwpowerlaw(p = c(p[1], y), x$time)
      
      fit_tau <- tryCatch(
        maxSANN(
          profile,
          start = c(x_min + 0.05, x_min + d + 0.1),
          constraints = list(ineqA = A_2, ineqB = d_2),
          control = list(iterlim = 500)
        )$estimate, error = function (t) return(1)
      )
      
      fit_tau
    }
    
    PMLE[, i] <- fit_tau
  }
  
  BIAS <- apply(PMLE - p[-1], 1, mean)
  SE <- apply(PMLE - p[-1], 1, function(x) sd(x))
  MSE <- apply(PMLE - p[-1], 1, function(x) mean(x^2))
  
  MC_BS <- SE / sqrt(Rep)
  MC_MS <- 
    sqrt(apply((PMLE - p[-1])^2 - MSE, 1,
               function(x) sum(x^2))/ (Rep * (Rep - 1)))
  
  list(Estimate = data.frame("n" = paste(n), BIAS, MSE),
       MC_SD = data.frame("n" = paste(n), MC_BS, MC_MS))
}

# Scenario 2 --------------------------------------------------------------

number_of_cores <- 8
clusters <- parallel::makeCluster(number_of_cores)
doParallel::registerDoParallel(clusters)
registerDoParallel(cores = 8)

ptm <- proc.time() #start

estimates <- 
  foreach(n_grid = c(750, 1000, 1250), .multicombine = TRUE, 
          .packages = c("maxLik", "numDeriv")) %dopar% 
  Table_tau(n_grid) 

proc.time() - ptm #final

estimates_Estimate <- list()
estimates_MC_SD <- list()

for (i in seq_along(estimates)) {
  estimate_df <- estimates[[i]]$Estimate
  mc_sd_df <- estimates[[i]]$MC_SD
  
  estimates_Estimate[[i]] <- estimate_df
  estimates_MC_SD[[i]] <- mc_sd_df
}

temp_1 <- do.call(rbind, estimates_Estimate)[
  c(seq(1, 6, by = 2), seq(1, 6, by = 2) + 1),]
temp_2 <- do.call(rbind, estimates_MC_SD)[
  c(seq(1, 6, by = 2), seq(1, 6, by = 2) + 1),]

print.xtable(xtable(temp_1, digits = 3), include.rownames = F)

print.xtable(xtable(temp_2, digits = 3), include.rownames = F)