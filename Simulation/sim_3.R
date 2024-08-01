require(numDeriv)
require(parallel)
require(doParallel)
#source("pwpowerlaw.R")
source("nbpwpl.R")
load("preliminar2.RData")

set.seed(2024)

FF <- function (n) {
  R <- 300
  Beta_real <- c(-0.74, 0.61,  0.00,
               0.67, 0.00,  0.00,
               -1.81, 0.00, 0.00,
               -5.71, 0.00, 4.81)
  p <- c(1, 1.5, 2)
  Alpha_real <-
    matrix(as.numeric(t(alpha_matrix(Beta_real, df[,-6], p))),
           nrow = 16, ncol = R)
  m <- 2
  
  MLE_Beta <- EE_Beta <-
    data.frame(matrix(ncol = R, nrow = length(Beta_real)))
  MLE <- EE <- data.frame(matrix(ncol = R, nrow = 16))
  i <- 1
  
  while (i <= R) {
    while (T) {
      D <- rnbpwplx(n, Beta_real, p = p, m = m)
      target <- function(w) -loglik_nbpwpl_reg_model(w, D, p)
      
      fit <- try(optim(par = runif(length(Beta_real), -0.5, 0.5), 
                       fn = target,
                       method = "CG",
                       hessian = T), silent = T)
      if (class(fit) != "try-error") break
    }
    
    H_temp <- solve((fit$hessian + t(fit$hessian)) / 2)
    
    if (all(diag(H_temp) > 0)) {
      ee_temp <- as.numeric(t(se_matrix(fit$par, H_temp, D, p)))
      if (sum(is.na(ee_temp)) == 0) {
        MLE_Beta[, i] <- fit$par
        EE_Beta[, i] <- sqrt(diag(H_temp))
        MLE[, i] <- as.numeric(t(alpha_matrix(fit$par, D, p)))
        EE[, i] <- ee_temp
        i <- i + 1
      }
    }
    cat("Iteration:", i, "\r")
  }
  
  li_Beta <- MLE_Beta - 1.96 * EE_Beta
  ls_Beta <- MLE_Beta + 1.96 * EE_Beta
  li <- MLE - 1.96 * EE
  ls <- MLE + 1.96 * EE
  
  BS_Beta <- apply(MLE_Beta - Beta_real, 1, mean)
  SE_Beta <- apply(MLE_Beta - Beta_real, 1, sd)
  MS_Beta <- apply(MLE_Beta - Beta_real, 1, function(x) mean(x^2))
  CP_Beta <- apply(li_Beta < Beta_real & ls_Beta > Beta_real, 1, mean)
  
  BS <- apply(MLE - Alpha_real, 1, mean)
  SE <- apply(MLE - Alpha_real, 1, sd)
  MS <- apply(MLE - Alpha_real, 1, function(x) mean(x^2))
  CP <- apply(li < Alpha_real & ls > Alpha_real, 1, mean)
  
  MC_BS_Beta <- SE_Beta / sqrt(R)
  MC_MS_Beta <-
    sqrt(apply((MLE_Beta - Beta_real) ^ 2 - matrix(MS_Beta, ncol = 3, nrow = length(MS_Beta)), 1,
               function(x) sum(x ^ 2)) / (R * (R - 1)))
  MC_CP_Beta <- sqrt(CP_Beta * (1 - CP_Beta) / R)
  
  MC_BS <- SE / sqrt(R)
  MC_MS <- 
    sqrt(apply((MLE - Alpha_real)^2 - matrix(MS, ncol = 3, nrow = length(MS)), 1,
               function(x) sum(x^2))/ (R * (R - 1)))
  MC_CP <- sqrt(CP * (1 - CP) / R)
  
  list(Beta = data.frame(BS_Beta, MS_Beta, CP_Beta, MC_BS_Beta, MC_MS_Beta, MC_CP_Beta),
       Alpha = data.frame(BS, MS, CP, MC_BS, MC_MS, MC_CP))
}

number_of_cores <- 20
clusters <- parallel::makeCluster(number_of_cores)
doParallel::registerDoParallel(clusters)
registerDoParallel(cores = 20)

ptm <- proc.time() #start
estimates <- foreach(z = c(750, 1000, 1250)) %dopar% { FF(z) }
proc.time() - ptm #final

estimates
