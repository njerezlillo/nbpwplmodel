#require(fastDummies)
require(numDeriv)

# q = -1 : bernoulli cure fraction
# q = 0 : poisson cure fraction
# q = -0.5 : negative binomial cure fraction
# q = 1 : geometric cure fraction

# Part 1 ------------------------------------------------------------------

dnbpwwei_base <- function(x, q, theta, alpha) {
  if (q == 0) {
    theta * dweibull(x, alpha[1], alpha[2]) * 
      exp(-theta * pweibull(x, alpha[1], alpha[2]))
  } else {
    theta * dweibull(x, alpha[1], alpha[2]) * 
      (1 + q * theta * pweibull(x, alpha[1], alpha[2]))^(-1/q - 1)
  }
}

snbpwwei_base <- function(x, q, theta, alpha) {
  if (q == 0) {
    exp(-theta * pweibull(x, alpha[1], alpha[2]))
  } else {
    (1 + q * theta * pweibull(x, alpha[1], alpha[2]))^(-1/q)
  }
}

dnbpwwei <- function(x, q, p0, alpha) {
  if (q == 0) {
    theta <- -log(p0)
    theta * dweibull(x, alpha[1], alpha[2]) * 
      exp(-theta * pweibull(x, alpha[1], alpha[2]))
  } else {
    theta <- (p0^(-q) - 1)/q
    theta * dweibull(x, alpha[1], alpha[2]) * 
      (1 + q * theta * pweibull(x, alpha[1], alpha[2]))^(-1/q - 1)
  }
}

snbpwwei <- function(x, q, p0, alpha) {
  if (q == 0) {
    theta <- -log(p0)
    exp(-theta * pweibull(x, alpha[1], alpha[2]))
  } else {
    theta <- (p0^(-q) - 1)/q
    (1 + q * theta * pweibull(x, alpha[1], alpha[2]))^(-1/q)
  }
}

# Part 2 ------------------------------------------------------------------

S_nbpwwei <- function(D, beta, q, alpha) {
  t <- D$time
  X <- D[,-(1:2)]
  ppwwei <- Vectorize(function(z) pweibull(z, alpha[1], alpha[2]))
  p0 <- plogis(as.numeric(as.matrix(X) %*% beta))
  if (q == 0) {
    theta <- -log(p0)
    out <- exp(-theta * ppwwei(t))
  } else {
    theta <- (p0^(-q) - 1)/q
    out <- (1 + q * theta * ppwwei(t))^(-1/q)
  }
  out
}

D_nbpwwei <- function(D, beta, q, alpha) {
  t <- D$time
  X <- D[,-(1:2)]
  dpwwei <- Vectorize(function(z) dweibull(z, alpha[1], alpha[2]))
  ppwwei <- Vectorize(function(z) pweibull(z, alpha[1], alpha[2]))
  p0 <- plogis(as.numeric(as.matrix(X) %*% beta))
  if (q == 0) {
    theta <- -log(p0)
    out <- theta * dpwwei(t) * exp(-theta * ppwwei(t))
  } else {
    theta <- (p0^(-q) - 1)/q
    out <- theta * dpwwei(t) * (1 + q * theta * ppwwei(t))^(-1/q - 1)
  }
  return(out)
}

loglik_nbpwwei <- function(arg, D, q) {
  n_beta <- ifelse(ncol(D) == 3, 1, ncol(D[,-(1:2)]))
  alpha <- arg[-((length(arg) - n_beta + 1):length(arg))] 
  beta <- arg[((length(arg) - n_beta + 1):length(arg))]
  
  D_obs <- D[D$status == 1,]
  D_cen <- D[D$status == 0,]
  
  aux1 <- D_nbpwwei(D_obs, beta, q, alpha)
  aux2 <- S_nbpwwei(D_cen, beta, q, alpha)
  sum(log(aux1)) + sum(log(aux2))
}

# Part 3 ------------------------------------------------------------------

preliminar_rnbpwwei <- function (n, alpha, p0, lambda) {
  b <- rbinom(n, 1, p0)
  t <- c() ; y <- c(); delta <- c()
  y[b == 1] <- 1.000e+54
  y[b == 0] <- rpwpowerlaw(sum(b == 0), alpha)
  
  cax <- runif(n, 0, lambda)
  
  t[y <= cax] <- y[y <= cax]
  delta[y <= cax] <- 1
  
  t[y > cax] <- cax[y > cax]
  delta[y > cax] <- 0
  
  out <- data.frame(time = t, status = delta)
  return(out)
}

rnbpwwei <- function (n, alpha, beta, q, tmax = 100) {
  n_beta <- length(beta)
  #X <- cbind(rep(1, n), 
  #           dummy_cols(sample(letters[1:n_beta], n, T),
  #remove_first_dummy = T)[,-1])
  X <- data.frame(i = rep(1, n))
  out <- data.frame(time = rep(NA, n), status = rep(NA, n))
  
  for (j in 1:n) {
    u <- runif(1)
    
    ff <- function (zz) {
      df <- cbind(data.frame(time = zz, s = 1), as.data.frame(t(X[j,])))
      u - as.numeric(S_nbpwwei(df, beta, q, alpha))
    }
    
    y <- try(uniroot(ff, c(1, 100))$root, silent = T)
    y <- ifelse(is.numeric(y), y, tmax)
    y_star <- runif(1, 1, tmax)
    temp <- min(y, y_star)
    out$time[j] <- temp
    out$status[j] <- ifelse(temp == tmax | temp == y_star, 0, 1)
  }
  
  cbind(out, X)
}

deltamethod <- function (mod, cov) {
  temp1 <- diag(rep(1, 5))
  pp <- function(z) plogis(z) * (1 - plogis(z))
  temp2 <- matrix(0, ncol = 4)
  temp3 <- matrix(pp(mod[5]), ncol = 1)
  j <- rbind(temp1, cbind(temp2, temp3))
  v <- j %*% cov %*% t(j)
  return(v)
}

deltamethod2 <- function (mod, cov) {
  myfun <- function (x) {
    c(x[1], x[2], x[3], x[4], plogis(x[4]))
  }
  j <- jacobian(myfun, mod)
  v <- j %*% cov %*% t(j)
  return(v)
}
