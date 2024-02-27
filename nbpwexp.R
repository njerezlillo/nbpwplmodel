source("/Users/nixonjerez/Library/CloudStorage/GoogleDrive-njerezlillo@gmail.com/Mi unidad/Manuscritos/Unification of piecewise models/pwexp.R")
require(fastDummies)
require(numDeriv)
require(Rfast)

# q = -1 : bernoulli cure fraction
# q = 0 : poisson cure fraction
# q = -0.5 : negative binomial cure fraction
# q = 1 : geometric cure fraction

# baseline functions ------------------------------------------------------

dnbpwexp_base <- function(x, q, theta, alpha, p) {
  if (q == 0) {
    theta * dpwexp(x, p, alpha) * 
      exp(-theta * ppwexp(x, p, alpha))
  } else {
    theta * dpwexp(x, p, alpha) * 
      (1 + q * theta * ppwexp(x, p, alpha))^(-1/q - 1)
  }
}

snbpwexp_base <- function(x, q, theta, alpha, p) {
  if (q == 0) {
    exp(-theta * ppwexp(x, p, alpha))
  } else {
    (1 + q * theta * ppwexp(x, p, alpha))^(-1/q)
  }
}

dnbpwexp <- function(x, q, p0, alpha, p) {
  if (q == 0) {
    theta <- -log(p0)
    theta * dpwexp(x, p, alpha) * 
      exp(-theta * ppwexp(x, p, alpha))
  } else {
    theta <- (p0^(-q) - 1)/q
    theta * dpwexp(x, p, alpha) * 
      (1 + q * theta * ppwexp(x, p, alpha))^(-1/q - 1)
  }
}

snbpwexp <- function(x, q, p0, alpha, p) {
  if (q == 0) {
    theta <- -log(p0)
    exp(-theta * ppwexp(x, p, alpha))
  } else {
    theta <- (p0^(-q) - 1)/q
    (1 + q * theta * ppwexp(x, p, alpha))^(-1/q)
  }
}

# log-likelihood function -------------------------------------------------

S_nbpwexp <- function(D, p0, q, alpha, p) {
  t <- D
  ppwpl <- Vectorize(function(z) ppwexp(z, p, alpha))
  if (q == 0) {
    theta <- -log(p0)
    out <- exp(-theta * ppwexp(t))
  } else {
    theta <- (p0^(-q) - 1)/q
    out <- (1 + q * theta * ppwpl(t))^(-1/q)
  }
  out
}

D_nbpwexp <- function(D, p0, q, alpha, p) {
  t <- D
  dpwpl <- Vectorize(function(z) dpwexp(z, p, alpha))
  ppwpl <- Vectorize(function(z) ppwexp(z, p, alpha))
  if (q == 0) {
    theta <- -log(p0)
    out <- theta * dpwpl(t) * exp(-theta * ppwexp(t))
  } else {
    theta <- (p0^(-q) - 1)/q
    out <- theta * dpwpl(t) * (1 + q * theta * ppwpl(t))^(-1/q - 1)
  }
  return(out)
}


loglik_nbpwexp <- function(arg, D, p, q) {
  alpha <- arg[1:length(p)] 
  beta <- arg[length(p) + 1]
  
  aux1 <- 1
  aux2 <- 1
  
  if(D[2] == T) {
    aux1 <- D_nbpwexp(D[1], beta, q, alpha, p)
  } else{
    aux2 <- S_nbpwexp(D[1], beta, q, alpha, p)
  }
  sum(log(aux1)) + sum(log(aux2))
}

# Regression model --------------------------------------------------------

alpha_matrix <- function(arg, df, p) {
  DF <- df[, 1:2]
  X <- df[, -(1:2)]
  
  n_beta <- ncol(X)
  n_link <- length(p) 
  n_group <- 2^(n_beta - 1)
  
  index_base <-
    Filter(function(x)
      any(x == 1),
      do.call("c", lapply(1:n_beta, function(i)
        combn(1:n_beta, i, FUN = list))))
  
  Alpha <- matrix(nrow = n_group, ncol = n_link + 1)
  colnames(Alpha) <- c(paste0("alpha", 1:n_link), "p0")
  rownames(Alpha) <- paste("group", 1:n_group)
  
  for (i in 1:n_group) {
    for (j in 0:(n_link - 2)) {
      index <- lapply(index_base, function(x) x + n_beta * j)
      Alpha[i, j + 1] <- 1.001 + exp(sum(arg[index[[i]]]))
    }
    
    index <- lapply(index_base, function(x) x + n_beta * (n_link - 1))
    Alpha[i, n_link] <- 3.001 + exp(sum(arg[index[[i]]]))
    
    index <- lapply(index_base, function(x) x + n_beta * n_link)
    Alpha[i, n_link + 1] <- plogis(sum(arg[index[[i]]]))
  }
  
  Alpha
}

theta_matrix <- function(arg, df, p) {
  X <- df[, -(1:2)]
  n_beta <- ncol(X)
  
  output <- matrix(arg, ncol = n_beta)
  colnames(output) <- paste("interval", 1:n_beta)
  temp <- paste0("beta", 0:(n_beta - 1))
  rownames(output) <- paste0(colnames(df[,-(1:2)]), 
                             "(", temp, ")")
  output
}

loglik_nbpwexp_reg_model_index <- function(arg, D, p, q = -1) {
  n_beta <- length(D[-(1:2)])
  n_link <- length(p) 
  n_group <- 2^(n_beta - 1)
  
  temp <- which(D[-(1:2)] == 1)
  
  Alpha <- vector(length = n_link + 1)
  
  for (k in 0:(n_link - 2)) {
    Alpha[k + 1] <- 1.001 + exp(sum(arg[temp + n_beta * k]))
  }
  
  Alpha[n_link] <- 3.001 + exp(sum(arg[temp + n_beta * (n_link - 1)]))
  Alpha[n_link + 1] <- plogis(sum(arg[temp + n_beta * (n_link)]))
  
  loglik_nbpwexp(Alpha, D[1:3], p, q)
}

loglik_nbpwexp_reg_model <- function(arg, D, p, q = -1) {
  l <- function(z) loglik_nbpwexp_reg_model_index(arg, z, p, q = -1)
  sum(apply(D, 1, l))
}

# Delta method ------------------------------------------------------------

deltamethod <- function (arg, cov, H, p) {
  n_beta <- length(cov)
  n_link <- length(p) 
  
  M <- matrix(0, nrow = n_link + 1, ncol = sum(arg != 0))
  index <- NULL
  
  for (i in 0:n_link) {
    index_new <- which(arg[1:n_beta + n_beta * i] != 0)
    aux <- length(index)
    zzz <- sum(arg[which(cov[index_new] == 1) + n_beta * i])
    if (i < n_link) {
      M[i + 1, (length(index) + 1):(length(index) + length(index_new))] <- 
        cov[index_new] * exp(zzz)
    } else {
      M[i + 1, (length(index) + 1):(length(index) + length(index_new))] <- 
        cov[index_new] * exp(zzz) / (1 + exp(zzz))^2
    }
    
    index <- c(index_new, index)
  }
  
  sqrt(diag(M %*% H %*% t(M)))
}

se_matrix <- function(arg, H, df, p) {
  n_beta <- ncol(df[, -(1:2)])
  n_link <- length(p) 
  n_group <- 2^(n_beta - 1)
  
  index_base <-
    Filter(function(x)
      any(x == 1),
      do.call("c", lapply(1:n_beta, function(i)
        combn(1:n_beta, i, FUN = list))))
  
  Alpha <- matrix(nrow = n_group, ncol = n_link + 1)
  colnames(Alpha) <- c(paste0("alpha", 1:n_link), "p0")
  rownames(Alpha) <- paste("group", 1:n_group)
  
  for (i in 1:n_group) {
    cov <- rep(0, n_beta)
    cov[index_base[[i]]] <- 1
    Alpha[i, ] <- deltamethod(arg, cov, H, p)
  }
  
  Alpha
}