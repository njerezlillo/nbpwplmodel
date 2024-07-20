source("/Users/nixonjerez/Library/CloudStorage/GoogleDrive-njerezlillo@gmail.com/Mi unidad/Manuscritos/Unification of piecewise models/pwexp.R")

#require(fastDummies)
require(numDeriv)

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
    out <- exp(-theta * ppwpl(t))
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
    out <- theta * dpwpl(t) * exp(-theta * ppwpl(t))
  } else {
    theta <- (p0^(-q) - 1)/q
    out <- theta * dpwpl(t) * (1 + q * theta * ppwpl(t))^(-1/q - 1)
  }
  return(out)
}

loglik_nbpwexp <- function(arg, D, p, q) {
  alpha <- arg[1:length(p)] 
  beta <- plogis(sum(arg[-(1:length(p))] * D[-(1:2)]))
  
  aux1 <- 1
  aux2 <- 1
  
  if(D[2] == T) {
    aux1 <- D_nbpwexp(D[1], beta, q, alpha, p)
  } else{
    aux2 <- S_nbpwexp(D[1], beta, q, alpha, p)
  }
  sum(log(aux1)) + sum(log(aux2))
}


