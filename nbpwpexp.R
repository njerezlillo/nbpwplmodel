source("/Users/nixonjerez/Library/CloudStorage/GoogleDrive-njerezlillo@gmail.com/Mi unidad/Manuscritos/Unification of piecewise models/pwexp.R")

#require(fastDummies)
require(numDeriv)

# q = -1 : bernoulli cure fraction
# q = 0 : poisson cure fraction
# q = -0.5 : negative binomial cure fraction
# q = 1 : geometric cure fraction

# baseline functions ------------------------------------------------------

dnbpwpexp_base <- function(x, q, theta, alpha, lambda, p) {
  if (q == 0) {
    theta * lambda * dpwexp(x, p, alpha) * 
      ppwexp(x, p, alpha)^(lambda - 1) * 
      exp(-theta * ppwexp(x, p, alpha)^lambda)
  } else {
    theta * lambda * dpwexp(x, p, alpha) * 
      ppwexp(x, p, alpha)^(lambda - 1) * 
      (1 + q * theta * ppwexp(x, p, alpha)^lambda)^(-1/q - 1)
  }
}

snbpwpexp_base <- function(x, q, theta, alpha, lambda, p) {
  if (q == 0) {
    exp(-theta * ppwexp(x, p, alpha)^lambda)
  } else {
    (1 + q * theta * ppwexp(x, p, alpha)^lambda)^(-1/q)
  }
}

dnbpwpexp <- function(x, q, p0, alpha, lambda, p) {
  if (q == 0) {
    theta <- -log(p0)
    theta * lambda * dpwexp(x, p, alpha) * 
      ppwexp(x, p, alpha)^(lambda - 1) * 
      exp(-theta * ppwexp(x, p, alpha)^lambda)
  } else {
    theta <- (p0^(-q) - 1)/q
    theta * lambda * dpwexp(x, p, alpha) * 
      ppwexp(x, p, alpha)^(lambda - 1) * 
      (1 + q * theta * ppwexp(x, p, alpha)^lambda)^(-1/q - 1)
  }
}

snbpwpexp <- function(x, q, p0, alpha, lambda, p) {
  if (q == 0) {
    theta <- -log(p0)
    exp(-theta * ppwexp(x, p, alpha)^lambda)
  } else {
    theta <- (p0^(-q) - 1)/q
    (1 + q * theta * ppwexp(x, p, alpha)^lambda)^(-1/q)
  }
}

# log-likelihood function -------------------------------------------------

S_nbpwpexp <- function(D, p0, q, alpha, lambda, p) {
  t <- D
  ppwpl <- Vectorize(function(z) ppwexp(z, p, alpha))
  if (q == 0) {
    theta <- -log(p0)
    out <- exp(-theta * ppwpl(t)^lambda)
  } else {
    theta <- (p0^(-q) - 1)/q
    out <- (1 + q * theta * ppwpl(t)^lambda)^(-1/q)
  }
  out
}

D_nbpwpexp <- function(D, p0, q, alpha, lambda, p) {
  t <- D
  dpwpl <- Vectorize(function(z) dpwexp(z, p, alpha))
  ppwpl <- Vectorize(function(z) ppwexp(z, p, alpha))
  if (q == 0) {
    theta <- -log(p0)
    out <- theta * lambda * dpwpl(t) * ppwpl(t)^(lambda - 1) * 
      exp(-theta * ppwpl(t)^lambda)
  } else {
    theta <- (p0^(-q) - 1)/q
    out <- theta * lambda * dpwpl(t) * ppwpl(t)^(lambda - 1) * 
      (1 + q * theta * ppwpl(t)^lambda)^(-1/q - 1)
  }
  return(out)
}

loglik_nbpwpexp <- function(arg, D, p, q) {
  alpha <- arg[1:length(p)] 
  temp <- arg[-(1:length(p))]
  lambda <- temp[1]
  beta <- plogis(sum(temp[-1] * D[-(1:2)]))
  
  aux1 <- 1
  aux2 <- 1
  
  if(D[2] == T) {
    aux1 <- D_nbpwpexp(D[1], beta, q, alpha, lambda, p)
  } else{
    aux2 <- S_nbpwpexp(D[1], beta, q, alpha, lambda, p)
  }
  sum(log(aux1)) + sum(log(aux2))
}
