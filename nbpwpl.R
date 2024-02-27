source("/Users/nixonjerez/Library/CloudStorage/GoogleDrive-njerezlillo@gmail.com/Mi unidad/Manuscritos/Unification of piecewise models/pwpowerlaw.R")
require(numDeriv)

# q = -1 : bernoulli cure fraction
# q = 0 : poisson cure fraction
# q = -0.5 : negative binomial cure fraction
# q = 1 : geometric cure fraction

# baseline functions ------------------------------------------------------

dnbpwpl_base <- function(x, q, theta, alpha, p) {
  if (q == 0) {
    theta * dpwpowerlaw(x, p, alpha) * 
      exp(-theta * ppwpowerlaw(x, p, alpha))
  } else {
    theta * dpwpowerlaw(x, p, alpha) * 
      (1 + q * theta * ppwpowerlaw(x, p, alpha))^(-1/q - 1)
  }
}

snbpwpl_base <- function(x, q, theta, alpha, p) {
  if (q == 0) {
    exp(-theta * ppwpowerlaw(x, p, alpha))
  } else {
    (1 + q * theta * ppwpowerlaw(x, p, alpha))^(-1/q)
  }
}

dnbpwpl <- function(x, q, p0, alpha, p) {
  if (q == 0) {
    theta <- -log(p0)
    theta * dpwpowerlaw(x, p, alpha) * 
      exp(-theta * ppwpowerlaw(x, p, alpha))
  } else {
    theta <- (p0^(-q) - 1)/q
    theta * dpwpowerlaw(x, p, alpha) * 
      (1 + q * theta * ppwpowerlaw(x, p, alpha))^(-1/q - 1)
  }
}

snbpwpl <- function(x, q, p0, alpha, p) {
  if (q == 0) {
    theta <- -log(p0)
    exp(-theta * ppwpowerlaw(x, p, alpha))
  } else {
    theta <- (p0^(-q) - 1)/q
    (1 + q * theta * ppwpowerlaw(x, p, alpha))^(-1/q)
  }
}

# Sampling ----------------------------------------------------------------

rnbpwpl <- function (n, alpha, p, beta, q, tmax = 100) {
  n_beta <- length(beta)
  X <- data.frame(i = rep(1, n))
  out <- data.frame(time = rep(NA, n), status = rep(NA, n))
  
  for (j in 1:n) {
    u <- runif(1)
    
    ff <- function (zz) {
      df <- cbind(data.frame(time = zz, s = 1), as.data.frame(t(X[j,])))
      u - as.numeric(S_nbpwpl(df, beta, q, alpha, p))
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

rnbpwplx <- function (n, arg, p, m, q = -1, tmax = 100) {
  n_link <- length(p) 
  
  X <- cbind(rep(1, n), matrix(rbinom(n * m, 1, 0.5), ncol = m))
  colnames(X) <- paste0("V", 0:(ncol(X) - 1))
  out <- data.frame(time = rep(NA, n), status = rep(NA, n))
  
  for (j in 1:n) {
    u <- runif(1)
    
    Alpha <- vector(length = n_link)
    for (k in 0:(n_link - 2)) {
      Alpha[k + 1] <- 1.001 + exp(sum(arg[which(X[j,] == 1) + (m + 1) * k]))
    }
    Alpha[n_link] <- 3.001 + exp(sum(arg[which(X[j,] == 1) + (m + 1) * (n_link - 1)]))
    p0 <- plogis(sum(arg[which(X[j,] == 1) + (m + 1) * (n_link)]))
    
    ff <- function (zz) {
      u - as.numeric(S_nbpwpl(zz, p0, q, Alpha, p))
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

# log-likelihood function -------------------------------------------------

S_nbpwpl <- function(D, p0, q, alpha, p) {
  t <- D
  ppwpl <- Vectorize(function(z) ppwpowerlaw(z, p, alpha))
  if (q == 0) {
    theta <- -log(p0)
    out <- exp(-theta * ppwpl(t))
  } else {
    theta <- (p0^(-q) - 1)/q
    out <- (1 + q * theta * ppwpl(t))^(-1/q)
  }
  out
}

D_nbpwpl <- function(D, p0, q, alpha, p) {
  t <- D
  dpwpl <- Vectorize(function(z) dpwpowerlaw(z, p, alpha))
  ppwpl <- Vectorize(function(z) ppwpowerlaw(z, p, alpha))
  if (q == 0) {
    theta <- -log(p0)
    out <- theta * dpwpl(t) * exp(-theta * ppwpl(t))
  } else {
    theta <- (p0^(-q) - 1)/q
    out <- theta * dpwpl(t) * (1 + q * theta * ppwpl(t))^(-1/q - 1)
  }
  return(out)
}


loglik_nbpwpl <- function(arg, D, p, q) {
  alpha <- arg[1:length(p)] 
  beta <- arg[length(p) + 1]
  
  aux1 <- 1
  aux2 <- 1
  
  if(D[2] == T) {
    aux1 <- D_nbpwpl(D[1], beta, q, alpha, p)
  } else{
    aux2 <- S_nbpwpl(D[1], beta, q, alpha, p)
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

loglik_nbpwpl_reg_model_index <- function(arg, D, p, q = -1) {
  n_beta <- length(D[-(1:2)])
  n_link <- length(p) 
  n_group <- 2^(n_beta - 1)
  
  temp <- which(D[-(1:2)] == 1)
  
  Alpha <- vector(length = n_link + 1)
  
  for (k in 0:(n_link - 2)) {
    Alpha[k + 1] <- 1.001 + exp(sum(arg[temp + n_beta * k]))
  }
  
  Alpha[n_link] <- 3.001 + exp(sum(arg[temp + n_beta * (n_link - 1)])) # 3.001
  Alpha[n_link + 1] <- plogis(sum(arg[temp + n_beta * (n_link)]))
  
  loglik_nbpwpl(Alpha, D[1:3], p, q)
}

loglik_nbpwpl_reg_model <- function(arg, D, p, q = -1) {
  l <- function(z) loglik_nbpwpl_reg_model_index(arg, z, p, q = -1)
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

# Influence ---------------------------------------------------------------

### Case weight

f_alpha <- function(t, alpha, p) {
  C <- auxiliar_pwpowerlaw(p, alpha)
  j <- which(index_each_interval(t, p) == 1)
  out <- rep(0, length(p))
  
  if (j > 1) {
    for(i in 1:(j - 1)){
      out[i] <- (alpha[j] - 1)/p[j] * (t/p[j])^(-alpha[j]) * log(p[i]/p[i + 1]) * C[j]
    }
  }
  
  out[j] <- (t/p[j])^(-alpha[j]) * (1 - (alpha[j] - 1) * log(t/p[j]))/p[j] * C[j]
  
  return(as.data.frame(out))
}

F_alpha <- function(t, alpha, p) {
  C <- auxiliar_pwpowerlaw(p, alpha)
  j <- which(index_each_interval(t, p) == 1)
  out <- rep(0, length(p))
  
  if (j > 1) {
    for(i in 1:(j - 1)){
      out[i] <- (t/p[j])^(1 - alpha[j]) * log(p[i + 1]/p[i]) * C[j]
    }
  }
  
  out[j] <- (t/p[j])^(1 - alpha[j]) * log(t/p[j]) * C[j]
  
  return(as.data.frame(out))
}

Deltai_beta_cw <- function(DF, arg, p) {
  t <- as.numeric(DF[1])
  status <- as.numeric(DF[2])
  x <- as.numeric(DF[-(1:2)])
  n_beta <- length(x)
  n_link <- length(p) 
  
  index_base <-
    Filter(function(z)
      any(z == 1),
      do.call("c", lapply(1:n_beta, function(i)
        combn(1:n_beta, i, FUN = list))))
  
  v <- which(x == 1) 
  position <- which(sapply(index_base, function(z) identical(z, v)))
  
  alpha <- alpha_matrix(arg, DF, p)[position, 1:length(p)]
  p0 <- alpha_matrix(arg, DF, p)[position, (length(p) + 1)]
  
  hi <- p0 ^ (-q) - 1
  gi <- 1 + hi * ppwpowerlaw(t, p, alpha)
  
  D_f_alpha <- status * f_alpha(t, alpha, p) / dpwpowerlaw(t, p, alpha) -
    ((status + q^(-1)) * hi * F_alpha(t, alpha, p))/gi

  out <- NULL
  
  for (i in 0:(n_link - 1)) {
    index_new <- 1:n_beta + (n_link + 1) * i 
    temp <- x[which(arg[index_new] != 0)]
    out <- c(out, temp * exp(sum(arg[index_new] * x)) * D_f_alpha[i + 1, 1])
  }
  
  out
}

Deltai_gamma_cw <- function(DF, arg, p) {
  t <- as.numeric(DF[1])
  status <- as.numeric(DF[2])
  x <- as.numeric(DF[-(1:2)])
  n_beta <- length(x)
  n_link <- length(p) 
  
  index_base <-
    Filter(function(z)
      any(z == 1),
      do.call("c", lapply(1:n_beta, function(i)
        combn(1:n_beta, i, FUN = list))))
  
  v <- which(x == 1) 
  position <- which(sapply(index_base, function(z) identical(z, v)))
  
  alpha <- alpha_matrix(arg, DF, p)[position, 1:length(p)]
  p0 <- alpha_matrix(arg, DF, p)[position, (length(p) + 1)]
  
  hi <- p0 ^ (-q) - 1
  gi <- 1 + hi * ppwpowerlaw(t, p, alpha)
  
  D_f_p0 <- -q * x[which(arg[1:n_beta + (n_link + 1) * n_link] != 0)] * 
    p0 ^ (-q) * (1 - p0) * (status/hi - (status + q^(-1)) * ppwpowerlaw(t, p, alpha)/gi)
  
  D_f_p0
}

### Response

f_wi <- function(t, alpha, p) {
  C <- auxiliar_pwpowerlaw(p, alpha)
  j <- max(which(p < t))
  -(alpha[j]*(alpha[j] - 1)/p[j]^2) * (t/p[j])^(-1 - alpha[j]) * C[j]
}

f_wi_alpha <- function(t, alpha, p) {
  C <- auxiliar_pwpowerlaw(p, alpha)
  j <- which(index_each_interval(t, p) == 1)
  out <- rep(0, length(p))
  
  if (j > 1) {
    for(i in 1:(j - 1)){
      out[i] <- alpha[j] * (alpha[j] - 1)/p[j]^2 * 
        (t/p[j])^(-1 - alpha[j]) * log(p[i + 1]/p[i]) * C[j]
    }
  }
  
  out[j] <- (1 - 2 * alpha[j] + (alpha[j]^2 - alpha[j])*log(t/p[j])) *
    t^(-1 - alpha[j])/ (p[j]^(1 - alpha[j])) * C[j]
  
  return(as.data.frame(out))
}

F_wi_alpha <- function(t, alpha, p) {
  f_alpha(t, alpha, p)
}

F_wi <- function(t, alpha, p) {
  dpwpowerlaw(t, p, alpha)
}

Deltai_beta_rs <- function(DF, arg, p) {
  t <- as.numeric(DF[1])
  status <- as.numeric(DF[2])
  x <- as.numeric(DF[-(1:2)])
  n_beta <- length(x)
  n_link <- length(p) 
  
  index_base <-
    Filter(function(z)
      any(z == 1),
      do.call("c", lapply(1:n_beta, function(i)
        combn(1:n_beta, i, FUN = list))))
  
  v <- which(x == 1) 
  position <- which(sapply(index_base, function(z) identical(z, v)))
  
  alpha <- alpha_matrix(arg, DF, p)[position, 1:length(p)]
  p0 <- alpha_matrix(arg, DF, p)[position, (length(p) + 1)]
  
  hi <- p0 ^ (-q) - 1
  gi <- 1 + hi * ppwpowerlaw(t, p, alpha)
  
  D_f_alpha <- status / dpwpowerlaw(t, p, alpha)^2 * 
    (f_wi_alpha(t, alpha, p) * dpwpowerlaw(t, p, alpha) - 
       f_alpha(t, alpha, p) * f_wi(t, alpha, p)) - (status + q^(-1)) * hi / gi^2 * 
    (gi * F_wi_alpha(t, alpha, p) - hi * F_wi(t, alpha, p) * F_alpha(t, alpha, p))
  
  out <- NULL
  
  for (i in 0:(n_link - 1)) {
    index_new <- 1:n_beta + (n_link + 1) * i 
    temp <- x[which(arg[index_new] != 0)]
    out <- c(out, temp * exp(sum(arg[index_new] * x)) * D_f_alpha[i + 1, 1])
  }
  
  out
  
}

Deltai_gamma_rs <- function(DF, arg, p) {
  t <- as.numeric(DF[1])
  status <- as.numeric(DF[2])
  x <- as.numeric(DF[-(1:2)])
  n_beta <- length(x)
  n_link <- length(p) 
  
  index_base <-
    Filter(function(z)
      any(z == 1),
      do.call("c", lapply(1:n_beta, function(i)
        combn(1:n_beta, i, FUN = list))))
  
  v <- which(x == 1) 
  position <- which(sapply(index_base, function(z) identical(z, v)))
  
  alpha <- alpha_matrix(arg, DF, p)[position, 1:length(p)]
  p0 <- alpha_matrix(arg, DF, p)[position, (length(p) + 1)]
  
  hi <- p0 ^ (-q) - 1
  gi <- 1 + hi * ppwpowerlaw(t, p, alpha)
  
  D_f_p0 <- q * p0 ^ (-q) * (status + q^(-1)) *
    x[which(arg[1:n_beta + (n_link + 1) * n_link] != 0)] * 
    F_wi(t, alpha, p) * (1 - p0)/ gi^2
  
  D_f_p0
}

### Out

Influence <- function(df, arg, p, hes, perturbation) {
  n <- nrow(df)
  
  if (perturbation == "cw") {
    Deltai <- 
      function(i) { c(Deltai_beta_cw(df[i,], arg, p),
                        Deltai_gamma_cw(df[i,], arg, p)) }
  } else if (perturbation == "rs") {
    Deltai <- 
      function(i) { c(Deltai_beta_rs(df[i,], arg, p),
                        Deltai_gamma_rs(df[i,], arg, p)) }
  }
  
  D <- as.matrix(do.call(cbind, lapply(1:n, Deltai)))
  C <- 2 * abs(t(D) %*% hes %*% D)
  C <- (C + t(C)) / 2
  Bi <- diag(C)/sum(diag(C))
  
  eigen_temp <- eigen(C, symmetric = T)
  
  valores <- eigen_temp$values
  vectores <- eigen_temp$vectors
  
  s <- valores[1] * vectores[, 1] %*% t(vectores[, 1])
  
  for(i in 2:n) {
    temp <- valores[i] * vectores[, i] %*% t(vectores[, i])
    s <- temp + s
  }
  
  Ui <- apply(eigen_temp$values * eigen_temp$vectors^2, 2, sum)
  dmax <- eigen_temp$vectors[, which.max(eigen_temp$values)]
  return(list(dmax = dmax, Bi = Bi, Ui = Ui))
}

