trunc_norm_gamma <- function(a, b) {
  norm_const <- b^a / gamma(a)
  return(pgamma(1, a, b) / norm_const)
}

slobe_update_gamma <- function(gamma, c, beta, p, sigma) {
  W <- diag(c * gamma + 1 - gamma)
  Wbeta <- W %*% beta
  ord_Wbeta <- order(Wbeta, decreasing=TRUE)
  
  gamma_new <- rep(0, p)
  
  for(j in 1:p) {
    exp1 <- exp(-c / sigma * abs(beta[j]) * lambda[ord_Wbeta[j]])
    exp2 <- exp(-1 / sigma * abs(beta[j]) * lambda[ord_Wbeta[j]])
    
    gamma_new[j] <- theta * c * exp1 / ((1 - theta) * exp2 + theta * c * exp1)
  }
  
  return(gamma_new)
}

slobe_update_theta <- function(gamma, a_prior, b_prior, p) {
  return ((a_prior + sum(gamma == 1)) / (a_prior + b_prior + p))
}

slobe_update_c <- function(gamma, c, beta, sigma) {
  W <- diag(c * gamma + 1 - gamma)
  Wbeta <- W %*% beta
  ord_Wbeta <- order(Wbeta, decreasing=TRUE)
  a_prime <- 1 + sum(gamma == 1)
  b_prime <- 1 / sigma * sum(abs(beta) * lambda[ord_Wbeta] * (gamma == 1))
  
  return(trunc_norm_gamma(a_prime+1, b_prime) / trunc_norm_gamma(a_prime, b_prime))
}
  
slobe_approx_Xmis <- function(X) {
  return(X)
}

