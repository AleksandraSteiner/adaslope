library(glmnet)
library(SLOPE)
library(nlshrink)
library(MASS)

trunc_norm_gamma <- function(a, b) {
  norm_const <- b^a / gamma(a)
  return(pgamma(1, a, b) / norm_const)
}

slope_admm <- function(X, y, lambda, p, rho, max_iter=500, tol=1e-6) {
  
  M <- t(X) %*% X
  M <- M + diag(p) * rho
  M <- solve(M)
  MXtY <- M %*% (t(X) %*% Y)
  lam_seq_rho <- lambda / rho
  
  i <- 0
  x <- rep(0, p)
  z <- rep(0, p)
  u <- rep(0, p)
  z_new <- rep(0, p)
  z_new_arma <- rep(0, p)
  x_plus_u <- rep(0, p)
  dual_feas = 0.0; primal_feas = 0.0
  
  while(i < max_iter) {
    x <- MXtY + M %*% (rho * (z - u))
    x_plus_u <- x+u
    z_new <- prox_sorted_L1(x_plus_u, lam_seq_rho)
    z_new_arma <- z_new
    
    u <- u + (x - z_new)
    
    dual_feas = norm(rho * (z_new - z))
    primal_feas = norm(z_new - x)
    
    z = z_new
    
    if(primal_feas < tol & dual_feas < tol) {
      i = max_iter
    }
  }
  
  return(z)
}

SLOBE <- function(X, y, lambda,
                  beta_start=NA,
                  sigma_known=NA, sigma_init=NA,
                  a_prior=1, b_prior=1,
                  max_iter=500, tol=1e-5,
                  verbose=FALSE)
{
  p <- ncol(X)
  n <- length(y)
  
  gamma <- rep(0, p)
  i <- 1
  converged <- FALSE
  
  # Initialize beta
  if(is.na(beta_start)){ # initialization not given
    # LASSO cv
    objstart = cv.glmnet(X, y)
    beta = coef(objstart,s = "lambda.1se")
    beta = beta[2:(p+1)]
  } else{beta = beta_start} # initialization given
  
  # Initialize sigma
  if(is.na(sigma_known)){ # real value of sigma is unknown
    if(is.na(sigma_init)){ # initialization not given
      #sigma = sd(y - X.sim %*% beta) 
      sigma = sqrt(sum((y - X %*% beta)^2)/(n-1))
    } else{sigma <- sigma_init} # initialization given
  } else{sigma <- sigma_known} # real value of sigma is known
  
  # Initialize theta
  theta <- (sum(beta != 0)+a_prior)/(p+b_prior+a_prior);
  
  # Initialize c
  lambda_sigma_inv <- lambda / sigma
  lambda_sigma <- lambda * sigma
  c <- min((sum(abs(beta)>0)+1)/sum(abs(beta[beta != 0]))/lambda_sigma_inv[p],1)
  if(!is.finite(c)) c <- 1
  
  # Initialize mu, Sigma
  mu = apply(X,2,mean)
  Sigma = linshrink_cov(X)
  
  # Initialize parameter history vectors
  sigmas <- c(sigma)
  
  # Main loop
  while(i < max_iter) {
    
    mu <- apply(X, 2, mean)
    Sigma <- linshrink_cov(X)
    
    if(verbose) {
      print(sprintf('Iteration %s', i))
    }
    
    # Update gammas
    
    W <- diag(c * gamma + 1 - gamma)
    Wbeta <- W %*% beta
    ord_Wbeta <- order(Wbeta, decreasing=TRUE)
    
    gamma_new <- rep(0, p)
    
    for(j in 1:p) {
      exp1 <- exp(-c / sigma * abs(beta[j]) * lambda[ord_Wbeta[j]])
      exp2 <- exp(-1 / sigma * abs(beta[j]) * lambda[ord_Wbeta[j]])
      
      gamma_new[j] <- theta * c * exp1 / ((1 - theta) * exp2 + theta * c * exp1)
    }
    
    # Update theta
    
    theta_new <- (a_prior + sum(gamma == 1)) / (a_prior + b_prior + p)
    
    # Update c
    
    Wbeta <- W %*% beta
    ord_Wbeta <- order(Wbeta, decreasing=TRUE)
    a_prime <- 1 + sum(gamma == 1)
    b_prime <- 1 / sigma * sum(abs(beta) * lambda[ord_Wbeta] * (gamma == 1))
    
    c_new <- trunc_norm_gamma(a_prime+1, b_prime) / trunc_norm_gamma(a_prime, b_prime)
    
    # Generate missing observations
    ## TODO
    
    # Get new parameters by maximizing likelihood
    
    # Maximize for beta
    beta_new <- SLOPE_solver(X, y, lambda=lambda_sigma)$x
      
    # Maximize for sigma
    rk <- p-rank(abs(beta),ties.method="max")+1;
    RSS<- sum((y-X%*%beta)^2)
    sum_lamwbeta <- sum(lambda[rk]*abs(beta_new))
    sigma_new <- (sqrt(sum_lamwbeta^2+4*n*RSS)+sum_lamwbeta)/(2*n)
    
    # Maximize for mu
    mu_new <- mu # Remains the same if no missing data
    
    # Maximize for Sigma
    Sigma_new <- Sigma # Remains the same if no missing data
    
    # Check convergence condition
    err <- t(beta_new - beta) %*% (beta_new - beta)
    if(err < tol) {
      i <- max_iter
      converged <- TRUE
    }
    
    # Update parameters, parameter history and increase iter count
    lambda_sigma <- lambda * sigma
    gamma <- gamma_new
    theta <- theta_new
    c <- c_new
    
    beta <- beta_new
    sigma <- sigma_new
    mu <- mu_new
    Sigma <- Sigma_new
    
    sigmas <- c(sigmas, sigma)
    i <- i+1
  }
  
  gamma <- which(beta != 0)
  
  return(list('beta'=beta, 'selected'=gamma, 'sigmas'=sigmas))
}
