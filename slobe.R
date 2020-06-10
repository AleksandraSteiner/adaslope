library(glmnet)
library(SLOPE)
library(nlshrink)
library(MASS)
source("auxiliary_slobe.R")

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
                  a=1, b=1,
                  beta.start=NA,
                  maxit=500, tol_em=1e-5,
                  impute='PCA',
                  sigma.known=NA, sigma.init=NA,
                  print_iter=FALSE)
{
  
  initialization_list <- get_initial_parameters(X, y,
                                                beta.start,
                                                sigma.known,
                                                sigma.init,
                                                lambda,
                                                calculate_missing_cols(X),
                                                a, b)

  beta = initialization_list[[1]]
  gamma = initialization_list[[2]] 
  sigma = initialization_list[[3]]
  theta = initialization_list[[4]]
  c = initialization_list[[5]]
  Sigma = initialization_list[[6]]
  rank = initialization_list[[7]]
  
  lambdas = calculate_lambdas(lambda, sigma)
  lambda_sigma = lambdas[[1]]
  lambda_sigma_inv = lambdas[[2]]
  
  betas <- matrix(NA, nrow=maxit, ncol=ncol(X))
  betas[1,] <- beta
  
  p <- ncol(X)
  n <- length(y)
  
  i <- 1
  converged <- FALSE
  
  # Initialize parameter history vectors
  sigmas <- c(sigma)
  
  # Main loop
  while(i < maxit && !converged) {
    
    lambda_sigma <- lambda * sigma
    mu <- apply(X, 2, mean)
    Sigma <- linshrink_cov(X)
    
    if(print_iter) {
      print(sprintf('Iteration %s', i))
    }
    
    gamma_new <- slobe_update_gamma(gamma, c, beta, p, sigma, lambda, theta)

    theta_new <- slobe_update_theta(gamma, a, b, p)
    
    c_new <- slobe_update_c(gamma, c, beta, sigma, lambda)
    
    # Generate missing observations
    ## TODO
    X <- slobe_approx_Xmis(X)
    
    # Get new parameters by maximizing likelihood
    
    # Maximize for beta
    beta_new <- SLOPE_solver(X, y, lambda=lambda_sigma)$x
    #beta_new <- estimate_beta_ML(gamma, c, X, y, lambda_sigma)
    
    sigma_new <- estimate_sigma_ML(X, y, beta, lambda,
                                   lambda_sigma, sigma.known, sigma,
                                   gamma, c)

    mu_new <- mu # Remains the same if no missing data

    Sigma_new <- Sigma # Remains the same if no missing data
    
    # Check convergence condition
    err <- t(beta_new - beta) %*% (beta_new - beta)
    if(err < tol_em) {
      converged <- TRUE
    }
    
    # Update parameters, parameter history and increase iter count
    
    
    beta <- beta_new
    betas[i+1,] <- beta_new
    sigma <- sigma_new
    mu <- mu_new
    Sigma <- Sigma_new
    
    lambda_sigma <- lambda * sigma
    gamma <- gamma_new
    theta <- theta_new
    c <- c_new
    
    sigmas <- c(sigmas, sigma)
    i <- i+1
  }
  
  gamma <- which(beta != 0)
  
  res = list('beta'=beta, 'selected'=gamma, 'sigmas'=sigmas, 'betas'=betas[1:i,])
  class(res) <- c('abslope')
  
  return(res)
}


X = matrix(rnorm(1000), nrow=100)
b = c(sample(-5:5, 5), rep(0, 5))
y = X %*% b + rnorm(100, 0, 0.1)
A <- SLOBE(X, y, lambda=seq(5, 1, length.out=10), print_iter=TRUE)
seq(5, 1, length.out=10)
A
class(A)
summary(A)
A$betas
plot(A)
