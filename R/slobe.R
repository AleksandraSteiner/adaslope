#' Create an adaSlope model with simplified algorithm
#' 
#' @param X Design matrix
#' @param y Response vecor
#' @param lambda Vector of coefficient L1 penalties
#' @param a Prior for coefficient calculation
#' @param b Prior for coefficient calculation
#' @param maxit
#' @param tol_em
#' @param impute
#' @param sigma.known
#' @param sigma.init
#' @param print_iter
#' 
#' @return A SLOPE model
#' 
#' @examples
#'X = matrix(rnorm(1000), nrow=100)
#'b = c(sample(-5:5, 5), rep(0, 5))
#'y = X %*% b + rnorm(100, 0, 0.1)
#'A <- SLOBE(X, y, lambda=seq(10, 5, length.out=10))
#'@export
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
    #beta_new <- SLOPE_solver(X, y, lambda=lambda_sigma)$x
    beta_new <- estimate_beta_ML(gamma, c, X, y, lambda_sigma)
    
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


#X = matrix(rnorm(1000), nrow=100)
#b = c(sample(-5:5, 5), rep(0, 5))
#y = X %*% b + rnorm(100, 0, 0.1)
#A <- SLOBE(X, y, lambda=seq(10, 5, length.out=10))
#seq(5, 1, length.out=10)
#A
#class(A)
#summary(A)
#A$betas
#plot(A)
