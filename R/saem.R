#' Create an adaSlope model with full algorithm
#' 
#' @param X Design matrix
#' @param y Response vector
#' @param lambda Vector of coefficient L1 penalties
#' @param a Prior for coefficient calculation
#' @param b Prior for coefficient calculation
#' @param beta.start Initial value for beta
#' @param maxit Maximum number of iterations
#' @param tol_em The tolerance number to stop SAEM
#' @param impute Imputation method
#' @param sigma.known Logical value that tells whether sigma is known
#' @param sigma.init Value of known sigma
#' @param print_iter Logical value that tells whether the estimated parameters will be printed in each iteration of SAEM
#' @param scale Logical value that tells whether scaling will be performed in each iteration
#' 
#' @return A SLOPE model
#' 
#' @examples
#'X = matrix(rnorm(1000), nrow=100)
#'b = c(sample(-5:5, 5), rep(0, 5))
#'y = X %*% b + rnorm(100, 0, 0.1)
#'A <- ABSLOPE(X, y,lambda=seq(10, 5, length.out=10), a=1, b=1)
#'@export
ABSLOPE = function(X, y,
                   lambda, #add seed
                   a, b,
                   beta.start = NA,
                   maxit = 300,
                   print_iter = FALSE,
                   tol_em = 1e-6,
                   impute = 'PCA',
                   sigma.known = NA,
                   sigma.init = NA,
                   scale = FALSE,
                   method_na = 'lineq') {
  
  #firstly checking missingness (will be done later)
  #...
  
  init_list = get_initial_parameters(X, y,
                                     beta.start,
                                     sigma.known,
                                     sigma.init,
                                     lambda,
                                     calculate_missing_cols(X),
                                     a, b)
  
  est_cache_result = iterate_saem_algorithm(cstop = 1,
                                            tol_em,
                                            init_list,
                                            create_estimations_cache(X,
                                                                     maxit,
                                                                     init_list),
                                            maxit,
                                            print_iter,
                                            X, y,
                                            lambda,
                                            sigma.known,
                                            scale,
                                            a, b)
  
  beta = as.data.frame(est_cache_result[[1]][1])[, maxit]
  betas = t(as.data.frame(est_cache_result[[1]][1])[, 1:maxit])
  gamma = which(as.data.frame(est_cache_result[[1]][2])[, maxit] > 0)
  sigmas = as.data.frame(est_cache_result[[1]][3])[, 1:maxit]
  
  res = list('beta'=beta, 'selected'=gamma, 'sigmas'=sigmas, 'betas'=betas)
  class(res) <- c('abslope')
  
  return(res)
}

#small test
lambda = rep(0.5, 10)
a = 3
b = 4
beta.start = 1:3
X = matrix(rnorm(30), nrow=10)
beta.start = c(sample(-5:5, 2), 0)
y = X %*% beta.start + rnorm(10) #important to add error (otherwise init sigma equals to 0!)


#ABSLOPE(X, y, lambda,
#        a, b, beta.start,
#        maxit = 300, print_iter = FALSE,
#        tol_em = 1e-6, impute = 'PCA',
#        sigma.known = NA, sigma.init = NA, 
#        scale = FALSE, method_na = 'lineq')
