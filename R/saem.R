#' Create an adaSlope model with full algorithm
#' 
#' @param X Design matrix
#' @param y Response vector
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
#'A <- ABSLOBE(X, y, lambda=seq(10, 5, length.out=10))
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
  
  return(beta = est_cache_result[1][, maxit],
         estimations_cache_result[1], 
         beta.new = est_cache_result[1][, maxit] * 
           (rowMeans(est_cache_result[2][ , -(1:20)], na.rm = TRUE) > 1/2))
}

#small test
X = matrix(1:9, 3, 3)
y = c(3, 6, 28)
lambda = rep(0.5, 10)
a = 3
b = 4
beta.start = 1:3

#Marcin data
X = matrix(rnorm(1000), nrow=100)
beta.start = c(sample(-5:5, 5), rep(0, 5))
y = X %*% beta.start + 1

#ABSLOPE(X, y, lambda,
#        a, b, beta.start,
#        maxit = 300, print_iter = FALSE,
#        tol_em = 1e-6, impute = 'PCA',
#        sigma.known = NA, sigma.init = NA, 
#        scale = FALSE, method_na = 'lineq')