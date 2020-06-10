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
  
  list(beta = as.data.frame(est_cache_result[[1]][1])[, maxit],
       beta.new = as.data.frame(est_cache_result[[1]][1])[, maxit] *
         (rowMeans(as.data.frame(est_cache_result[[1]][2])[ , -(1:20)], na.rm = TRUE) > 1/2))
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
