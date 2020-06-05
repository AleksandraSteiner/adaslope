library(SLOPE)
library(truncdist)
library(nlshrink)
library(MASS)
library(glmnet)
library(missMDA)
library(mice)

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

  init_list = initialization_functions(X, y,
                                       beta.start, 
                                       sigma.known, 
                                       sigma.init, 
                                       lambda,
                                       calculate_missing_cols(X),
                                       a, b)
  
  est_cache_result = iterate_algorithm(t = 0, 
                                       cstop = 1, 
                                       tol_em,
                                       init_list, 
                                       create_estimations_cache(X, 
                                                                maxit, 
                                                                init_list),
                                       maxit,
                                       print_iter,
                                       X,
                                       lambda,
                                       sigma.known)
  
  return(beta = est_cache_result[1][, maxit],
         estimations_cache_result[1], 
         beta.new = est_cache_result[1][, maxit] * 
           (rowMeans(est_cache_result[2][ , - (1 : 20)], na.rm = TRUE) > 1/2))
}

#small test
X = matrix(1:9, 3, 3)
y = c(3, 6, 28)
lambda = rep(0.5, 10)
a = 3
b = 4
beta.start = 1:3

ABSLOPE(X, y, lambda, #add seed
         a, b, beta.start,
         maxit = 300, print_iter = FALSE,
         tol_em = 1e-6, impute = 'PCA',
         sigma.known = NA, sigma.init = NA, scale = FALSE,
         method_na = 'lineq')
