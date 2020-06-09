#checking missingnes functions (will be done later)
#...

check_missingness = function(X) {
  as.matrix(is.na(X))
}

detect_missing_col = function(X) {
  miss_presence = check_missingness(X)
  if (sum(miss_presence) > 0) 
    (1:ncol(miss_presence))[apply(miss_presence, 2, sum) > 0]
  else
    NULL
}

calculate_missing_cols = function(X) {
  length(detect_missing_col(X))
}

initialize_beta = function(beta.start, X, y) {
  p = ncol(X)
  if (is.na(beta.start[1])) {
    objstart = cv.glmnet(X, y)
    coeff = coef(objstart, s = "lambda.1se")
    coeff[2:(p + 1)]
  } 
  else beta.start
}

initialize_sigma = function(sigma.known, sigma.init, 
                            beta, X, y) {
  n = nrow(X)
  if (is.na(sigma.known)) {
    if (is.na(sigma.init)) {
      sqrt(sum((y - X %*% beta) ^ 2) / (n - 1))
    } 
    else sigma.init
  } 
  else sigma.known 
}

initialize_theta = function(beta, a, b, X) {
  p = ncol(X)
  (sum(beta != 0) + a) / (a + b + p)
}

initialize_rank = function(beta, X) {
  p = ncol(X)
  p - rank(abs(beta), ties.method = "max") + 1
}

calculate_lambdas = function(lambda, sigma) {
  list(lambda * sigma, lambda / sigma)
}

initialize_c = function(beta, lambda_sigma_inv, X) {
  p = ncol(X)
  c = min((sum(abs(beta) > 0) + 1) / 
            (sum(abs(beta[beta != 0])) * lambda_sigma_inv[p]), 1)
  if (!is.finite(c)) c = 1
  c
}

initialize_gamma = function(beta) {
  (abs(beta) > 0) * 1
}

initialize_mu_BSigma = function(X) {
  if (calculate_missing_cols(X) != 0) {
    mu = apply(X, 2, mean)
    Big_Sigma = linshrink_cov(X)
  }
  else
    mu = Big_Sigma = NULL
  list(mu, Big_Sigma)
}

get_initial_parameters = function(X, y, 
                                  beta.start, 
                                  sigma.known, 
                                  sigma.init, 
                                  lambda, 
                                  missingcols,
                                  a, b) {
  beta = initialize_beta(beta.start, X, y)
  
  sigma = initialize_sigma(sigma.known, sigma.init, 
                           beta, X, y)
  
  theta = initialize_theta(beta, a, b, X)
  
  rk = initialize_rank(beta, X)
  
  lambdas = calculate_lambdas(lambda, sigma)
  lambda_sigma = lambdas[[1]]
  lambda_sigma_inv = lambdas[[2]]
  
  c = initialize_c(beta, lambda_sigma_inv, X)
  
  gamma = initialize_gamma(beta)
  
  mu.Big_Sigma = initialize_mu_BSigma(X)
  mu = mu.Big_Sigma[1]
  Big_Sigma = mu.Big_Sigma[2]
  
  list(beta, gamma, sigma, theta, c, Big_Sigma, rk, mu)
}

create_estimations_cache = function(X, maxit, initialization_list) {
  betas_matrix = matrix(NA, ncol(X), (maxit + 1))
  gammas_matrix = matrix(NA, ncol(X), (maxit + 1))
  sigmas_vector = matrix(NA, 1, (maxit + 1))
  theta_vector = matrix(NA, 1, (maxit + 1))
  c_vector = matrix(NA, 1, (maxit + 1))
  betas_matrix[, 1] = initialization_list[[1]]
  gammas_matrix[, 1] = initialization_list[[2]]
  sigmas_vector[, 1] = initialization_list[[3]]
  theta_vector[, 1] = initialization_list[[4]]
  c_vector[, 1] = initialization_list[[5]]
  list(betas_matrix, 
       gammas_matrix, 
       sigmas_vector,
       theta_vector, 
       c_vector)
} 

update_step_size = function(t) {
  ifelse(t <= 20, 1, 1 / (t - 20))
}

generate_binomial_prob = function(beta, theta, c, lambda_sigma_inv, rank) {
  product = abs(beta) * lambda_sigma_inv[rank]
  1 / (1 + ((1 - theta) / (theta * c)) * exp((c - 1) * product))
}

generate_gamma_t = function(p, binomial_prob) {
  rbinom(p, 1, binomial_prob)
}

generate_c_t = function(beta, gamma, lambda_sigma_inv, rank) {
  a_gamma = 1 + sum(gamma == 1)
  b_gamma = sum(abs(beta) * lambda_sigma_inv[rank] * (gamma == 1))
  if (a_gamma > 1) {
    if (b_gamma > 0)
      rtrunc(1, "gamma", 0, 1, shape = a_gamma, rate = b_gamma)
    else
      rbeta(1, shape1 = a_gamma, shape2 = 1)
  }
  else
    runif(1, 0, 1)
}

scale_mean.w = function(V, weight) {
  res = sum(V * weight, na.rm = TRUE) / sum(weight[!is.na(V)])
}
scale_std.w = function(V, weight) {
  res = sqrt(sum(V ^ 2 * weight, na.rm = TRUE) / sum(weight[!is.na(V)]))
}

scale_X_mis = function(X_mis, scale, row.w, std.w, mean.w) {
  n = nrow(X)
  if (scale) {
    X_mis = t(t(X_mis) * std.w)
    X_mis = t(t(X_mis) + mean.w)
    mean.w = apply(X_mis, 2, scale_mean.w, row.w)
    X_mis = t(t(X_mis) - mean.w)
    std.w = apply(X_mis, 2, scale_std.w, row.w) * sqrt(n)
    X_mis = t(t(X_mis) / std.w)
  }
}

#TODO: this function is still too composite
generate_X = function(Big_Sigma, mu, beta, sigma, 
                      y, X, row.w, std.w, mean.w) {
  missingcols = calculate_missing_cols(X)
  if (missingcols != 0) {
    Big_Sigma_inv = solve(Big_Sigma)
    tau = sqrt(diag(Big_Sigma_inv + beta ^ 2 / sigma ^ 2))
    
    sapply(1:length(y), function(i) {
      mis = which(is.na(X[i, ])) 
      if (length(mis) > 0) {
        X_obs = X[i, -mis] 
        beta_mis = beta[mis]
        beta_obs = beta[-mis]
        m_i = (Big_Sigma_inv %*% mu)[mis]
        u_i = Big_Sigma_inv[mis, -mis] %*% X_obs
        r = (y[i] - X_obs %*% beta_obs)[1, 1]
        tau_i = tau[mis]
        
        el_1 = ((r * beta_mis) / sigma ^ 2 + m_i - u_i) / tau_i
        el_2 = ((beta_mis %*% t(beta_mis)) / sigma ^ 2 + Big_Sigma_inv[mis, mis]) /
          (tau_i %*% t(tau_i))
        diag(el_2) = 1 
        mu_tilde = solve(el_2, el_1)
        B_inv = solve(el_2)
        Z = mvrnorm(n = 1, mu = mu_tilde, Big_Sigma = B_inv)
        X[i, mis] = Z / tau
      }
    })
    if (scale) scale_X_mis(X, scale, row.w, std.w, mean.w)
  }
  X
}

create_old_list = function(beta, sigma, theta, mu, 
                           Big_Sigma, missingcols) {
  beta.old = beta 
  sigma.old = sigma
  theta.old = theta
  if (missingcols != 0) {
    mu.old = mu
    Big_Sigma.old = Big_Sigma
    list(beta.old, sigma.old, theta.old, mu.old, Big_Sigma.old)
  }
  else
    list(beta.old, sigma.old, theta.old, NULL, NULL)
}

#PSOBCZYK function
slope_admm <- function(A, b, z, u, lambda_seq, rho,
                       max_iter = 100, tol_infeas = 1e-3,
                       verbose = FALSE) {
  M <- solve(crossprod(A) + diag(rho, ncol(A)))
  MtAb <- M %*% crossprod(A,b)
  lambda_seq_rho <- lambda_seq/rho
  z_new <- NULL
  for(iter in 1:max_iter){ #just until we do not choose some reasonable convergence criterion
    
    x <- MtAb + crossprod(M, (rho*(z - u)))
    z_new <- SLOPE::prox_sorted_L1(x = as.vector(x + u), lambda = lambda_seq_rho)
    u <- u + x - z_new
    
    dual_feasibility <- norm(rho*(z_new-z), type = "2")
    primal_feasibility <- norm((z_new - x), type = "2")
    
    z <- z_new
    
    if(verbose)
      message(sprintf("Iter %i\nprimal: %f\ndual: %f\n",
                      iter, primal_feasibility, dual_feasibility))
    
    if(dual_feasibility < tol_infeas & primal_feasibility < tol_infeas){
      break;
    }
  }
  return(list(x = x, z = z, u = u,
              primal_feasibility = primal_feasibility,
              dual_feasibility = dual_feasibility))
}

estimate_beta_ML = function(gamma, c, X, y, lambda_sigma) {
  p = ncol(X)
  W = gamma * c + (rep(1, p) - gamma)
  revW = 1 / W
  Xtemp = sweep(X, 2, revW, '*');
  z = slope_admm(Xtemp, y, rep(0, p), rep(0, p), lambda_seq = lambda_sigma, 1)$z
  revW * z
}

estimate_sigma_ML = function(X, y, beta, 
                             lambda, lambda_sigma,
                             sigma.known, sigma,
                             gamma, c) {
  if (is.na(sigma.known)) {
    n = nrow(X)
    p = ncol(X)
    W = gamma * c + (rep(1, p) - gamma)
    revW = 1 / W
    Xtemp = sweep(X, 2, revW, '*');
    z = slope_admm(Xtemp, y, rep(0, p), rep(0, p), lambda_seq = lambda_sigma, 1)$z
    rk = p - rank(abs(z), ties.method = "max") + 1
    RSS = sum((y - X %*% beta) ^ 2)
    sum_lamwbeta = sum(lambda[rk] * abs(z))
    sigma = (sqrt(sum_lamwbeta ^ 2 + 4 * n * RSS) + sum_lamwbeta) / (2 * n)
  }
  sigma
}

estimate_theta_ML = function(gamma, a, b, X) {
  p = ncol(X)
  (sum(gamma == 1) + a) / (a + b + p)
}

estimate_mu_Big_Sigma = function(X, mu.old, Big_Sigma.old) {
  missingcols = calculate_missing_cols(X)
  if (missingcols != 0) {
    list(apply(X, 2, mean), linshrink_cov(X))
  }
  else
    list(mu.old, Big_Sigma.old)
}

update_params = function(beta.old, beta, eta,
                         sigma.old, sigma,
                         theta.old, theta,
                         mu.old, mu,
                         Big_Sigma.old, Big_Sigma,
                         X) {
  missingcols = calculate_missing_cols(X)
  if (eta != 1) {
    beta = beta.old + eta * (beta - beta.old)  
    sigma = sigma.old + eta * (sigma - sigma.old)
    theta = theta.old + eta * (theta - theta.old)
    if (missingcols != 0) { 
      mu = mu.old + eta*(mu - mu.old)
      Big_Sigma = Big_Sigma.old + eta * (Big_Sigma - Big_Sigma.old)
    }
  }
  list(beta, sigma, theta, mu, Big_Sigma)
  
}

scale_X = function(X, scale) {
  n = nrow(X)
  row.w = rep(1, nrow(X)) / nrow(X)
  mean.w = apply(X, 2, scale_mean.w, row.w)
  X.sim = t(t(X) - mean.w)
  std.w = apply(X.sim, 2, scale_std.w, row.w)*sqrt(n)
  if (scale) {
    X.sim = t(t(X.sim) / std.w)
  }
  else X.sim = X
  list(X.sim, row.w, mean.w, std.w)
}

iterate_saem_algorithm = function(cstop = 1, 
                                  tol_em,
                                  initialization_list, 
                                  estimations_cache, 
                                  maxit,
                                  print_iter,
                                  X, y,
                                  lambda,
                                  sigma.known,
                                  scale,
                                  a, b) {
  t = 0
  beta = initialization_list[[1]] 
  gamma = initialization_list[[2]] 
  sigma = initialization_list[[3]]
  theta = initialization_list[[4]]
  c = initialization_list[[5]]
  Big_Sigma = initialization_list[[6]]
  rank = initialization_list[[7]]
  mu = initialization_list[[8]]
  eta = 0
  lambdas = calculate_lambdas(lambda, sigma)
  lambda_sigma = lambdas[[1]]
  lambda_sigma_inv = lambdas[[2]]
  scaling = scale_X(X, scale)
  X = scaling[[1]]
  row.w = scaling[[2]]
  mean.w = scaling[[3]]
  std.w = scaling[[4]]
  while (cstop > tol_em & t < maxit | t < 20) { #why exactly this condition
    t = t + 1
    
    if ((print_iter == TRUE) & (t %% 100 == 0)) {
      cat(sprintf("Iteration no."), t, "\n",
          "beta = ", beta, ", 
                sigma = ", sigma, "\n",
          "Distance from last iter = ", cstop)
    }
    
    eta = update_step_size(t)
    binomial_prob =  generate_binomial_prob(beta, theta, c, lambda_sigma_inv, rank)
    gamma = generate_gamma_t(ncol(X), binomial_prob)
    c = generate_c_t(beta, gamma, lambda_sigma_inv, rank)
    X = generate_X(Big_Sigma, mu, beta, sigma, 
                   y, X, row.w, std.w, mean.w)
    old_list = create_old_list(beta, sigma, theta, 
                               mu, Big_Sigma, 
                               calculate_missing_cols(X))
    beta.old = old_list[1]
    sigma.old = old_list[2]
    theta.old = old_list[3]
    mu.old = old_list[4]
    Big_Sigma.old = old_list[5]
    
    beta = estimate_beta_ML(gamma, c, X, y, lambda_sigma)
    
    cstop = sum((beta - beta.old) ^ 2)
    
    sigma = estimate_sigma_ML(X, y, beta, 
                              lambda, lambda_sigma,
                              sigma.known, sigma)
    
    lambdas = calculate_lambdas(lambda, sigma)
    lambda_sigma = lambdas[[1]]
    lambda_sigma_inv = lambdas[[2]]
    
    theta = estimate_theta_ML(gamma, a, b, X)
    
    mu.Big_Sigma = estimate_mu_Big_Sigma(X, mu.old, Big_Sigma.old)
    mu = mu.Big_Sigma[1]
    Big_Sigma = mu.Big_Sigma[2]
    
    updated_params = update_params(beta.old, beta, eta,
                                   sigma.old, sigma,
                                   theta.old, theta,
                                   mu.old, mu,
                                   Big_Sigma.old, Big_Sigma,
                                   X)
    beta = updated_params[[1]]
    sigma = updated_params[[2]]
    theta = updated_params[[3]]
    mu = updated_params[[4]]
    Big_Sigma = updated_params[[5]]
    
    estimations_cache[1][, t] = beta
    estimations_cache[2][, t] = sigma
    estimations_cache[3][, t] = gamma
    estimations_cache[4][, t] = theta
    estimations_cache[5][, t] = c
  }
  list(estimations_cache, X)
}
