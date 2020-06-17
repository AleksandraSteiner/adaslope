X_mis0 <- matrix(rnorm(100), nrow=5)
X_mis1 <- matrix(rnorm(100), nrow=5)
X_mis3 <- matrix(rnorm(100), nrow=5)
X_misall <- matrix(rnorm(100), nrow=5)

X_mis1[1,1] <- NA
X_mis3[1,1:3] <- NA
X_misall[1,] <- NA

test_that("Missing columns detection", {
  expect_equal(calculate_missing_cols(X_mis0), 0)
  expect_equal(calculate_missing_cols(X_mis1), 1)
  expect_equal(calculate_missing_cols(X_mis3), 3)
  expect_equal(calculate_missing_cols(X_misall), ncol(X_misall))
})

X <- matrix(rnorm(200), ncol=10)
b <- c(5, -5, 4, -4, 2, rep(0, 5))
y <- X %*% b + rnorm(20, sd=0.1)

beta_testing <- initialize_beta(NA, X, y)
sigma_testing <- initialize_sigma(NA, NA, b, X, y)

test_that("Initializations", {
  expect_equal(initialize_beta(1:10, X, y), 1:10)
  expect_equal(length(beta_testing), 10)
  expect_lt(sum(abs(beta_testing[6:10])), 2)
  
  expect_equal(initialize_sigma(1.5, NA, b, X, y), 1.5)
  expect_equal(initialize_sigma(NA, 1.4, b, X, y), 1.4)
  expect_lt(abs(sigma_testing), 0.5)
  
  expect_equal(initialize_theta(b, 1, 1, X), 6/12)
  expect_equal(initialize_rank(b, X), c(1,1,3,3,5,6,6,6,6,6))
  
  expect_equal(initialize_c(b, rep(1,10), X), 0.3)
  expect_equal(initialize_c(b, 1/(1:10), X), 1)
  
  expect_null(initialize_mu_BSigma(X)[[1]])
  expect_null(initialize_mu_BSigma(X)[[2]])
})

beta_ML <- estimate_beta_ML(rep(1:10), 0.3, X, y, rep(1,10))
sigma_ML <- estimate_sigma_ML(X, y, b, rep(1, 10), rep(1, 10), NA, NA, rep(1:10), 1)
theta_ML <- estimate_theta_ML(rep(1,10), 1, 1, X)
test_that("ML maximizing", {
  expect_lt(sum(abs(beta_ML - b)), 0.6)
  expect_lt(abs(sigma_ML - 1), 0.1)
  expect_lt(abs(theta_ML - 1), 0.2)
})

A <- ABSLOPE(X, y, lambda=(10:1)/3, a=1, b=1, maxit=20)
test_that("Model check", {
  expect_lt(sum(abs(A$beta - b)), 0.5)
  expect_true(all(A$selected ==c(1,2,3,4,5)))
})

