test_that("Gamma generating", {
  expect_lt(abs(trunc_norm_gamma(1,1) - 0.63212), 0.01)
  expect_lt(abs(trunc_norm_gamma(1,2) - 0.43233), 0.01)
  expect_lt(abs(trunc_norm_gamma(2,1) - 0.26424), 0.01)
})

X <- matrix(rnorm(200), ncol=10)
b <- c(5, -5, 4, -4, 2, rep(0, 5))
y <- X %*% b + rnorm(20, sd=0.1)

new_gamma <- slobe_update_gamma(rep(1, 10), 0.3, b, 10, 0.1, rep(1,10), 0.1)
new_c <- slobe_update_c(rep(1,10), 0.3, b, 0.1, rep(1,10))
new_theta <- slobe_update_theta(rep(1:10), 1, 1, 10)

test_that("Updating params", {
  expect_lt(abs(sum(new_gamma[1:5]) - 5), 0.1)
  expect_lt(sum(abs(new_gamma[6:10])), 0.5)
  expect_lt(abs(new_c - 0.055), 0.01)
  expect_lt(abs(new_theta - 1/6), 0.01)
})


X = matrix(rnorm(1000), nrow=100)
b = c(sample(-5:5, 5), rep(0, 5))
y = X %*% b + rnorm(100, 0, 0.1)
A <- SLOBE(X, y, lambda=seq(10, 5, length.out=10))

test_that("Model check", {
  expect_lt(sum(abs(A$beta - b)), 0.5)
  expect_true(all(A$selected == which(b != 0)))
})