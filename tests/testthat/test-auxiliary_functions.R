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


