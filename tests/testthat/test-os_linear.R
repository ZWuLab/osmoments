# ==============================================================================
# File: test-os_linear.R
# Purpose: The linear layer: the quantile weights against quantile() for all
#          nine types, the moments of linear combinations against the entrywise
#          expansion, and Lloyd's BLUE against the reference engine, an earlier
#          implementation of the same method for the normal
#          (fixtures/03-frozen-blue.csv), the generalized least squares formula
#          and its structural identities.

test_that("os_quantile_weights() reproduces quantile() on sorted samples for all nine types", {
  set.seed(20261013)
  p <- c(0, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1)
  worst <- 0
  for (n in c(1:12, 100)) {
    x <- sort(rnorm(n))
    for (type in 1:9) {
      qw <- os_quantile_weights(p, n, type)
      expect_true(all(qw$ranks >= 1 & qw$ranks <= n))
      expect_equal(unname(colSums(qw$W)), rep(1, length(p)))
      got <- colSums(qw$W * x[qw$ranks])
      worst <- max(worst, abs(got - quantile(x, p, type = type, names = FALSE)))
    }
  }
  expect_lt(worst, 1e-12)
})

test_that("os_quantile_weights(): the even-sample median, the endpoints, the names, the errors", {
  qw <- os_quantile_weights(0.5, 10)
  expect_identical(qw$ranks, c(5L, 6L))
  expect_equal(unname(qw$W[, 1]), c(0.5, 0.5))
  expect_identical(colnames(qw$W), "50%")
  expect_identical(rownames(qw$W), c("r5", "r6"))
  qw <- os_quantile_weights(c(0, 1), 7, type = 6)
  expect_identical(qw$ranks, c(1L, 7L))
  expect_equal(unname(qw$W), diag(2))
  qw <- os_quantile_weights(c(0.25, 0.5, 0.75), 9)           # exact ranks at n = 4Q + 1
  expect_identical(qw$ranks, c(3L, 5L, 7L))
  expect_equal(unname(qw$W), diag(3))
  expect_error(os_quantile_weights(1.5, 10), "'p'")
  expect_error(os_quantile_weights(0.5, 0), "'n'")
  expect_error(os_quantile_weights(0.5, 10, type = 10), "'type'")
})

test_that("os_lin_moments() equals the entrywise expansion and handles several statistics", {
  rk <- c(1, 7, 13, 19, 25)
  m <- os_moments(25, rk)
  w <- c(0.10, 0.20, 0.30, 0.25, 0.15)
  lm <- os_lin_moments(w, m)
  M <- m$B + tcrossprod(m$alpha)                                   # E[X_r X_s]
  expect_lt(abs(lm$var[[1]] - (sum(outer(w, w) * M) - sum(w * m$alpha)^2)), 1e-12)
  expect_lt(abs(lm$mean[[1]] - sum(w * m$alpha)), 1e-14)
  expect_equal(dim(lm$cov), c(1L, 1L))
  W <- cbind(mean = rep(0.2, 5), range = c(-1, 0, 0, 0, 1))
  lm2 <- os_lin_moments(W, m)
  expect_identical(names(lm2$mean), c("mean", "range"))
  expect_lt(abs(lm2$mean[["mean"]]), 1e-10)                        # symmetric ranks
  expect_lt(abs(lm2$cov["mean", "range"]), 1e-10)                  # cov of a symmetric and an antisymmetric combination
  expect_equal(lm2$cov["range", "range"], m$B[1, 1] + m$B[5, 5] - 2 * m$B[1, 5])
  # the even-sample median of a normal sample of size 10
  qw <- os_quantile_weights(0.5, 10)
  m10 <- os_moments(10, qw$ranks)
  med <- os_lin_moments(qw$W, m10)
  expect_lt(abs(med$mean[[1]]), 1e-10)
  expect_equal(unname(med$var[[1]]), (m10$B[1, 1] + m10$B[2, 2] + 2 * m10$B[1, 2]) / 4)
  expect_error(os_lin_moments(c(1, 2), m), "one weight per rank")
  expect_error(os_lin_moments(w, os_moments(25, rk, cov = FALSE)), "cov = FALSE")
  expect_error(os_lin_moments(w, list()), "os_moments")
})

test_that("os_blue() reproduces the reference engine's BLUE for the five-number scenarios", {
  fx <- read.csv(test_path("fixtures", "03-frozen-blue.csv"))
  worst_w <- 0; worst_v <- 0
  for (n in unique(fx$n)) for (sc in unique(fx$scenario)) {
    f <- fx[fx$n == n & fx$scenario == sc, ]
    b <- os_blue(os_moments(n, f$rank))
    worst_w <- max(worst_w, abs(b$w_location - f$w_mu), abs(b$w_scale - f$w_sig))
    worst_v <- max(worst_v, abs(b$var_location - f$var_mean[1]), abs(b$var_scale - f$var_sd[1]),
                   abs(b$cov - f$cov[1]))
  }
  expect_lt(worst_w, 1e-10)
  expect_lt(worst_v, 1e-10)
})

test_that("os_blue() is the generalized least squares estimator, with its identities", {
  for (n in c(25, 41)) for (rk in list(c(1, (n + 1) / 2, n), c(1, 5, 9, (n + 1) / 2, n - 2, n))) {
    m <- os_moments(n, rk)
    b <- os_blue(m)
    A <- cbind(1, m$alpha); Om <- solve(m$B)
    G <- solve(crossprod(A, Om %*% A), crossprod(A, Om))             # rows: the two weight vectors
    expect_lt(max(abs(b$w_location - G[1, ])), 1e-12)
    expect_lt(max(abs(b$w_scale - G[2, ])), 1e-12)
    expect_lt(abs(sum(b$w_location) - 1), 1e-12)                    # unbiased for location
    expect_lt(abs(sum(b$w_location * m$alpha)), 1e-12)
    expect_lt(abs(sum(b$w_scale)), 1e-12)                            # unbiased for scale
    expect_lt(abs(sum(b$w_scale * m$alpha) - 1), 1e-12)
    V <- solve(crossprod(A, Om %*% A))
    expect_lt(abs(b$var_location - V[1, 1]), 1e-12)
    expect_lt(abs(b$var_scale - V[2, 2]), 1e-12)
    expect_lt(abs(b$cov - V[1, 2]), 1e-12)
    expect_lt(abs(b$var_location - os_lin_moments(b$w_location, m)$var[[1]]), 1e-12)
  }
  # three symmetric ranks: the scale BLUE is the range over its expectation, and the
  # location and scale estimators are uncorrelated
  m <- os_moments(41, c(1, 21, 41))
  b <- os_blue(m)
  expect_lt(max(abs(b$w_scale - c(-1, 0, 1) / (2 * m$alpha[[3]]))), 1e-7)
  expect_lt(abs(b$cov), 1e-8)
})

test_that("os_blue() does not depend on location and scale, and refuses what it cannot do", {
  z <- os_blue(os_moments(13, c(1, 4, 7, 10, 13)))
  s <- os_blue(os_moments(13, c(1, 4, 7, 10, 13), location = -2, scale = 3))
  expect_lt(max(abs(z$w_location - s$w_location)), 1e-12)
  expect_lt(max(abs(z$w_scale - s$w_scale)), 1e-12)
  expect_lt(abs(z$var_scale - s$var_scale), 1e-12)
  # at a location of 1e15 the subtraction (alpha - location)/scale lost the means' digits
  # (the scale weights' expectation came out as 0.994); the standard member's moments
  # kept in the object give the same weights as at location 0
  big <- os_blue(os_moments(13, c(1, 4, 7, 10, 13), location = 1e15))
  expect_identical(big$w_location, z$w_location)
  expect_identical(big$w_scale, z$w_scale)
  e <- os_blue(os_moments(10, c(1, 10), "exp"))                     # the scale is 1/rate
  expect_lt(abs(sum(e$w_scale * (os_moments(10, c(1, 10), "exp")$alpha)) - 1), 1e-12)
  expect_error(os_blue(os_moments(10, 5)), "two distinct ranks")
  expect_error(os_blue(os_moments(10, c(1, 10), cov = FALSE)), "cov = FALSE")
  expect_error(os_blue(42), "os_moments")
})

test_that("os_blue() keeps the unbiasedness conditions when the means are nearly collinear", {
  # Two adjacent uniform ranks at n = 1e8: 1 and alpha differ by 1e-8 in relative terms,
  # and the determinant formula returned weights summing to 0.555 instead of 1. The
  # exact weights follow from the two conditions: (n, -(n - 1)) and (-(n + 1), n + 1).
  n <- 1e8
  m <- os_moments(n, c(n - 1, n), "unif")
  expect_silent(b <- os_blue(m))
  a <- unname(m$alpha)
  expect_lt(max(abs(c(sum(b$w_location) - 1, sum(b$w_location * a),
                      sum(b$w_scale), sum(b$w_scale * a) - 1))), 1e-6)
  expect_lt(max(abs(unname(b$w_location) / c(n, -(n - 1)) - 1)), 1e-6)
  expect_lt(max(abs(unname(b$w_scale) / c(-(n + 1), n + 1) - 1)), 1e-6)
})

test_that("the quantile labels are those of quantile(), also at 100 probabilities", {
  for (p in list(c(0.1, 0.25, 0.5), seq(0, 1, length.out = 100)))
    expect_identical(colnames(os_quantile_weights(p, 10)$W), names(quantile(1:10, p)))
})

test_that("os_lin_moments() requires finite weights", {
  m <- os_moments(5, c(1, 5))
  expect_error(os_lin_moments(c(Inf, 1), m), "finite")
  expect_error(os_lin_moments(c(NA, 1), m), "finite")
})

test_that("os_lin_moments() keeps a location-free statistic exact at a large location", {
  # the range of two uniform observations at location 1e16 came out as 0 instead of 1/3
  # when the means were multiplied after the location was added
  m <- os_moments(2, dist = "unif", location = 1e16)
  r <- os_lin_moments(c(-1, 1), m)
  expect_equal(r$mean, 1 / 3)
  expect_equal(r$var, 1 / 18)
  W <- cbind(range = c(-1, 1), sum = c(1, 1))
  r2 <- os_lin_moments(W, m)
  expect_equal(unname(r2$mean), c(1 / 3, 2e16 + 1))
})

test_that("os_quantile_weights() at the largest integer n", {
  n <- .Machine$integer.max
  for (type in c(2L, 6L, 7L)) {
    expect_silent(qw <- os_quantile_weights(1, n, type = type))
    expect_identical(qw$ranks, n)
    expect_equal(unname(qw$W[, 1]), 1)
  }
  expect_error(os_quantile_weights(0.5, n + 1), "at most")
})
