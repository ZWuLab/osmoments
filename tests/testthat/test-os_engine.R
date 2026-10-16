# ==============================================================================
# File: test-os_engine.R
# Purpose: The package's normal moments against the reference engine, an
#          earlier implementation of the same method for the normal written
#          separately before the package, on 39 sample sizes from 5 to 1001 at
#          the five-number ranks (fixtures/01-frozen-grid.csv), the
#          user-distribution path against the built-in normal, and the
#          cancellation-free difference.

test_that("the normal at the five-number ranks agrees with the reference engine on the grid", {
  skip_on_cran()                     # about 15 s
  fx <- read.csv(test_path("fixtures", "01-frozen-grid.csv"))
  expect_equal(length(unique(fx$n)), 39L)
  worst_alpha <- 0; worst_B <- 0
  for (n in unique(fx$n)) {
    f <- fx[fx$n == n, ]
    rk <- f$rank_i[f$what == "alpha"]
    m <- os_moments(n, rk, "norm")
    worst_alpha <- max(worst_alpha, abs(m$alpha - f$value[f$what == "alpha"]))
    fb <- f[f$what == "B", ]
    worst_B <- max(worst_B, abs(m$B[cbind(fb$i, fb$j)] - fb$value))
  }
  expect_lt(worst_alpha, 1e-10)
  expect_lt(worst_B, 1e-10)
})

test_that("a normal supplied through os_dist() reproduces the built-in normal", {
  user <- os_dist(pnorm, dnorm, qnorm, symmetric = TRUE, name = "user-normal")
  a <- os_moments(20, c(1, 5, 10, 20), user)
  b <- os_moments(20, c(1, 5, 10, 20), "norm")
  expect_lt(max(abs(a$alpha - b$alpha)), 1e-12)
  expect_lt(max(abs(a$B - b$B)), 1e-12)
  expect_identical(a$route, "integrate")
})

test_that("the cancellation-free difference agrees with the direct difference where that is accurate", {
  cd <- osmoments:::.os_centered(osmoments:::.os_family("norm"))
  y <- -1.2; z <- c(-0.9, -0.3, 0.4, 2.5)
  expect_equal(osmoments:::.os_logdiff(y, z, cd), log(pnorm(z) - pnorm(y)), tolerance = 1e-13)
  y <- 0.7; z <- c(0.8, 1.5, 3)
  expect_equal(osmoments:::.os_logdiff(y, z, cd), log(pnorm(z) - pnorm(y)), tolerance = 1e-13)
  # deep in the upper tail the direct difference is zero in doubles and the tail form is not
  expect_true(is.finite(osmoments:::.os_logdiff(40, 41, cd)))
  expect_identical(osmoments:::.os_logdiff(0.5, 0.5, cd), -Inf)
})

test_that("a divergent second moment stops even when the tail index is not declared", {
  # t with 1.9 degrees of freedom has no variance. Declared without its tail index, the
  # tail pieces of E[X^2] came back as -11.3 each, flagged probably divergent with an
  # error estimate below the acceptance threshold, and gave a variance of -19.
  t19 <- os_dist(function(q, lower.tail = TRUE, log.p = FALSE)
                   pt(q, 1.9, lower.tail = lower.tail, log.p = log.p),
                 function(x, log = FALSE) dt(x, 1.9, log = log),
                 function(p) qt(p, 1.9), symmetric = TRUE)
  expect_error(os_moments(1, dist = t19), "probably divergent")
})

test_that("a divergence flag on a bounded piece is spurious and the piece is accepted", {
  # the beta(1/2, 1/2) has an infinite density at both endpoints; integrate() flagged an
  # inner piece near 1 as probably divergent with an error estimate of 1e-11, and the
  # rule that stopped on any such flag refused the covariance of two observations, whose
  # exact value is 4/pi^4 (X_(1) X_(2) = X_1 X_2 and the beta moments)
  b <- os_dist(function(q, lower.tail = TRUE, log.p = FALSE)
                 pbeta(q, 0.5, 0.5, lower.tail = lower.tail, log.p = log.p),
               function(x, log = FALSE) dbeta(x, 0.5, 0.5, log = log),
               function(p) qbeta(p, 0.5, 0.5), support = c(0, 1), symmetric = TRUE)
  m <- os_moments(2, dist = b)
  expect_lt(abs(m$B[1, 2] - 4 / pi^4), 1e-9)
})
