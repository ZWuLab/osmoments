# ==============================================================================
# File: test-os_moments.R
# Purpose: The interface of os_moments(): argument validation, cov = FALSE, the
#          methods and the route, the symmetry halving, the location-scale
#          transform and its limit, the hand-derived normal closed forms at
#          n = 2 and 3, and print().

test_that("arguments are validated", {
  expect_error(os_moments(0), "'n'")
  expect_error(os_moments(2.5), "'n'")
  expect_error(os_moments(5, c(3, 1)), "increasing")
  expect_error(os_moments(5, c(1, 1)), "increasing")
  expect_error(os_moments(5, 6), "1, ..., n")
  expect_error(os_moments(5, 1.5), "integers")
  expect_error(os_moments(5, scale = 0), "'scale'")
  expect_error(os_moments(5, scale = -1), "'scale'")
  expect_error(os_moments(5, location = NA), "'location'")
  expect_error(os_moments(5, method = "exact"))
  expect_error(os_moments(5, dist = "cauchy"), "os_dist")
  expect_error(os_moments(5, dist = list()), "os_dist")
  expect_error(os_moments(5, cov = NA), "'cov'")
  expect_error(os_moments(5, rel.tol = 0), "'rel.tol'")
})

test_that("cov = FALSE gives the means and the variances by single integrals only", {
  a <- os_moments(10, c(1, 5, 10), "norm", cov = FALSE)
  b <- os_moments(10, c(1, 5, 10), "norm")
  expect_null(a$B)
  expect_false(a$cov)
  expect_lt(max(abs(a$alpha - b$alpha)), 1e-12)
  expect_lt(max(abs(a$var - diag(b$B))), 1e-12)
  expect_lt(max(abs(b$var - diag(b$B))), 1e-15)
  expect_output(print(a), "cov = FALSE")
})

test_that("the methods and the route", {
  expect_error(os_moments(5, 1:2, "norm", method = "closed"), "no closed form")
  expect_error(os_moments(5, 1:2, "logis", method = "closed"), "covariance")
  expect_error(os_moments(5, 1, "gumbel", method = "closed"), "mean")
  expect_identical(os_moments(5, 1:5, "unif")$route, "closed")
  expect_identical(os_moments(5, 1:5, "exp")$route, "closed")
  expect_identical(os_moments(5, 1:5, "norm")$route, "integrate")
  expect_identical(os_moments(5, 1:5, "logis")$route, "mixed")
  expect_identical(os_moments(5, 1:5, "unif", method = "integrate")$route, "integrate")
  expect_identical(os_moments(5, 5, "gumbel", method = "closed")$route, "closed")
})

test_that("the symmetry halving gives the full computation at n = 13", {
  f <- osmoments:::.os_family("norm")
  full <- os_dist(f$p, f$d, f$q, symmetric = FALSE, name = "norm-full")
  a <- os_moments(13, 1:13, full)
  b <- os_moments(13, 1:13, "norm")
  expect_lt(max(abs(a$B - b$B)), 1e-12)
  expect_lt(max(abs(a$alpha - b$alpha)), 1e-14)
  # a rank set that is not closed under reflection is computed in full
  d <- os_moments(13, c(1, 2, 7), "norm")
  expect_true(is.na(d$checks$reflection))
  expect_lt(max(abs(d$B - a$B[c(1, 2, 7), c(1, 2, 7)])), 1e-12)
})

test_that("location and scale transform the standard member, at a location of 1e5 and 1e8", {
  rk <- c(1, 5, 10, 20)
  cl <- os_moments(20, rk, "unif")                                  # closed forms
  m <- os_moments(20, rk, "unif", location = 1e5, scale = 1, method = "integrate")
  expect_lt(max(abs((m$alpha - 1e5) - cl$alpha)), 1e-9)
  expect_lt(max(abs(m$B - cl$B)), 1e-10)
  m <- os_moments(20, rk, "unif", location = 1e8, scale = 1, method = "integrate")
  expect_lt(max(abs((m$alpha - 1e8) - cl$alpha)), 1e-7)            # the double keeps 8 decimals of 1e8 + x
  expect_lt(max(abs(m$B - cl$B)), 1e-10)
  s <- os_moments(20, rk, "norm", location = -3, scale = 2.5)
  z <- os_moments(20, rk, "norm")
  expect_equal(unname(s$alpha), unname(-3 + 2.5 * z$alpha))
  expect_equal(unname(s$B), unname(2.5^2 * z$B))
  # a user distribution whose own location is 1e8 defeats double precision and stops
  shifted <- os_dist(function(q, lower.tail = TRUE, log.p = FALSE)
                       punif(q, 1e8, 1e8 + 1, lower.tail = lower.tail, log.p = log.p),
                     function(x, log = FALSE) dunif(x, 1e8, 1e8 + 1, log = log),
                     function(p) qunif(p, 1e8, 1e8 + 1),
                     support = c(1e8, 1e8 + 1), symmetric = TRUE, name = "shifted")
  expect_error(os_moments(20, rk, shifted), "integral failed")
})

test_that("the normal closed forms at n = 2 and 3, and n = 1 returns the parent", {
  # E[Z_(2:2)] = 1/sqrt(pi); Cov(Z_(1:2), Z_(2:2)) = 1/pi; E[Z_(3:3)] = 3/(2 sqrt(pi))
  m2 <- os_moments(2)
  expect_lt(max(abs(m2$alpha - c(-1, 1) / sqrt(pi))), 1e-9)
  expect_lt(max(abs(m2$B - matrix(c(1 - 1 / pi, 1 / pi, 1 / pi, 1 - 1 / pi), 2))), 1e-9)
  m3 <- os_moments(3)
  expect_lt(max(abs(m3$alpha - c(-1.5, 0, 1.5) / sqrt(pi))), 1e-9)
  m1 <- os_moments(1, 1, "norm")
  expect_lt(abs(m1$alpha[[1]]), 1e-12)
  expect_lt(abs(m1$var[[1]] - 1), 1e-10)
  expect_equal(dim(m1$B), c(1L, 1L))
})

test_that("the object and its print method", {
  m <- os_moments(5, c(1, 3, 5), "logis", location = 2, scale = 3)
  expect_s3_class(m, "osmoments")
  expect_identical(names(m$alpha), c("r1", "r3", "r5"))
  expect_identical(dimnames(m$B), list(c("r1", "r3", "r5"), c("r1", "r3", "r5")))
  expect_identical(m$n, 5L)
  expect_identical(m$ranks, c(1L, 3L, 5L))
  expect_true(is.numeric(m$quad_err) && m$quad_err >= 0)
  expect_identical(names(m$standard$alpha), c("r1", "r3", "r5"))
  expect_output(print(m), "logis \\(location 2, scale 3\\), n = 5, ranks 1 3 5")
  expect_output(print(m), "route: mixed")
})

test_that("a distribution supplied on any scale is split around its own peak", {
  # A normal supplied with standard deviation 0.001 to 1e6: with the split floor in
  # absolute units the mass of a narrow one fell between the integrator's nodes and the
  # variance came out as exactly zero, with every run-time check satisfied, and with the
  # floor multiplied by the spread once too often a wide one failed the same way (the
  # second review round). In standardized units the moments are the standard normal's
  # scaled by sd and sd^2.
  for (sd in c(0.001, 0.01, 100, 1e4, 1e6)) {
    narrow <- os_dist(function(q, lower.tail = TRUE, log.p = FALSE)
                        pnorm(q, sd = sd, lower.tail = lower.tail, log.p = log.p),
                      function(x, log = FALSE) dnorm(x, sd = sd, log = log),
                      function(p) qnorm(p, sd = sd), symmetric = TRUE)
    for (n in c(1L, 5L, 25L)) {
      rk <- unique(c(1L, (n + 1L) %/% 2L, n))
      expect_silent(m <- os_moments(n, rk, narrow, cov = n > 1))   # the checks are scale-free
      z <- os_moments(n, rk, cov = n > 1)
      expect_lt(max(abs(m$alpha - sd * z$alpha)) / sd, 1e-12)
      expect_lt(max(abs(m$var / (sd^2 * z$var) - 1)), 1e-11)
      if (n > 1) expect_lt(max(abs(m$B / (sd^2 * z$B) - 1)), 1e-11)
    }
  }
  # the beta(1e6, 1e6), whose standard deviation is 3.5e-4, against its variance formula,
  # to the accuracy of pbeta() itself at these shapes (2.5e-11 relative when written)
  conc <- os_dist(function(q, lower.tail = TRUE, log.p = FALSE)
                    pbeta(q, 1e6, 1e6, lower.tail = lower.tail, log.p = log.p),
                  function(x, log = FALSE) dbeta(x, 1e6, 1e6, log = log),
                  function(p) qbeta(p, 1e6, 1e6), support = c(0, 1), symmetric = TRUE)
  expect_lt(abs(os_moments(1, dist = conc)$var[[1]] / (1 / (4 * (2e6 + 1))) - 1), 1e-9)
})

test_that("a variance that comes out zero or negative stops with an error", {
  # The guard is defensive: with each moment about its own location no built-in case
  # reaches it (the uniform minimum at n = 1000001 by quadrature did, at -6.6e-11 against
  # 1e-12, when the moments were taken about the median), so the single integrals are
  # mocked to return zero, which gives a variance of exactly zero.
  # (the mass, pow = 0, is mocked to 1 so that the mass check passes and the guard is reached)
  local_mocked_bindings(.os_single = function(r, n, pow, ...) list(value = if (pow == 0) 1 else 0, abs.error = 0),
                        .package = "osmoments")
  expect_error(os_moments(5, 3, "norm"), "variance of X_\\(3\\) came out as 0")
})

test_that("the integration route keeps its digits at large n, against the closed forms", {
  # Each moment is taken about its own location c_r = Q(r/(n+1)); about the median, the
  # variance of the uniform minimum was off by 1.3e-4 at n = 1e4 and by 0.19 at n = 1e5,
  # and the covariance of the extremes by a factor of three at n = 1e4. The two extremes
  # are nearly independent (correlation of order 1/n), so their covariance keeps fewer
  # digits (8.4e-5 relative for the uniform at n = 1e5 when written) and is held to 1e-3
  # at n = 1e5.
  for (fam in c("unif", "exp")) for (n in c(10000L, 100000L, 1000000L)) {
    rk <- c(1L, round((n + 1) / 4), round((n + 1) / 2), n)
    cl <- os_moments(n, rk, fam)                                          # closed forms
    qu <- os_moments(n, rk, fam, method = "integrate", check = FALSE)
    expect_identical(qu$route, "integrate")
    expect_lt(max(abs(qu$var / cl$var - 1)), 1e-9)
    expect_lt(max(abs(qu$B[1, 2] / cl$B[1, 2] - 1), abs(qu$B[2, 3] / cl$B[2, 3] - 1)), 1e-9)
    expect_lt(abs(qu$B[1, 4] / cl$B[1, 4] - 1), if (n > 10000L) 1e-3 else 1e-9)
  }
})

test_that("a heavy-tailed parent's covariance keeps its far tail: the lognormal at n = 2", {
  # With sdlog = 2 the outer nodes reach y = 1e5, where an inner integral resolved at the
  # fixed absolute tolerance 1e-15 lost 3.5e-6 of the covariance 38.77; scaled by the
  # outer factor, 5.1e-8 is left (the exact value from X_(1) X_(2) = X_1 X_2).
  ln2 <- os_dist(function(q, lower.tail = TRUE, log.p = FALSE)
                   plnorm(q, sdlog = 2, lower.tail = lower.tail, log.p = log.p),
                 function(x, log = FALSE) dlnorm(x, sdlog = 2, log = log),
                 function(p) qlnorm(p, sdlog = 2), support = c(0, Inf))
  mu <- 2 * exp(2) * pnorm(c(-1, 1) * sqrt(2))                         # E[X_(1)], E[X_(2)]
  m <- suppressWarnings(os_moments(2, dist = ln2))
  expect_lt(max(abs(unname(m$alpha) - mu)), 1e-12)
  expect_lt(abs(m$B[1, 2] - (exp(4) - prod(mu))), 2e-7)
})

test_that("the object carries the standard member's moments", {
  z <- os_moments(9, c(1, 5, 9), "logis")
  m <- os_moments(9, c(1, 5, 9), "logis", location = 1e15, scale = 2)
  expect_identical(names(m$standard), c("alpha", "var", "B"))
  expect_identical(m$standard$alpha, z$alpha)                          # the same integrals
  expect_identical(m$standard$B, z$B)
  expect_equal(m$alpha, 1e15 + 2 * z$alpha)
  expect_equal(m$B, 4 * z$B)
  expect_null(os_moments(9, c(1, 5, 9), "logis", cov = FALSE)$standard$B)
})

test_that("the scale is applied one factor at a time, and an overflow stops", {
  z <- os_moments(1001, 501, "unif")
  m <- os_moments(1001, 501, "unif", scale = 1e155)                    # scale^2 overflows
  expect_true(is.finite(m$var[[1]]))
  expect_equal(m$var[[1]], (z$var[[1]] * 1e155) * 1e155, tolerance = 1e-14)
  expect_gt(os_moments(1, dist = "logis", scale = 1e-162)$var[[1]], 0) # scale^2 underflows
  expect_error(os_moments(1001, 501, "unif", scale = 1e200), "overflow")
})

test_that("the normal row-sum identity follows the functions, not the label", {
  # an exponential labeled "norm" gets the general identities only
  d <- os_dist(pexp, dexp, qexp, support = c(0, Inf), name = "norm")
  expect_silent(m <- os_moments(3, dist = d))
  expect_false("row_sums" %in% names(m$checks$identities))
  # R's normal functions under another label get the row sums
  z <- os_dist(pnorm, dnorm, qnorm, symmetric = TRUE, name = "z")
  expect_true("row_sums" %in% names(os_moments(3, dist = z)$checks$identities))
})

test_that("n is bounded by the integer range", {
  expect_error(os_moments(2^31, 1), "at most")
})

test_that("a parent with a second mode is refused by the mass check", {
  # 0.9 N(0, 0.01) + 0.1 N(10, 0.01): the mode at 10 fell between the nodes of the tail
  # piece at n = 1, and the moments came back as those of the first mode (mean 0.0014
  # instead of 1) with every identity satisfied, because the parent moments are computed
  # by the same quadrature. The density of X_(1) then integrates to 0.9.
  s <- 0.1
  log_add <- function(a, b) { hi <- pmax(a, b); out <- hi + log1p(exp(-abs(a - b))); out[is.infinite(hi) & hi < 0] <- -Inf; out }
  pp <- function(q, lower.tail = TRUE, log.p = FALSE) {
    z <- log_add(log(0.9) + pnorm(q, 0, s, lower.tail = lower.tail, log.p = TRUE),
                 log(0.1) + pnorm(q, 10, s, lower.tail = lower.tail, log.p = TRUE))
    if (log.p) z else exp(z)
  }
  dd <- function(x, log = FALSE) {
    z <- log_add(log(0.9) + dnorm(x, 0, s, log = TRUE), log(0.1) + dnorm(x, 10, s, log = TRUE))
    if (log) z else exp(z)
  }
  qq <- function(p) vapply(p, function(u) uniroot(function(x) pp(x) - u, c(-1, 11), tol = 1e-13)$root, numeric(1))
  mix <- os_dist(pp, dd, qq, name = "mixture")
  expect_error(os_moments(1, dist = mix), "integrates to 0.9")
})

test_that("the extremes of a parent with an infinite endpoint density: computed, then refused", {
  # the minimum of n observations from a beta(1/2, 1/2), whose density is infinite at 0,
  # concentrates against the endpoint on the scale n^-2. The reference is the quadrature
  # in the probability scale, Q(u) = sin(pi u/2)^2 and U_(1) ~ n (1 - u)^(n - 1), v = n u.
  # At n = 1e5 the engine of the morning returned a variance 99 percent too small without
  # a word (the second review round); now the extremes are computed to n = 100 and stop
  # with an error from n = 1000 on, where the mass check or a bounded piece near the
  # singular endpoint fails.
  b <- os_dist(function(q, lower.tail = TRUE, log.p = FALSE)
                 pbeta(q, 0.5, 0.5, lower.tail = lower.tail, log.p = log.p),
               function(x, log = FALSE) dbeta(x, 0.5, 0.5, log = log),
               function(p) qbeta(p, 0.5, 0.5), support = c(0, 1), symmetric = TRUE)
  for (n in c(10L, 100L)) {
    a <- vapply(1:2, function(k) integrate(function(v) (n^2 * sin(pi * v / (2 * n))^2)^k *
                                              exp((n - 1) * log1p(-v / n)), 0, n, rel.tol = 1e-12)$value,
                numeric(1))
    m <- os_moments(n, c(1L, n), dist = b, cov = FALSE)
    expect_lt(abs(m$alpha[[1]] / (a[1] / n^2) - 1), 1e-10)
    expect_lt(abs(m$var[[1]] / ((a[2] - a[1]^2) / n^4) - 1), 1e-10)
  }
  expect_error(os_moments(1000L, c(1L, 1000L), dist = b, cov = FALSE), "os_moments")
  expect_error(os_moments(100000L, c(1L, 100000L), dist = b, cov = FALSE), "os_moments")
})

test_that("the spread is applied one factor at a time", {
  # an exponential supplied with scale 2e-162: sp^2 underflows to 0 while sp * (sp * v)
  # reaches the subnormal range
  s <- 2e-162
  d <- os_dist(function(q, lower.tail = TRUE, log.p = FALSE) pexp(q, rate = 1 / s, lower.tail = lower.tail, log.p = log.p),
               function(x, log = FALSE) dexp(x, rate = 1 / s, log = log),
               function(p) qexp(p, rate = 1 / s), support = c(0, Inf))
  # (the second moment of this parent is subnormal, so the identity residuals are noise:
  # the run-time checks are off)
  expect_equal(os_moments(1, dist = d, check = FALSE)$var[[1]], s * s)
})
