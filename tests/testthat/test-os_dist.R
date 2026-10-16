# ==============================================================================
# File: test-os_dist.R
# Purpose: The distribution object: the validation of os_dist(), the built-in
#          families' centers and symmetry flags, and the Gumbel functions'
#          log-tails.

test_that("os_dist() validates its arguments and the functions' behavior at the median", {
  expect_error(os_dist(1, dnorm, qnorm), "functions")
  expect_error(os_dist(pnorm, dnorm, qnorm, support = c(1, -1)), "support")
  expect_error(os_dist(pnorm, dnorm, qnorm, symmetric = NA), "symmetric")
  expect_error(os_dist(pnorm, dnorm, qnorm, tail = c(0, 1)), "tail")
  expect_error(os_dist(pnorm, dnorm, qnorm, name = 1), "name")
  expect_error(os_dist(pnorm, dnorm, qexp), "disagree")                     # p and q of different laws
  expect_error(os_dist(function(q, lower.tail = TRUE, log.p = FALSE) pnorm(q), dnorm, qnorm),
               "lower.tail")                                                # ignores both arguments
  expect_error(os_dist(function(q, lower.tail = TRUE, log.p = FALSE) pnorm(q, lower.tail = lower.tail),
                       dnorm, qnorm), "log.p")                              # ignores log.p
  expect_error(os_dist(pnorm, function(x, log = FALSE) dnorm(x), qnorm), "log = TRUE")
  expect_error(os_dist(pnorm, dnorm, qnorm, support = c(1, 2)), "outside")
  d <- os_dist(pnorm, dnorm, qnorm, symmetric = TRUE, name = "z")
  expect_s3_class(d, "os_dist")
  expect_equal(d$center, 0)
  expect_identical(d$tail, c(Inf, Inf))
  expect_null(d$closed$mean)
})

test_that("the built-in families have the right centers, supports and symmetry", {
  f <- osmoments:::.os_family
  expect_equal(f("norm")$center, 0); expect_true(f("norm")$symmetric)
  expect_equal(f("unif")$center, 0.5); expect_true(f("unif")$symmetric)
  expect_identical(f("unif")$support, c(0, 1))
  expect_equal(f("exp")$center, log(2)); expect_false(f("exp")$symmetric)
  expect_identical(f("exp")$support, c(0, Inf))
  expect_equal(f("logis")$center, 0); expect_true(f("logis")$symmetric)
  expect_equal(f("gumbel")$center, -log(log(2))); expect_false(f("gumbel")$symmetric)
  expect_error(f("cauchy"))
})

test_that("the Gumbel functions keep their log-tails finite and accurate", {
  p <- osmoments:::.pgumbel; d <- osmoments:::.dgumbel; q <- osmoments:::.qgumbel
  x <- seq(-5, 700, by = 0.5)
  expect_true(all(is.finite(p(x, log.p = TRUE))))
  expect_true(all(is.finite(p(x, lower.tail = FALSE, log.p = TRUE))))
  # F + S = 1 in the center
  xc <- seq(-3, 3, by = 0.25)
  expect_lt(max(abs(p(xc) + p(xc, lower.tail = FALSE) - 1)), 1e-15)
  # the two branches of log S agree where both are accurate, and log S ~ -x far out
  xb <- c(29, 29.5, 29.9, 30.1, 30.5, 31)
  expect_lt(max(abs(p(xb, lower.tail = FALSE, log.p = TRUE) - (-xb + log1p(-exp(-xb) / 2)))), 1e-14)
  xf <- c(40, 100, 500)
  expect_lt(max(abs(p(xf, lower.tail = FALSE, log.p = TRUE) + xf)), 1e-15)
  # the density integrates to one and q inverts p
  expect_lt(abs(integrate(d, -Inf, Inf, rel.tol = 1e-12)$value - 1), 1e-10)
  expect_lt(max(abs(q(p(xc)) - xc)), 1e-12)
  expect_lt(max(abs(d(xc, log = TRUE) - log(d(xc)))), 1e-14)
})

test_that("a gamma supplied through os_dist() reproduces its exact rational moments at n = 5 and 10", {
  # The gamma is not a built-in family, so this exercises the user-distribution path end
  # to end against an external exact value. For the gamma of shape 2,
  # F(x) = 1 - (1 + x) exp(-x), so every moment of an order statistic is a finite sum of
  # rationals: expand F(x)^(k-1) and then (1 + x)^m by the binomial theorem and integrate
  # x^(s+1+i) exp(-c x) termwise, which gives (s+1+i)!/c^(s+2+i). The fixture holds those
  # sums evaluated in exact rational arithmetic (fixtures/README.md); the printed table
  # of Gupta (1960, Technometrics 2, Table I) is wrong at n = 10 from rank 4 on, so the
  # exact values and not the table are the reference here.
  ex <- read.csv(test_path("fixtures", "04-gamma2-exact.csv"))
  gamma2 <- os_dist(p = function(q, lower.tail = TRUE, log.p = FALSE)
                      pgamma(q, shape = 2, lower.tail = lower.tail, log.p = log.p),
                    d = function(x, log = FALSE) dgamma(x, shape = 2, log = log),
                    q = function(p) qgamma(p, shape = 2),
                    support = c(0, Inf), name = "gamma(shape = 2)")
  for (n in c(5L, 10L)) {
    m <- os_moments(n, dist = gamma2, cov = FALSE)
    e <- ex[ex$n == n, ]
    expect_identical(m$route, "integrate")
    expect_lt(max(abs(unname(m$alpha) - e$mean)), 1e-10)
    expect_lt(max(abs(unname(m$var) - e$var)), 1e-10)
  }
})

test_that("a function without the required arguments, or an NA, gives an error that says what to fix", {
  expect_error(os_dist(function(q) pnorm(q), dnorm, qnorm),
               "p\\(q, lower.tail = TRUE, log.p = FALSE\\)", class = "error")
  expect_error(os_dist(pnorm, function(x) dnorm(x), qnorm), "d\\(x, log = FALSE\\)")
  expect_error(os_dist(function(q, lower.tail = TRUE, log.p = FALSE) NA_real_, dnorm, qnorm),
               "disagree")
  # a gamma declared on the whole line: its log-tail is -Inf below zero
  g <- os_dist(function(q, lower.tail = TRUE, log.p = FALSE)
                 pgamma(q, 2, lower.tail = lower.tail, log.p = log.p),
               function(x, log = FALSE) dgamma(x, 2, log = log), function(p) qgamma(p, 2))
  expect_error(os_moments(5, dist = g), "support")
})

test_that("a density that disagrees with the distribution function is refused", {
  # 1.001 dnorm passed every check at the median and gave a variance of 1.001
  expect_error(os_dist(pnorm, function(x, log = FALSE) if (log) log(1.001) + dnorm(x, log = TRUE)
                                                        else 1.001 * dnorm(x), qnorm),
               "integrates to")
})
