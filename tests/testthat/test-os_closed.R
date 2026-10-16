# ==============================================================================
# File: test-os_closed.R
# Purpose: The quadrature against the closed forms: the uniform and the
#          exponential at every pair, the logistic single moments, the Gumbel
#          maximum, and the n = 2 anchors.

five_ranks <- function(n) unique(pmin(n, c(1, round(n / 4), round(n / 2), round(3 * n / 4), n)))

test_that("uniform: quadrature equals the beta and Dirichlet closed forms, every pair, to n = 1001", {
  for (n in c(20, 100, 1001)) {
    rk <- five_ranks(n)
    q <- os_moments(n, rk, "unif", method = "integrate")
    cl <- os_moments(n, rk, "unif", method = "closed")
    expect_identical(cl$route, "closed")
    expect_lt(max(abs(q$alpha - cl$alpha)), 1e-10)
    expect_lt(max(abs(q$B - cl$B)), 1e-10)
    expect_equal(unname(cl$alpha), rk / (n + 1))
    expect_equal(unname(cl$B), outer(rk, rk, function(r, s) pmin(r, s) * (n - pmax(r, s) + 1)) /
                   ((n + 1)^2 * (n + 2)))
  }
})

test_that("exponential: quadrature equals Renyi's closed forms, every pair, to n = 1001", {
  for (n in c(20, 100, 1001)) {
    rk <- five_ranks(n)
    q <- os_moments(n, rk, "exp", method = "integrate")
    cl <- os_moments(n, rk, "exp", method = "closed")
    expect_lt(max(abs(q$alpha - cl$alpha)), 1e-10)
    expect_lt(max(abs(q$B - cl$B)), 1e-10)
    a <- sapply(rk, function(r) sum(1 / (n - seq_len(r) + 1)))
    v <- sapply(rk, function(r) sum(1 / (n - seq_len(r) + 1)^2))
    expect_equal(unname(cl$alpha), a)
    expect_equal(unname(cl$B), outer(seq_along(rk), seq_along(rk), function(i, j) v[pmin(i, j)]))
  }
})

test_that("logistic: means and variances equal the digamma and trigamma forms", {
  for (n in c(20, 100)) {
    rk <- c(1, round(n / 4), round(n / 2), n)
    q <- os_moments(n, rk, "logis", method = "integrate")
    expect_lt(max(abs(q$alpha - (digamma(rk) - digamma(n - rk + 1)))), 1e-10)
    expect_lt(max(abs(q$var - (trigamma(rk) + trigamma(n - rk + 1)))), 1e-10)
    cl <- os_moments(n, rk, "logis")
    expect_identical(cl$route, "mixed")
    expect_equal(unname(cl$alpha), digamma(rk) - digamma(n - rk + 1))
    expect_lt(max(abs(cl$B - q$B)), 1e-10)
  }
})

test_that("Gumbel: the maximum has mean gamma + log n and variance pi^2/6", {
  for (n in c(20, 100)) {
    q <- os_moments(n, c(1, n), "gumbel", method = "integrate")
    expect_lt(abs(q$alpha[[2]] - (-digamma(1) + log(n))), 1e-10)
    expect_lt(abs(q$var[[2]] - pi^2 / 6), 1e-10)
    cl <- os_moments(n, c(1, n), "gumbel")
    expect_identical(cl$route, "mixed")
    expect_equal(cl$alpha[[2]], -digamma(1) + log(n))
    expect_lt(max(abs(cl$B - q$B)), 1e-10)
  }
})

test_that("n = 2 anchors: X_(1) X_(2) = X_1 X_2 gives Cov = (E X)^2 - mu_{1:2} mu_{2:2}", {
  g <- os_moments(2, 1:2, "logis")
  expect_lt(abs(g$B[1, 2] - 1), 1e-10)             # (E X)^2 = 0, mu_{1:2} = -1, mu_{2:2} = 1
  g <- os_moments(2, 1:2, "gumbel")
  expect_lt(abs(g$B[1, 2] - log(2)^2), 1e-10)      # gamma^2 - (gamma - log 2)(gamma + log 2)
  g <- os_moments(2, 1:2, "norm")
  expect_lt(abs(g$B[1, 2] - 1 / pi), 1e-9)         # mu_{2:2} = 1/sqrt(pi)
})
