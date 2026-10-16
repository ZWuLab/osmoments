# ==============================================================================
# File: test-os_checks.R
# Purpose: The moment condition and the run-time checks: the Cauchy median
#          against the quantile-function route, the refusals, the permutation
#          identities on complete rank sets, the unit row sums for the normal,
#          and the reflection symmetry computed without halving.

cauchy <- os_dist(pcauchy, dcauchy, qcauchy, symmetric = TRUE, tail = c(1, 1), name = "cauchy")

test_that("the Cauchy median of five has the variance of the quantile-function route", {
  # E[X_(3:5)^2] = int_0^1 tan(pi (u - 1/2))^2 30 u^2 (1 - u)^2 du, mean 0 by symmetry
  ref <- integrate(function(u) tan(pi * (u - 0.5))^2 * 30 * u^2 * (1 - u)^2, 0, 1,
                   rel.tol = 1e-12)$value
  g <- os_moments(5, 3, cauchy)
  expect_lt(abs(g$var[[1]] - ref), 1e-9)
  expect_lt(abs(g$alpha[[1]]), 1e-9)
})

test_that("the moment condition refuses ranks whose moments do not exist, naming the rank", {
  expect_error(os_moments(5, 1:5, cauchy), "X_\\(1\\)")
  expect_error(os_moments(5, c(3, 5), cauchy), "X_\\(5\\)")
  t15 <- os_dist(function(q, lower.tail = TRUE, log.p = FALSE) pt(q, 1.5, lower.tail = lower.tail, log.p = log.p),
                 function(x, log = FALSE) dt(x, 1.5, log = log), function(p) qt(p, 1.5),
                 symmetric = TRUE, tail = c(1.5, 1.5), name = "t1.5")
  expect_error(os_moments(5, 1:5, t15), "X_\\(1\\)")
  expect_silent(m <- os_moments(5, 2:4, t15, check = FALSE))     # 1.5 * 2 = 3 > 2
  t3 <- os_dist(function(q, lower.tail = TRUE, log.p = FALSE) pt(q, 3, lower.tail = lower.tail, log.p = log.p),
                function(x, log = FALSE) dt(x, 3, log = log), function(p) qt(p, 3),
                symmetric = TRUE, tail = c(3, 3), name = "t3")
  m <- os_moments(5, 1:5, t3)                                        # 3 * 1 = 3 > 2
  expect_lt(max(m$checks$identities), 1e-8)
})

test_that("the permutation identities hold on complete rank sets for every family, by quadrature", {
  for (fam in c("norm", "exp", "logis", "gumbel", "unif")) for (n in c(5L, 9L, 13L)) {
    m <- os_moments(n, 1:n, fam, method = "integrate")
    expect_lt(max(m$checks$identities), 1e-9, label = paste(fam, n))
    # the same, explicitly: sum alpha = n E[X], sum E[X^2] = n E[X^2], sum B = n Var(X)
    p1 <- os_moments(1, 1, fam, method = "integrate")
    expect_lt(abs(sum(m$alpha) - n * p1$alpha[[1]]), 1e-9)
    expect_lt(abs(sum(m$var + m$alpha^2) - n * (p1$var[[1]] + p1$alpha[[1]]^2)), 1e-9)
    expect_lt(abs(sum(m$B) - n * p1$var[[1]]), 1e-9)
  }
  # the unit row sums are an identity for the normal and not for the exponential
  m <- os_moments(9, 1:9, "norm")
  expect_lt(max(abs(rowSums(m$B) - 1)), 1e-9)
  expect_true("row_sums" %in% names(m$checks$identities))
  e <- os_moments(2, 1:2, "exp")
  expect_equal(unname(rowSums(e$B)), c(1 / 2, 3 / 2))
  expect_false("row_sums" %in% names(e$checks$identities))
})

test_that("reflection symmetry holds without being imposed, for the normal and the logistic at n = 25", {
  rk <- c(1, 7, 13, 19, 25)
  for (fam in c("norm", "logis")) {
    f <- osmoments:::.os_family(fam)
    full <- os_dist(f$p, f$d, f$q, symmetric = FALSE, name = fam)    # no halving, no closed forms
    m <- os_moments(25, rk, full, method = "integrate")
    expect_true(is.na(m$checks$reflection))
    M <- m$B + tcrossprod(m$alpha)                                   # E[X_r X_s]
    expect_lt(max(abs(m$alpha + rev(m$alpha))), 1e-8)
    expect_lt(max(abs(M - M[5:1, 5:1])), 1e-8)
    h <- os_moments(25, rk, fam, method = "integrate")               # halved by symmetry
    expect_lt(max(abs(h$B - m$B)), 1e-10)
    expect_lt(h$checks$reflection, 1e-8)
  }
})
