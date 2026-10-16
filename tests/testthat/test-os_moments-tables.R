# ==============================================================================
# File: test-os_moments-tables.R
# Purpose: The printed tables: Tietjen, Kahaner and Beckman (1976) at n = 13 and
#          on the five-number grid to n = 49 (fixtures/02-tietjen-grid.csv), and
#          Harter (1960) for the mean and variance of the normal range.
#
#          Tietjen orders the sample descending (x_1 > ... > x_N) and the package
#          ascending. The five-number rank set is invariant under r -> N + 1 - r,
#          so the covariance matrix at the ascending ranks equals the matrix at
#          Tietjen's descending indices with Cov(x_i, x_j) = Cov(x_{N+1-j}, x_{N+1-i})
#          (their relation (1)), and the means change sign (their relation (2)).
#          The n = 13 entries below are transcribed from the OSTI scan of report
#          LA-6482-MS (doi:10.2172/7216273), Table I rows "13 i j" and Table II
#          rows "13 i". Tietjen
#          states seven decimals as the worst-case accuracy (p. 3), so 1e-7 is
#          what the table certifies.
#
#          Harter (1960), Ann. Math. Statist. 31, Table 2 (pp. 1130-1131), the
#          columns "Mean, E(W)" and "Variance, sigma_W^2" of the range of n
#          standard normal variables at n = 4m + 1 from 5 to 97, transcribed by
#          eye from the Project Euclid scan (doi:10.1214/aoms/1177705684).
#          Harter states accuracy to a unit in the last place.

test_that("Tietjen's table at n = 13, ranks 1, 4, 7, 10, 13, within its stated accuracy", {
  N <- 13L; idx <- c(1L, 4L, 7L, 10L, 13L)
  tj <- c("1,1" = 0.3152053842, "1,4" = 0.0808649736, "1,7" = 0.0468833088,
          "1,10" = 0.0309322744, "1,13" = 0.0184348220, "4,4" = 0.1330111819,
          "4,7" = 0.0780173339, "4,10" = 0.0517328051, "7,7" = 0.1167989950)
  cov13 <- function(i, j) {
    a <- min(i, j); b <- max(i, j); key <- paste(a, b, sep = ",")
    if (!key %in% names(tj)) key <- paste(N + 1L - b, N + 1L - a, sep = ",")
    unname(tj[key])
  }
  B_tj <- outer(idx, idx, Vectorize(cov13))
  m_tj <- c(1.6679901770, 1.1640771940, 0.8498346324, 0.6028500882,
            0.3883271210, 0.1905236911, 0)
  alpha_desc <- c(m_tj, -rev(m_tj[1:6]))
  alpha_tj <- rev(alpha_desc)[idx]
  m <- os_moments(N, idx, "norm")
  expect_lt(max(abs(m$alpha - alpha_tj)), 1e-7)
  expect_lt(max(abs(m$B - B_tj)), 1e-7)
})

test_that("Tietjen's grid to n = 49: agreement where the printed entry is correct, disagreement where it is not", {
  tj <- read.csv(test_path("fixtures", "02-tietjen-grid.csv"))
  expect_equal(nrow(tj), 240L)
  worst_pass <- 0; n_fail_rows <- 0
  for (N in unique(tj$n)) {
    Q <- (N - 1L) %/% 4L; rk <- c(1L, Q + 1L, 2L * Q + 1L, 3L * Q + 1L, N)
    m <- os_moments(N, rk, "norm")
    t <- tj[tj$n == N, ]
    for (k in seq_len(nrow(t))) {
      if (t$doubtful[k]) next                                   # the illegible cell at n = 49
      got <- if (t$table[k] == "II") m$alpha[[match(t$rank_i[k], rk)]]
             else m$B[match(t$rank_i[k], rk), match(t$rank_j[k], rk)]
      d <- abs(got - t$tietjen[k])
      if (t$pass[k]) worst_pass <- max(worst_pass, d)
      else { n_fail_rows <- n_fail_rows + 1L; expect_gt(d, 1e-7) }  # the table's error, not ours
    }
  }
  expect_lt(worst_pass, 1e-7)
  expect_equal(n_fail_rows, sum(!tj$pass & !tj$doubtful))
})

test_that("Harter's mean and variance of the normal range at n = 5 to 97", {
  harter <- read.table(header = TRUE, text = "
n   mean_W        var_W
5   2.3259289473  0.7466376009
9   2.9700263244  0.6525962151
13  3.3359803541  0.5935411244
17  3.5878839618  0.5536130572
21  3.7783358298  0.5244270274
25  3.9306292195  0.5018883188
29  4.0570442921  0.4837838184
33  4.1648166719  0.4688083838
37  4.2585540507  0.4561382080
41  4.3413643695  0.4452248562
45  4.4154391280  0.4356873457
49  4.4823791942  0.4272516817
53  4.5433881668  0.4197150671
57  4.5993928961  0.4129237256
61  4.6511217036  0.4067586207
65  4.6991569040  0.4011259526
69  4.7439711082  0.3959506527
73  4.7859529522  0.3911718179
77  4.8254257340  0.3867394380
81  4.8626611740  0.3826120024
85  4.8978897515  0.3787547214
89  4.9313085861  0.3751381825
93  4.9630875311  0.3717373222
97  4.9933739426  0.3685306286
")
  eng <- t(vapply(harter$n, function(n) {
    m <- os_moments(n, c(1L, n), "norm")
    c(unname(m$alpha[2] - m$alpha[1]), unname(m$B[1, 1] + m$B[2, 2] - 2 * m$B[1, 2]))
  }, numeric(2)))
  expect_lt(max(abs(eng[, 1] - harter$mean_W)), 2e-8)
  expect_lt(max(abs(eng[, 2] - harter$var_W)), 2e-8)
})
