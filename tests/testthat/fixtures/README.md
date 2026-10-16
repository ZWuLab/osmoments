
# tests/testthat/fixtures/

Reference values the tests compare against. Every file was written by a script and none is
edited by hand. The scripts are kept in the development tree under `dev/` (see
`dev/README.md` there), which is not part of the distributed package, so the Content column
states where each file's values come from.

| File | Written by | Content |
| --- | --- | --- |
| `01-frozen-grid.csv` | `dev/01-frozen-grid-fixture.R` | The means (`what = "alpha"`) and the upper triangle of the covariance matrix (`what = "B"`) of the standard normal order statistics at the five-number ranks, for 39 sample sizes from 5 to 1001, computed by the reference engine, an earlier implementation of the same method for the normal, written separately before the package. Columns: `n`, `what`, `i`, `j` (positions in the rank set), `rank_i`, `rank_j`, `value`, `quad_err`, `seconds` |
| `02-tietjen-grid.csv` | `dev/02-tietjen-fixture.R` | The printed variances, covariances (`table = "I"`) and means (`table = "II"`) of Tietjen, Kahaner and Beckman (1976) at the five-number ranks for $n = 5$ to $49$, in the ascending convention, with `pass` (the printed entry is correct to $10^{-7}$) and `doubtful` (the one illegible cell) |
| `03-frozen-blue.csv` | `dev/03-frozen-blue-fixture.R` | Lloyd's BLUE weights and the variance factors of the location and scale estimators from the reference engine, for the scenarios S1 (minimum, median, maximum), S2 (the three quartiles) and S3 (the five-number summary) at $n = 25$ and $41$. Columns: `n`, `scenario`, `rank`, `w_mu`, `w_sig`, `var_mean`, `var_sd`, `cov` |
| `04-gamma2-exact.csv` | `dev/04-gamma2-exact-fixture.R`, which copies the output of the exact evaluation it names | The means and variances of the order statistics of the standardized gamma with shape 2 ($f(x) = x e^{-x}$) at every rank for $n = 5$ and $10$, computed as finite sums of rationals in exact arithmetic (the binomial expansion of $F(x)^{k-1}$ and $(1+x)^m$, with $F(x) = 1 - (1+x)e^{-x}$, integrated termwise) and printed to 17 significant digits. Columns: `n`, `k`, `mean`, `var` |
