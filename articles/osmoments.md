# Exact moments of order statistics with osmoments

## What this package computes

Order statistics are the sorted values of a sample,
$`X_{(1)} \le \dots \le X_{(n)}`$. Many statistics are built from them,
for example, a sample median or quantile, the range, a trimmed mean, the
best linear unbiased estimators of location and scale, or a mean and a
standard deviation imputed from a reported five-number summary. The mean
and the standard error of such a statistic require the means, variances
and covariances of the order statistics it uses. Without the exact
values, one falls back on simulation or on large-sample approximations,
whose error is largest at the extreme ranks and in small samples. The
tables of the covariances of normal order statistics printed in the
literature stop at $`n = 50`$. The software we found gives the means, or
the means and variances, of single order statistics, or approximates the
covariances.

The `osmoments` package computes these moments exactly, for any ranks
and any sample size, by numerical integration of their defining
integrals, or by closed-form expressions where the distribution has
them. It covers the normal, uniform, exponential, logistic and Gumbel
distributions, and any continuous distribution the user supplies through
its distribution, density and quantile functions. Building on the
moments, three more functions are provided.
[`os_quantile_weights()`](https://zwulab.github.io/osmoments/reference/os_quantile_weights.md)
gives the weights of a sample quantile under any of the nine definitions
of R’s [`quantile()`](https://rdrr.io/r/stats/quantile.html).
[`os_lin_moments()`](https://zwulab.github.io/osmoments/reference/os_lin_moments.md)
gives the means and covariances of any linear combinations of order
statistics.
[`os_blue()`](https://zwulab.github.io/osmoments/reference/os_blue.md)
gives the best linear unbiased estimators (BLUE) of location and scale
of Lloyd (1952).

The package has two requirements. The distribution must be continuous,
and a moment is computed only if it exists, which the package checks
rank by rank from the tail behavior of the distribution (the `tail`
argument of
[`os_dist()`](https://zwulab.github.io/osmoments/reference/os_dist.md)).
Here “exact” means the finite-sample integral evaluated numerically, as
opposed to a large-sample or series approximation. The integrator is
asked for a relative tolerance of $`10^{-11}`$ on a single integral and
$`10^{-9}`$ on the outer integral of a covariance. A tolerance is the
accuracy asked of the integrator. It is not a bound on the error of the
result (see “What the result holds”).

The package’s tests, which `R CMD check` runs, compare the quadrature
with values known from other sources: closed forms, printed tables and
exact identities. All the tolerances below are absolute differences,
unless they say otherwise. The tests check the quadrature against the
closed forms of the uniform and the exponential at the five ranks of a
five-number summary and all their pairs, for $`n = 20`$, 100 and 1001,
within $`10^{-10}`$. For the normal, they check the ten-decimal tables
of Tietjen, Kahaner and Beckman (1976) within $`10^{-7}`$, the accuracy
those tables certify, on the entries printed correctly (Appendix B says
which are not). They also check the mean and variance of the range
tabulated by Harter (1960), within $`2 \times 10^{-8}`$ up to
$`n = 97`$. On a complete set of ranks, the order statistics are a
permutation of the sample, so their moments satisfy identities such as
“the means sum to $`n`$ times the parent mean”. The tests assert these
identities within $`10^{-9}`$. The largest residual, the difference
between the two sides of an identity, that we measured on the full
matrices of the five families up to $`n = 100`$ was
$`1.1 \times 10^{-12}`$. Wu and Tong (2026) describe the computation and
these checks in full, and compare the package with other software.

The call is one function:

``` r

m <- os_moments(
  n = 20,                   # the sample size
  ranks = c(1, 5, 10, 20),  # the order statistics to compute (default: all)
  dist = "norm"             # "norm", "unif", "exp", "logis" or "gumbel",
                            # or a distribution built with os_dist()
)
m
#> Exact moments of order statistics: norm (location 0, scale 1), n = 20, ranks 1 5 10 20
#> alpha (means):
#>          r1          r5         r10         r20 
#> -1.86747500 -0.74538300 -0.06199629  1.86747500 
#> B (covariance matrix):
#>             r1         r5        r10        r20
#> r1  0.27569660 0.05705664 0.02965625 0.01020472
#> r5  0.05705664 0.09399600 0.04933743 0.01707114
#> r10 0.02965625 0.04933743 0.07694744 0.02688388
#> r20 0.01020472 0.01707114 0.02688388 0.27569660
#> route: integrate; quadrature error estimate: 1.5e-11
#> checks: reflection n/a; identities n/a
```

The print shows the means of the four order statistics of a standard
normal sample of 20 (`alpha`) and their covariance matrix (`B`), with
the variances on its diagonal. “What the result holds” describes the
other fields. The last line, `checks`, shows `n/a` here. The two
run-time checks of Appendix B need all the ranks 1 to $`n`$, or a rank
set that contains $`n + 1 - r`$ whenever it contains $`r`$, such as 1,
5, 16 and 20 at $`n = 20`$.

Use the family name when the distribution is one of the five. For
another distribution, build it once with
[`os_dist()`](https://zwulab.github.io/osmoments/reference/os_dist.md)
(see “A user-supplied distribution”). For the means and the variances at
many ranks without the covariances, pass `cov = FALSE`, which skips the
double integrals. The help page
[`?os_moments`](https://zwulab.github.io/osmoments/reference/os_moments.md)
lists every argument.

## Computing cost

This section helps the user size a computation before asking for a full
matrix. Each mean and each variance is one integral over the real line.
Each covariance is a double integral, evaluated as nested quadrature,
and costs far more. The chunk below prints the times, in milliseconds,
on the machine that built this vignette.

``` r

t_single <- system.time(os_moments(100, ranks = c(1, 50, 100),
                                   cov = FALSE))[["elapsed"]]
t_pairs  <- system.time(os_moments(100, ranks = c(1, 50, 100)))[["elapsed"]]
t_full10 <- system.time(os_moments(10))[["elapsed"]]
t_full20 <- system.time(os_moments(20))[["elapsed"]]
round(1000 * c(three_means_and_variances = t_single,
               with_three_covariances = t_pairs,
               full_matrix_n10 = t_full10, full_matrix_n20 = t_full20))
#> three_means_and_variances    with_three_covariances           full_matrix_n10 
#>                         2                       186                       957 
#>           full_matrix_n20 
#>                      4100
```

The computing cost of one integral hardly depends on $`n`$. On an Apple
M3 Pro, for the standard normal and $`n`$ from 5 to 1001, a single
integral took 0.19 to 0.25 ms and a double integral 36 to 76 ms. A call
at the five ranks of a five-number summary, with their ten covariances,
took 0.6 to 0.8 seconds for the uniform and the exponential at
$`n = 10^4`$, $`10^5`$ and $`10^6`$. So a few ranks cost well under a
second, whatever the sample size.

A full matrix is different. At sample size $`n`$ it needs $`n(n-1)/2`$
covariances, or about half as many for a symmetric family, so its time
grows with the square of $`n`$. By that count, the full matrix at
$`n = 100`$ costs about 26 times the one at $`n = 20`$ above. On the
machine that built this vignette, that predicts about 107 seconds. On
the M3 Pro, the full normal matrix took 2.4 minutes at $`n = 100`$ on
one core, and 70 minutes at $`n = 1001`$ on 12 cores. For $`n`$ in the
hundreds, run one call of
[`os_moments()`](https://zwulab.github.io/osmoments/reference/os_moments.md)
per pair of ranks in parallel, rather than one call for the whole
matrix. The package has no helper for this. Each call returns the 2 by 2
covariance matrix of its pair.

## Installation and the families

This section says how to install the package and what the five built-in
families are. Install it from GitHub, with the vignette:

``` r

# install.packages("remotes")
remotes::install_github("ZWuLab/osmoments", build_vignettes = TRUE)
```

To cite the package in a publication, use `citation("osmoments")`, which
gives the paper of Wu and Tong (2026).

The package imports nothing beyond base R and `stats`. It computes the
moments of each family’s standard member $`Z`$, the member with location
0 and scale 1, then shifts and scales them by the `location` and `scale`
arguments of
[`os_moments()`](https://zwulab.github.io/osmoments/reference/os_moments.md),
so $`X = \mathrm{location} + \mathrm{scale} \cdot Z`$. The means
(`alpha`) become
$`\mathrm{location} + \mathrm{scale} \cdot E[Z_{(r)}]`$, and the
covariance matrix (`B`) is multiplied by $`\mathrm{scale}^2`$.

| `dist` | Standard member $`Z`$ | What `location` and `scale` mean | Closed forms returned by default |
|----|----|----|----|
| `"norm"` | $`N(0, 1)`$ | the mean and the standard deviation | none; all moments by quadrature |
| `"unif"` | uniform on $`[0, 1]`$ | the lower end and the length of the support | every mean, variance and covariance |
| `"exp"` | exponential with rate 1 | the lower end of the support, and the mean $`1/\mathrm{rate}`$ | every mean, variance and covariance |
| `"logis"` | standard logistic | the location and scale parameters | the means and the variances |
| `"gumbel"` | Gumbel, maximum form, $`F(x) = \exp(-e^{-x})`$ | the location and scale parameters | the mean and variance of the maximum |

By default the package uses a closed form wherever the family has one.
For 20 exponential waiting times with mean 2, the means of the minimum,
the tenth order statistic and the maximum come from closed forms:

``` r

# the order statistics of 20 exponential waiting times with mean 2 (rate 1/2)
e <- os_moments(20, ranks = c(1, 10, 20), dist = "exp", scale = 2)
e$alpha
#>       r1      r10      r20 
#> 0.100000 1.337543 7.195479
e$route
#> [1] "closed"
```

`route` says where the numbers came from. `"closed"` means the closed
forms, and `"integrate"` the quadrature. `"mixed"` means that the family
has closed forms for some of the requested quantities and not for
others. `method = "integrate"` forces the quadrature everywhere, which
is how the package’s tests check the quadrature against the closed
forms.

## What the result holds

This section names the fields of the result that a user reads first, and
says what the error estimate is and is not.
[`os_moments()`](https://zwulab.github.io/osmoments/reference/os_moments.md)
returns an object of class `"osmoments"`.

| Field | Content |
|----|----|
| `alpha` | the means $`E[X_{(r)}]`$, named by rank |
| `var` | the variances |
| `B` | the covariance matrix, with the variances on its diagonal (`NULL` when `cov = FALSE`) |
| `route` | `"closed"`, `"integrate"` or `"mixed"` |
| `quad_err` | the sum of the error estimates that [`integrate()`](https://rdrr.io/r/stats/integrate.html) reported for the single integrals and for the outer integrals of the covariances |
| `standard` | a list with the `alpha`, `var` and `B` of the standard member, before `location` and `scale` are applied, which [`os_blue()`](https://zwulab.github.io/osmoments/reference/os_blue.md) reads |
| `checks` | the residuals of the run-time checks (Appendix B), `NA` where a check did not apply |

`quad_err` is an estimate, not a bound, for two reasons. The inner
integral of a covariance does not pass its error estimate to the outer
one. And the error of a mean is not propagated into the variance that
subtracts its square, or into the covariances that subtract its product
with another mean. Treat `quad_err` as a diagnostic. No certified error
bound is stored. The evidence of accuracy is the agreement with closed
forms, printed tables and identities, summarized above and in Appendix
B.

That evidence has a limit worth knowing. The covariances with closed
forms belong to the uniform and the exponential, which are both bounded
below, so they do not test a double integral over an unbounded lower
tail. For the normal beyond $`n = 50`$, where no printed table exists,
the checks are the identities on complete rank sets, which test sums of
moments rather than each moment.

## A user-supplied distribution

This section shows how to supply a distribution that is not built in,
and the four points that matter for such a distribution.
[`os_dist()`](https://zwulab.github.io/osmoments/reference/os_dist.md)
builds the object that `dist` accepts, from the three functions of a
distribution in R’s convention. Here is Student’s $`t`$ with 5 degrees
of freedom:

``` r

t5 <- os_dist(
  p = function(q, lower.tail = TRUE, log.p = FALSE)
    pt(q, df = 5, lower.tail = lower.tail, log.p = log.p),
  d = function(x, log = FALSE) dt(x, df = 5, log = log),
  q = function(p) qt(p, df = 5),
  symmetric = TRUE,     # the density is symmetric about its median
  tail = c(5, 5),       # the density decays like |x|^(-5-1) in both tails
  name = "t5"
)
os_moments(5, dist = t5)
#> Exact moments of order statistics: t5 (location 0, scale 1), n = 5, ranks 1 2 3 4 5
#> alpha (means):
#>            r1            r2            r3            r4            r5 
#> -1.450975e+00 -5.580065e-01 -8.649563e-18  5.580065e-01  1.450975e+00 
#> B (covariance matrix):
#>           r1        r2        r3        r4        r5
#> r1 1.1478770 0.3802510 0.2288315 0.1665120 0.1398297
#> r2 0.3802510 0.4270062 0.2648562 0.1959689 0.1665120
#> r3 0.2288315 0.2648562 0.3501660 0.2648562 0.2288315
#> r4 0.1665120 0.1959689 0.2648562 0.4270062 0.3802510
#> r5 0.1398297 0.1665120 0.2288315 0.3802510 1.1478770
#> route: integrate; quadrature error estimate: 2.17e-10
#> checks: reflection 1.61e-17; identities sum_alpha 2.94e-17, sum_sq 8.88e-16, sum_B 3.29e-14
```

The mean of the middle order statistic prints as a tiny number rather
than 0. It is zero up to rounding, as it must be for a symmetric parent.
Every variance exists here, because $`\nu r > 2`$ and
$`\nu (n - r + 1) > 2`$ at every rank when $`\nu = 5`$ and $`n = 5`$.
The `checks` line gives the residuals of the reflection check and of the
three permutation identities of Appendix B.

1.  **`p` must honor `log.p` and `lower.tail`, and `d` must honor
    `log`.** The package evaluates the integrands on the log scale in
    both tails, through `p(x, log.p = TRUE)` and
    `p(x, lower.tail = FALSE, log.p = TRUE)`. The log scale keeps powers
    like $`F(x)^{r-1}`$ from underflowing to zero at large $`n`$. A `p`
    that ignores these arguments, or that computes `log(p(x))` from a
    probability that has already underflowed, loses the tails.
    [`os_dist()`](https://zwulab.github.io/osmoments/reference/os_dist.md)
    checks that the arguments are honored at the median and the lower
    quartile, and that `d` integrates to one half between the quartiles.
    [`os_moments()`](https://zwulab.github.io/osmoments/reference/os_moments.md)
    checks that the log-tails and the log-density are finite at a few
    points around the requested ranks. Neither check certifies accuracy,
    so prefer the `p*`, `d*` and `q*` functions of `stats`, or of a
    package that implements them with care. A distribution whose support
    has a finite end gives its endpoints in `support`, for example
    `support = c(0, Inf)` for a gamma, a Weibull or a lognormal.
2.  **The tail indices decide which moments exist.** Suppose the density
    decays like $`|x|^{-\nu_L - 1}`$ in the lower tail and like
    $`x^{-\nu_U - 1}`$ in the upper tail. Then the order statistic
    $`X_{(r)}`$ has a finite second moment if and only if
    $`\nu_L r > 2`$ and $`\nu_U (n - r + 1) > 2`$. Declare
    $`(\nu_L, \nu_U)`$ in `tail`, with `Inf` for a tail that is bounded
    or decays faster than any power. A rank whose moment does not exist
    stops the computation, with an error that names the rank. The
    default `tail = c(Inf, Inf)` declares that every moment exists, so
    for a heavy-tailed distribution, omitting `tail` removes the check.
    A moment that does not exist then usually stops the computation with
    an integration error, but nothing guarantees it. So declare the tail
    indices of a heavy-tailed distribution.
3.  **Supply the standard member.** Pass the location and the scale to
    [`os_moments()`](https://zwulab.github.io/osmoments/reference/os_moments.md)
    rather than building them into the functions. The integrals are
    taken about the median of the standard member. A standard member
    with a large location of its own, say $`10^8`$ with unit scale,
    makes $`x - Q(1/2)`$ a difference of two large numbers, which loses
    precision, and the computation stops with an error. The `location`
    argument itself may be as large as needed: the tests compute the
    range of two uniform observations at `location = 1e16` exactly. The
    scale of the standard member does not matter, because the
    integration works in units of its spread (Appendix A). The tests
    check that a normal supplied with a standard deviation from 0.001 to
    $`10^6`$ agrees with the scaled standard normal within $`10^{-11}`$
    in relative terms.
4.  **Some parents stop with an error rather than return a number.** The
    integration range is split around the peak of each order statistic
    and follows the tails of the parent. The split does not follow a
    second mode. For a mixture of two narrow normals, the mode that
    holds a tenth of the mass fell between the quadrature points. The
    package therefore integrates the density of every order statistic it
    computes by quadrature, and stops, naming the rank, when that mass
    is not 1. It also stops when a piece of the range cannot be
    resolved. The package stops in this way for a tail as heavy as a
    lognormal with `sdlog = 3` or a Weibull with shape 0.3, whose
    moments exist but whose integrands put their mass far out. It also
    happens for the extremes of a parent whose density is infinite at an
    endpoint, once they concentrate against it: the minimum of a
    beta(1/2, 1/2) sample is computed up to $`n = 100`$ and refused from
    $`n = 1000`$ on. These checks catch some failures and certify none.

The Cauchy distribution shows the second point. The parent has no mean,
but the median of five Cauchy observations has a variance. The ranks 1,
2, 4 and 5 have no second moment. At rank 2, for example,
$`\nu_L r = 2`$, which is not above 2.
[`os_moments()`](https://zwulab.github.io/osmoments/reference/os_moments.md)
stops at the first such rank:

``` r

cauchy <- os_dist(pcauchy, dcauchy, qcauchy, symmetric = TRUE,
                  tail = c(1, 1), name = "cauchy")
os_moments(5, ranks = 3, dist = cauchy)$var
#>       r3 
#> 1.221253
os_moments(5, ranks = 1:5, dist = cauchy)
#> Error:
#> ! os_moments(): the moment of order 2 of X_(1) does not exist: the lower tail index 1 times the rank 1 is not above 2
```

## Sample quantiles, linear combinations and the BLUE

This section shows the three functions that turn the moments of order
statistics into the moments of the statistics computed from them. A
sample quantile under any of the nine definitions of
[`quantile()`](https://rdrr.io/r/stats/quantile.html) (Hyndman and Fan,
1996) is a fixed combination of two adjacent order statistics, so its
moments follow from theirs. `os_quantile_weights(p, n, type)` gives the
ranks and the weights, with R’s default `type = 7` unless another of the
nine definitions is named.
[`os_lin_moments()`](https://zwulab.github.io/osmoments/reference/os_lin_moments.md)
turns them into moments. It returns the means in `mean`, the variances
in `var`, and the covariance matrix of the statistics in `cov`. As a
first case, take the median of an even-sized sample, which averages two
order statistics. How does its variance compare with that of the median
of the next odd size?

``` r

qw <- os_quantile_weights(p = 0.5, n = 10)   # type 7, R's default
qw$W   # one column per probability, rows named by rank
#>    50%
#> r5 0.5
#> r6 0.5
m10 <- os_moments(10, ranks = qw$ranks)
med10 <- os_lin_moments(qw$W, m10)
m11 <- os_moments(11, ranks = 6)
c(sd_median_n10 = sqrt(med10$var[[1]]), sd_median_n11 = sqrt(m11$var[[1]]))
#> sd_median_n10 sd_median_n11 
#>     0.3719226     0.3703545
```

For a standard normal sample of size 10, the median is the average of
the fifth and sixth order statistics. The chunk prints its standard
deviation next to that of the median of 11 observations, which is a
single order statistic. The standard deviation of the median of ten is
1.0042 times that of the median of eleven.

[`os_lin_moments()`](https://zwulab.github.io/osmoments/reference/os_lin_moments.md)
takes any weights. With one column per statistic, it returns their means
and their covariance matrix, so the range and the midrange of the same
sample come together:

``` r

m <- os_moments(20, ranks = c(1, 20), dist = "logis")
W <- cbind(range = c(-1, 1), midrange = c(0.5, 0.5))
os_lin_moments(W, m)[c("mean", "cov")]
#> $mean
#>    range midrange 
#> 7.095479 0.000000 
#> 
#> $cov
#>             range  midrange
#> range    3.284619 0.0000000
#> midrange 0.000000 0.8750502
```

`cov` is the 2 by 2 covariance matrix of the range and the midrange.
Their covariance is half the difference of the variances of the maximum
and the minimum, which is zero for a symmetric parent such as the
logistic.

The BLUE of location and scale of Lloyd (1952) use the moments of the
available order statistics. Suppose a study reports only its five-number
summary: the minimum, the three quartiles and the maximum of $`n = 25`$
observations. If the data are normal,
[`os_blue()`](https://zwulab.github.io/osmoments/reference/os_blue.md)
gives the weights of the linear unbiased estimators of the mean and the
standard deviation with the smallest variance, in `w_location` and
`w_scale`. It also gives the variances of those estimators per unit of
$`\sigma^2`$, in `var_location` and `var_scale`. The chunk prints the
weights and the standard errors in units of $`\sigma`$:

``` r

five <- c(1, 7, 13, 19, 25)   # the five-number ranks when n = 25 = 4 * 6 + 1
b <- os_blue(os_moments(25, ranks = five))
round(b$w_location, 4)
#>     r1     r7    r13    r19    r25 
#> 0.0821 0.2857 0.2644 0.2857 0.0821
round(b$w_scale, 4)
#>      r1      r7     r13     r19     r25 
#> -0.1713 -0.2566  0.0000  0.2566  0.1713
c(se_of_mean_estimate = sqrt(b$var_location),
  se_of_sd_estimate = sqrt(b$var_scale))
#> se_of_mean_estimate   se_of_sd_estimate 
#>           0.2068675           0.1585884
```

The standard error of the estimated mean is 0.2069 $`\sigma`$, against
$`\sigma/\sqrt{25} = 0.2\,\sigma`$ for the sample mean of all 25
observations. So the standard error from the five-number summary is 3.4%
larger here. The scale that
[`os_blue()`](https://zwulab.github.io/osmoments/reference/os_blue.md)
estimates is the family’s scale parameter: the standard deviation for
the normal, the scale parameter for the logistic or the Gumbel, and the
length of the support for the uniform. The same call works for any
family and any subset of two or more ranks whose covariance matrix is
invertible. For example, here is the BLUE of a logistic location and
scale from a complete sample of nine:

``` r

bl <- os_blue(os_moments(9, dist = "logis"))
round(bl$w_location, 4)
#>     r1     r2     r3     r4     r5     r6     r7     r8     r9 
#> 0.0320 0.0901 0.1335 0.1600 0.1689 0.1600 0.1335 0.0901 0.0320
round(c(var_location = bl$var_location, var_scale = bl$var_scale,
        cov = bl$cov), 5)
#> var_location    var_scale          cov 
#>      0.34247      0.08627      0.00000
```

The covariance of the two estimators prints as zero, because the
standard logistic is symmetric about zero and the rank set is complete.
So the location and scale estimators are uncorrelated.

## Appendix A: how the integrals are evaluated

This appendix states the integrals and the four numerical devices that
evaluate them, for the reader who wants to know what “exact” rests on.
The density of $`X_{(r)}`$ and the joint density of $`X_{(r)}`$ and
$`X_{(s)}`$, $`r < s`$, are
``` math
\begin{aligned}
f_{r:n}(x) &= \frac{n!}{(r-1)!\,(n-r)!} F(x)^{r-1} \{1 - F(x)\}^{n-r} f(x), \\
f_{r,s:n}(y, z) &= \frac{n!}{(r-1)!\,(s-r-1)!\,(n-s)!} F(y)^{r-1} \{F(z) - F(y)\}^{s-r-1}
\{1 - F(z)\}^{n-s} f(y) f(z)
\end{aligned}
```
for $`y < z`$, and the moments are their integrals against $`x`$,
$`x^2`$ and $`yz`$ (David and Nagaraja, 2003). The package evaluates
them with [`integrate()`](https://rdrr.io/r/stats/integrate.html) and
four devices.

1.  **Everything on the log scale.** The factorials overflow double
    precision ($`n!`$ alone does from $`n = 171`$ on), and the powers of
    the tail probabilities underflow to zero at large $`n`$. So the
    package works with logarithms. The log of the coefficient comes from
    [`lbeta()`](https://rdrr.io/r/base/Special.html), the log powers of
    the tails from `p(x, log.p = TRUE)` and
    `p(x, lower.tail = FALSE, log.p = TRUE)`, and the log density from
    `d(x, log = TRUE)`. The sum is exponentiated once. Working on the
    log scale is the device of Royston (1982) for the means of normal
    order statistics, carried to the product moments.
2.  **No cancellation in $`F(z) - F(y)`$.** When both points lie above
    the median, $`F(y)`$ and $`F(z)`$ are close to one, and their
    difference loses its relative precision. The package therefore forms
    the difference on the side where the probabilities are small. It
    computes $`F(z)\{1 - F(y)/F(z)\}`$ when both points lie below the
    median, and $`S(y)\{1 - S(z)/S(y)\}`$ with $`S = 1 - F`$ when both
    lie above it. The logarithm of the difference then keeps its
    relative precision.
3.  **A split around the peak.** The density of $`X_{(r)}`$ concentrates
    near $`Q(r/(n+1))`$, in a region whose width shrinks like
    $`n^{-1/2}`$. An adaptive rule started on the whole line can miss
    that region, and return a value near zero with a small claimed
    error. The package therefore splits the range at the peak plus and
    minus ten large-sample standard deviations of the order statistic,
    with a floor of $`4/\sqrt{n+1}`$ in the standardized units of
    device 4. The peak then sits inside a finite piece. The inner
    integral of a covariance is split the same way, around the peak of
    $`X_{(s)}`$. Its outer integral is split at the split points of both
    order statistics, because the inner quadrature changes its pieces
    where the outer variable crosses them.
4.  **Integration in standardized units, each moment about its own
    location.** The integrals are taken in
    $`t = \{x - Q(1/2)\}/\tilde\sigma`$, where $`Q`$ is the quantile
    function of the standard member and $`\tilde\sigma`$ is its
    interquartile range divided by that of the standard normal. So
    $`\tilde\sigma = 1`$ for the normal. The scaling puts the split
    floor and the integrator’s absolute tolerances in the distribution’s
    own units, so a distribution supplied on any scale is computed as
    accurately as on scale one. The moments of $`X_{(r)}`$ are taken
    about its own location $`c_r`$, the point $`Q(r/(n+1))`$ in these
    units, which is also a split point. Without this centering, the
    variance of an order statistic far from the median would be the
    small difference of two moments of size $`c_r^2`$, and it would lose
    its precision at large $`n`$. Taken about $`c_r`$, the variances of
    the uniform and the exponential order statistics we tested keep
    about ten significant digits up to $`n = 10^6`$, and every
    correlation is reproduced within about $`10^{-10}`$ in absolute
    terms. The one quantity that keeps fewer digits is the covariance of
    two nearly independent order statistics, such as the two extremes,
    whose correlation is of order $`1/n`$. Its relative error is the
    error of its correlation divided by the correlation, up to the
    relative errors of the two variances. The help page of
    [`os_moments()`](https://zwulab.github.io/osmoments/reference/os_moments.md)
    states the measured cases. Location and scale are applied after the
    integration.

The tolerances decide when
[`integrate()`](https://rdrr.io/r/stats/integrate.html) stops. The
single integrals and the inner integrals run at a relative tolerance of
$`10^{-11}`$, and the outer integrals at $`10^{-9}`$. The inner
tolerance is a hundred times tighter, so that the noise of the inner
quadrature stays below what the outer rule resolves. The absolute
tolerances and the acceptance thresholds are multiplied by
$`\min(1, w_r)`$, where $`w_r`$ is ten large-sample standard deviations
of the order statistic. The reason is that a narrow order statistic has
small moments: the variance of the minimum of a million uniform
observations is $`7 \times 10^{-12}`$ in standardized units, below a
fixed tolerance that suits the normal. The inner integral’s absolute
tolerance is also divided by $`\max(1, |y - c_r|)`$, the outer factor
that multiplies it, so that the far tail of a heavy-tailed parent is not
lost. With this, a lognormal with `sdlog = 2` is accurate to about
$`10^{-9}`$ in relative terms, rather than $`10^{-11}`$.

[`integrate()`](https://rdrr.io/r/stats/integrate.html) may report a
roundoff or subdivision problem on a piece of the range. The package
accepts such a piece when its error estimate is below $`10^{-10}`$,
multiplied by the same factor $`\min(1, w_r)`$ as the tolerances.
Otherwise it stops the computation with an error. A piece on an infinite
range that [`integrate()`](https://rdrr.io/r/stats/integrate.html)
reports as probably divergent is accepted only when its value is below
that scaled $`10^{-10}`$ as well.

The closed forms used by default are the beta and Dirichlet moments of
the uniform, the exponential’s representation as a sum of independent
exponential spacings, the digamma and trigamma forms of the logistic
means and variances, and the Gumbel maximum, whose law is a shifted
Gumbel.

## Appendix B: the checks behind the numbers

This appendix lists the checks that run inside the package and the
checks that its tests assert. Their purpose is that a quadrature failure
shows up as a residual or a failed test, rather than as a wrong number.
A variance that comes out zero or negative always stops the computation
with an error. When `check = TRUE`, the default, two checks run if their
conditions are met, and their residuals are stored in `checks`. A
residual above $`10^{-6}`$ raises a warning, and so does a covariance
matrix with a negative eigenvalue beyond rounding. These checks catch
some failures of the quadrature and certify none, because a failure that
every integral shares can satisfy them.

- **Reflection symmetry**, for a symmetric distribution and a set of
  ranks that is closed under $`r \mapsto n + 1 - r`$. The means are
  antisymmetric and the second moments symmetric about the median. The
  integrators do not impose this, so a violation signals a quadrature
  failure. For such rank sets, the covariances are computed for half of
  the pairs and reflected.
- **The permutation identities**, when every rank is requested. The
  order statistics are then a permutation of the sample, so
  $`\sum_r \alpha_r = n E[X]`$, $`\sum_r E[X^2_{(r)}] = n E[X^2]`$ and
  $`\sum_{r,s} B_{rs} = n\,\mathrm{Var}(X)`$, for any parent with a
  finite second moment. For the normal, every row of $`B`$ also sums to
  $`\sigma^2`$, because the sample mean is independent of every
  deviation from it. That row identity does not hold for other parents
  in general. The exponential at $`n = 2`$ has row sums $`1/2`$ and
  $`3/2`$. A symmetric parent has unit row sums at $`n = 2`$ by
  symmetry, but the uniform fails at $`n = 3`$.

The package’s tests assert the following when `R CMD check` runs them.
The comparison with the earlier implementation, the second-last item, is
skipped on CRAN to save time.

- the uniform and exponential closed forms at the five ranks of a
  five-number summary and at all their pairs, for $`n = 20`$, $`100`$
  and $`1001`$, within $`10^{-10}`$.
- the logistic means and variances and the Gumbel maximum within
  $`10^{-10}`$, and the $`n = 2`$ covariances of the logistic and the
  Gumbel, which follow from $`X_{(1)} X_{(2)} = X_1 X_2`$.
- the variance of the Cauchy median of five against the integral on the
  uniform scale, within $`10^{-9}`$.
- the permutation identities for all five families by quadrature at
  $`n = 5`$, $`9`$ and $`13`$, within $`10^{-9}`$, and the unit row sums
  for the normal.
- the tables of Tietjen, Kahaner and Beckman (1976) at the five-number
  ranks for $`n = 5`$ to $`49`$, which is 240 printed entries. One of
  them is illegible in the available scan and is skipped. The package
  agrees with 206 of the others within the table’s certified accuracy of
  $`10^{-7}`$. It differs by more than $`10^{-7}`$ from 33 quartile
  entries at $`n = 29`$ to $`49`$, by up to $`1.4 \times 10^{-4}`$. We
  evaluated those 33 entries once more in the probability scale, through
  the beta density of the uniform order statistic, which shares nothing
  with the package’s integrand. That evaluation agrees with the package
  within $`4 \times 10^{-16}`$, so we take the 33 entries to be errors
  of the printed table.
- Harter’s (1960) mean and variance of the normal range at 24 sample
  sizes from 5 to 97, within $`2 \times 10^{-8}`$.
- a gamma distribution with shape 2, supplied through
  [`os_dist()`](https://zwulab.github.io/osmoments/reference/os_dist.md),
  against its exact means and variances at every rank for $`n = 5`$ and
  $`10`$, which are finite sums of rationals, within $`10^{-10}`$.
- the quadrature against the uniform and exponential closed forms at the
  minimum, the lower quartile, the median and the maximum for
  $`n = 10^4`$, $`10^5`$ and $`10^6`$. The variances and the covariances
  of neighboring ranks agree within $`10^{-9}`$ in relative terms. The
  covariance of the two extremes agrees within $`10^{-9}`$ at
  $`n = 10^4`$ and within $`10^{-3}`$ at $`n = 10^5`$ and $`10^6`$.
- a normal supplied with standard deviations from $`10^{-3}`$ to
  $`10^6`$ against the standard normal’s moments scaled, within
  $`10^{-11}`$ in relative terms.
- the minimum of 10 and of 100 observations from a beta(1/2, 1/2)
  against a quadrature in the probability scale, within $`10^{-10}`$ in
  relative terms, and the minimum of a thousand and of $`10^5`$ refused
  with an error.
- the covariance of the two order statistics of two beta(1/2, 1/2)
  observations against its exact value $`4/\pi^4`$, within $`10^{-9}`$.
- a mixture of two narrow normals refused by the mass check, and the
  range of two uniform observations at `location = 1e16` computed
  exactly.
- the covariance of the two order statistics of two lognormal
  observations with `sdlog = 2` against its exact value, within
  $`2 \times 10^{-7}`$ on a value of 38.8.
- the normal at the five-number ranks against an earlier implementation
  of the same method for the normal, written separately before the
  package, at 39 sample sizes from 5 to 1001, within $`10^{-10}`$.
- the quantile weights against
  [`quantile()`](https://rdrr.io/r/stats/quantile.html) for all nine
  types, and the BLUE against its generalized least squares definition.

## References

David, H. A. and Nagaraja, H. N. (2003). *Order Statistics*, third
edition. John Wiley and Sons, Hoboken. <doi:10.1002/0471722162>.

Harter, H. L. (1960). Tables of range and studentized range. *The Annals
of Mathematical Statistics*, 31(4), 1122-1147.
<doi:10.1214/aoms/1177705684>.

Hyndman, R. J. and Fan, Y. (1996). Sample quantiles in statistical
packages. *The American Statistician*, 50(4), 361-365.
<doi:10.1080/00031305.1996.10473566>.

Lloyd, E. H. (1952). Least-squares estimation of location and scale
parameters using order statistics. *Biometrika*, 39(1-2), 88-95.
<doi:10.1093/biomet/39.1-2.88>.

Royston, J. P. (1982). Algorithm AS 177: Expected normal order
statistics (exact and approximate). *Applied Statistics*, 31(2),
161-165. <doi:10.2307/2347982>.

Tietjen, G. L., Kahaner, D. K. and Beckman, R. J. (1976). Variances and
covariances of the normal order statistics for sample sizes 2 to 50.
Report LA-6482-MS, Los Alamos Scientific Laboratory, Los Alamos, New
Mexico.

Wu, Z. and Tong, T. (2026). osmoments: Exact means, variances and
covariances of order statistics in R. Manuscript.
