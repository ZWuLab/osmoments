# Define a distribution for the moments of its order statistics

Builds the distribution object that
[`os_moments`](https://zwulab.github.io/osmoments/reference/os_moments.md)
takes in its `dist` argument, from the distribution function, the
density and the quantile function of a continuous distribution. The
object describes the standard member of a location-scale family:
[`os_moments()`](https://zwulab.github.io/osmoments/reference/os_moments.md)
computes the moments for this member and applies its own `location` and
`scale` arguments afterwards. A user who wants a \\t_5\\ distribution
with location 3 and scale 2 supplies the \\t_5\\ functions here and
passes `location = 3, scale = 2` to
[`os_moments()`](https://zwulab.github.io/osmoments/reference/os_moments.md).

## Usage

``` r
os_dist(
  p,
  d,
  q,
  support = c(-Inf, Inf),
  symmetric = FALSE,
  tail = c(Inf, Inf),
  name = "user"
)
```

## Arguments

- p:

  The distribution function, with the signature of R's `p*` functions:
  `p(q, lower.tail = TRUE, log.p = FALSE)`, vectorized in `q`. The
  `log.p` and `lower.tail` arguments must be honored, because the
  integrands are evaluated on the log scale in both tails.

- d:

  The density, `d(x, log = FALSE)`, vectorized in `x`, with `log`
  honored.

- q:

  The quantile function, `q(p)`, vectorized in `p`.

- support:

  Numeric vector of length two, the endpoints of the support, finite or
  infinite. They are the limits of integration.

- symmetric:

  Logical. `TRUE` if the density is symmetric about its median, that is
  \\Q(1 - u) = 2 Q(1/2) - Q(u)\\. A symmetric distribution lets
  [`os_moments()`](https://zwulab.github.io/osmoments/reference/os_moments.md)
  halve the number of double integrals of a full matrix and check the
  reflection symmetry of its result.

- tail:

  Numeric vector of length two, the tail indices \\(\nu_L, \nu_U)\\: the
  density decays like \\\|x\|^{-\nu_L - 1}\\ in the lower tail and like
  \\x^{-\nu_U - 1}\\ in the upper tail. `Inf` for a tail that decays
  faster than any power, or for a bounded tail. The indices decide which
  moments exist (see Details).

- name:

  A label for printing.

## Value

An object of class `"os_dist"`: a list with the components `name`, `p`,
`d`, `q`, `support`, `symmetric`, `tail`, `center` (the median
\\Q(1/2)\\) and `closed` (the closed-form moments the distribution has,
none for a user distribution).

## Details

**The moment condition.** For a parent whose density decays like
\\\|x\|^{-\nu - 1}\\, the density of \\X\_{(r)}\\ decays like
\\\|x\|^{-\nu r - 1}\\ below and like \\x^{-\nu (n - r + 1) - 1}\\
above, so \\E\|X\_{(r)}\|^k \< \infty\\ if and only if \\\nu_L r \> k\\
and \\\nu_U (n - r + 1) \> k\\.
[`os_moments()`](https://zwulab.github.io/osmoments/reference/os_moments.md)
applies this condition with \\k = 2\\ to every requested rank and
refuses a moment that does not exist. A Cauchy parent, `tail = c(1, 1)`,
has a median of five with a finite variance and extremes without a mean.

**Accuracy requirement.** The log-tails `p(x, log.p = TRUE)` and
`p(x, lower.tail = FALSE, log.p = TRUE)` must keep relative precision
deep in both tails, as R's base distributions do. A function that
computes `log(p(x))` from an underflowing probability loses the tail,
and the integrals with it. `os_dist()` checks that the three functions
honor their arguments at the median and are vectorized, and that `d`
integrates to \\1/2\\ between the quartiles, where `p` rises by \\1/2\\.
[`os_moments()`](https://zwulab.github.io/osmoments/reference/os_moments.md)
checks that the log-tails and the log-density are finite at the median
and at the split points of the smallest and the largest requested rank,
and that the density of every order statistic it integrates has mass 1.
None of these certifies accuracy.

**Location limit.** The integrals are taken about the median of the
standard member, and the integration range is split in units of its
spread (its interquartile range), so the scale on which the standard
member is supplied does not matter: the package's tests check a normal
supplied with standard deviations from \\10^{-3}\\ to \\10^6\\. Its
location does: a standard member whose own location is large, say
\\10^8\\ with unit scale, loses the digits of \\x - Q(1/2)\\ in double
precision and the computation stops with an error. Supply the standard
member and pass the location to
[`os_moments()`](https://zwulab.github.io/osmoments/reference/os_moments.md)
instead, which also keeps
[`os_blue()`](https://zwulab.github.io/osmoments/reference/os_blue.md)
exact.

**Common parents.** Each function wraps the corresponding `stats`
function with its shape fixed and passes `lower.tail`, `log.p` and `log`
by name, as in the examples. The declarations that differ between
parents are these.

|  |  |  |
|----|----|----|
| Parent | `support` | `tail` |
| Student \\t\_\nu\\ (`pt`, `dt`, `qt`), symmetric | `c(-Inf, Inf)` | `c(nu, nu)` |
| gamma (`pgamma`, `dgamma`, `qgamma`) | `c(0, Inf)` | `c(Inf, Inf)` |
| Weibull (`pweibull`, `dweibull`, `qweibull`) | `c(0, Inf)` | `c(Inf, Inf)` |
| lognormal (`plnorm`, `dlnorm`, `qlnorm`) | `c(0, Inf)` | `c(Inf, Inf)` |
| beta (`pbeta`, `dbeta`, `qbeta`) | `c(0, 1)` | `c(Inf, Inf)` |

The lognormal tail decays faster than any power, so every moment exists.
A distribution declared on `c(-Inf, Inf)` whose density is zero below
some point stops with an error that its log-tail is not finite, which is
the sign that `support` is missing.

## See also

[`os_moments`](https://zwulab.github.io/osmoments/reference/os_moments.md)

## Examples

``` r
# Student t with 5 degrees of freedom: tail index 5 in both tails
t5 <- os_dist(function(q, lower.tail = TRUE, log.p = FALSE)
                stats::pt(q, df = 5, lower.tail = lower.tail, log.p = log.p),
              function(x, log = FALSE) stats::dt(x, df = 5, log = log),
              function(p) stats::qt(p, df = 5),
              symmetric = TRUE, tail = c(5, 5), name = "t5")
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

# The Cauchy: the median of five has a variance, the extremes have no mean
cauchy <- os_dist(stats::pcauchy, stats::dcauchy, stats::qcauchy,
                  symmetric = TRUE, tail = c(1, 1), name = "cauchy")
os_moments(5, ranks = 3, dist = cauchy)
#> Exact moments of order statistics: cauchy (location 0, scale 1), n = 5, ranks 3
#> alpha (means):
#>            r3 
#> -1.028762e-17 
#> B (covariance matrix):
#>          r3
#> r3 1.221253
#> route: integrate; quadrature error estimate: 3.89e-14
#> checks: reflection 1.39e-17; identities n/a

# A gamma with shape 2: a semi-bounded support, every moment finite
gamma2 <- os_dist(function(q, lower.tail = TRUE, log.p = FALSE)
                    stats::pgamma(q, shape = 2, lower.tail = lower.tail, log.p = log.p),
                  function(x, log = FALSE) stats::dgamma(x, shape = 2, log = log),
                  function(p) stats::qgamma(p, shape = 2),
                  support = c(0, Inf), name = "gamma2")
os_moments(10, ranks = c(1, 5, 10), dist = gamma2, scale = 3)
#> Exact moments of order statistics: gamma2 (location 0, scale 3), n = 10, ranks 1 5 10
#> alpha (means):
#>        r1        r5       r10 
#>  1.398065  4.725263 13.868870 
#> B (covariance matrix):
#>            r1        r5        r10
#> r1  0.6842539 0.3962219  0.2986211
#> r5  0.3962219 2.0064850  1.5292810
#> r10 0.2986211 1.5292810 20.0825000
#> route: integrate; quadrature error estimate: 1.11e-10
#> checks: reflection n/a; identities n/a
```
