# Exact moments of order statistics

Computes the means \\\alpha_i = E\[X\_{(r_i)}\]\\, the variances and the
covariance matrix \\B\_{ij} = \mathrm{Cov}(X\_{(r_i)}, X\_{(r_j)})\\ of
the order statistics \\X\_{(1)} \le \dots \le X\_{(n)}\\ of an i.i.d.
sample of size `n` from a continuous distribution, at the ranks `ranks`,
by numerical evaluation of the exact integrals, or from a closed form
where the family has one.

## Usage

``` r
os_moments(
  n,
  ranks = seq_len(n),
  dist = "norm",
  location = 0,
  scale = 1,
  cov = TRUE,
  method = c("auto", "integrate", "closed"),
  rel.tol = 1e-11,
  check = TRUE
)
```

## Arguments

- n:

  The sample size, an integer \\\ge 1\\.

- ranks:

  The ranks, strictly increasing integers in \\1, \dots, n\\. The
  default is every rank, which gives the full covariance matrix.

- dist:

  The distribution: one of `"norm"`, `"unif"`, `"exp"`, `"logis"`,
  `"gumbel"`, or an object made by
  [`os_dist`](https://zwulab.github.io/osmoments/reference/os_dist.md).
  See Details for the standard members.

- location, scale:

  The location and scale applied to the standard member: \\X =
  \mathrm{location} + \mathrm{scale} \cdot Z\\, with `scale > 0`. The
  means become \\\mathrm{location} + \mathrm{scale}\\\alpha\\ and the
  covariances \\\mathrm{scale}^2 B\\.

- cov:

  `TRUE` for the full covariance matrix, `FALSE` for the means and the
  variances only, which need single integrals and no double integral.

- method:

  `"auto"` takes the closed form where the family has one and the
  quadrature elsewhere, `"integrate"` forces the quadrature, `"closed"`
  requires a closed form for every requested quantity and stops
  otherwise.

- rel.tol:

  The relative tolerance of the quadrature. The single integrals and the
  inner integral of a product moment run at `rel.tol`, the outer
  integral at `100 * rel.tol`.

- check:

  Whether to run the run-time checks of the Details and warn when a
  residual exceeds \\10^{-6}\\.

## Value

An object of class `"osmoments"`, a list with the components `n`,
`ranks`, `alpha` (the means, named by rank), `var` (the variances), `B`
(the covariance matrix, `NULL` when `cov = FALSE`), `quad_err` (the
summed error estimates of the quadrature, see Details), `standard` (a
list with the `alpha`, `var` and `B` of the standard member, before
`location` and `scale` are applied, which
[`os_blue`](https://zwulab.github.io/osmoments/reference/os_blue.md)
reads), `dist` (the `os_dist` object), `location`, `scale`, `cov`,
`method`, `route` (`"closed"`, `"integrate"` or `"mixed"`), `checks`
(the residuals of the run-time checks, `NA` where a check did not apply)
and `call`.

## Details

**The integrals.** With cdf \\F\\, density \\f\\ and \\y \< z\\,
\$\$f\_{r:n}(x) = \frac{n!}{(r-1)!\\(n-r)!} F(x)^{r-1} \\1 -
F(x)\\^{n-r} f(x),\$\$ \$\$f\_{r,s:n}(y, z) =
\frac{n!}{(r-1)!\\(s-r-1)!\\(n-s)!} F(y)^{r-1} \\F(z) - F(y)\\^{s-r-1}
\\1 - F(z)\\^{n-s} f(y) f(z),\$\$ and the moments are \$\$\mu\_{r:n} =
\int x f\_{r:n}(x)\\dx, \qquad \mu^{(2)}\_{r:n} = \int x^2
f\_{r:n}(x)\\dx, \qquad \mu\_{r,s:n} = \int\\\\\int\_{y \< z} y z
f\_{r,s:n}(y, z)\\dy\\dz,\$\$ with \\B\_{ij} = \mu\_{r_i, r_j:n} -
\alpha_i \alpha_j\\. The integrals are evaluated by adaptive quadrature
([`integrate`](https://rdrr.io/r/stats/integrate.html)) with four
devices: the integrands on the log scale, so that the powers of the
tails never underflow; \\\log\\F(z) - F(y)\\\\ formed on the smaller
tail, so that nothing cancels; the range split at the spike of the order
statistic's density, at \\c_r \pm h_r\\ with \\c_r = Q(p_r)\\, \\p_r =
r/(n+1)\\ and \\h_r = 10 s_r / f(c_r) + 4\tilde\sigma/\sqrt{n+1}\\,
\\s_r = \sqrt{p_r(1-p_r)/(n+2)}\\, so that the adaptive rule sees the
narrow peak at large \\n\\ at the end of a piece rather than between its
nodes, the outer range of a product moment being split at the points of
both ranks; and integration in standardized units of the standard
member, \\t = \\x - Q(1/2)\\/\tilde\sigma\\, where \\\tilde\sigma\\ is
the spread of the distribution, its interquartile range divided by that
of the standard normal (so \\\tilde\sigma = 1\\ for the normal), with
each moment taken about its own order statistic's location \\c_r\\. The
scaling puts the split floor and the integrator's absolute tolerances in
the distribution's own units, so that a distribution supplied on any
scale is computed as accurately as on the scale one. The centering at
\\c_r\\ keeps a variance from being the small difference of two moments
of size \\(c_r - c)^2\\: taken about the median \\c\\, the variance of
the uniform minimum kept seven digits at \\n = 1000\\ and none at \\n =
10^5\\.

**The standard members and the closed forms.** `"norm"` is \\N(0, 1)\\;
`"unif"` is the uniform on \\\[0, 1\]\\, so that `location` is the lower
end and `scale` the length; `"exp"` has unit rate, so that `scale` is
\\1/\mathrm{rate}\\; `"logis"` is the standard logistic; `"gumbel"` is
the standard Gumbel in the maximum form, \\F(x) = \exp(-e^{-x})\\. The
uniform and the exponential have closed forms for every moment, the
logistic for the means and variances, the Gumbel for the mean and
variance of the maximum. With `method = "auto"` these are returned and
the rest is quadrature.

**The moment condition.** A moment is computed only if it exists: with
the tail indices \\(\nu_L, \nu_U)\\ of the distribution,
\\E\|X\_{(r)}\|^2 \< \infty\\ if and only if \\\nu_L r \> 2\\ and
\\\nu_U (n - r + 1) \> 2\\, and a rank that fails stops the computation
with an error. The built-in families always pass.

**The error estimate.** `quad_err` is the sum of the error estimates
that [`integrate()`](https://rdrr.io/r/stats/integrate.html) reports for
the single integrals and for the outer integrals of the product moments,
in the units of the standard member before `location` and `scale` are
applied. It is an estimate, not a bound: the inner integral's estimate
is not carried into the outer one, and the errors of the means are not
propagated into the covariances. The accuracy evidence is the agreement
with the closed forms, the identities and the printed tables, documented
in the package's tests.

**Accuracy at large n.** With `method = "integrate"`, the means and the
variances of the uniform and the exponential at the minimum, the
quartiles, the median and the maximum agree with their closed forms to
within \\10^{-9}\\ in relative terms at \\n = 10^4\\, \\10^5\\ and
\\10^6\\, and so does the correlation of every pair of them, in absolute
terms. Up to the relative errors of the two variances, the relative
error of a covariance is the error of its correlation divided by the
correlation, so the covariance of two nearly independent order
statistics keeps fewer digits: the two extremes, whose correlation is of
order \\1/n\\, keep five to ten significant digits of their covariance
at those sizes, depending on the parent and on \\n\\. The package's
tests assert the variances and the covariances of adjacent ranks at
\\10^{-9}\\ and the extremes' covariance at \\10^{-3}\\ beyond \\n =
10^4\\. A parent as heavy-tailed as the lognormal with log-SD 2 is
accurate to about \\10^{-9}\\ in relative terms rather than
\\10^{-11}\\, because the far tail of the inner integral of a covariance
is resolved at a fixed tolerance. These are the cases that were
measured, not a bound for every parent.

**What stops with an error.** The quadrature follows the peak of each
order statistic's density and the tails of the parent. It does not
follow a second mode: for a mixture of two narrow normals, the mode that
holds a tenth of the mass fell between the nodes, and the moments came
back as those of the other mode. The package therefore integrates the
density of every order statistic it computes by quadrature and stops
when the mass is not 1 within \\10^{-7}\\, naming the rank and the mass
it found. It also stops when a piece of the range cannot be resolved,
which happens for a moment that does not exist when the tail indices
were not declared, for a tail as heavy as a lognormal with log-SD 3 or a
Weibull with shape 0.3, for the extremes of a parent whose density is
infinite at an endpoint of its support once they concentrate against
that endpoint (the minimum of a beta(1/2, 1/2) is computed to \\n =
100\\ and refused from \\n = 1000\\ on), and for a variance that comes
out zero or negative. A wrong number is never returned knowingly, but
these checks detect some failures and certify none.

**The run-time checks.** For a symmetric distribution and a rank set
closed under \\r \mapsto n + 1 - r\\, the means are antisymmetric and
the second moments symmetric about the median, which the integrators do
not impose; the product moments are then computed for the pairs with
\\r + s \le n + 1\\ only and reflected. For a complete rank set, the
order statistics are a permutation of the sample, so \\\sum_r \alpha_r =
n E\[X\]\\, \\\sum_r \mu^{(2)}\_{r:n} = n E\[X^2\]\\ and \\\sum\_{r,s}
B\_{rs} = n\\\mathrm{Var}(X)\\ for any parent with a finite second
moment, and for the normal parent (recognized by its functions `pnorm`,
`dnorm` and `qnorm`, not by its name) every row of \\B\\ sums to
\\\sigma^2\\ as well. The residuals are stored in `checks`, in the
standardized units of the integration (the standard member centered at
its median and scaled by its spread, so that the threshold means the
same on every scale), and a residual above \\10^{-6}\\ raises a warning,
as does a covariance matrix with a negative eigenvalue beyond rounding.
The checks detect some failures of the quadrature and do not certify its
accuracy: a failure that every integral shares, such as a peak missed in
the same way at every rank, can satisfy them.

**Cost.** A single integral takes a fraction of a millisecond and a
double integral a few hundredths of a second. When the covariances are
not closed forms, a full matrix needs one double integral per pair,
\\n(n-1)/2\\ of them, halved for a symmetric family.

## References

Wu, Z. and Tong, T. (2026). osmoments: Exact means, variances and
covariances of order statistics in R. Manuscript. The computation, its
verification and the comparison with other software.

Royston, J. P. (1982). Algorithm AS 177: Expected normal order
statistics (exact and approximate). *Applied Statistics*, 31(2),
161-165.

## See also

[`os_dist`](https://zwulab.github.io/osmoments/reference/os_dist.md)

## Examples

``` r
m <- os_moments(5)            # all ranks of n = 5, standard normal
m$alpha
#>            r1            r2            r3            r4            r5 
#> -1.162964e+00 -4.950190e-01  5.413025e-29  4.950190e-01  1.162964e+00 
m$B
#>            r1        r2        r3        r4         r5
#> r1 0.44753407 0.2243310 0.1481477 0.1057720 0.07421527
#> r2 0.22433096 0.3115190 0.2084354 0.1499427 0.10577198
#> r3 0.14814773 0.2084354 0.2868337 0.2084354 0.14814773
#> r4 0.10577198 0.1499427 0.2084354 0.3115190 0.22433096
#> r5 0.07421527 0.1057720 0.1481477 0.2243310 0.44753407
os_moments(20, ranks = c(1, 5, 10, 20), dist = "exp", scale = 2)
#> Exact moments of order statistics: exp (location 0, scale 2), n = 20, ranks 1 5 10 20
#> alpha (means):
#>        r1        r5       r10       r20 
#> 0.1000000 0.5590213 1.3375430 7.1954790 
#> B (covariance matrix):
#>       r1         r5        r10        r20
#> r1  0.01 0.01000000 0.01000000 0.01000000
#> r5  0.01 0.06289184 0.06289184 0.06289184
#> r10 0.01 0.06289184 0.18558210 0.18558210
#> r20 0.01 0.06289184 0.18558210 6.38465300
#> route: closed; quadrature error estimate: 0
#> checks: reflection n/a; identities n/a
os_moments(100, ranks = c(1, 100), dist = "gumbel", cov = FALSE)
#> Exact moments of order statistics: gumbel (location 0, scale 1), n = 100, ranks 1 100
#> alpha (means):
#>        r1      r100 
#> -1.618284  5.182386 
#> var (variances; cov = FALSE):
#>         r1       r100 
#> 0.05442057 1.64493400 
#> route: mixed; quadrature error estimate: 1.1e-13
#> checks: reflection n/a; identities n/a
```
