# Lloyd's best linear unbiased estimators of location and scale

For a location-scale family with standard member \\Z\\, the order
statistics at the ranks of `mom` satisfy \\X\_{(r_i)} = \mu + \sigma
Z\_{(r_i)}\\, with \\E\[Z\_{(r_i)}\] = \alpha_i\\ and
\\\mathrm{Cov}(Z\_{(r_i)}, Z\_{(r_j)}) = B\_{ij}\\. Lloyd's estimators
are the generalized least squares estimators of \\(\mu, \sigma)\\ in
this linear model: \\\hat\mu = \bm w\_\mu'\bm X\\ and \\\hat\sigma = \bm
w\_\sigma'\bm X\\ with the weights of the Details. They are unbiased and
have the smallest variance among linear unbiased estimators from these
ranks.

## Usage

``` r
os_blue(mom)
```

## Arguments

- mom:

  An object of class `"osmoments"` computed with `cov = TRUE` at two or
  more ranks. The weights are computed from its `standard` component,
  the moments of the standard member, so they do not depend on the
  `location` and `scale` it was computed with.

## Value

A list with `w_location` and `w_scale` (the weight vectors, named by
rank), `var_location`, `var_scale` and `cov` (the variances and the
covariance of the two estimators per unit of the squared scale
parameter, so that \\\mathrm{Var}(\hat\mu) = \sigma^2 \cdot\\
`var_location`), `ranks`, `n` and `dist`.

## Details

With \\\bm\Omega = \bm B^{-1}\\, \\a\_{11} = \mathbf 1'\bm\Omega\mathbf
1\\, \\a\_{1\alpha} = \mathbf 1'\bm\Omega\bm\alpha\\,
\\a\_{\alpha\alpha} = \bm\alpha'\bm\Omega\bm\alpha\\ and \\\Delta =
a\_{11} a\_{\alpha\alpha} - a\_{1\alpha}^2\\, \$\$\bm w\_\mu' =
(a\_{\alpha\alpha}\mathbf 1' - a\_{1\alpha}\bm\alpha')\bm\Omega/\Delta,
\qquad \bm w\_\sigma' = (a\_{11}\bm\alpha' - a\_{1\alpha}\mathbf
1')\bm\Omega/\Delta,\$\$ \$\$\mathrm{Var}(\hat\mu) = \sigma^2
a\_{\alpha\alpha}/\Delta, \qquad \mathrm{Var}(\hat\sigma) = \sigma^2
a\_{11}/\Delta, \qquad \mathrm{Cov}(\hat\mu, \hat\sigma) = -\sigma^2
a\_{1\alpha}/\Delta.\$\$

**Computation.** The formulas above define the estimators. The function
computes them by whitening the design \\\[\mathbf 1, \bm\alpha\]\\ with
the Cholesky factor of \\\bm B\\ and solving the least-squares problem
by QR, which keeps its accuracy when \\\mathbf 1\\ and \\\bm\alpha\\ are
nearly collinear, as for two close ranks at a large \\n\\, where
\\\Delta\\ itself loses its digits. It then checks the four unbiasedness
conditions, \\\bm w\_\mu'\mathbf 1 = 1\\, \\\bm w\_\mu'\bm\alpha = 0\\,
\\\bm w\_\sigma'\mathbf 1 = 0\\ and \\\bm w\_\sigma'\bm\alpha = 1\\, and
warns when one fails by more than \\10^{-6}\\.

**Conditions.** The distribution of `mom` is treated as the standard
member of a location-scale family, which every
[`os_dist()`](https://zwulab.github.io/osmoments/reference/os_dist.md)
object is by construction; the required moments exist, which
[`os_moments()`](https://zwulab.github.io/osmoments/reference/os_moments.md)
has checked; at least two distinct ranks are present, so that the design
\\\[\mathbf 1, \bm\alpha\]\\ has rank two; and \\\bm B\\ is positive
definite. The scale estimated is the family's scale parameter: the
standard deviation for the normal, the scale of the logistic or Gumbel,
and the length of the support for the uniform.

## References

Lloyd, E. H. (1952). Least-squares estimation of location and scale
parameters using order statistics. *Biometrika*, 39(1-2), 88-95.

## See also

[`os_moments`](https://zwulab.github.io/osmoments/reference/os_moments.md),
[`os_lin_moments`](https://zwulab.github.io/osmoments/reference/os_lin_moments.md)

## Examples

``` r
m <- os_moments(9, ranks = c(1, 3, 5, 7, 9))
b <- os_blue(m)
b$w_location                     # sums to one
#>        r1        r3        r5        r7        r9 
#> 0.1457682 0.2382168 0.2320301 0.2382168 0.1457682 
sqrt(b$var_location)             # SD of the location estimate per unit sigma
#> [1] 0.3375654
```
