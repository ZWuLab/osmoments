# osmoments

The means, variances and covariances of order statistics, the sorted
values of a sample, are what the standard error of a sample quantile,
the variance of a range, or the best linear estimator of a mean from a
five-number summary are built from. Printed tables of the covariances of
normal order statistics stop at n = 50 in the sources we found, and the
software we found gives the moments of single order statistics or
approximates the covariances of the normal by a series. osmoments
computes the exact means, variances and covariances of order statistics,
at any ranks and any sample size, by numerical evaluation of the
defining integrals. Five families are built in, the normal, uniform,
exponential, logistic and Gumbel, and any other continuous distribution
enters through its distribution, density and quantile functions. On top
of the moments it gives those of sample quantiles and of any linear
combination of order statistics, and Lloyd’s best linear unbiased
estimators of location and scale.

The package has two requirements: the distribution must be continuous,
and a moment is computed only if it exists, which the package checks
rank by rank from the declared tail indices. The integrals are evaluated
on the log scale with the integration range split around the peak of
each order statistic’s density, so a few ranks at a sample size in the
thousands cost about as much as at n = 20, while a full covariance
matrix grows with the square of n. Where a family has closed forms, they
are returned instead of the quadrature.

The package’s tests check the quadrature against four kinds of evidence.
The closed forms of the uniform and the exponential are reproduced at
the minimum, quartiles, median and maximum and all their pairs for n =
20, 100 and 1001 within 1e-10. The ten-decimal normal tables of Tietjen,
Kahaner and Beckman (1976) are reproduced within their stated accuracy
of 1e-7 on the entries that are printed correctly. Harter’s (1960)
moments of the normal range are reproduced within 2e-8. The identities
that hold on a complete set of ranks are satisfied within 1e-9. The
vignette states each check and its limits.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("ZWuLab/osmoments", build_vignettes = TRUE)
```

## Example

``` r

library(osmoments)

# the means and the covariance matrix of the minimum, the quartiles and the maximum
# of 25 standard normal observations
m <- os_moments(n = 25, ranks = c(1, 7, 13, 19, 25), dist = "norm")
m$alpha
m$B

# the same ranks for a logistic distribution with location 10 and scale 2
os_moments(25, ranks = c(1, 7, 13, 19, 25), dist = "logis", location = 10, scale = 2)

# the standard deviation of the median of an even-sized normal sample
qw <- os_quantile_weights(p = 0.5, n = 10)
sqrt(os_lin_moments(qw$W, os_moments(10, ranks = qw$ranks))$var)

# Lloyd's BLUE of the mean and the standard deviation from a five-number summary
os_blue(m)
```

A distribution that is not built in enters through
[`os_dist()`](https://zwulab.github.io/osmoments/reference/os_dist.md),
with its distribution, density and quantile functions and the tail
indices that decide which moments exist:

``` r

t5 <- os_dist(
  p = function(q, lower.tail = TRUE, log.p = FALSE) pt(q, 5, lower.tail = lower.tail, log.p = log.p),
  d = function(x, log = FALSE) dt(x, 5, log = log),
  q = function(p) qt(p, 5),
  symmetric = TRUE, tail = c(5, 5), name = "t5"
)
os_moments(5, dist = t5)
```

The vignette walks through the families, a user-supplied distribution,
the linear layer, what the result holds and how the integrals are
evaluated. It is at <https://zwulab.github.io/osmoments/> and, after the
installation above,
[`vignette("osmoments")`](https://zwulab.github.io/osmoments/articles/osmoments.md).

## Citation

The paper that describes the package, its verification and the
comparison with other software is

Wu, Z. and Tong, T. (2026). osmoments: Exact means, variances and
covariances of order statistics in R. Manuscript.

`citation("osmoments")` prints this reference.

## References

Harter, H. L. (1960). Tables of range and studentized range. The Annals
of Mathematical Statistics, 31(4), 1122-1147.

Tietjen, G. L., Kahaner, D. K. and Beckman, R. J. (1976). Variances and
covariances of the normal order statistics for sample sizes 2 to 50.
Report LA-6482-MS, Los Alamos Scientific Laboratory.

## License

GPL (\>= 3).

> Developed with AI assistance; see
> [AI-USE.md](https://zwulab.github.io/osmoments/AI-USE.md).
