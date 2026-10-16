# osmoments: Exact Moments of Order Statistics

Computes the exact means, variances and covariances of the order
statistics \\X\_{(1)} \le \dots \le X\_{(n)}\\ of an i.i.d. sample from
a continuous distribution, at any ranks and any sample size, by
numerical evaluation of the defining integrals, with closed forms where
the distribution has them. Five families are built in (normal, uniform,
exponential, logistic and Gumbel), and any other continuous distribution
enters through its distribution, density and quantile functions via
[`os_dist`](https://zwulab.github.io/osmoments/reference/os_dist.md).

## Details

The main function is
[`os_moments`](https://zwulab.github.io/osmoments/reference/os_moments.md),
which returns the mean vector and the covariance matrix of the order
statistics at the requested ranks. On top of it,
[`os_quantile_weights`](https://zwulab.github.io/osmoments/reference/os_quantile_weights.md),
[`os_lin_moments`](https://zwulab.github.io/osmoments/reference/os_lin_moments.md)
and [`os_blue`](https://zwulab.github.io/osmoments/reference/os_blue.md)
give the moments of sample quantiles and of any linear combination of
order statistics, and Lloyd's best linear unbiased estimators of
location and scale. The vignette walks through all of it:
[`vignette("osmoments")`](https://zwulab.github.io/osmoments/articles/osmoments.md).
The paper of Wu and Tong (2026) describes the computation, its
verification against closed forms, identities and printed tables, and
the comparison with other software; `citation("osmoments")` gives the
reference.

## References

Wu, Z. and Tong, T. (2026). osmoments: Exact means, variances and
covariances of order statistics in R. Manuscript.

## See also

Useful links:

- <https://github.com/ZWuLab/osmoments>

- <https://zwulab.github.io/osmoments/>

- Report bugs at <https://github.com/ZWuLab/osmoments/issues>

## Author

**Maintainer**: Zheyang Wu <zheyangwu@wpi.edu>
