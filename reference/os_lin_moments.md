# Moments of linear combinations of order statistics

For a statistic \\T = \sum_i w_i X\_{(r_i)}\\ over the ranks of an
`osmoments` object, \\E\[T\] = \bm w'\bm\alpha\\ and \\\mathrm{Var}(T) =
\bm w'\bm B\bm w\\, and for several statistics with the weights in the
columns of a matrix \\\bm W\\, the means are \\\bm W'\bm\alpha\\ and the
covariance matrix is \\\bm W'\bm B\bm W\\.

## Usage

``` r
os_lin_moments(w, mom)
```

## Arguments

- w:

  A numeric vector with one weight per rank of `mom`, or a matrix with
  one such column per statistic.

- mom:

  An object of class `"osmoments"` computed with `cov = TRUE`, whose
  ranks are those the weights refer to, in order.

## Value

A list with `mean` (one value per statistic), `var` (the variances),
`cov` (the covariance matrix of the statistics, one by one for a single
statistic) and `ranks`. The statistics are named by the column names of
`w` when it has them.

## See also

[`os_moments`](https://zwulab.github.io/osmoments/reference/os_moments.md),
[`os_quantile_weights`](https://zwulab.github.io/osmoments/reference/os_quantile_weights.md)

## Examples

``` r
# the even-sample median of a normal sample of size 10
qw <- os_quantile_weights(0.5, n = 10)
m <- os_moments(10, ranks = qw$ranks)
os_lin_moments(qw$W, m)
#> $mean
#>          50% 
#> 2.775558e-17 
#> 
#> $var
#>       50% 
#> 0.1383264 
#> 
#> $cov
#>           50%
#> 50% 0.1383264
#> 
#> $ranks
#> [1] 5 6
#> 
# the range of an exponential sample of size 20
m <- os_moments(20, ranks = c(1, 20), dist = "exp")
os_lin_moments(c(-1, 1), m)
#> $mean
#> [1] 3.54774
#> 
#> $var
#> [1] 1.593663
#> 
#> $cov
#>          [,1]
#> [1,] 1.593663
#> 
#> $ranks
#> [1]  1 20
#> 
```
