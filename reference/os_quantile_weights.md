# Ranks and weights of the sample quantile definitions

A sample quantile under any of the nine definitions of
[`quantile`](https://rdrr.io/r/stats/quantile.html) is a fixed linear
combination of two adjacent order statistics, \\(1 - h) X\_{(j)} + h
X\_{(j+1)}\\, with the ranks clipped to \\1, \dots, n\\. This function
returns, for each probability in `p`, the ranks and the weights, so that
the moments of the sample quantile follow from the moments of the order
statistics through
[`os_lin_moments`](https://zwulab.github.io/osmoments/reference/os_lin_moments.md).
The even-sample median, for example, is \\(X\_{(n/2)} +
X\_{(n/2+1)})/2\\ under type 7.

## Usage

``` r
os_quantile_weights(p, n, type = 7)
```

## Arguments

- p:

  The probabilities, in \\\[0, 1\]\\.

- n:

  The sample size, an integer \\\ge 1\\.

- type:

  The definition, an integer from 1 to 9 as in
  [`quantile`](https://rdrr.io/r/stats/quantile.html); 7 is R's default.

## Value

A list with `ranks` (the sorted union of the ranks that carry a weight,
for every probability), `W` (a matrix with one row per rank and one
column per probability, each column summing to one, with the rows named
by rank and the columns as
[`quantile()`](https://rdrr.io/r/stats/quantile.html) names them), `n`,
`p` and `type`.

## Details

The pair \\(j, h)\\ is computed by the same arithmetic as
`stats:::quantile.default`. For type 7, \\j = \lfloor 1 + (n-1)p
\rfloor\\ and \\h = 1 + (n-1)p - j\\. For types 4 to 9 except 7, with
\\(a, b)\\ equal to \\(0, 1)\\, \\(1/2, 1/2)\\, \\(0, 0)\\, \\(1/3,
1/3)\\ and \\(3/8, 3/8)\\ for types 4, 5, 6, 8 and 9, \\m = a + p(n +
1 - a - b)\\, \\j = \lfloor m \rfloor\\ (with a fuzz of four machine
epsilons) and \\h = m - j\\. For types 1 to 3, \\m = np\\ (types 1 and
2) or \\m = np - 1/2\\ (type 3), \\j = \lfloor m \rfloor\\, and \\h\\ is
1 if \\m \> j\\ and 0 otherwise (type 1), 1 if \\m \> j\\ and 1/2
otherwise (type 2), and 0 if \\m = j\\ with \\j\\ even and 1 otherwise
(type 3). The ranks \\j\\ and \\j + 1\\ are then clipped to \\\[1,
n\]\\, and a rank whose weight is zero is dropped. On a sorted sample
the weighted sum reproduces
[`quantile()`](https://rdrr.io/r/stats/quantile.html) exactly, which the
package's tests check for all nine types.

## References

Hyndman, R. J. and Fan, Y. (1996). Sample quantiles in statistical
packages. *The American Statistician*, 50(4), 361-365.

## See also

[`os_lin_moments`](https://zwulab.github.io/osmoments/reference/os_lin_moments.md),
[`quantile`](https://rdrr.io/r/stats/quantile.html)

## Examples

``` r
qw <- os_quantile_weights(c(0.25, 0.5, 0.75), n = 10)
qw$ranks
#> [1] 3 4 5 6 7 8
qw$W
#>     25% 50%  75%
#> r3 0.75 0.0 0.00
#> r4 0.25 0.0 0.00
#> r5 0.00 0.5 0.00
#> r6 0.00 0.5 0.00
#> r7 0.00 0.0 0.25
#> r8 0.00 0.0 0.75
m <- os_moments(10, ranks = qw$ranks)
os_lin_moments(qw$W, m)          # moments of the three quartiles of a normal sample
#> $mean
#>           25%           50%           75% 
#> -5.859855e-01  2.775558e-17  5.859855e-01 
#> 
#> $var
#>       25%       50%       75% 
#> 0.1584864 0.1383264 0.1584864 
#> 
#> $cov
#>            25%       50%        75%
#> 25% 0.15848637 0.1030256 0.06910974
#> 50% 0.10302557 0.1383264 0.10302557
#> 75% 0.06910974 0.1030256 0.15848637
#> 
#> $ranks
#> [1] 3 4 5 6 7 8
#> 
```
