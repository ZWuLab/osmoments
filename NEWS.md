# osmoments 0.1.0

* First release.
* `os_moments()` computes the exact means, variances and covariances of order
  statistics at any ranks and any sample size, for the normal, uniform,
  exponential, logistic and Gumbel distributions and for a user-supplied
  distribution built with `os_dist()`. Closed forms are returned where the
  family has them (the uniform, the exponential, the logistic single moments,
  the Gumbel maximum) and the exact integrals are evaluated by adaptive
  quadrature elsewhere.
* `os_quantile_weights()` gives the ranks and weights of the sample quantiles
  under the nine definitions of `quantile()`, `os_lin_moments()` the moments
  of any linear combination of order statistics, and `os_blue()` Lloyd's best
  linear unbiased estimators of location and scale.
