# ==============================================================================
# File: osmoments-package.R
# Purpose: Package-level documentation of osmoments and the imports from stats.

#' osmoments: Exact Moments of Order Statistics
#'
#' Computes the exact means, variances and covariances of the order statistics
#' \eqn{X_{(1)} \le \dots \le X_{(n)}} of an i.i.d. sample from a continuous
#' distribution, at any ranks and any sample size, by numerical evaluation of
#' the defining integrals, with closed forms where the distribution has them.
#' Five families are built in (normal, uniform, exponential, logistic and
#' Gumbel), and any other continuous distribution enters through its
#' distribution, density and quantile functions via \code{\link{os_dist}}.
#'
#' The main function is \code{\link{os_moments}}, which returns the mean vector
#' and the covariance matrix of the order statistics at the requested ranks.
#' On top of it, \code{\link{os_quantile_weights}}, \code{\link{os_lin_moments}}
#' and \code{\link{os_blue}} give the moments of sample quantiles and of any
#' linear combination of order statistics, and Lloyd's best linear unbiased
#' estimators of location and scale. The vignette walks through all of it:
#' \code{vignette("osmoments")}. The paper of Wu and Tong (2026) describes the
#' computation, its verification against closed forms, identities and printed
#' tables, and the comparison with other software; \code{citation("osmoments")}
#' gives the reference.
#'
#' @references
#' Wu, Z. and Tong, T. (2026). osmoments: Exact means, variances and covariances
#' of order statistics in R. Manuscript.
#'
#' @keywords internal
#' @importFrom stats integrate pnorm dnorm qnorm punif dunif qunif pexp dexp qexp
#'   plogis dlogis qlogis
"_PACKAGE"
