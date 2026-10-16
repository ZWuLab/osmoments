# ==============================================================================
# File: os_closed.R
# Purpose: The closed-form moments of order statistics that the built-in
#          families have, in the raw units of the standard member: the uniform
#          (beta and Dirichlet moments), the exponential (Renyi's
#          representation), the logistic single moments (digamma and trigamma)
#          and the Gumbel maximum. Each hook returns the mean, the variance or
#          the covariance, or NA where the family has no closed form for the
#          requested rank. os_moments() uses them when method = "auto" or
#          "closed" and checks them against the quadrature in the tests.
#
# INTERNAL HELPERS:
#   - .os_closed_unif_mean(), _var(), _cov()
#   - .os_closed_exp_mean(), _var(), _cov()
#   - .os_closed_logis_mean(), _var()
#   - .os_closed_gumbel_mean(), _var()

#################### INTERNAL HELPER FUNCTIONS ####################

# ---- uniform on [0, 1]: U_(r) ~ Beta(r, n - r + 1) ----------------------------
# E[U_(r)] = r/(n+1); Var = r(n-r+1)/{(n+1)^2 (n+2)};
# Cov(U_(r), U_(s)) = r(n-s+1)/{(n+1)^2 (n+2)} for r <= s, from the Dirichlet law
# of the spacings (E[U_(r) U_(s)] = r(s+1)/{(n+1)(n+2)}).

#' @keywords internal
#' @noRd
.os_closed_unif_mean <- function(r, n) r / (n + 1)

#' @keywords internal
#' @noRd
.os_closed_unif_var <- function(r, n) r * (n - r + 1) / ((n + 1)^2 * (n + 2))

#' @keywords internal
#' @noRd
.os_closed_unif_cov <- function(r, s, n) r * (n - s + 1) / ((n + 1)^2 * (n + 2))

# ---- exponential with unit rate: X_(r) = sum_{i <= r} E_i/(n - i + 1) ----------
# with independent unit exponentials E_i (Renyi's representation), so
# E[X_(r)] = sum 1/(n-i+1), Var = sum 1/(n-i+1)^2, and Cov(X_(r), X_(s)) = Var(X_(r))
# for r <= s because X_(s) - X_(r) is independent of X_(r).

#' @keywords internal
#' @noRd
.os_closed_exp_mean <- function(r, n) sum(1 / (n - seq_len(r) + 1))

#' @keywords internal
#' @noRd
.os_closed_exp_var <- function(r, n) sum(1 / (n - seq_len(r) + 1)^2)

#' @keywords internal
#' @noRd
.os_closed_exp_cov <- function(r, s, n) .os_closed_exp_var(r, n)

# ---- standard logistic: X_(r) = log U_(r) - log(1 - U_(r)), U_(r) beta -------
# The log-moments of a beta variable are digamma and trigamma values, so
# E[X_(r)] = psi(r) - psi(n-r+1) and Var = psi'(r) + psi'(n-r+1). No closed
# covariance is used.

#' @keywords internal
#' @noRd
.os_closed_logis_mean <- function(r, n) digamma(r) - digamma(n - r + 1)

#' @keywords internal
#' @noRd
.os_closed_logis_var <- function(r, n) trigamma(r) + trigamma(n - r + 1)

# ---- standard Gumbel (maximum form): F^n(x) = F(x - log n) ---------------------
# so the maximum is Gumbel with location log n: E[X_(n)] = gamma + log n and
# Var[X_(n)] = pi^2/6, with gamma = -digamma(1) the Euler constant. Other ranks
# have no closed form here (NA).

#' @keywords internal
#' @noRd
.os_closed_gumbel_mean <- function(r, n) if (r == n) -digamma(1) + log(n) else NA_real_

#' @keywords internal
#' @noRd
.os_closed_gumbel_var <- function(r, n) if (r == n) pi^2 / 6 else NA_real_
