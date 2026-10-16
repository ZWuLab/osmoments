# ==============================================================================
# File: os_linear.R
# Purpose: The linear layer of osmoments over the moments of os_moments(): the
#          ranks and weights of the nine sample quantile definitions of R's
#          quantile(), the mean and covariance of any linear combination of
#          order statistics, and Lloyd's best linear unbiased estimators of
#          location and scale from any rank subset that meets their conditions.
#
# EXPORTED FUNCTIONS:
#   - os_quantile_weights()  Ranks and weights of quantile() types 1 to 9
#   - os_lin_moments()       E[W' X] and Cov(W' X) for one or several weight vectors
#   - os_blue()              Lloyd's BLUE of location and scale: weights and variance factors
#
# INTERNAL HELPERS:
#   - .os_quantile_jh()      The (j, h) of each definition, as quantile.default computes them
#   - .os_pct_names()        The "25%" style names of quantile()
#   - .os_check_mom()        An osmoments object with its covariance matrix

#################### EXPORTED MAIN FUNCTIONS ####################

#' Ranks and weights of the sample quantile definitions
#'
#' A sample quantile under any of the nine definitions of
#' \code{\link[stats]{quantile}} is a fixed linear combination of two adjacent
#' order statistics, \eqn{(1 - h) X_{(j)} + h X_{(j+1)}}, with the ranks clipped
#' to \eqn{1, \dots, n}. This function returns, for each probability in
#' \code{p}, the ranks and the weights, so that the moments of the sample
#' quantile follow from the moments of the order statistics through
#' \code{\link{os_lin_moments}}. The even-sample median, for example, is
#' \eqn{(X_{(n/2)} + X_{(n/2+1)})/2} under type 7.
#'
#' @param p The probabilities, in \eqn{[0, 1]}.
#' @param n The sample size, an integer \eqn{\ge 1}.
#' @param type The definition, an integer from 1 to 9 as in
#'   \code{\link[stats]{quantile}}; 7 is R's default.
#'
#' @return A list with \code{ranks} (the sorted union of the ranks that carry a
#'   weight, for every probability), \code{W} (a matrix with one row per rank
#'   and one column per probability, each column summing to one, with the
#'   rows named by rank and the columns as \code{quantile()} names them),
#'   \code{n}, \code{p} and \code{type}.
#'
#' @details
#' The pair \eqn{(j, h)} is computed by the same arithmetic as
#' \code{stats:::quantile.default}. For type 7, \eqn{j = \lfloor 1 + (n-1)p
#' \rfloor} and \eqn{h = 1 + (n-1)p - j}. For types 4 to 9 except 7, with
#' \eqn{(a, b)} equal to \eqn{(0, 1)}, \eqn{(1/2, 1/2)}, \eqn{(0, 0)},
#' \eqn{(1/3, 1/3)} and \eqn{(3/8, 3/8)} for types 4, 5, 6, 8 and 9,
#' \eqn{m = a + p(n + 1 - a - b)}, \eqn{j = \lfloor m \rfloor} (with a fuzz of
#' four machine epsilons) and \eqn{h = m - j}. For types 1 to 3, \eqn{m = np}
#' (types 1 and 2) or \eqn{m = np - 1/2} (type 3), \eqn{j = \lfloor m \rfloor},
#' and \eqn{h} is 1 if \eqn{m > j} and 0 otherwise (type 1), 1 if \eqn{m > j}
#' and 1/2 otherwise (type 2), and 0 if \eqn{m = j} with \eqn{j} even and 1
#' otherwise (type 3). The ranks \eqn{j} and \eqn{j + 1} are then clipped to
#' \eqn{[1, n]}, and a rank whose weight is zero is dropped. On a sorted
#' sample the weighted sum reproduces \code{quantile()} exactly, which the
#' package's tests check for all nine types.
#'
#' @examples
#' qw <- os_quantile_weights(c(0.25, 0.5, 0.75), n = 10)
#' qw$ranks
#' qw$W
#' m <- os_moments(10, ranks = qw$ranks)
#' os_lin_moments(qw$W, m)          # moments of the three quartiles of a normal sample
#' @references
#' Hyndman, R. J. and Fan, Y. (1996). Sample quantiles in statistical packages.
#' \emph{The American Statistician}, 50(4), 361-365.
#' @seealso \code{\link{os_lin_moments}}, \code{\link[stats]{quantile}}
#' @export
os_quantile_weights <- function(p, n, type = 7) {
  fail <- function(...) stop("os_quantile_weights(): ", ..., call. = FALSE)
  if (!is.numeric(p) || length(p) < 1L || anyNA(p) || any(p < 0) || any(p > 1))
    fail("'p' must be probabilities in [0, 1]")
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n != round(n) || n < 1)
    fail("'n' must be one integer >= 1")
  if (n > .Machine$integer.max)
    fail("'n' must be at most ", .Machine$integer.max, ", the largest integer R stores")
  if (!is.numeric(type) || length(type) != 1L || !(type %in% 1:9))
    fail("'type' must be an integer from 1 to 9")
  n <- as.integer(n); type <- as.integer(type)
  jh <- .os_quantile_jh(p, n, type)        # j in doubles: j + 1 overflows the integer range at n = 2^31 - 1
  lo <- pmin(pmax(jh$j, 1), n)             # X_(j), clipped to the sample
  hi <- pmin(pmax(jh$j + 1, 1), n)         # X_(j+1), clipped
  w_lo <- 1 - jh$h; w_hi <- jh$h
  used <- as.integer(sort(unique(c(lo[w_lo > 0], hi[w_hi > 0]))))
  W <- matrix(0, length(used), length(p),
              dimnames = list(paste0("r", used), .os_pct_names(p, type)))
  for (i in seq_along(p)) {
    if (w_lo[i] > 0) W[match(lo[i], used), i] <- W[match(lo[i], used), i] + w_lo[i]
    if (w_hi[i] > 0) W[match(hi[i], used), i] <- W[match(hi[i], used), i] + w_hi[i]
  }
  list(ranks = used, W = W, n = n, p = p, type = type)
}

#' Moments of linear combinations of order statistics
#'
#' For a statistic \eqn{T = \sum_i w_i X_{(r_i)}} over the ranks of an
#' \code{osmoments} object, \eqn{E[T] = \bm w'\bm\alpha} and
#' \eqn{\mathrm{Var}(T) = \bm w'\bm B\bm w}, and for several statistics with
#' the weights in the columns of a matrix \eqn{\bm W}, the means are
#' \eqn{\bm W'\bm\alpha} and the covariance matrix is \eqn{\bm W'\bm B\bm W}.
#'
#' @param w A numeric vector with one weight per rank of \code{mom}, or a
#'   matrix with one such column per statistic.
#' @param mom An object of class \code{"osmoments"} computed with
#'   \code{cov = TRUE}, whose ranks are those the weights refer to, in order.
#'
#' @return A list with \code{mean} (one value per statistic), \code{var} (the
#'   variances), \code{cov} (the covariance matrix of the statistics, one by
#'   one for a single statistic) and \code{ranks}. The statistics are named by
#'   the column names of \code{w} when it has them.
#'
#' @examples
#' # the even-sample median of a normal sample of size 10
#' qw <- os_quantile_weights(0.5, n = 10)
#' m <- os_moments(10, ranks = qw$ranks)
#' os_lin_moments(qw$W, m)
#' # the range of an exponential sample of size 20
#' m <- os_moments(20, ranks = c(1, 20), dist = "exp")
#' os_lin_moments(c(-1, 1), m)
#' @seealso \code{\link{os_moments}}, \code{\link{os_quantile_weights}}
#' @export
os_lin_moments <- function(w, mom) {
  .os_check_mom(mom, "os_lin_moments()")
  if (!is.numeric(w) || !all(is.finite(w)))
    stop("os_lin_moments(): 'w' must be finite numbers", call. = FALSE)
  W <- if (is.matrix(w)) w else matrix(w, ncol = 1L)
  k <- length(mom$ranks)
  if (nrow(W) != k)
    stop("os_lin_moments(): 'w' must have one weight per rank of 'mom' (", k, " ranks: ",
         paste(mom$ranks, collapse = ", "), ")", call. = FALSE)
  # from the standard member's moments, so that a large location does not take the digits
  # of a location-free statistic (the range at location 1e16 came out as 0 instead of 1/3)
  mean <- mom$location * colSums(W) + mom$scale * drop(crossprod(W, mom$standard$alpha))
  cov <- mom$scale * (mom$scale * crossprod(W, mom$standard$B %*% W))
  nm <- colnames(W)
  if (!is.null(nm)) { names(mean) <- nm; dimnames(cov) <- list(nm, nm) } else names(mean) <- NULL
  list(mean = mean, var = diag(cov), cov = cov, ranks = mom$ranks)
}

#' Lloyd's best linear unbiased estimators of location and scale
#'
#' For a location-scale family with standard member \eqn{Z}, the order
#' statistics at the ranks of \code{mom} satisfy \eqn{X_{(r_i)} = \mu + \sigma
#' Z_{(r_i)}}, with \eqn{E[Z_{(r_i)}] = \alpha_i} and
#' \eqn{\mathrm{Cov}(Z_{(r_i)}, Z_{(r_j)}) = B_{ij}}. Lloyd's estimators are the
#' generalized least squares estimators of \eqn{(\mu, \sigma)} in this linear
#' model: \eqn{\hat\mu = \bm w_\mu'\bm X} and \eqn{\hat\sigma = \bm w_\sigma'\bm X}
#' with the weights of the Details. They are unbiased and have the smallest
#' variance among linear unbiased estimators from these ranks.
#'
#' @param mom An object of class \code{"osmoments"} computed with
#'   \code{cov = TRUE} at two or more ranks. The weights are computed from
#'   its \code{standard} component, the moments of the standard member, so
#'   they do not depend on the \code{location} and \code{scale} it was
#'   computed with.
#'
#' @return A list with \code{w_location} and \code{w_scale} (the weight vectors,
#'   named by rank), \code{var_location}, \code{var_scale} and \code{cov} (the
#'   variances and the covariance of the two estimators per unit of the squared
#'   scale parameter, so that \eqn{\mathrm{Var}(\hat\mu) = \sigma^2 \cdot}
#'   \code{var_location}), \code{ranks}, \code{n} and \code{dist}.
#'
#' @details
#' With \eqn{\bm\Omega = \bm B^{-1}}, \eqn{a_{11} = \mathbf 1'\bm\Omega\mathbf 1},
#' \eqn{a_{1\alpha} = \mathbf 1'\bm\Omega\bm\alpha},
#' \eqn{a_{\alpha\alpha} = \bm\alpha'\bm\Omega\bm\alpha} and
#' \eqn{\Delta = a_{11} a_{\alpha\alpha} - a_{1\alpha}^2},
#' \deqn{\bm w_\mu' = (a_{\alpha\alpha}\mathbf 1' - a_{1\alpha}\bm\alpha')\bm\Omega/\Delta,
#'   \qquad \bm w_\sigma' = (a_{11}\bm\alpha' - a_{1\alpha}\mathbf 1')\bm\Omega/\Delta,}
#' \deqn{\mathrm{Var}(\hat\mu) = \sigma^2 a_{\alpha\alpha}/\Delta, \qquad
#'   \mathrm{Var}(\hat\sigma) = \sigma^2 a_{11}/\Delta, \qquad
#'   \mathrm{Cov}(\hat\mu, \hat\sigma) = -\sigma^2 a_{1\alpha}/\Delta.}
#'
#' \strong{Computation.} The formulas above define the estimators. The function
#' computes them by whitening the design \eqn{[\mathbf 1, \bm\alpha]} with the
#' Cholesky factor of \eqn{\bm B} and solving the least-squares problem by QR,
#' which keeps its accuracy when \eqn{\mathbf 1} and \eqn{\bm\alpha} are nearly
#' collinear, as for two close ranks at a large \eqn{n}, where \eqn{\Delta}
#' itself loses its digits. It then checks the four unbiasedness conditions,
#' \eqn{\bm w_\mu'\mathbf 1 = 1}, \eqn{\bm w_\mu'\bm\alpha = 0},
#' \eqn{\bm w_\sigma'\mathbf 1 = 0} and \eqn{\bm w_\sigma'\bm\alpha = 1}, and
#' warns when one fails by more than \eqn{10^{-6}}.
#'
#' \strong{Conditions.} The distribution of \code{mom} is treated as the
#' standard member of a location-scale family, which every \code{os_dist()}
#' object is by construction; the required moments exist, which
#' \code{os_moments()} has checked; at least two distinct ranks are present, so
#' that the design \eqn{[\mathbf 1, \bm\alpha]} has rank two; and \eqn{\bm B} is
#' positive definite. The scale estimated is the family's scale parameter: the
#' standard deviation for the normal, the scale of the logistic or Gumbel, and
#' the length of the support for the uniform.
#'
#' @examples
#' m <- os_moments(9, ranks = c(1, 3, 5, 7, 9))
#' b <- os_blue(m)
#' b$w_location                     # sums to one
#' sqrt(b$var_location)             # SD of the location estimate per unit sigma
#' @references
#' Lloyd, E. H. (1952). Least-squares estimation of location and scale
#' parameters using order statistics. \emph{Biometrika}, 39(1-2), 88-95.
#' @seealso \code{\link{os_moments}}, \code{\link{os_lin_moments}}
#' @export
os_blue <- function(mom) {
  .os_check_mom(mom, "os_blue()")
  k <- length(mom$ranks)
  if (k < 2L)
    stop("os_blue(): at least two distinct ranks are needed to estimate location and scale",
         call. = FALSE)
  alpha0 <- mom$standard$alpha                          # the standard member's moments,
  B0 <- mom$standard$B                                  # kept exact whatever the location
  # Whiten by the Cholesky factor, B0 = R'R, so that the design R'^{-1} X has
  # uncorrelated unit-variance errors, and solve that least-squares problem by QR.
  # This avoids forming Delta = a11 aaa - a1a^2, which loses its digits when the
  # columns 1 and alpha are nearly collinear (two close ranks at a large n).
  Rc <- tryCatch(chol(B0), error = function(e)
    stop("os_blue(): the covariance matrix is not positive definite, so the estimators ",
         "are not defined: ", conditionMessage(e), call. = FALSE))
  # The alpha column is centered and scaled first: X = [1, z], z = (alpha0 - abar)/sa,
  # whose coefficients are g1 = mu + sigma abar and g2 = sigma sa.
  abar <- mean(alpha0); sa <- sqrt(mean((alpha0 - abar)^2))
  if (!is.finite(sa) || sa == 0)
    stop("os_blue(): the means at the ranks are all equal, so location and scale cannot ",
         "be separated", call. = FALSE)
  Xt <- backsolve(Rc, cbind(1, (alpha0 - abar) / sa), transpose = TRUE)   # R'^{-1} X
  qx <- qr(Xt)
  if (qx$rank < 2L)
    stop("os_blue(): the design [1, alpha] is not of full rank", call. = FALSE)
  Rx_inv <- backsolve(qr.R(qx), diag(2))                # (X'B0^{-1}X)^{-1} = Rx_inv Rx_inv'
  G <- backsolve(Rc, qr.Q(qx) %*% t(Rx_inv))            # k x 2: the weights of (g1, g2)
  Tm <- rbind(c(1, -abar / sa), c(0, 1 / sa))           # (mu, sigma) = Tm (g1, g2)
  Wt <- G %*% t(Tm)                                     # k x 2: the weights of (mu, sigma)
  M <- Tm %*% Rx_inv
  V <- tcrossprod(M)                                    # Cov(mu_hat, sigma_hat) / sigma^2
  w_mu <- Wt[, 1]; w_sig <- Wt[, 2]
  # The four unbiasedness conditions; rounding of alpha limits them when the weights are large
  resid <- max(abs(c(sum(w_mu) - 1, sum(w_mu * alpha0), sum(w_sig), sum(w_sig * alpha0) - 1)))
  if (!is.finite(resid) || resid > 1e-6)
    warning("os_blue(): the weights meet the unbiasedness conditions only to ",
            signif(resid, 3), "; the means at these ranks are too close together for the ",
            "estimators to be determined in double precision", call. = FALSE)
  names(w_mu) <- names(w_sig) <- names(mom$alpha)
  list(w_location = w_mu, w_scale = w_sig,
       var_location = V[1, 1], var_scale = V[2, 2], cov = V[1, 2],
       ranks = mom$ranks, n = mom$n, dist = mom$dist$name)
}

#################### INTERNAL HELPER FUNCTIONS ####################

#' The (j, h) of each quantile definition, as stats:::quantile.default computes them
#' @keywords internal
#' @noRd
.os_quantile_jh <- function(p, n, type) {
  if (type == 7L) {
    index <- 1 + max(n - 1, 0) * p
    j <- floor(index)
    h <- index - j
  } else if (type <= 3L) {
    nppm <- if (type == 3L) n * p - 0.5 else n * p
    j <- floor(nppm)
    h <- switch(type,
                as.numeric(nppm > j),                               # type 1
                (as.numeric(nppm > j) + 1) / 2,                     # type 2
                as.numeric((nppm != j) | ((j %% 2) == 1)))          # type 3
  } else {
    ab <- switch(type - 3L, c(0, 1), c(0.5, 0.5), c(0, 0), NULL, c(1 / 3, 1 / 3), c(3 / 8, 3 / 8))
    fuzz <- 4 * .Machine$double.eps
    nppm <- ab[1] + p * (n + 1 - ab[1] - ab[2])
    j <- floor(nppm + fuzz)
    h <- nppm - j
    h[abs(h) < fuzz] <- 0
  }
  list(j = j, h = h)
}

#' The "25%" style names of quantile(), taken from quantile() itself
#'
#' The names depend on the probabilities only (and on how many there are, since
#' quantile() changes its formatter at 100 probabilities), so a two-point sample
#' gives them, in the form of the R version in use.
#' @keywords internal
#' @noRd
.os_pct_names <- function(p, type) names(stats::quantile(c(0, 1), probs = p, type = type))

#' An osmoments object with its covariance matrix, or an error
#' @keywords internal
#' @noRd
.os_check_mom <- function(mom, who) {
  if (!inherits(mom, "osmoments"))
    stop(who, ": 'mom' must be an object returned by os_moments()", call. = FALSE)
  if (is.null(mom$B))
    stop(who, ": 'mom' was computed with cov = FALSE; the covariance matrix is needed",
         call. = FALSE)
  invisible(TRUE)
}
