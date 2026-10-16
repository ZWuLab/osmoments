# ==============================================================================
# File: os_moments.R
# Purpose: The main function of osmoments: os_moments() returns the exact mean
#          vector and covariance matrix of the order statistics of a continuous
#          distribution at the requested ranks, from the closed forms where the
#          family has them and from the integrators of os_engine.R otherwise,
#          with the location-scale transform, the symmetry halving and the
#          run-time checks, and keeps the moments of the standard member in the
#          object for the linear layer. Also the print method of the returned object.
#
# EXPORTED FUNCTIONS:
#   - os_moments()        Exact alpha, variances and B at the requested ranks
#   - print.osmoments()   Print method
#
# INTERNAL HELPERS:
#   - .os_validate_args() Argument checks of os_moments()
#   - .os_mass_check()    The density of X_(r) integrates to 1, or the quadrature missed a mode
#   - .os_is_normal()     Whether the distribution's functions are R's normal functions
#   - .os_closed_value()  A closed-form hook's value, or NA when absent or not wanted

#################### EXPORTED MAIN FUNCTIONS ####################

#' Exact moments of order statistics
#'
#' Computes the means \eqn{\alpha_i = E[X_{(r_i)}]}, the variances and the
#' covariance matrix \eqn{B_{ij} = \mathrm{Cov}(X_{(r_i)}, X_{(r_j)})} of the
#' order statistics \eqn{X_{(1)} \le \dots \le X_{(n)}} of an i.i.d. sample of
#' size \code{n} from a continuous distribution, at the ranks \code{ranks}, by
#' numerical evaluation of the exact integrals, or from a closed form where the
#' family has one.
#'
#' @param n The sample size, an integer \eqn{\ge 1}.
#' @param ranks The ranks, strictly increasing integers in \eqn{1, \dots, n}.
#'   The default is every rank, which gives the full covariance matrix.
#' @param dist The distribution: one of \code{"norm"}, \code{"unif"},
#'   \code{"exp"}, \code{"logis"}, \code{"gumbel"}, or an object made by
#'   \code{\link{os_dist}}. See Details for the standard members.
#' @param location,scale The location and scale applied to the standard member:
#'   \eqn{X = \mathrm{location} + \mathrm{scale} \cdot Z}, with \code{scale > 0}.
#'   The means become \eqn{\mathrm{location} + \mathrm{scale}\,\alpha} and the
#'   covariances \eqn{\mathrm{scale}^2 B}.
#' @param cov \code{TRUE} for the full covariance matrix, \code{FALSE} for the
#'   means and the variances only, which need single integrals and no double
#'   integral.
#' @param method \code{"auto"} takes the closed form where the family has one
#'   and the quadrature elsewhere, \code{"integrate"} forces the quadrature,
#'   \code{"closed"} requires a closed form for every requested quantity and
#'   stops otherwise.
#' @param rel.tol The relative tolerance of the quadrature. The single integrals
#'   and the inner integral of a product moment run at \code{rel.tol}, the
#'   outer integral at \code{100 * rel.tol}.
#' @param check Whether to run the run-time checks of the Details and warn when
#'   a residual exceeds \eqn{10^{-6}}.
#'
#' @return An object of class \code{"osmoments"}, a list with the components
#'   \code{n}, \code{ranks}, \code{alpha} (the means, named by rank),
#'   \code{var} (the variances), \code{B} (the covariance matrix, \code{NULL}
#'   when \code{cov = FALSE}), \code{quad_err} (the summed error estimates of
#'   the quadrature, see Details), \code{standard} (a list with the
#'   \code{alpha}, \code{var} and \code{B} of the standard member, before
#'   \code{location} and \code{scale} are applied, which \code{\link{os_blue}}
#'   reads), \code{dist} (the \code{os_dist} object),
#'   \code{location}, \code{scale}, \code{cov}, \code{method}, \code{route}
#'   (\code{"closed"}, \code{"integrate"} or \code{"mixed"}), \code{checks}
#'   (the residuals of the run-time checks, \code{NA} where a check did not
#'   apply) and \code{call}.
#'
#' @details
#' \strong{The integrals.} With cdf \eqn{F}, density \eqn{f} and
#' \eqn{y < z},
#' \deqn{f_{r:n}(x) = \frac{n!}{(r-1)!\,(n-r)!} F(x)^{r-1} \{1 - F(x)\}^{n-r} f(x),}
#' \deqn{f_{r,s:n}(y, z) = \frac{n!}{(r-1)!\,(s-r-1)!\,(n-s)!} F(y)^{r-1}
#'   \{F(z) - F(y)\}^{s-r-1} \{1 - F(z)\}^{n-s} f(y) f(z),}
#' and the moments are
#' \deqn{\mu_{r:n} = \int x f_{r:n}(x)\,dx, \qquad
#'   \mu^{(2)}_{r:n} = \int x^2 f_{r:n}(x)\,dx, \qquad
#'   \mu_{r,s:n} = \int\!\!\int_{y < z} y z f_{r,s:n}(y, z)\,dy\,dz,}
#' with \eqn{B_{ij} = \mu_{r_i, r_j:n} - \alpha_i \alpha_j}. The integrals are
#' evaluated by adaptive quadrature (\code{\link[stats]{integrate}}) with four
#' devices: the integrands on the log scale, so that the powers of the tails
#' never underflow; \eqn{\log\{F(z) - F(y)\}} formed on the smaller tail, so
#' that nothing cancels; the range split at the spike of the order statistic's
#' density, at \eqn{c_r \pm h_r} with \eqn{c_r = Q(p_r)}, \eqn{p_r = r/(n+1)}
#' and \eqn{h_r = 10 s_r / f(c_r) + 4\tilde\sigma/\sqrt{n+1}},
#' \eqn{s_r = \sqrt{p_r(1-p_r)/(n+2)}}, so that the adaptive rule sees the
#' narrow peak at large \eqn{n} at the end of a piece rather than between its
#' nodes, the outer range of a product moment being split at the points of
#' both ranks; and integration in standardized units of the standard member,
#' \eqn{t = \{x - Q(1/2)\}/\tilde\sigma}, where \eqn{\tilde\sigma} is the
#' spread of the distribution, its interquartile range divided by that of the
#' standard normal (so \eqn{\tilde\sigma = 1} for the normal), with each
#' moment taken about its own order statistic's location \eqn{c_r}. The
#' scaling puts the split floor and the integrator's absolute tolerances in
#' the distribution's own units, so that a distribution supplied on any scale
#' is computed as accurately as on the scale one. The centering at \eqn{c_r}
#' keeps a variance from being the small difference of two moments of size
#' \eqn{(c_r - c)^2}: taken about the median \eqn{c}, the variance of the
#' uniform minimum kept seven digits at \eqn{n = 1000} and none at
#' \eqn{n = 10^5}.
#'
#' \strong{The standard members and the closed forms.} \code{"norm"} is
#' \eqn{N(0, 1)}; \code{"unif"} is the uniform on \eqn{[0, 1]}, so that
#' \code{location} is the lower end and \code{scale} the length;
#' \code{"exp"} has unit rate, so that \code{scale} is \eqn{1/\mathrm{rate}};
#' \code{"logis"} is the standard logistic; \code{"gumbel"} is the standard
#' Gumbel in the maximum form, \eqn{F(x) = \exp(-e^{-x})}. The uniform and the
#' exponential have closed forms for every moment, the logistic for the means
#' and variances, the Gumbel for the mean and variance of the maximum. With
#' \code{method = "auto"} these are returned and the rest is quadrature.
#'
#' \strong{The moment condition.} A moment is computed only if it exists:
#' with the tail indices \eqn{(\nu_L, \nu_U)} of the distribution,
#' \eqn{E|X_{(r)}|^2 < \infty} if and only if \eqn{\nu_L r > 2} and
#' \eqn{\nu_U (n - r + 1) > 2}, and a rank that fails stops the computation
#' with an error. The built-in families always pass.
#'
#' \strong{The error estimate.} \code{quad_err} is the sum of the error
#' estimates that \code{integrate()} reports for the single integrals and for
#' the outer integrals of the product moments, in the units of the standard
#' member before \code{location} and \code{scale} are applied. It is an
#' estimate, not a bound: the inner integral's estimate is not carried into the
#' outer one, and the errors of the means are not propagated into the
#' covariances. The accuracy evidence is the agreement with the closed forms,
#' the identities and the printed tables, documented in the package's tests.
#'
#' \strong{Accuracy at large n.} With \code{method = "integrate"}, the means
#' and the variances of the uniform and the exponential at the minimum, the
#' quartiles, the median and the maximum agree with their closed forms to
#' within \eqn{10^{-9}} in relative terms at \eqn{n = 10^4}, \eqn{10^5} and
#' \eqn{10^6}, and so does the correlation of every pair of them, in absolute
#' terms. Up to the relative errors of the two variances, the relative error
#' of a covariance is the error of its correlation divided by the
#' correlation, so the covariance of two nearly independent
#' order statistics keeps fewer digits: the two extremes, whose correlation is
#' of order \eqn{1/n}, keep five to ten significant digits of their covariance
#' at those sizes, depending on the parent and on \eqn{n}. The package's tests
#' assert the variances and the covariances of adjacent ranks at \eqn{10^{-9}}
#' and the extremes' covariance at \eqn{10^{-3}} beyond \eqn{n = 10^4}. A
#' parent as heavy-tailed as the lognormal with log-SD 2 is accurate to about
#' \eqn{10^{-9}} in relative terms rather than \eqn{10^{-11}}, because the far
#' tail of the inner integral of a covariance is resolved at a fixed tolerance.
#' These are the cases that were measured, not a bound for every parent.
#'
#' \strong{What stops with an error.} The quadrature follows the peak of each
#' order statistic's density and the tails of the parent. It does not follow a
#' second mode: for a mixture of two narrow normals, the mode that holds a
#' tenth of the mass fell between the nodes, and the moments came back as those
#' of the other mode. The package therefore integrates the density of every
#' order statistic it computes by quadrature and stops when the mass is not 1
#' within \eqn{10^{-7}}, naming the rank and the mass it found. It also stops
#' when a piece of the range cannot be resolved, which happens for a moment
#' that does not exist when the tail indices were not declared, for a tail as
#' heavy as a lognormal with log-SD 3 or a Weibull with shape 0.3, for the
#' extremes of a parent whose density is infinite at an endpoint of its support
#' once they concentrate against that endpoint (the minimum of a beta(1/2, 1/2)
#' is computed to \eqn{n = 100} and refused from \eqn{n = 1000} on), and for a
#' variance that comes out zero or negative. A wrong number is never returned
#' knowingly, but these checks detect some failures and certify none.
#'
#' \strong{The run-time checks.} For a symmetric distribution and a rank set
#' closed under \eqn{r \mapsto n + 1 - r}, the means are antisymmetric and the
#' second moments symmetric about the median, which the integrators do not
#' impose; the product moments are then computed for the pairs with
#' \eqn{r + s \le n + 1} only and reflected. For a complete rank set, the order
#' statistics are a permutation of the sample, so \eqn{\sum_r \alpha_r = n E[X]},
#' \eqn{\sum_r \mu^{(2)}_{r:n} = n E[X^2]} and \eqn{\sum_{r,s} B_{rs} = n\,\mathrm{Var}(X)}
#' for any parent with a finite second moment, and for the normal parent
#' (recognized by its functions \code{pnorm}, \code{dnorm} and \code{qnorm}, not
#' by its name) every row of \eqn{B} sums to \eqn{\sigma^2} as well. The
#' residuals are stored in \code{checks}, in the standardized units of the
#' integration (the standard member centered at its median and scaled by its
#' spread, so that the threshold means the same on every scale), and a
#' residual above \eqn{10^{-6}}
#' raises a warning, as does a covariance matrix with a negative eigenvalue
#' beyond rounding. The checks detect some failures of the quadrature and do not
#' certify its accuracy: a failure that every integral shares, such as a peak
#' missed in the same way at every rank, can satisfy them.
#'
#' \strong{Cost.} A single integral takes a fraction of a millisecond and a
#' double integral a few hundredths of a second. When the covariances are not
#' closed forms, a full matrix needs one double integral per pair,
#' \eqn{n(n-1)/2} of them, halved for a symmetric family.
#'
#' @examples
#' m <- os_moments(5)            # all ranks of n = 5, standard normal
#' m$alpha
#' m$B
#' os_moments(20, ranks = c(1, 5, 10, 20), dist = "exp", scale = 2)
#' os_moments(100, ranks = c(1, 100), dist = "gumbel", cov = FALSE)
#' @references
#' Wu, Z. and Tong, T. (2026). osmoments: Exact means, variances and covariances
#' of order statistics in R. Manuscript. The computation, its verification and
#' the comparison with other software.
#'
#' Royston, J. P. (1982). Algorithm AS 177: Expected normal order statistics
#' (exact and approximate). \emph{Applied Statistics}, 31(2), 161-165.
#' @seealso \code{\link{os_dist}}
#' @export
os_moments <- function(n, ranks = seq_len(n), dist = "norm", location = 0, scale = 1,
                       cov = TRUE, method = c("auto", "integrate", "closed"),
                       rel.tol = 1e-11, check = TRUE) {
  method <- match.arg(method)
  .os_validate_args(n, ranks, location, scale, cov, rel.tol, check)
  n <- as.integer(n); ranks <- as.integer(ranks); k <- length(ranks)
  dist <- .os_resolve_dist(dist)
  ctrl <- .os_ctrl(rel.tol)
  .os_moment_condition(ranks, n, dist$tail, k = 2)
  cd <- .os_centered(dist)
  .os_tail_probe(ranks, n, cd, ctrl)
  c0 <- dist$center
  sp <- cd$spread     # the integrals are in units of the spread; sp and sp^2 bring them back
  cr <- cd$q(ranks / (n + 1))   # each moment is taken about its own location, in those units
  m1 <- numeric(k)              # E[X_(r)] - c_r in those units, the first centered moment
  cl <- dist$closed
  err <- 0
  used_closed <- logical(0)                   # one entry per quantity: closed form or not
  need_quad <- function(what, r, s = NULL) {  # stop under method = "closed"
    if (method == "closed")
      stop("os_moments(): no closed form for the ", what, " of rank ",
           if (is.null(s)) r else paste(r, "and", s), " of the \"", dist$name,
           "\" family; use method = \"auto\" or \"integrate\"", call. = FALSE)
  }

  # ---- means and variances, in centered units of the standard member ----
  a_c <- numeric(k); v <- numeric(k)
  for (i in seq_len(k)) {
    r <- ranks[i]
    m <- .os_closed_value(cl$mean, method, r, n)
    if (is.na(m)) {
      need_quad("mean", r)
      .os_mass_check(r, n, cd, ctrl, cr[i])
      o <- .os_single(r, n, 1, cd, ctrl, cr[i]); a_c[i] <- sp * (cr[i] + o$value)
      err <- err + sp * o$abs.error
      used_closed <- c(used_closed, FALSE)
    } else {
      a_c[i] <- m - c0; used_closed <- c(used_closed, TRUE)
    }
    m1[i] <- a_c[i] / sp - cr[i]
    vv <- .os_closed_value(cl$var, method, r, n)
    if (is.na(vv)) {
      need_quad("variance", r)
      o <- .os_single(r, n, 2, cd, ctrl, cr[i]); v[i] <- sp * (sp * (o$value - m1[i]^2))
      err <- err + sp * (sp * o$abs.error)
      used_closed <- c(used_closed, FALSE)
    } else {
      v[i] <- vv; used_closed <- c(used_closed, TRUE)
    }
  }

  # ---- a variance is positive; one that is not was lost by the quadrature ----
  bad <- which(!is.finite(v) | v <= 0)
  if (length(bad))
    stop("os_moments(): the variance of X_(", ranks[bad[1]], ") came out as ",
         signif(v[bad[1]], 3), ", so the quadrature could not resolve it. This happens when ",
         "the order statistic is concentrated on a scale far below the spread of the ",
         "distribution: an extreme rank of a bounded or semi-bounded distribution at a ",
         "large n on the integration route (method = \"auto\" uses the closed form where one ",
         "exists), or a distribution supplied on a scale far from one (supply the standard ",
         "member and pass the scale to os_moments())", call. = FALSE)

  # ---- covariances: the pairs with r + s <= n + 1 when the parent is symmetric ----
  # and the rank set is closed under reflection, the mirror pairs filled afterwards;
  # the rank sums are formed in doubles, so that they cannot overflow the integer range
  reflect <- dist$symmetric && all(as.numeric(ranks) + rev(ranks) == n + 1)
  B_c <- NULL
  if (cov && k >= 2L) {
    B_c <- matrix(NA_real_, k, k); diag(B_c) <- v
    for (i in seq_len(k - 1L)) for (j in seq.int(i + 1L, k)) {
      r <- ranks[i]; s <- ranks[j]
      if (reflect && as.numeric(r) + s > n + 1) next
      cv <- .os_closed_value(cl$cov, method, r, s, n)
      if (is.na(cv)) {
        need_quad("covariance", r, s)
        o <- .os_prod(r, s, n, cd, ctrl, cr[i], cr[j]); cv <- sp * (sp * (o$value - m1[i] * m1[j]))
        err <- err + sp * (sp * o$abs.error)
        used_closed <- c(used_closed, FALSE)
      } else {
        used_closed <- c(used_closed, TRUE)
      }
      B_c[i, j] <- B_c[j, i] <- cv
    }
    if (reflect) for (i in seq_len(k - 1L)) for (j in seq.int(i + 1L, k))
      if (as.numeric(ranks[i]) + ranks[j] > n + 1)  # Cov(X_(r), X_(s)) = Cov(X_(n+1-s), X_(n+1-r))
        B_c[i, j] <- B_c[j, i] <- B_c[k + 1L - j, k + 1L - i]
  } else if (cov) {
    B_c <- matrix(v, 1L, 1L)
  }
  route <- if (method == "integrate" || !any(used_closed)) "integrate"
           else if (all(used_closed)) "closed" else "mixed"

  # ---- run-time checks, in the standardized units of the engine, so that the warning
  # threshold means the same on every scale (a normal supplied with standard deviation
  # 1e6 has residuals of 1e-3 in its own units, 1e-9 in standardized ones) ----
  checks <- list(reflection = NA_real_, identities = NA)
  if (check) {
    a_s <- a_c / sp; v_s <- v / (sp * sp)
    B_s <- if (is.null(B_c)) NULL else B_c / (sp * sp)
    if (reflect) checks$reflection <- .os_reflection_check(a_s, v_s + a_s^2)
    if (k == n) {
      pm <- .os_parent_moments(dist, cd, ctrl)
      pm_s <- list(m1 = pm$m1 / sp, m2 = pm$m2 / (sp * sp))
      checks$identities <- .os_identity_check(a_s, v_s, B_s, n, pm_s, .os_is_normal(dist))
    }
    resid <- c(checks$reflection, checks$identities)
    resid <- resid[!is.na(resid)]
    if (length(resid) && max(resid) > 1e-6)
      warning("os_moments(): a run-time check has residual ", signif(max(resid), 3),
              " at n = ", n, " for \"", dist$name, "\"; the quadrature is suspect",
              call. = FALSE)
    # a covariance matrix is positive semidefinite; rounding may leave a negative
    # eigenvalue of relative size near the double epsilon, a failed quadrature more
    if (!is.null(B_c) && k >= 2L) {
      ev <- eigen(B_c, symmetric = TRUE, only.values = TRUE)$values
      if (min(ev) < -sqrt(.Machine$double.eps) * max(ev))
        warning("os_moments(): the covariance matrix has the negative eigenvalue ",
                signif(min(ev), 3), " (largest ", signif(max(ev), 3), ") at n = ", n,
                " for \"", dist$name, "\"; the quadrature is suspect", call. = FALSE)
    }
  }

  # ---- add the center back, apply location and scale, name by rank ----
  # The scale multiplies the variances and covariances one factor at a time, so that
  # scale^2 cannot overflow or underflow where the products themselves are representable.
  nm <- paste0("r", ranks)
  alpha0 <- a_c + c0; names(alpha0) <- nm             # the standard member's moments
  var0 <- v; names(var0) <- nm
  B0 <- if (is.null(B_c)) NULL else { Bm <- B_c; dimnames(Bm) <- list(nm, nm); Bm }
  alpha <- location + scale * alpha0
  var <- scale * (scale * var0)
  B <- if (is.null(B0)) NULL else scale * (scale * B0)
  if (!all(is.finite(alpha)) || !all(is.finite(var)) || (!is.null(B) && !all(is.finite(B))))
    stop("os_moments(): the moments overflow double precision at location ", format(location),
         " and scale ", format(scale), call. = FALSE)
  structure(list(n = n, ranks = ranks, alpha = alpha, var = var, B = B, quad_err = err,
                 standard = list(alpha = alpha0, var = var0, B = B0),
                 dist = dist, location = location, scale = scale, cov = cov,
                 method = method, route = route, checks = checks, call = match.call()),
            class = "osmoments")
}

#' Print an osmoments object
#'
#' @param x An object of class \code{"osmoments"}.
#' @param digits The number of significant digits.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @export
print.osmoments <- function(x, digits = 7L, ...) {
  cat("Exact moments of order statistics: ", x$dist$name,
      " (location ", format(x$location), ", scale ", format(x$scale), "), n = ", x$n,
      ", ranks ", paste(x$ranks, collapse = " "), "\n", sep = "")
  cat("alpha (means):\n"); print(signif(x$alpha, digits))
  if (is.null(x$B)) {
    cat("var (variances; cov = FALSE):\n"); print(signif(x$var, digits))
  } else {
    cat("B (covariance matrix):\n"); print(signif(x$B, digits))
  }
  cat("route: ", x$route, "; quadrature error estimate: ", format(signif(x$quad_err, 3)),
      "\n", sep = "")
  ck <- x$checks
  cat("checks: reflection ",
      if (is.na(ck$reflection)) "n/a" else format(signif(ck$reflection, 3)),
      "; identities ",
      if (length(ck$identities) == 1L && is.na(ck$identities)) "n/a"
      else paste(names(ck$identities), format(signif(ck$identities, 3)), collapse = ", "),
      "\n", sep = "")
  invisible(x)
}

#################### INTERNAL HELPER FUNCTIONS ####################

#' Argument checks of os_moments()
#' @keywords internal
#' @noRd
.os_validate_args <- function(n, ranks, location, scale, cov, rel.tol, check) {
  fail <- function(...) stop("os_moments(): ", ..., call. = FALSE)
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n != round(n) || n < 1)
    fail("'n' must be one integer >= 1")
  if (n > .Machine$integer.max)
    fail("'n' must be at most ", .Machine$integer.max, ", the largest integer R stores")
  if (!is.numeric(ranks) || length(ranks) < 1L || anyNA(ranks) || any(ranks != round(ranks)))
    fail("'ranks' must be integers")
  if (any(ranks < 1) || any(ranks > n))
    fail("'ranks' must lie in 1, ..., n = ", n)
  if (length(ranks) > 1L && any(diff(ranks) <= 0))
    fail("'ranks' must be strictly increasing")
  if (!is.numeric(location) || length(location) != 1L || !is.finite(location))
    fail("'location' must be one finite number")
  if (!is.numeric(scale) || length(scale) != 1L || !is.finite(scale) || scale <= 0)
    fail("'scale' must be one positive number")
  if (!is.logical(cov) || length(cov) != 1L || is.na(cov))
    fail("'cov' must be TRUE or FALSE")
  if (!is.numeric(rel.tol) || length(rel.tol) != 1L || !is.finite(rel.tol) || rel.tol <= 0)
    fail("'rel.tol' must be one positive number")
  if (!is.logical(check) || length(check) != 1L || is.na(check))
    fail("'check' must be TRUE or FALSE")
  invisible(TRUE)
}

#' The mass of the density of X_(r), against 1
#'
#' Device (3) splits the range at the peak of each order statistic and the
#' parent's tails, and a second mode of the parent can fall between the nodes of
#' a piece: for the mixture 0.9 N(0, 0.01) + 0.1 N(10, 0.01) the mass at 10 was
#' missed entirely at n = 1, and the mean came back as 0.0014 instead of 1 with
#' every identity satisfied, because the parent moments are computed by the same
#' quadrature. The zeroth moment, which costs one single integral per rank, is
#' the one quantity with a known value, so it is checked before the moments are
#' trusted.
#' @keywords internal
#' @noRd
.os_mass_check <- function(r, n, cd, ctrl, shift) {
  mass <- .os_single(r, n, 0, cd, ctrl, shift)$value
  if (!is.finite(mass) || abs(mass - 1) > 1e-7)
    stop("os_moments(): the density of X_(", r, ") integrates to ", signif(mass, 6),
         " instead of 1, so the quadrature missed part of its mass. The integration range ",
         "is split around the peak of each order statistic's density, which does not follow ",
         "a parent with more than one mode, or with mass far from its median on the scale ",
         "of its spread", call. = FALSE)
  invisible(mass)
}

#' Whether a distribution object is the normal, by its functions
#'
#' The unit row sums of B hold for the normal parent only, so the identity check
#' applies them when the three functions are R's own normal functions, and never
#' on the strength of the label, which is the user's to choose.
#' @keywords internal
#' @noRd
.os_is_normal <- function(dist) {
  identical(dist$p, stats::pnorm) && identical(dist$d, stats::dnorm) &&
    identical(dist$q, stats::qnorm)
}

#' A closed-form hook's value, or NA when the hook is absent, not wanted, or
#' returns NA for this rank
#' @keywords internal
#' @noRd
.os_closed_value <- function(fn, method, ...) {
  if (method == "integrate" || is.null(fn)) return(NA_real_)
  v <- fn(...)
  if (is.null(v) || length(v) != 1L || is.na(v)) NA_real_ else as.numeric(v)
}
