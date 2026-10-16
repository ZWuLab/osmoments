# ==============================================================================
# File: os_engine.R
# Purpose: The internal integrators of osmoments: the exact single and product
#          moments of order statistics by adaptive quadrature of the defining
#          integrals, with four numerical devices: (1) log-scale integrands
#          through log.p and lower.tail, (2) the cancellation-free
#          log{F(z) - F(y)} evaluated on the smaller tail, (3) the range split at
#          the spike of the order statistic's density, (4) integration in the
#          standardized units t = {x - Q(1/2)}/spread, each moment taken about
#          its own order statistic's location c_r = Q(r/(n+1)), so that
#          E[X^2] - E[X]^2 does not cancel at a large location or at a large n.
#          Everything here works in those units; os_moments() brings the
#          moments back and applies the location and scale.
#
# INTERNAL HELPERS:
#   - .os_ctrl()            The quadrature control constants from rel.tol
#   - .os_centered()        The log-tails, log-density and quantile function in centered units
#   - .os_breaks()          Device (3): the split points of rank r, c_r +/- h_r
#   - .os_width()           c_r and the delta-method half-width w_r of rank r
#   - .os_tolscale()        min(1, w): the scale of the absolute tolerances of rank r
#   - .os_integrate_split() integrate() over the pieces, values and error estimates summed
#   - .os_logdiff()         Device (2): log{F(z) - F(y)} on the smaller tail
#   - .os_xlog()            k * log(v) with the k = 0 guard against 0 * (-Inf)
#   - .os_single()          E[X_(r)^pow], pow = 1 or 2, one integral
#   - .os_prod()            E[X_(r) X_(s)], r < s, nested quadrature

#################### INTERNAL HELPER FUNCTIONS ####################

#' The quadrature control constants
#'
#' The constants of the reference engine, an earlier implementation of the same
#' method for the normal, written separately before the package: the single
#' integrals and the inner integral of
#' a product moment at rel.tol (1e-11 by default), the outer integral at 100
#' times rel.tol, so that the inner quadrature noise sits two orders below what
#' the outer adaptive rule chases; the absolute tolerances 1e-13 (single), 1e-15
#' (inner) and 1e-12 (outer); 400 subdivisions; a non-"OK" piece accepted when
#' its error estimate is below 1e-10; the half-width k = 10 of the split in
#' units of the large-sample standard deviation of the uniform order statistic.
#' @keywords internal
#' @noRd
.os_ctrl <- function(rel.tol = 1e-11) {
  list(rel.tol = rel.tol, rel.tol.outer = 100 * rel.tol,
       abs.tol.single = 1e-13, abs.tol.inner = 1e-15, abs.tol.outer = 1e-12,
       subdivisions = 400L, accept.err = 1e-10, k = 10)
}

#' The distribution in standardized units t = {x - Q(1/2)} / spread
#'
#' Device (4). The spread is the standard member's interquartile range divided
#' by that of the standard normal, 2 qnorm(3/4) = 1.349, so that it is exactly 1
#' for the normal, for which t = x. Returns closures lF(t) = log F, lS(t) =
#' log(1 - F), ld(t) = log of the density of t (the density of x times the
#' spread), pf(t) = F and q(u) = {Q(u) - c} / spread, the support endpoints in
#' the same units, and the center c = Q(1/2) and the spread, which os_moments()
#' uses to bring the results back to the standard member's units. In these
#' units the median is 0, the quartiles are 1.349 apart (at -0.674 and 0.674 for
#' a symmetric parent, as for the normal), and the split floor and the integrator's absolute tolerances, which
#' are fixed numbers, mean the same thing for every distribution, whatever the
#' scale on which it is supplied. (With the floor and the tolerances in absolute
#' units, a normal supplied with standard deviation 0.01 came back with variance
#' exactly zero, its whole mass between the nodes of a central piece 290 of its
#' standard deviations wide, and at 0.001 the second moments, of order 1e-6,
#' were integrated to seven digits only.)
#' @keywords internal
#' @noRd
.os_centered <- function(dist) {
  c0 <- dist$center
  p <- dist$p; d <- dist$d; q <- dist$q
  sp <- (q(0.75) - q(0.25)) / (2 * qnorm(0.75))
  list(lF = function(t) p(c0 + sp * t, lower.tail = TRUE, log.p = TRUE),
       lS = function(t) p(c0 + sp * t, lower.tail = FALSE, log.p = TRUE),
       ld = function(t) d(c0 + sp * t, log = TRUE) + log(sp),
       pf = function(t) p(c0 + sp * t),
       q  = function(u) (q(u) - c0) / sp,
       lower = (dist$support[1] - c0) / sp, upper = (dist$support[2] - c0) / sp,
       center = c0, spread = sp)
}

#' The split points of rank r (device 3)
#'
#' Why this exists: the density of X_(r) concentrates around c_r = Q(p_r),
#' p_r = r/(n+1), in a region whose width shrinks like n^{-1/2}, while the
#' adaptive rule starts from a fixed set of nodes on the whole (transformed)
#' range. When no initial node lands on the spike, the rule returns a value near
#' zero with a small claimed error (observed with the reference engine for
#' the pair (1, 251) at n = 1001). Splitting the range at c_r - h_r and c_r + h_r puts the spike
#' inside one finite piece whose endpoints the rule sees, and leaves monotone
#' tail pieces. The half-width is the reference engine's,
#'   h_r = 10 s_r / f(c_r) + 4 / sqrt(n + 1),
#' in the standardized units of .os_centered(), with s_r = sqrt{p_r (1 - p_r)/(n + 2)}
#' the standard deviation of the uniform order statistic U_(r), so that
#' s_r / f(c_r) is the large-sample standard deviation of X_(r) by the delta
#' method, generous by the factor 10 for the skewed densities of the extreme
#' ranks, and floored at 4/sqrt(n + 1) spreads of the distribution, which is
#' 4/sqrt(n + 1) in the normal's own units. The delta-method half-width
#' w_r = k s_r / f(c_r), returned by .os_width(), also sets the scale of the
#' absolute tolerances of the integrals of rank r (see .os_single()). A split at
#' c_r +/- w_r as well, when w_r is below the floor, was tried and dropped: it
#' drives the bisection onto an endpoint where the density is infinite (the
#' beta(1/2, 1/2)), and the scaled tolerances alone recover the sharp extremes.
#' A rule in probability units, Q(p_r +/- 10 s_r), was
#' tried first and rejected: for an extreme rank it gives at most one break, so
#' the spike sits in a semi-infinite piece, and the rule then converges to a
#' value off by 1e-10 (the normal at n = 97) to 1e-7 (the logistic at n = 9),
#' which a tighter tolerance does not repair. The breaks are in centered units,
#' and a break outside the support is dropped. If f(c_r) is zero or not finite,
#' only the floor is used.
#' @keywords internal
#' @noRd
.os_breaks <- function(r, n, cd, k) {
  cw <- .os_width(r, n, cd, k)
  h <- cw$w + 4 / sqrt(n + 1)               # the floor is in the standardized units
  b <- c(cw$c - h, cw$c + h)
  b[b > cd$lower & b < cd$upper]
}

#' The location c_r and the delta-method half-width w_r = k s_r / f(c_r) of rank r
#'
#' In standardized units; w_r is 0 where f(c_r) is zero or not finite. Besides
#' the split of .os_breaks(), w_r sets the scale of the absolute tolerances of
#' the integrals of rank r through .os_tolscale(): the moments of a sharp order
#' statistic are small (the variance of the minimum of 10^6 uniform
#' observations is 7e-12 in these units), and the fixed absolute tolerances of
#' the reference engine, set for the normal's scale one, were met before
#' such a moment was resolved (the minimum of 10^6 uniform or exponential
#' observations by quadrature was off by 1e-4 to 1e-3). Each absolute tolerance
#' and acceptance threshold is multiplied by min(1, w_r)^pow for a single
#' integral, min(1, w_s) for the inner integral of a product moment and
#' min(1, w_r) min(1, w_s) for its outer integral, so they are never looser than
#' the engine's constants and tighten with the order statistic; the normal, whose
#' w_r is of order one at n in the thousands, moves by 2e-14.
#' @keywords internal
#' @noRd
.os_width <- function(r, n, cd, k) {
  p <- r / (n + 1)
  s <- sqrt(p * (1 - p) / (n + 2))
  c0 <- cd$q(p)
  w <- k * s / exp(cd$ld(c0))               # k large-sample standard deviations of X_(r)
  list(c = c0, w = if (is.finite(w)) w else 0)
}

#' @keywords internal
#' @noRd
.os_tolscale <- function(w) if (w > 0) min(1, w) else 1

#' integrate() over the pieces cut by the breaks, values and estimates summed
#'
#' Pieces entirely in a tail are monotone-decaying and safe for the adaptive
#' rule. A piece whose convergence message is not "OK" (for example "roundoff
#' error was detected" on a near-zero tail piece, where the relative tolerance is
#' unattainable) is accepted if its reported error estimate is below
#' ctrl$accept.err, and stops the computation otherwise, naming the piece. A
#' piece reported as "probably divergent" is treated by where it lies. On a
#' finite piece the flag is spurious for a probability density, since any power
#' of t times an integrable density is integrable on a bounded interval, and it
#' arises at an integrable singularity of the density (the beta(1/2, 1/2) at its
#' endpoints), so the piece is accepted on its error estimate like any other.
#' On an infinite piece it is accepted only when its value is below
#' ctrl$accept.err as well: a negligible tail piece of a convergent integral can
#' carry the flag (a t parent with 1.5 degrees of freedom gave one of size
#' 1e-15), while a divergent second moment returns a finite value of ordinary
#' size with a small error estimate (a t parent with 1.9 degrees of freedom
#' declared without its tail index gave -11.3 on each tail piece). The error
#' names the piece and the usual causes.
#' @keywords internal
#' @noRd
.os_integrate_split <- function(f, lower, upper, breaks, rel.tol, abs.tol, ctrl, ctx,
                                accept = ctrl$accept.err) {
  breaks <- sort(breaks[breaks > lower & breaks < upper])
  bounds <- c(lower, breaks, upper)
  val <- 0; err <- 0
  for (i in seq_len(length(bounds) - 1L)) {
    out <- tryCatch(integrate(f, bounds[i], bounds[i + 1L], rel.tol = rel.tol, abs.tol = abs.tol,
                              subdivisions = ctrl$subdivisions, stop.on.error = FALSE),
                    error = function(e)
      stop("os_moments(): the integral failed (", ctx, "), piece [", signif(bounds[i], 4), ", ",
           signif(bounds[i + 1L], 4), "] in standardized units: ", conditionMessage(e),
           ". A non-finite value of the integrand means that the density is infinite or ",
           "undefined at a point the quadrature reached: an order statistic within about ",
           "1e-13 of a finite endpoint of the support, relative to the spread, cannot be ",
           "placed in double precision", call. = FALSE))
    divergent <- out$message == "the integral is probably divergent"
    infinite_piece <- !is.finite(bounds[i]) || !is.finite(bounds[i + 1L])
    if (out$message != "OK" &&
        (out$abs.error > accept ||
         (divergent && infinite_piece && abs(out$value) > accept)))
      stop("os_moments(): the integral failed (", ctx, "), piece [",
           signif(bounds[i], 4), ", ", signif(bounds[i + 1L], 4), "] in standardized units: ",
           out$message, " (error estimate ", signif(out$abs.error, 3), "). ",
           if (infinite_piece)
             paste0("On a tail piece this usually means that the moment does not exist ",
                    "(declare the tail indices in os_dist()), or that the tail is so heavy ",
                    "that the mass of the integrand lies far out, as for a lognormal with ",
                    "log-SD 3, which the quadrature does not resolve")
           else
             paste0("On a bounded piece this points at a density that is zero, infinite ",
                    "or not smooth inside the piece"),
           call. = FALSE)
    val <- val + out$value; err <- err + out$abs.error
  }
  list(value = val, abs.error = err)
}

#' log{F(z) - F(y)} for y <= z without cancellation (device 2)
#'
#' Why this exists: the joint density carries {F(z) - F(y)}^(s-r-1), and the
#' log-scale integrand takes (s-r-1) log of the difference, which therefore
#' needs relative precision even when y and z sit deep in the same tail, where
#' F(y) and F(z) are both near 0 or both near 1. The difference is formed on
#' the side where the numbers are small: below the median (z <= 0 in centered
#' units) as log F(z) + log{1 - F(y)/F(z)}, above it (y >= 0) as
#' log S(y) + log{1 - S(z)/S(y)} with S = 1 - F, and directly when the pair
#' straddles the median, where the difference is of ordinary size. A pair whose
#' two log-tails both underflow to -Inf gives NaN, which is set to -Inf: the
#' difference is zero in doubles. y is a scalar (the outer variable) and z a
#' vector (the inner nodes). The caller may pass the values it already holds,
#' log F(y), log S(y), F(y) and log S(z), so that each is evaluated once per
#' node; the arithmetic is the same either way.
#' @keywords internal
#' @noRd
.os_logdiff <- function(y, z, cd, lFy = cd$lF(y), lSy = cd$lS(y), pfy = cd$pf(y),
                        lSz = NULL) {
  out <- numeric(length(z))
  a <- z <= 0
  b <- !a & (y >= 0)
  m <- !a & !b
  if (any(a)) { lFz <- cd$lF(z[a]); out[a] <- lFz + log1p(-exp(lFy - lFz)) }
  if (any(b)) out[b] <- lSy + log1p(-exp((if (is.null(lSz)) cd$lS(z[b]) else lSz[b]) - lSy))
  if (any(m)) out[m] <- log(cd$pf(z[m]) - pfy)
  out[is.nan(out)] <- -Inf
  out
}

#' k * log(v), returning 0 when k = 0 (guards 0 * -Inf = NaN)
#' @keywords internal
#' @noRd
.os_xlog <- function(k, logv) if (k == 0) 0 else k * logv

#' E[(X_(r) - shift)^pow] in the engine's units by one exact integral
#'
#' Density: f_{r:n}(t) = C F(t)^(r-1) {1 - F(t)}^(n-r) f(t) with
#' C = n!/{(r-1)! (n-r)!}, assembled on the log scale (device 1): lbeta for the
#' coefficient (the difference of three lgamma values lost 1e-9 of a central
#' rank's covariance at n = 1e6), the log-tails for the powers, the log-density,
#' exponentiated once. With pow = 0 the integral is the mass of the density,
#' which os_moments() checks against 1.
#' The moment is taken about `shift`, which os_moments() sets to the order
#' statistic's own location c_r = Q(r/(n+1)) (device 4): with the moments taken
#' about the median, a variance was the small difference of two moments of size
#' (c_r - c)^2, and its relative error grew like that size over the variance,
#' n^2 for the uniform minimum (7 digits left at n = 1000, none at
#' n = 1e5). About c_r the first moment is small and nothing cancels. For pow = 1,
#' the only power at which the integrand changes sign, the shift is added to the
#' split points, so that (t - shift) times the density is one-signed on each
#' piece: a piece straddling the sign change has a near-zero value that the
#' relative tolerance cannot resolve, and integrate() then flags roundoff. For
#' pow = 0 and 2 the split would only drive the bisection onto a singular
#' endpoint. The absolute tolerance and the acceptance threshold are scaled by
#' the order statistic's width to the power pow (.os_width()).
#' @keywords internal
#' @noRd
.os_single <- function(r, n, pow, cd, ctrl, shift = 0) {
  lC <- -lbeta(r, n - r + 1)                 # log n!/{(r-1)! (n-r)!}, without the cancellation
                                             # of three large lgamma values at n = 1e6
  f <- function(t) {
    lw <- lC + .os_xlog(r - 1, cd$lF(t)) + .os_xlog(n - r, cd$lS(t)) + cd$ld(t)
    (t - shift)^pow * exp(lw)
  }
  sc <- .os_tolscale(.os_width(r, n, cd, ctrl$k)$w)^pow   # the moment's own scale, at most 1
  br <- .os_breaks(r, n, cd, ctrl$k)
  if (pow == 1) br <- c(br, shift)
  .os_integrate_split(f, cd$lower, cd$upper, br, ctrl$rel.tol, ctrl$abs.tol.single * sc, ctrl,
                      ctx = paste0("E[X^", pow, "], r = ", r, ", n = ", n),
                      accept = ctrl$accept.err * sc)
}

#' E[(X_(r) - sr)(X_(s) - ss)] in the engine's units, r < s, by nested exact quadrature
#'
#' Joint density for y < z:
#'   f(y, z) = C F(y)^(r-1) {F(z) - F(y)}^(s-r-1) {1 - F(z)}^(n-s) f(y) f(z),
#'   C = n!/{(r-1)! (s-r-1)! (n-s)!}.
#' The outer integral runs over y in the support and the inner over z in
#' (y, upper), split at the breaks of rank s. The outer range is split at the
#' breaks of both ranks: at those of rank r for the spike of y's own density, and
#' at those of rank s because the inner quadrature changes its piece structure
#' where y crosses them, a feature of the outer integrand that the adaptive rule
#' otherwise does not see. With the breaks of rank r alone, the pair (9, 20) at
#' n = 70 of the normal converged to a value 5.5e-10 off with an error estimate
#' of 1e-12, which the permutation identity caught and an integration in the
#' probability scale confirmed; with both sets of breaks the value is recovered
#' to 3e-14. The middle factor uses .os_logdiff() and is skipped when s = r + 1
#' (exponent 0). The inner error estimate is not carried into the outer one: the
#' returned estimate is the outer integral's own. The tails at the outer point y
#' are evaluated once per inner integral, and log S(z) once per inner node, which
#' saves about a tenth of the time and leaves every value unchanged. The moment
#' is taken about the two locations sr and ss of the order statistics (device 4,
#' see .os_single()), which are added to the split points of their ranks, so
#' that the covariance E[(Y - c_r)(Z - c_s)] - m_r m_s is formed from small
#' quantities. The inner absolute tolerance is divided by max(1, |y - sr|),
#' the outer factor the inner value is multiplied by: at a fixed 1e-15 a
#' lognormal with log-SD 2, whose outer nodes reach y = 1e5, lost 3.5e-6 of a
#' covariance of 38.8 at n = 2, and the scaling leaves 5e-8 at no measurable
#' cost (the outer nodes of the ordinary families have |y| below 10). The
#' inner absolute tolerance and acceptance threshold are scaled by the width of
#' rank s, the outer ones by the widths of both ranks (.os_width()).
#' @keywords internal
#' @noRd
.os_prod <- function(r, s, n, cd, ctrl, sr = 0, ss = 0) {
  lC <- -lbeta(r, n - r + 1) - lbeta(s - r, n - s + 1)   # log n!/{(r-1)! (s-r-1)! (n-s)!}
  m <- s - r - 1
  br_r <- c(.os_breaks(r, n, cd, ctrl$k), sr)
  br_s <- c(.os_breaks(s, n, cd, ctrl$k), ss)
  sc_r <- .os_tolscale(.os_width(r, n, cd, ctrl$k)$w)   # the two order statistics' widths
  sc_s <- .os_tolscale(.os_width(s, n, cd, ctrl$k)$w)
  ctx_in <- paste0("inner, r = ", r, ", s = ", s, ", n = ", n)
  inner <- function(y) {                      # y scalar: int_y^upper z f(y, z) dz
    lFy <- cd$lF(y); lSy <- cd$lS(y); pfy <- cd$pf(y)
    ly <- lC + .os_xlog(r - 1, lFy) + cd$ld(y)
    g <- function(z) {
      lSz <- cd$lS(z)
      lmid <- if (m == 0) 0 else m * .os_logdiff(y, z, cd, lFy, lSy, pfy, lSz)
      (z - ss) * exp(ly + lmid + .os_xlog(n - s, lSz) + cd$ld(z))
    }
    .os_integrate_split(g, y, cd$upper, br_s, ctrl$rel.tol,
                        ctrl$abs.tol.inner * sc_s / max(1, abs(y - sr)), ctrl, ctx = ctx_in,
                        accept = ctrl$accept.err * sc_s)$value
  }
  f_outer <- function(y) (y - sr) * vapply(y, inner, numeric(1))
  .os_integrate_split(f_outer, cd$lower, cd$upper, sort(unique(c(br_r, br_s))), ctrl$rel.tol.outer,
                      ctrl$abs.tol.outer * sc_r * sc_s, ctrl,
                      ctx = paste0("outer, r = ", r, ", s = ", s, ", n = ", n),
                      accept = ctrl$accept.err * sc_r * sc_s)
}
