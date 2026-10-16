# ==============================================================================
# File: os_checks.R
# Purpose: The conditions and run-time checks of osmoments: the moment
#          condition by rank from the tail indices, the finiteness probe of a
#          distribution's log-tails at the split points, the reflection check
#          for a symmetric parent, and the permutation identities on a complete
#          rank set.
#
# INTERNAL HELPERS:
#   - .os_moment_condition()  E|X_(r)|^k < Inf iff nu_L r > k and nu_U (n-r+1) > k
#   - .os_tail_probe()        log F, log S, log f finite at the split points and the center
#   - .os_reflection_check()  alpha_i = -alpha_{k+1-i}, E[X_(r_i)^2] = E[X_(r_{k+1-i})^2]
#   - .os_parent_moments()    E[X] and E[X^2] of the standard member, centered
#   - .os_identity_check()    the permutation identities, and the unit row sums for the normal

#################### INTERNAL HELPER FUNCTIONS ####################

#' The moment condition, rank by rank
#'
#' For a parent whose density decays like |x|^(-nu_L - 1) below and x^(-nu_U - 1)
#' above, the density of X_(r) decays like |x|^(-nu_L r - 1) below and
#' x^(-nu_U (n-r+1) - 1) above, so E|X_(r)|^k is finite iff nu_L r > k and
#' nu_U (n - r + 1) > k. Stops with the rank and the failing bound.
#' @keywords internal
#' @noRd
.os_moment_condition <- function(ranks, n, tail, k = 2) {
  for (r in ranks) {
    if (!(tail[1] * r > k))
      stop("os_moments(): the moment of order ", k, " of X_(", r, ") does not exist: ",
           "the lower tail index ", tail[1], " times the rank ", r, " is not above ", k,
           call. = FALSE)
    if (!(tail[2] * (n - r + 1) > k))
      stop("os_moments(): the moment of order ", k, " of X_(", r, ") does not exist: ",
           "the upper tail index ", tail[2], " times n - r + 1 = ", n - r + 1,
           " is not above ", k, call. = FALSE)
  }
  invisible(TRUE)
}

#' Probe the log-tails and the log-density where the integrands are evaluated
#'
#' Evaluates log F, log S and log f of the centered standard member at the split
#' points of the smallest and the largest requested rank and at the center, and
#' stops if any value is not finite. Catches a user function that ignores log.p
#' or lower.tail, or whose tail underflows where the integrand still matters. It
#' does not certify accuracy.
#' @keywords internal
#' @noRd
.os_tail_probe <- function(ranks, n, cd, ctrl) {
  pts <- unique(c(0, .os_breaks(min(ranks), n, cd, ctrl$k), .os_breaks(max(ranks), n, cd, ctrl$k)))
  xs <- cd$center + cd$spread * pts            # the same points in the standard member's units
  vals <- list(`p(x, log.p = TRUE)` = cd$lF(pts),
               `p(x, lower.tail = FALSE, log.p = TRUE)` = cd$lS(pts),
               `d(x, log = TRUE)` = cd$ld(pts))
  for (nm in names(vals)) {
    bad <- !is.finite(vals[[nm]])
    if (any(bad))
      stop("os_moments(): ", nm, " is not finite at x = ",
           paste(signif(xs[bad], 6), collapse = ", "),
           "; the distribution's log-tails must stay finite ",
           "where the order statistics have their mass. Check that 'support' in os_dist() ",
           "gives the endpoints of a bounded or semi-bounded distribution, and that 'p' ",
           "honors log.p and lower.tail", call. = FALSE)
  }
  invisible(TRUE)
}

#' The reflection residual of a symmetric parent on a reflection-closed rank set
#'
#' In centered units, alpha_i = -alpha_{k+1-i} and E[X_(r_i)^2] = E[X_(r_{k+1-i})^2].
#' The integrators do not impose these, so a violation signals a quadrature
#' failure. The reference engine, an earlier implementation of the same method
#' for the normal, checks its results in the same way.
#' @keywords internal
#' @noRd
.os_reflection_check <- function(a_c, sq_c) {
  max(abs(a_c + rev(a_c)), abs(sq_c - rev(sq_c)))
}

#' E[X] and E[X^2] of the standard member, in centered units
#'
#' From the closed forms at n = 1 when the family has them, and from the
#' engine's own n = 1 integrals otherwise, so that the identities need no
#' closed form of the parent.
#' @keywords internal
#' @noRd
.os_parent_moments <- function(dist, cd, ctrl) {
  cl <- dist$closed
  m <- if (is.null(cl$mean)) NA_real_ else cl$mean(1L, 1L)
  v <- if (is.null(cl$var)) NA_real_ else cl$var(1L, 1L)
  if (is.na(m) || is.na(v)) .os_mass_check(1L, 1L, cd, ctrl, 0)
  if (is.na(m)) m <- cd$spread * .os_single(1L, 1L, 1, cd, ctrl)$value + dist$center
  m1 <- m - dist$center
  if (is.na(v)) {
    m2 <- cd$spread * (cd$spread * .os_single(1L, 1L, 2, cd, ctrl)$value)
  } else {
    m2 <- v + m1^2
  }
  list(m1 = m1, m2 = m2)
}

#' The permutation identities on a complete rank set
#'
#' With every rank requested, the order statistics are a permutation of the
#' sample, so for any parent with a finite second moment
#'   sum_r alpha_r = n E[X], sum_r E[X_(r)^2] = n E[X^2], sum_{r,s} B_rs = n Var(X).
#' For the normal, where the sample mean is independent of every deviation
#' X_(r) - mean, each row of B sums to sigma^2 as well. Returns the absolute
#' residuals, in centered units of the standard member.
#' @keywords internal
#' @noRd
.os_identity_check <- function(a_c, v, B_c, n, pm, is_norm) {
  out <- c(sum_alpha = abs(sum(a_c) - n * pm$m1),
           sum_sq = abs(sum(v + a_c^2) - n * pm$m2))
  if (!is.null(B_c)) {
    out <- c(out, sum_B = abs(sum(B_c) - n * (pm$m2 - pm$m1^2)))
    if (is_norm) out <- c(out, row_sums = max(abs(rowSums(B_c) - 1)))
  }
  out
}
