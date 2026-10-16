# ==============================================================================
# File: os_dist.R
# Purpose: The distribution object of osmoments: the public constructor
#          os_dist() for a user-supplied distribution, the internal constructor
#          of the five built-in families, the validator, and the internal
#          Gumbel functions (maximum form) with accurate log-tails.
#
# EXPORTED FUNCTIONS:
#   - os_dist()           Build a distribution object from (p, d, q)
#
# INTERNAL HELPERS:
#   - .os_family()        The built-in families: norm, unif, exp, logis, gumbel
#   - .os_validate_dist() Check the functions, the support, the tails, the center
#   - .os_resolve_dist()  A family name or an os_dist object -> the object
#   - .pgumbel(), .dgumbel(), .qgumbel()  Gumbel (max form), F(x) = exp(-e^{-x})

#################### EXPORTED MAIN FUNCTIONS ####################

#' Define a distribution for the moments of its order statistics
#'
#' Builds the distribution object that \code{\link{os_moments}} takes in its
#' \code{dist} argument, from the distribution function, the density and the
#' quantile function of a continuous distribution. The object describes the
#' standard member of a location-scale family: \code{os_moments()} computes the
#' moments for this member and applies its own \code{location} and
#' \code{scale} arguments afterwards. A user who wants a \eqn{t_5} distribution
#' with location 3 and scale 2 supplies the \eqn{t_5} functions here and passes
#' \code{location = 3, scale = 2} to \code{os_moments()}.
#'
#' @param p The distribution function, with the signature of R's
#'   \code{p*} functions: \code{p(q, lower.tail = TRUE, log.p = FALSE)},
#'   vectorized in \code{q}. The \code{log.p} and \code{lower.tail} arguments
#'   must be honored, because the integrands are evaluated on the log scale in
#'   both tails.
#' @param d The density, \code{d(x, log = FALSE)}, vectorized in \code{x}, with
#'   \code{log} honored.
#' @param q The quantile function, \code{q(p)}, vectorized in \code{p}.
#' @param support Numeric vector of length two, the endpoints of the support,
#'   finite or infinite. They are the limits of integration.
#' @param symmetric Logical. \code{TRUE} if the density is symmetric about its
#'   median, that is \eqn{Q(1 - u) = 2 Q(1/2) - Q(u)}. A symmetric distribution
#'   lets \code{os_moments()} halve the number of double integrals of a full
#'   matrix and check the reflection symmetry of its result.
#' @param tail Numeric vector of length two, the tail indices
#'   \eqn{(\nu_L, \nu_U)}: the density decays like \eqn{|x|^{-\nu_L - 1}} in the
#'   lower tail and like \eqn{x^{-\nu_U - 1}} in the upper tail. \code{Inf}
#'   for a tail that decays faster than any power, or for a bounded tail. The
#'   indices decide which moments exist (see Details).
#' @param name A label for printing.
#'
#' @return An object of class \code{"os_dist"}: a list with the components
#'   \code{name}, \code{p}, \code{d}, \code{q}, \code{support},
#'   \code{symmetric}, \code{tail}, \code{center} (the median \eqn{Q(1/2)}) and
#'   \code{closed} (the closed-form moments the distribution has, none for a
#'   user distribution).
#'
#' @details
#' \strong{The moment condition.} For a parent whose density decays like
#' \eqn{|x|^{-\nu - 1}}, the density of \eqn{X_{(r)}} decays like
#' \eqn{|x|^{-\nu r - 1}} below and like \eqn{x^{-\nu (n - r + 1) - 1}} above,
#' so \eqn{E|X_{(r)}|^k < \infty} if and only if \eqn{\nu_L r > k} and
#' \eqn{\nu_U (n - r + 1) > k}. \code{os_moments()} applies this condition with
#' \eqn{k = 2} to every requested rank and refuses a moment that does not exist.
#' A Cauchy parent, \code{tail = c(1, 1)}, has a median of five with a finite
#' variance and extremes without a mean.
#'
#' \strong{Accuracy requirement.} The log-tails \code{p(x, log.p = TRUE)} and
#' \code{p(x, lower.tail = FALSE, log.p = TRUE)} must keep relative precision
#' deep in both tails, as R's base distributions do. A function that computes
#' \code{log(p(x))} from an underflowing probability loses the tail, and the
#' integrals with it. \code{os_dist()} checks that the three functions honor
#' their arguments at the median and are vectorized, and that \code{d}
#' integrates to \eqn{1/2} between the quartiles, where \code{p} rises by
#' \eqn{1/2}. \code{os_moments()} checks that the log-tails and the log-density
#' are finite at the median and at the split points of the smallest and the
#' largest requested rank, and that the density of every order statistic it
#' integrates has mass 1. None of these certifies accuracy.
#'
#' \strong{Location limit.} The integrals are taken about the median of the
#' standard member, and the integration range is split in units of its spread
#' (its interquartile range), so the scale on which the standard member is
#' supplied does not matter: the package's tests check a normal supplied with
#' standard deviations from \eqn{10^{-3}} to \eqn{10^6}. Its location does: a standard member whose own
#' location is large, say \eqn{10^8} with unit scale, loses the digits of
#' \eqn{x - Q(1/2)} in double precision and the computation stops with an
#' error. Supply the standard member and pass the location to
#' \code{os_moments()} instead, which also keeps \code{os_blue()} exact.
#'
#' \strong{Common parents.} Each function wraps the corresponding
#' \code{stats} function with its shape fixed and passes \code{lower.tail},
#' \code{log.p} and \code{log} by name, as in the examples. The declarations
#' that differ between parents are these.
#' \tabular{lll}{
#'   Parent \tab \code{support} \tab \code{tail} \cr
#'   Student \eqn{t_\nu} (\code{pt}, \code{dt}, \code{qt}), symmetric \tab \code{c(-Inf, Inf)} \tab \code{c(nu, nu)} \cr
#'   gamma (\code{pgamma}, \code{dgamma}, \code{qgamma}) \tab \code{c(0, Inf)} \tab \code{c(Inf, Inf)} \cr
#'   Weibull (\code{pweibull}, \code{dweibull}, \code{qweibull}) \tab \code{c(0, Inf)} \tab \code{c(Inf, Inf)} \cr
#'   lognormal (\code{plnorm}, \code{dlnorm}, \code{qlnorm}) \tab \code{c(0, Inf)} \tab \code{c(Inf, Inf)} \cr
#'   beta (\code{pbeta}, \code{dbeta}, \code{qbeta}) \tab \code{c(0, 1)} \tab \code{c(Inf, Inf)}
#' }
#' The lognormal tail decays faster than any power, so every moment exists.
#' A distribution declared on \code{c(-Inf, Inf)} whose density is zero below
#' some point stops with an error that its log-tail is not finite, which is
#' the sign that \code{support} is missing.
#'
#' @examples
#' # Student t with 5 degrees of freedom: tail index 5 in both tails
#' t5 <- os_dist(function(q, lower.tail = TRUE, log.p = FALSE)
#'                 stats::pt(q, df = 5, lower.tail = lower.tail, log.p = log.p),
#'               function(x, log = FALSE) stats::dt(x, df = 5, log = log),
#'               function(p) stats::qt(p, df = 5),
#'               symmetric = TRUE, tail = c(5, 5), name = "t5")
#' os_moments(5, dist = t5)
#'
#' # The Cauchy: the median of five has a variance, the extremes have no mean
#' cauchy <- os_dist(stats::pcauchy, stats::dcauchy, stats::qcauchy,
#'                   symmetric = TRUE, tail = c(1, 1), name = "cauchy")
#' os_moments(5, ranks = 3, dist = cauchy)
#'
#' # A gamma with shape 2: a semi-bounded support, every moment finite
#' gamma2 <- os_dist(function(q, lower.tail = TRUE, log.p = FALSE)
#'                     stats::pgamma(q, shape = 2, lower.tail = lower.tail, log.p = log.p),
#'                   function(x, log = FALSE) stats::dgamma(x, shape = 2, log = log),
#'                   function(p) stats::qgamma(p, shape = 2),
#'                   support = c(0, Inf), name = "gamma2")
#' os_moments(10, ranks = c(1, 5, 10), dist = gamma2, scale = 3)
#' @seealso \code{\link{os_moments}}
#' @export
os_dist <- function(p, d, q, support = c(-Inf, Inf), symmetric = FALSE,
                    tail = c(Inf, Inf), name = "user") {
  obj <- structure(list(name = name, p = p, d = d, q = q, support = support,
                        symmetric = symmetric, tail = tail, center = NA_real_,
                        closed = list(mean = NULL, var = NULL, cov = NULL)),
                   class = "os_dist")
  .os_validate_dist(obj)
}

#################### INTERNAL HELPER FUNCTIONS ####################

#' The five built-in families as os_dist objects
#'
#' Each family is its standard member: N(0, 1); the uniform on [0, 1], whose
#' location is the lower end and whose scale is the length; the unit-rate
#' exponential, whose scale is 1/rate; the standard logistic; the standard
#' Gumbel in the maximum form, F(x) = exp(-exp(-x)). The closed-form hooks
#' return the mean, the variance and the covariance in the raw units of the
#' standard member, or NA where the family has no closed form for that rank
#' (the Gumbel below the maximum).
#' @keywords internal
#' @noRd
.os_family <- function(name) {
  name <- match.arg(name, c("norm", "unif", "exp", "logis", "gumbel"))
  obj <- switch(name,
    norm = os_dist(pnorm, dnorm, qnorm, support = c(-Inf, Inf), symmetric = TRUE,
                   name = "norm"),
    unif = os_dist(punif, dunif, qunif, support = c(0, 1), symmetric = TRUE,
                   name = "unif"),
    exp = os_dist(pexp, dexp, qexp, support = c(0, Inf), symmetric = FALSE,
                  name = "exp"),
    logis = os_dist(plogis, dlogis, qlogis, support = c(-Inf, Inf), symmetric = TRUE,
                    name = "logis"),
    gumbel = os_dist(.pgumbel, .dgumbel, .qgumbel, support = c(-Inf, Inf),
                     symmetric = FALSE, name = "gumbel"))
  obj$closed <- switch(name,
    norm = list(mean = NULL, var = NULL, cov = NULL),
    unif = list(mean = .os_closed_unif_mean, var = .os_closed_unif_var,
                cov = .os_closed_unif_cov),
    exp = list(mean = .os_closed_exp_mean, var = .os_closed_exp_var,
               cov = .os_closed_exp_cov),
    logis = list(mean = .os_closed_logis_mean, var = .os_closed_logis_var, cov = NULL),
    gumbel = list(mean = .os_closed_gumbel_mean, var = .os_closed_gumbel_var, cov = NULL))
  obj
}

#' A family name or an os_dist object, resolved to the object
#' @keywords internal
#' @noRd
.os_resolve_dist <- function(dist) {
  if (inherits(dist, "os_dist")) return(dist)
  families <- c("norm", "unif", "exp", "logis", "gumbel")
  if (is.character(dist) && length(dist) == 1L && dist %in% families) return(.os_family(dist))
  stop("os_moments(): 'dist' must be one of ", paste0("\"", families, "\"", collapse = ", "),
       ", or an object made by os_dist()", call. = FALSE)
}

#' Validate a distribution object and compute its center
#'
#' Checks the types, the support, the tail indices and that the three functions
#' honor their arguments at the median and are vectorized. Returns the object
#' with its center Q(1/2) filled in.
#' @keywords internal
#' @noRd
.os_validate_dist <- function(obj) {
  fail <- function(...) stop("os_dist(): ", ..., call. = FALSE)
  if (!is.function(obj$p) || !is.function(obj$d) || !is.function(obj$q))
    fail("'p', 'd' and 'q' must be functions")
  if (!is.numeric(obj$support) || length(obj$support) != 2L || anyNA(obj$support) ||
      obj$support[1] >= obj$support[2])
    fail("'support' must be two ordered endpoints, finite or infinite")
  if (!is.logical(obj$symmetric) || length(obj$symmetric) != 1L || is.na(obj$symmetric))
    fail("'symmetric' must be TRUE or FALSE")
  if (!is.numeric(obj$tail) || length(obj$tail) != 2L || anyNA(obj$tail) || any(obj$tail <= 0))
    fail("'tail' must be two positive tail indices (Inf allowed)")
  if (!is.character(obj$name) || length(obj$name) != 1L)
    fail("'name' must be a single string")
  # Each call names the signature it needs when it fails, and each value is checked for
  # being one finite number before it is compared, so that a missing argument or an NA
  # gives an error that says what to fix.
  sig_p <- "p(q, lower.tail = TRUE, log.p = FALSE)"
  sig_d <- "d(x, log = FALSE)"
  sig_q <- "q(p)"
  call_fn <- function(what, sig, expr) tryCatch(expr, error = function(e)
    fail(what, " failed (", conditionMessage(e), "); the function must have the form ", sig))
  near <- function(v, target) is.numeric(v) && length(v) == 1L && is.finite(v) &&
    abs(v - target) <= tol
  show <- function(v) if (is.numeric(v) && length(v) == 1L) signif(v, 8) else deparse(v)
  tol <- 1e-8
  # the median, and the three functions evaluated there
  cen <- call_fn("q(0.5)", sig_q, obj$q(0.5))
  if (!is.numeric(cen) || length(cen) != 1L || !is.finite(cen))
    fail("q(0.5) must return one finite number")
  if (cen < obj$support[1] || cen > obj$support[2])
    fail("the median q(0.5) = ", signif(cen, 6), " lies outside 'support'")
  pm <- call_fn("p(q(0.5))", sig_p, obj$p(cen))
  if (!near(pm, 0.5))
    fail("p(q(0.5)) must be 0.5, got ", show(pm), ": 'p' and 'q' disagree")
  # the two tails differ at the lower quartile, where the median cannot tell them apart
  x1 <- call_fn("q(0.25)", sig_q, obj$q(0.25))
  if (!is.numeric(x1) || length(x1) != 1L || !is.finite(x1))
    fail("q(0.25) must return one finite number")
  if (!near(call_fn("p(x)", sig_p, obj$p(x1)), 0.25))
    fail("p(q(0.25)) must be 0.25: 'p' and 'q' disagree")
  if (!near(call_fn("p(x, lower.tail = FALSE)", sig_p, obj$p(x1, lower.tail = FALSE)), 0.75))
    fail("p(x, lower.tail = FALSE) does not return the upper tail 1 - F(x)")
  if (!near(call_fn("p(x, log.p = TRUE)", sig_p, obj$p(x1, log.p = TRUE)), log(0.25)))
    fail("p(x, log.p = TRUE) does not return the log probability log F(x)")
  if (!near(call_fn("p(x, lower.tail = FALSE, log.p = TRUE)", sig_p,
                    obj$p(x1, lower.tail = FALSE, log.p = TRUE)), log(0.75)))
    fail("p(x, lower.tail = FALSE, log.p = TRUE) does not return the log upper tail")
  dm <- call_fn("d(x)", sig_d, obj$d(cen))
  dl <- call_fn("d(x, log = TRUE)", sig_d, obj$d(cen, log = TRUE))
  if (!is.numeric(dm) || length(dm) != 1L || !is.finite(dm) || dm <= 0)
    fail("d() must be finite and positive at the median")
  if (!near(dl, log(dm)))
    fail("d(x, log = TRUE) does not return the log density at the median")
  # the density against the distribution function: d integrates to F(q3) - F(q1) = 1/2
  # between the quartiles (a density with the wrong constant or Jacobian, 1.001 dnorm,
  # passed every check at the median and gave a variance of 1.001)
  x3 <- call_fn("q(0.75)", sig_q, obj$q(0.75))
  mass <- tryCatch(integrate(function(x) obj$d(x), x1, x3, rel.tol = 1e-8)$value,
                   error = function(e) fail("d() could not be integrated between the quartiles: ",
                                            conditionMessage(e)))
  if (!is.finite(mass) || abs(mass - 0.5) > 1e-6)
    fail("d() integrates to ", signif(mass, 7), " between the quartiles, where F rises by 0.5: ",
         "'d' and 'p' disagree (a normalizing constant or a Jacobian is off)")
  # vectorization
  xx <- c(cen, cen, cen)
  if (length(call_fn("p(x)", sig_p, obj$p(xx))) != 3L ||
      length(call_fn("d(x)", sig_d, obj$d(xx))) != 3L)
    fail("'p' and 'd' must be vectorized in their first argument")
  qq <- call_fn("q(p)", sig_q, obj$q(c(0.25, 0.75)))
  if (length(qq) != 2L || !all(is.finite(qq)) || qq[1] >= qq[2])
    fail("'q' must be vectorized and increasing in p")
  obj$center <- cen
  obj
}

#' The standard Gumbel distribution in the maximum form, F(x) = exp(-exp(-x))
#'
#' log F = -exp(-x) directly. log S = log(1 - F) as log(-expm1(log F)) up to
#' x = 30, and as the series -x + log(1 - exp(-x)/2) beyond, where exp(-x) is
#' near the double epsilon and -expm1(log F) loses its relative precision. These
#' are the two tail repairs of the feasibility check of the literature report.
#' @keywords internal
#' @noRd
.pgumbel <- function(q, lower.tail = TRUE, log.p = FALSE) {
  lF <- -exp(-q)
  if (lower.tail) {
    v <- lF
  } else {
    v <- log(-expm1(lF))
    big <- q > 30
    v[big] <- -q[big] + log1p(-exp(-q[big]) / 2)
  }
  if (log.p) v else exp(v)
}

#' @keywords internal
#' @noRd
.dgumbel <- function(x, log = FALSE) {
  v <- -x - exp(-x)
  if (log) v else exp(v)
}

#' @keywords internal
#' @noRd
.qgumbel <- function(p) -log(-log(p))
