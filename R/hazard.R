#' Graph-agnostic hazard functions for BOCPD
#'
#' Hazard functions return the changepoint probability at time \code{t}
#' conditional on current run length(s) \code{run_length}.
#'
#'
#' @name hazard
NULL



#' Constant hazard function
#'
#' @param p Scalar changepoint probability in `[0,1]`.
#'
#' @return A hazard function \code{h(run_length, t, ...)} returning a numeric
#'   vector of probabilities (same length as \code{run_length}).
#'
#' @examples
#' h <- hazard_constant(0.01)
#' h(0:5, t = 10)
#'
#' @export
hazard_constant <- function(p) {
  if (!is.numeric(p) || length(p) != 1L || !is.finite(p)) {
    stop("`p` must be a single finite numeric value.", call. = FALSE)
  }
  p <- as.numeric(p)
  if (p < 0 || p > 1) stop("`p` must be in [0,1].", call. = FALSE)

  function(run_length, t, ...) {
    rep(p, length(run_length))
  }
}

#' Geometric hazard (constant changepoint probability)
#'
#' This is an alias for \code{hazard_constant()}, included for clarity because
#' the induced segment length distribution is geometric.
#'
#' @param p Scalar changepoint probability in `[0,1]`.
#'
#' @return A hazard function.
#'
#' @export
hazard_geometric <- function(p) {
  hazard_constant(p)
}



#' Piecewise-constant hazard over time
#'
#' The hazard depends on time \code{t} only via user-supplied cutpoints.
#'
#' @param times Integer vector of time cutpoints (strictly increasing).
#'   Hazard \code{p[i]} is used for \code{t <= times[i]} and the last value is used
#'   after the final cutpoint.
#' @param p Numeric vector of probabilities in `[0,1]` of length \code{length(times)+1}.
#'
#' @return A hazard function.
#'
#' @examples
#' h <- hazard_piecewise_time(times = c(50, 100), p = c(0.01, 0.05, 0.02))
#' h(0:3, t = 25)
#' h(0:3, t = 80)
#'
#' @export
hazard_piecewise_time <- function(times, p) {
  if (!is.numeric(times) || any(!is.finite(times))) {
    stop("`times` must be a numeric/integer vector of finite values.", call. = FALSE)
  }
  times <- as.integer(times)
  if (length(times) > 0L && any(diff(times) <= 0L)) {
    stop("`times` must be strictly increasing.", call. = FALSE)
  }

  if (!is.numeric(p) || any(!is.finite(p))) {
    stop("`p` must be a numeric vector of finite values.", call. = FALSE)
  }
  p <- as.numeric(p)
  if (length(p) != length(times) + 1L) {
    stop("`p` must have length length(times) + 1.", call. = FALSE)
  }
  if (any(p < 0) || any(p > 1)) stop("All entries of `p` must be in [0,1].", call. = FALSE)

  function(run_length, t, ...) {
    if (!is.numeric(t) || length(t) != 1L || !is.finite(t)) {
      stop("`t` must be a single finite number.", call. = FALSE)
    }
    t <- as.integer(t)
    idx <- if (length(times) == 0L) 1L else (findInterval(t, times) + 1L)
    rep(p[[idx]], length(run_length))
  }
}



#' Logistic hazard increasing with run length
#'
#' A simple parametric hazard that increases with run length:
#' \deqn{h(r) = \mathrm{logit}^{-1}(a + b r).}
#'
#' @param a Intercept on logit scale.
#' @param b Slope on logit scale.
#'
#' @return A hazard function depending on \code{run_length} only.
#'
#' @examples
#' h <- hazard_logistic_run_length(a = -6, b = 0.05)
#' h(0:10, t = 1)
#'
#' @export
hazard_logistic_run_length <- function(a, b) {
  if (!is.numeric(a) || length(a) != 1L || !is.finite(a)) {
    stop("`a` must be a single finite number.", call. = FALSE)
  }
  if (!is.numeric(b) || length(b) != 1L || !is.finite(b)) {
    stop("`b` must be a single finite number.", call. = FALSE)
  }
  a <- as.numeric(a)
  b <- as.numeric(b)

  function(run_length, t, ...) {
    r <- as.numeric(run_length)
    eta <- a + b * r
    # stable inverse-logit
    out <- ifelse(
      eta >= 0,
      1 / (1 + exp(-eta)),
      exp(eta) / (1 + exp(eta))
    )
    # enforce bounds for numerical noise
    pmin(pmax(out, 0), 1)
  }
}



#' Discrete hazard from a segment-length PMF
#'
#' Converts a discrete segment length distribution into a hazard.
#' For a segment length L with PMF \code{pmf[k] = P(L = k)} for k=1..K,
#' the hazard at run length r (0-indexed) is:
#' \deqn{h(r) = P(L = r+1 \mid L \ge r+1) = pmf[r+1] / S[r+1],}
#' where \eqn{S[j] = P(L \ge j)}.
#'
#' @param pmf Numeric vector of nonnegative values (will be normalized).
#'   Interpreted as probabilities for lengths 1..K.
#' @param tail_hazard Hazard used for run lengths r >= K (default 1, i.e. force a changepoint).
#'
#' @return A hazard function depending on \code{run_length}.
#'
#' @examples
#' # Length is uniform on {1,2,3}
#' h <- hazard_from_pmf(c(1, 1, 1))
#' h(0:5, t = 1)
#'
#' @export
hazard_from_pmf <- function(pmf, tail_hazard = 1) {
  if (!is.numeric(pmf) || length(pmf) < 1L || any(!is.finite(pmf))) {
    stop("`pmf` must be a non-empty numeric vector of finite values.", call. = FALSE)
  }
  pmf <- as.numeric(pmf)
  if (any(pmf < 0)) stop("`pmf` entries must be nonnegative.", call. = FALSE)
  s <- sum(pmf)
  if (s <= 0) stop("`pmf` must have positive sum.", call. = FALSE)
  pmf <- pmf / s

  if (!is.numeric(tail_hazard) || length(tail_hazard) != 1L || !is.finite(tail_hazard)) {
    stop("`tail_hazard` must be a single finite number.", call. = FALSE)
  }
  tail_hazard <- as.numeric(tail_hazard)
  if (tail_hazard < 0 || tail_hazard > 1) stop("`tail_hazard` must be in [0,1].", call. = FALSE)

  K <- length(pmf)
  # Survival S[j] = P(L >= j), j=1..K
  surv <- rev(cumsum(rev(pmf)))  # length K, positive

  function(run_length, t, ...) {
    r <- as.integer(run_length)
    out <- numeric(length(r))

    # For 0 <= r <= K-1, use pmf[r+1] / surv[r+1]
    inside <- which(r >= 0L & r <= (K - 1L))
    if (length(inside) > 0L) {
      j <- r[inside] + 1L
      out[inside] <- pmf[j] / surv[j]
    }

    # For r >= K, use tail_hazard
    tail <- which(r >= K)
    if (length(tail) > 0L) out[tail] <- tail_hazard

    # For negative run lengths (shouldn't happen), return NA
    neg <- which(r < 0L)
    if (length(neg) > 0L) out[neg] <- NA_real_

    pmin(pmax(out, 0), 1)
  }
}
