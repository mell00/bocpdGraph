#' Run-length recursion utilities (BOCPD)
#'
#' Internal helpers for maintaining and normalizing the run-length posterior
#' in log-space. These are graph-agnostic and do not depend on any likelihood.
#'
#' @name run_length
NULL

#' Stable log-sum-exp
#'
#' @param x Numeric vector (may contain -Inf).
#'
#' @return A single number: log(sum(exp(x))) computed stably.
#'
#' @keywords internal
log_sum_exp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}




#' Normalize log-probabilities
#'
#' @param logp Numeric vector of unnormalized log-probabilities.
#'
#' @return Numeric vector of normalized log-probabilities (log-scale),
#'   such that sum(exp(out)) == 1 (up to floating point error).
#'
#' @keywords internal
normalize_log_probs <- function(logp) {
  logp - log_sum_exp(logp)
}




#' Coerce hazard specification to a function
#'
#' Accepts either a function hazard or a scalar constant hazard probability.
#'
#' @param hazard Either:
#' \itemize{
#'   \item a function \code{hazard(run_length, t, ...)} returning probabilities in [0,1], or
#'   \item a numeric scalar in [0,1] indicating a constant hazard.
#' }
#'
#' @return A function \code{f(run_length, t, ...)}.
#'
#' @keywords internal
as_hazard_fun <- function(hazard) {
  if (is.function(hazard)) return(hazard)

  if (is.numeric(hazard) && length(hazard) == 1L && is.finite(hazard)) {
    h0 <- as.numeric(hazard)
    if (h0 < 0 || h0 > 1) {
      stop("Constant hazard must be in [0,1].", call. = FALSE)
    }
    return(function(run_length, t, ...) rep(h0, length(run_length)))
  }

  stop("`hazard` must be a function or a finite numeric scalar in [0,1].", call. = FALSE)
}

#' Validate and truncate run-length support
#'
#' @param log_R Log run-length probabilities for r = 0..Rmax (numeric vector).
#' @param stats List of sufficient statistics aligned with \code{log_R}.
#' @param max_run_length Nonnegative integer truncation level.
#'
#' @return A list with truncated \code{log_R} and \code{stats}.
#'
#' @keywords internal
truncate_run_length <- function(log_R, stats, max_run_length) {
  if (!is.numeric(log_R) || length(log_R) < 1L) {
    stop("`log_R` must be a non-empty numeric vector.", call. = FALSE)
  }
  if (!is.list(stats) || length(stats) != length(log_R)) {
    stop("`stats` must be a list with the same length as `log_R`.", call. = FALSE)
  }
  if (!is.numeric(max_run_length) || length(max_run_length) != 1L || max_run_length < 0) {
    stop("`max_run_length` must be a single nonnegative integer.", call. = FALSE)
  }
  max_run_length <- as.integer(max_run_length)

  keep <- min(length(log_R), max_run_length + 1L)
  list(
    log_R = as.numeric(log_R)[seq_len(keep)],
    stats = stats[seq_len(keep)]
  )
}




#' One-step run-length update given log-likelihood and hazard
#'
#' Performs BOCPD run-length recursion update *only*:
#' \itemize{
#'   \item growth: r -> r+1 with probability (1-h)
#'   \item changepoint: r -> 0 with probability h
#' }
#'
#' NOTE: user supplies likelihoods.
#'
#' @param log_R_prev Log posterior for run-lengths r=0..Rmax at time t-1.
#' @param loglik Numeric vector of length \code{length(log_R_prev)} giving
#'   log predictive densities for the new observation y_t under each run-length.
#' @param loglik_prior Single number: log predictive density under the prior/reset
#'   statistics used for changepoints.
#' @param hazard Hazard specification (function or scalar in [0,1]).
#' @param t Integer time index after the update (i.e., current time).
#' @param max_run_length Truncation level for run-length support.
#' @param hazard_args Named list forwarded to hazard if hazard is a function.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{log_R}: normalized log run-length posterior at time t
#'   \item \code{R}: same on probability scale
#' }
#'
#' @keywords internal
update_run_length <- function(
    log_R_prev,
    loglik,
    loglik_prior,
    hazard,
    t,
    max_run_length,
    hazard_args = list()
) {
  if (!is.numeric(log_R_prev) || length(log_R_prev) < 1L) {
    stop("`log_R_prev` must be a non-empty numeric vector.", call. = FALSE)
  }
  if (!is.numeric(loglik) || length(loglik) != length(log_R_prev)) {
    stop("`loglik` must be numeric with the same length as `log_R_prev`.", call. = FALSE)
  }
  if (!is.numeric(loglik_prior) || length(loglik_prior) != 1L || !is.finite(loglik_prior)) {
    stop("`loglik_prior` must be a single finite number.", call. = FALSE)
  }
  if (!is.numeric(t) || length(t) != 1L || t < 1) {
    stop("`t` must be a single integer >= 1.", call. = FALSE)
  }
  t <- as.integer(t)

  if (!is.numeric(max_run_length) || length(max_run_length) != 1L || max_run_length < 0) {
    stop("`max_run_length` must be a single nonnegative integer.", call. = FALSE)
  }
  max_run_length <- as.integer(max_run_length)

  hz_fun <- as_hazard_fun(hazard)

  # Truncate if needed
  Rmax_prev <- length(log_R_prev) - 1L
  if (Rmax_prev > max_run_length) {
    log_R_prev <- log_R_prev[seq_len(max_run_length + 1L)]
    loglik <- loglik[seq_len(max_run_length + 1L)]
    Rmax_prev <- max_run_length
  }

  r_vec <- 0L:Rmax_prev
  h <- do.call(hz_fun, c(list(run_length = r_vec, t = t), hazard_args))
  h <- as.numeric(h)
  if (length(h) != length(r_vec)) stop("Hazard function returned wrong length.", call. = FALSE)
  if (any(!is.finite(h)) || any(h < 0) || any(h > 1)) stop("Hazard must be in [0,1].", call. = FALSE)

  # Allocate new support r=0..max_run_length
  log_R_new <- rep(-Inf, max_run_length + 1L)

  # Growth: r -> r+1
  grow_max_r <- min(Rmax_prev, max_run_length - 1L)
  if (grow_max_r >= 0L) {
    idx_r <- 0L:grow_max_r
    log_R_new[idx_r + 2L] <- log_R_prev[idx_r + 1L] + log1p(-h[idx_r + 1L]) + loglik[idx_r + 1L]
  }

  # for changepoint, r -> 0 uses prior predictive
  log_cp_terms <- log_R_prev + log(h) + loglik_prior
  log_R_new[1L] <- log_sum_exp(log_cp_terms)

  # Normalize
  log_R_new <- normalize_log_probs(log_R_new)

  list(
    log_R = log_R_new,
    R = exp(log_R_new)
  )
}
