#' Bayesian Online Changepoint Detection (BOCPD)
#'
#' Implements the Bayesian online changepoint detection recursion
#' (run-length message passing) described by Adams and MacKay (2007).
#' The algorithm maintains a posterior distribution over the current run length
#' and updates it online using a hazard function and a predictive model.
#'
#' @param data A vector or list of observations in time order. If a vector is
#'   supplied it will be treated as a sequence of scalar observations.
#' @param model A model object with methods \code{model_init()} and
#'   \code{model_update()}, and internal components \code{pred_loglik} and
#'   \code{prior}. See Details.
#' @param hazard A function \code{hazard(run_length, t, ...)} returning the hazard
#'   probabilities for run lengths \code{run_length}. For a constant hazard use
#'   \code{hazard_constant()}.
#' @param control Optional list of controls:
#'   \itemize{
#'     \item \code{r_max} integer; maximum run length retained (default: Inf).
#'     \item \code{prune_eps} numeric; drop tail mass below this threshold after
#'           normalization (default: 0).
#'     \item \code{return_rl_matrix} logical; if TRUE, returns a dense matrix of
#'           run-length posteriors (default: FALSE).
#'   }
#'
#' @details
#' The recursion uses only two transitions per time step: growth (\eqn{r \to r+1})
#' and changepoint (\eqn{r \to 0}). The predictive log-likelihood is evaluated
#' for each run-length hypothesis before posterior sufficient statistics are updated.
#'
#' Model requirements:
#' \itemize{
#'   \item \code{model_init(model)} returns a state list with fields \code{nu} (numeric)
#'         and \code{chi} (list), each aligned with the current run-length support.
#'   \item \code{model_update(model, state, x_t)} returns the *next* state, with
#'         run-length-shifted sufficient statistics (including prepending the prior).
#'   \item \code{model$pred_loglik(nu, chi, x_t)} returns log predictive density for one hypothesis.
#'   \item \code{model$prior} contains \code{nu} and \code{chi} to prepend on reset.
#' }
#'
#' @return If \code{return_rl_matrix = FALSE}, a list with:
#' \itemize{
#'   \item \code{rl} list of numeric vectors; \code{rl[[t]]} is \eqn{P(r_t | x_{1:t})}.
#'   \item \code{state} final model state (useful for streaming).
#' }
#' If \code{return_rl_matrix = TRUE}, also includes \code{rl_matrix}, a dense matrix
#' with rows indexed by time and columns by run length (truncated to \code{r_max}).
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(50)
#' model <- gaussian_iid_model(mu0 = 0, sigma = 1)
#' h <- hazard_constant(1/100)
#' fit <- bocpd(x, model, h, control = list(r_max = 200, prune_eps = 1e-12))
#' length(fit$rl)
#'
#' @export
bocpd <- function(data, model, hazard, control = list()) {

  if (!is.function(hazard)) stop("`hazard` must be a function.", call. = FALSE)

  ctrl <- modifyList(
    list(r_max = Inf, prune_eps = 0, return_rl_matrix = FALSE),
    control
  )

  if (!is.list(data)) data <- as.list(data)

  # init run-length posterior at t=0: P(r0=0)=1
  log_rl <- 0
  state <- model_init(model)

  rl_list <- vector("list", length(data))

  rl_mat <- NULL
  if (isTRUE(ctrl$return_rl_matrix)) {
    Tn <- length(data)
    Rn <- min(Tn, if (is.finite(ctrl$r_max)) as.integer(ctrl$r_max) else Tn)
    rl_mat <- matrix(0, nrow = Tn, ncol = Rn + 1)
  }

  for (t in seq_along(data)) {
    step <- bocpd_step(
      log_rl = log_rl,
      state  = state,
      x_t    = data[[t]],
      model  = model,
      hazard = hazard,
      t      = t,
      r_max  = ctrl$r_max,
      prune_eps = ctrl$prune_eps
    )

    rl_list[[t]] <- step$rl

    if (!is.null(rl_mat)) {
      Rn <- ncol(rl_mat) - 1L
      fill <- min(length(step$rl), Rn + 1L)
      rl_mat[t, seq_len(fill)] <- step$rl[seq_len(fill)]
    }

    log_rl <- step$log_rl
    state  <- step$state
  }

  out <- list(rl = rl_list, state = state)
  if (!is.null(rl_mat)) out$rl_matrix <- rl_mat
  out
}

#' Estimate changepoint locations from BOCPD output
#'
#' Provides simple changepoint extraction heuristics from the run-length
#' posterior returned by \code{\link{bocpd}}.
#'
#' The changepoint selection heuristics implemented here are inspired by
#' the post-processing procedures described in Bretz (2021), which proposes
#' identifying changepoints either via elevated posterior mass at
#' \eqn{r_t = 0} or via sharp drops in the maximum a posteriori (MAP)
#' run length. These heuristics operate after the BOCPD recursion and do
#' not alter the underlying Bayesian filtering algorithm of
#' Adams and MacKay (2007).
#'
#' Two methods are available:
#' \itemize{
#'   \item \code{"cp_prob"}: selects times where \eqn{P(r_t=0 \mid x_{1:t})}
#'         exceeds \code{threshold}.
#'   \item \code{"map_drop"}: selects times where the MAP run length decreases
#'         by at least \code{threshold} compared to the previous time step.
#' }
#'
#' @param fit Output from \code{\link{bocpd}} (a list with element \code{rl}),
#'   or a list with element \code{rl_matrix}.
#' @param method Character; one of \code{"cp_prob"} or \code{"map_drop"}.
#' @param threshold Numeric threshold. Interpretation depends on \code{method}:
#'   \itemize{
#'     \item \code{"cp_prob"}: probability in \eqn{[0,1]} (default \code{0.5}).
#'     \item \code{"map_drop"}: non-negative drop in MAP run length (default \code{10}).
#'   }
#' @param min_sep Integer; minimum separation (in time points) between returned
#'   changepoints. When multiple candidates occur within \code{min_sep}, the
#'   strongest one (largest score) is kept. Default is \code{1} (no filtering).
#'
#' @return Integer vector of time indices (1-based) flagged as changepoints.
#'
#' @examples
#' set.seed(1)
#' x <- c(rnorm(50, 0), rnorm(50, 3))
#' model <- gaussian_iid_model(mu0 = 0, sigma = 1)
#' h <- hazard_constant(1/100)
#' fit <- bocpd(x, model, h, control = list(r_max = 200, prune_eps = 1e-12))
#' bocpd_changepoints(fit, method = "cp_prob", threshold = 0.3)
#'
#' @export
bocpd_changepoints <- function(fit,
                               method = c("cp_prob", "map_drop"),
                               threshold = NULL,
                               min_sep = 1L) {

  method <- match.arg(method)
  min_sep <- as.integer(min_sep)
  if (min_sep < 1L) stop("`min_sep` must be >= 1.", call. = FALSE)

  if (!is.null(fit$rl)) {
    rl_list <- fit$rl
  } else if (!is.null(fit$rl_matrix)) {
    rl_mat <- fit$rl_matrix
    rl_list <- lapply(seq_len(nrow(rl_mat)), function(i) {
      v <- rl_mat[i, ]
      k <- max(which(v > 0), 1L)
      v[seq_len(k)]
    })
  } else {
    stop("`fit` must contain `rl` or `rl_matrix`.", call. = FALSE)
  }

  Tn <- length(rl_list)
  if (Tn < 1L) return(integer(0))

  if (method == "cp_prob") {

    if (is.null(threshold)) threshold <- 0.5
    if (!is.numeric(threshold) || length(threshold) != 1L ||
        threshold < 0 || threshold > 1) {
      stop("For method='cp_prob', `threshold` must be a single number in [0,1].",
           call. = FALSE)
    }

    score <- vapply(rl_list, function(v) v[1L], numeric(1))
    idx <- which(score >= threshold)

  } else {

    if (is.null(threshold)) threshold <- 10
    if (!is.numeric(threshold) || length(threshold) != 1L || threshold < 0) {
      stop("For method='map_drop', `threshold` must be a single non-negative number.",
           call. = FALSE)
    }

    map_rl <- vapply(rl_list, function(v) as.integer(which.max(v) - 1L), integer(1))
    score <- c(0, pmax(0, map_rl[-Tn] - map_rl[-1L]))  # previous - current when drop
    idx <- which(score >= threshold)
  }

  if (length(idx) == 0L) return(integer(0))

  if (min_sep > 1L) {
    ord <- idx[order(-score[idx], idx)]
    keep <- logical(Tn)
    out <- integer(0)
    for (i in ord) {
      lo <- max(1L, i - min_sep + 1L)
      hi <- min(Tn, i + min_sep - 1L)
      if (!any(keep[lo:hi])) {
        keep[i] <- TRUE
        out <- c(out, i)
      }
    }
    idx <- sort(out)
  }

  idx
}

