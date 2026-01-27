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
#' @param hazard A function \code{hazard(r)} returning the hazard probability for
#'   run lengths \code{r} (non-negative integers). For a constant hazard use
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
#' model <- gaussian_mean_model(mu0 = 0, kappa0 = 1, alpha0 = 1, beta0 = 1)
#' h <- hazard_constant(100)
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

  # normalize input to list for consistent indexing
  if (!is.list(data)) data <- as.list(data)

  # init run-length posterior at t=0: P(r0=0)=1
  log_rl <- 0  # log(1)
  state <- model_init(model)

  rl_list <- vector("list", length(data))

  # optional dense storage
  rl_mat <- NULL
  if (isTRUE(ctrl$return_rl_matrix)) {
    # r_t <= t, and also <= r_max
    Tn <- length(data)
    Rn <- min(Tn, if (is.finite(ctrl$r_max)) as.integer(ctrl$r_max) else Tn)
    rl_mat <- matrix(0, nrow = Tn, ncol = Rn + 1) # columns: r=0..Rn
  }

  for (t in seq_along(data)) {

    x_t <- data[[t]]

    # current support is r = 0..(length(log_rl)-1)
    r_prev <- seq_len(length(log_rl)) - 1L

    # predictive loglik for each hypothesis (aligned with r_prev)
    log_pred <- vapply(
      seq_along(state$nu),
      function(i) model$pred_loglik(state$nu[i], state$chi[[i]], x_t),
      numeric(1)
    )

    # hazard evaluated at r_prev+1 per paper convention (H(rt-1+1))
    hvals <- hazard(r_prev + 1L)
    if (any(!is.finite(hvals)) || any(hvals < 0) || any(hvals > 1)) {
      stop("`hazard` must return probabilities in [0,1].", call. = FALSE)
    }

    # message passing in log space
    # growth: r_t = r_prev + 1
    log_growth <- log_rl + log1p(-hvals) + log_pred

    # changepoint: r_t = 0 is logsumexp over all previous hypotheses
    log_cp <- .logsumexp(log_rl + log(hvals) + log_pred)

    # new joint over r_t: (0, 1..)
    log_rl_new <- c(log_cp, log_growth)

    # truncate by r_max (keep r=0..r_max)
    if (is.finite(ctrl$r_max)) {
      keep_len <- min(length(log_rl_new), as.integer(ctrl$r_max) + 1L)
      log_rl_new <- log_rl_new[seq_len(keep_len)]
    }

    # normalize
    log_Z <- .logsumexp(log_rl_new)
    log_rl_new <- log_rl_new - log_Z
    rl_new <- exp(log_rl_new)

    # optional pruning (after normalization)
    if (is.numeric(ctrl$prune_eps) && ctrl$prune_eps > 0) {
      keep <- rl_new > ctrl$prune_eps
      # always keep r=0
      keep[1L] <- TRUE
      rl_new <- rl_new[keep]
      rl_new <- rl_new / sum(rl_new)
      log_rl_new <- log(rl_new)
      # NOTE: state pruning is handled by updating state first, then pruning consistently below.
      prune_keep <- keep
    } else {
      prune_keep <- rep(TRUE, length(rl_new))
    }

    # update sufficient stats for next step (prepend prior for r=0, shift others)
    state_new <- model_update(model, state, x_t)

    # align state length to rl_new length, then apply pruning mask
    if (length(state_new$nu) < length(log_rl_new)) {
      stop("Model state shorter than run-length support after update.", call. = FALSE)
    }
    state_new$nu  <- state_new$nu[seq_len(length(log_rl_new))]
    state_new$chi <- state_new$chi[seq_len(length(log_rl_new))]
    state_new$nu  <- state_new$nu[prune_keep]
    state_new$chi <- state_new$chi[prune_keep]

    # save and roll forward
    rl_list[[t]] <- rl_new
    if (!is.null(rl_mat)) {
      Rn <- ncol(rl_mat) - 1L
      fill <- min(length(rl_new), Rn + 1L)
      rl_mat[t, seq_len(fill)] <- rl_new[seq_len(fill)]
    }

    log_rl <- log(rl_new)
    state <- state_new
  }

  out <- list(rl = rl_list, state = state)
  if (!is.null(rl_mat)) out$rl_matrix <- rl_mat
  out
}

# stable log-sum-exp for numeric vectors
.logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}
