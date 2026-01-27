# internal BOCPD engine (not exported)

.logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

.normalize_log_probs <- function(log_w) {
  log_Z <- .logsumexp(log_w)

  if (!is.finite(log_Z)) {
    p <- numeric(length(log_w))
    p[1L] <- 1
    return(list(p = p, log_p = log(p)))
  }

  log_p <- log_w - log_Z
  p <- exp(log_p)

  if (any(!is.finite(p))) {
    p <- numeric(length(log_w))
    p[1L] <- 1
    return(list(p = p, log_p = log(p)))
  }

  s <- sum(p)
  if (!is.finite(s) || s <= 0) {
    p <- numeric(length(log_w))
    p[1L] <- 1
    return(list(p = p, log_p = log(p)))
  }

  p <- p / s
  list(p = p, log_p = log(p))
}

# one BOCPD step (message passing + trunc/prune + state update)
bocpd_step <- function(log_rl, state, x_t, model, hazard, t,
                       r_max = Inf, prune_eps = 0) {

  # current support r = 0..(len-1)
  r_prev <- seq_len(length(log_rl)) - 1L

  # predictive loglik for each hypothesis (aligned with r_prev)
  log_pred <- vapply(
    seq_along(state$nu),
    function(i) model$pred_loglik(state$nu[i], state$chi[[i]], x_t),
    numeric(1)
  )

  # hazard evaluated at r_prev+1 per paper convention (H(rt-1+1))
  hvals <- hazard(run_length = r_prev + 1L, t = t)
  if (any(!is.finite(hvals)) || any(hvals < 0) || any(hvals > 1)) {
    stop("`hazard` must return probabilities in [0,1].", call. = FALSE)
  }

  # message passing in log space
  log_growth <- log_rl + log1p(-hvals) + log_pred
  log_cp <- .logsumexp(log_rl + log(hvals) + log_pred)

  log_rl_new <- c(log_cp, log_growth)

  # truncate by r_max (keep r=0..r_max)
  if (is.finite(r_max)) {
    keep_len <- min(length(log_rl_new), as.integer(r_max) + 1L)
    log_rl_new <- log_rl_new[seq_len(keep_len)]
  }

  norm <- .normalize_log_probs(log_rl_new)
  rl_new <- norm$p
  log_rl_new <- norm$log_p

  # optional pruning (after normalization)
  if (is.numeric(prune_eps) && prune_eps > 0) {
    keep <- rl_new > prune_eps
    keep[1L] <- TRUE

    rl_new <- rl_new[keep]
    rl_new <- rl_new / sum(rl_new)
    log_rl_new <- log(rl_new)

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

  list(
    rl = rl_new,
    log_rl = log_rl_new,
    state = state_new
  )
}
