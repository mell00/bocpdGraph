#' Graph IID Gaussian model (nodewise, no coupling)
#'
#' Minimal runnable graph model for examples and vignettes. Runs the scalar
#' Gaussian IID model independently at each node.
#'
#' @param mu0 Prior mean.
#' @param sigma Known standard deviation.
#' @return An object of class \code{"graph_gaussian_iid"}.
#' @export
graph_gaussian_iid_model <- function(mu0 = 0, sigma = 1) {
  base <- gaussian_iid_model(mu0 = mu0, sigma = sigma)
  structure(
    list(base_model = base, state = NULL, n_nodes = NULL, max_run_length = NULL),
    class = "graph_gaussian_iid"
  )
}

#' @export
graph_model_init.graph_gaussian_iid <- function(model, graph, data, control = list(), ...) {
  n_nodes <- nrow(data)
  base <- model$base_model

  max_run_length <- if (!is.null(control$r_max)) as.integer(control$r_max) else 200L
  if (!is.finite(max_run_length) || max_run_length < 1L) {
    stop("control$r_max must be a finite integer >= 1 for graph_gaussian_iid_model().", call. = FALSE)
  }

  prior_nu  <- base$prior$nu
  prior_chi <- base$prior$chi

  # preallocate each node state to length (r_max + 1) at the prior
  states <- vector("list", n_nodes)
  for (i in seq_len(n_nodes)) {
    states[[i]] <- list(
      nu  = rep(prior_nu, max_run_length + 1L),
      chi = replicate(max_run_length + 1L, prior_chi, simplify = FALSE)
    )
  }

  model$state <- states
  model$n_nodes <- n_nodes
  model$max_run_length <- max_run_length
  model
}

#' @export
graph_pred_loglik.graph_gaussian_iid <- function(model, y_t, node, run_length,
                                                 neighbors = integer(), t = NULL, ...) {
  st <- model$state[[node]]
  base <- model$base_model

  # only evaluate up to what we have, then pad with prior
  L <- length(run_length)
  have <- min(L, length(st$nu), length(st$chi))

  out <- numeric(L)

  if (have > 0L) {
    out[seq_len(have)] <- vapply(seq_len(have), function(k) {
      base$pred_loglik(nu = st$nu[k], chi = st$chi[[k]], x_t = y_t)
    }, numeric(1))
  }

  # if run_length asks for more than current state (shouldn't happen after init),
  # fall back to prior predictive for the remaining entries
  if (have < L) {
    prior <- base$prior
    out[(have + 1L):L] <- base$pred_loglik(nu = prior$nu, chi = prior$chi, x_t = y_t)
  }

  out
}

#' @export
graph_model_update.graph_gaussian_iid <- function(model, y_t, node, run_length,
                                                  neighbors = integer(), t = NULL, ...) {
  base <- model$base_model
  st <- model$state[[node]]
  r_max <- model$max_run_length

  prior_nu  <- base$prior$nu
  prior_chi <- base$prior$chi

  # update each hypothesis (k corresponds to run length r=k-1 in storage)
  chi_new <- lapply(seq_along(st$chi), function(k) {
    mu_k <- st$chi[[k]]$mu
    nu_k <- st$nu[k]
    list(mu = (mu_k * nu_k + y_t) / (nu_k + 1))
  })

  # shift run lengths forward and prepend prior for r=0
  nu_next  <- c(prior_nu, st$nu + 1)
  chi_next <- c(list(prior_chi), chi_new)

  # truncate to 0..r_max (length r_max+1)
  keep <- seq_len(r_max + 1L)
  st$nu  <- nu_next[keep]
  st$chi <- chi_next[keep]

  model$state[[node]] <- st
  model
}

