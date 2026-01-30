#' BOCPD on Graph-Structured Multivariate Series
#'
#' Runs Bayesian online changepoint detection (BOCPD) for data observed on
#' nodes of a graph. The BOCPD recursion is applied per node; graph structure
#' enters only through the graph-aware model's predictive log-likelihood.
#'
#' @param data Node-by-time observations. Accepts:
#'   \itemize{
#'     \item a numeric matrix with rows = nodes and cols = time,
#'     \item a list of length n_nodes where each element is a numeric vector of length T.
#'   }
#' @param graph Graph structure. Accepts:
#'   \itemize{
#'     \item a square adjacency matrix (0/1 or weighted),
#'     \item an \code{igraph} object (if igraph is installed).
#'   }
#' @param model A graph-aware model object. Must be compatible with:
#'   \code{graph_model_init()}, \code{graph_pred_loglik()}, \code{graph_model_update()}.
#' @param hazard A hazard function used by \code{update_run_length()}.
#' @param control Optional list of controls. Uses \code{r_max} (default 200),
#'   \code{store_R} (default TRUE), \code{verbose} (default FALSE).
#' @param ... Additional arguments forwarded to the hazard via \code{hazard_args}.
#'
#' @return An object of class \code{"bocpdGraph"}.
#' @export
bocpdGraph <- function(data, graph, model, hazard, control = list(), ...) {
  control <- .bocpdGraph_control_default(control)

  Y <- .bocpdGraph_coerce_data(data)
  n_nodes <- nrow(Y)
  T <- ncol(Y)

  G <- .bocpdGraph_coerce_graph(graph, n_nodes)
  neighbors <- .bocpdGraph_neighbors(G)

  max_run_length <- as.integer(control$r_max)

  # Initialize the graph model (graph + data + controls)
  model <- graph_model_init(model = model, graph = G, data = Y, control = control)

  # Initialize per-node BOCPD run-length state
  states <- vector("list", n_nodes)
  for (i in seq_len(n_nodes)) {
    states[[i]] <- .bocpdGraph_init_node_state(max_run_length)
  }

  # Storage (optional)
  store_R <- isTRUE(control$store_R)
  R_store <- NULL
  if (store_R) {
    R_store <- array(NA_real_, dim = c(n_nodes, T, max_run_length + 1L))
  }

  hazard_args <- list(...)

  # Main recursion
  for (t in seq_len(T)) {
    if (isTRUE(control$verbose) && (t %% 50L == 0L)) {
      message("bocpdGraph: t = ", t, "/", T)
    }

    for (i in seq_len(n_nodes)) {
      y_it <- Y[i, t]

      # Optional NA behavior: carry forward unchanged
      if (is.na(y_it)) {
        if (store_R) R_store[i, t, ] <- states[[i]]$R
        next
      }

      r <- states[[i]]$run_length

      # Predictive log-likelihood vector over run lengths (length = max_run_length + 1)
      loglik <- graph_pred_loglik(
        model = model,
        y_t = y_it,
        node = i,
        run_length = r,
        neighbors = neighbors[[i]],
        t = t
      )

      if (!is.numeric(loglik) || length(loglik) != (max_run_length + 1L)) {
        stop("graph_pred_loglik() must return numeric vector of length max_run_length + 1.")
      }

      # Prior predictive for changepoint (r -> 0). Minimal consistent choice:
      # use the r=0 predictive already returned.
      loglik_prior <- loglik[1L]

      res <- update_run_length(
        log_R_prev = states[[i]]$log_R,
        loglik = loglik,
        loglik_prior = loglik_prior,
        hazard = hazard,
        t = t,
        max_run_length = max_run_length,
        hazard_args = hazard_args
      )

      if (is.null(res$log_R)) stop("update_run_length() must return res$log_R.")
      if (!is.numeric(res$log_R) || length(res$log_R) != (max_run_length + 1L)) {
        stop("res$log_R must be numeric of length max_run_length + 1.")
      }

      states[[i]]$log_R <- res$log_R

      # Keep probability-scale copy for storage/convenience (stable normalization)
      lr <- states[[i]]$log_R
      lr <- lr - max(lr)
      R <- exp(lr)
      states[[i]]$R <- R / sum(R)

      # Model update (graph-aware)
      model <- graph_model_update(
        model = model,
        y_t = y_it,
        node = i,
        run_length = r,
        neighbors = neighbors[[i]],
        t = t
      )

      if (store_R) {
        R_store[i, t, ] <- states[[i]]$R
      }
    }
  }

  out <- list(
    R = if (store_R) R_store else lapply(states, `[[`, "R"),
    model = model,
    graph = G,
    neighbors = neighbors,
    control = control,
    call = match.call()
  )
  class(out) <- "bocpdGraph"
  out
}

# ---- helpers (internal) ----

.bocpdGraph_control_default <- function(control) {
  if (is.null(control) || !is.list(control)) control <- list()
  if (is.null(control$r_max)) control$r_max <- 200L
  if (is.null(control$store_R)) control$store_R <- TRUE
  if (is.null(control$verbose)) control$verbose <- FALSE
  control$r_max <- as.integer(control$r_max)
  if (control$r_max < 1L) stop("control$r_max must be >= 1.")
  control
}

.bocpdGraph_coerce_data <- function(data) {
  if (is.matrix(data)) {
    if (!is.numeric(data)) stop("If data is a matrix it must be numeric.")
    if (nrow(data) < 1L || ncol(data) < 1L) stop("data matrix must be non-empty.")
    return(data)
  }

  if (is.list(data)) {
    if (length(data) < 1L) stop("data list must be non-empty.")
    lens <- vapply(data, length, integer(1))
    if (length(unique(lens)) != 1L) stop("All elements of data list must have the same length.")
    T <- lens[1]
    n_nodes <- length(data)
    Y <- matrix(NA_real_, nrow = n_nodes, ncol = T)
    for (i in seq_len(n_nodes)) {
      if (!is.numeric(data[[i]])) stop("Each element of data list must be numeric.")
      Y[i, ] <- data[[i]]
    }
    return(Y)
  }

  stop("data must be a numeric matrix (nodes x time) or a list of numeric vectors.")
}

.bocpdGraph_coerce_graph <- function(graph, n_nodes) {
  if (inherits(graph, "igraph")) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
      stop("igraph object provided but igraph is not installed.")
    }
    A <- igraph::as_adjacency_matrix(graph, sparse = FALSE)
    graph <- as.matrix(A)
  }

  if (!is.matrix(graph)) stop("graph must be an adjacency matrix or an igraph object.")
  if (nrow(graph) != ncol(graph)) stop("graph adjacency matrix must be square.")
  if (nrow(graph) != n_nodes) stop("graph size must match number of nodes in data.")
  if (!is.numeric(graph)) stop("graph adjacency matrix must be numeric (0/1 or weights).")

  graph
}

.bocpdGraph_neighbors <- function(A) {
  n <- nrow(A)
  neigh <- vector("list", n)
  for (i in seq_len(n)) {
    idx <- which(A[i, ] != 0)
    idx <- idx[idx != i]
    neigh[[i]] <- idx
  }
  neigh
}

.bocpdGraph_init_node_state <- function(max_run_length) {
  max_run_length <- as.integer(max_run_length)
  run_length <- 0:max_run_length

  # Minimal safe default: start with P(r=0)=1
  R <- c(1, rep(0, max_run_length))
  log_R <- log(pmax(R, .Machine$double.xmin))

  list(
    max_run_length = max_run_length,
    run_length = run_length,
    R = R,
    log_R = log_R
  )
}
