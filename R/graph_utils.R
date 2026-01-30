#' Graph utilities (internal)
#'
#' Internal helpers for validating and normalizing graph inputs used by
#' \code{\link{bocpdGraph}}.
#'
#' These functions avoid hard dependencies on graph packages. If an \pkg{igraph}
#' object is supplied, \pkg{igraph} must be available (typically in Suggests).
#' Sparse matrices from \pkg{Matrix} are supported without coercing to dense.
#'
#' @name graph_utils
#' @keywords internal
NULL

.is_sparse_matrix <- function(x) {
  inherits(x, "Matrix")
}

.matrix_ns <- function() {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix package required for sparse matrix operations.")
  }
  asNamespace("Matrix")
}

.igraph_ns <- function() {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("`graph` is an igraph object but the igraph package is not available.")
  }
  asNamespace("igraph")
}

# codetools-safe setter for Matrix::diag<-
.Matrix_diag_set <- function(A, value) {
  M <- .matrix_ns()
  diag_setter <- get("diag<-", envir = M)
  diag_setter(A, value = value)
}

#' Normalize a graph input to an adjacency matrix (dense or sparse)
#'
#' Accepts:
#' \itemize{
#'   \item numeric matrix: treated as adjacency/weight matrix
#'   \item \pkg{Matrix} matrix: sparse/dense Matrix classes
#'   \item list: adjacency list where \code{graph[[i]]} are neighbor indices
#'   \item \pkg{igraph} object (if available)
#' }
#'
#' Returns a numeric adjacency matrix with zeros on the diagonal.
#' For \pkg{Matrix} inputs, returns a \pkg{Matrix} object (typically sparse).
#'
#' @param graph Graph representation (matrix, Matrix, adjacency list, or igraph)
#' @param n_nodes Optional integer number of nodes (needed for adjacency lists)
#' @param allow_weights Logical; if FALSE, adjacency is coerced to \eqn{\{0,1\}}.
#' @param directed Logical; if FALSE, enforces symmetry by max(A, t(A))
#' @return Adjacency matrix (base matrix or Matrix).
#' @importFrom methods as
#' @keywords internal
as_adjacency_matrix <- function(graph,
                                n_nodes = NULL,
                                allow_weights = TRUE,
                                directed = FALSE) {
  if (is.data.frame(graph)) {
    stop("Unsupported `graph` type: data.frame. Provide a matrix, Matrix, adjacency list, or igraph object.")
  }

  # Base matrix
  if (is.matrix(graph) && !inherits(graph, "Matrix")) {
    A <- graph
    if (!is.numeric(A)) {
      A <- suppressWarnings(matrix(as.numeric(A), nrow = nrow(A), ncol = ncol(A)))
    }
    if (nrow(A) != ncol(A)) stop("`graph` adjacency matrix must be square.")
    A[is.na(A)] <- 0
    if (any(!is.finite(A))) stop("Adjacency matrix contains non-finite values.")
    if (any(A < 0)) stop("Adjacency/weight matrix must be nonnegative.")
    diag(A) <- 0
    if (!directed) A <- pmax(A, t(A))
    if (!allow_weights) A <- (A != 0) * 1
    return(A)
  }

  # Matrix package classes (sparse or dense)
  if (.is_sparse_matrix(graph)) {
    M <- .matrix_ns()
    A <- graph
    if (nrow(A) != ncol(A)) stop("`graph` adjacency matrix must be square.")

    # Validate numeric + finite + nonnegative using summary (works broadly across Matrix classes)
    s <- M$summary(A) # i, j, x for nonzeros
    if (nrow(s) && anyNA(s$x)) {
      stop("Adjacency matrix contains NA values in nonzero entries.")
    }
    if (nrow(s) && any(!is.finite(s$x))) stop("Adjacency matrix contains non-finite values.")
    if (nrow(s) && any(s$x < 0)) stop("Adjacency/weight matrix must be nonnegative.")

    # No self-loops (codetools-safe)
    A <- .Matrix_diag_set(A, 0)

    # Symmetrize if undirected without densifying
    if (!directed) {
      A <- .sparse_pmax(A, M$t(A))
      A <- .Matrix_diag_set(A, 0)
    }

    if (!allow_weights) {
      # Coerce all nonzeros to 1 while preserving sparsity
      s <- M$summary(A)
      if (nrow(s)) {
        A <- M$sparseMatrix(
          i = s$i, j = s$j,
          x = rep(1, nrow(s)),
          dims = dim(A),
          giveCsparse = TRUE
        )
        A <- .Matrix_diag_set(A, 0)
      }
      A <- M$drop0(A)
    }

    return(A)
  }

  # Adjacency list
  if (is.list(graph) && !inherits(graph, "igraph")) {
    if (is.null(n_nodes)) {
      n_nodes <- length(graph)
      if (!is.numeric(n_nodes) || n_nodes <= 0) stop("Could not infer `n_nodes` from adjacency list.")
    }
    A <- matrix(0, nrow = n_nodes, ncol = n_nodes)
    for (i in seq_len(n_nodes)) {
      nbrs <- graph[[i]]
      if (is.logical(nbrs)) stop("Adjacency list entries must be integer indices, not logical.")
      if (length(nbrs) == 0) next
      if (!is.numeric(nbrs) && !is.integer(nbrs)) stop("Adjacency list entries must be integer indices.")
      nbrs <- as.integer(nbrs)
      nbrs <- nbrs[!is.na(nbrs)]
      if (any(nbrs < 1L | nbrs > n_nodes)) stop("Adjacency list contains out-of-range node indices.")
      A[i, nbrs] <- 1
    }
    diag(A) <- 0
    if (!directed) A <- pmax(A, t(A))
    if (!allow_weights) A <- (A != 0) * 1
    return(A)
  }

  # igraph
  if (inherits(graph, "igraph")) {
    .igraph_ns() # validates availability
    M <- .matrix_ns()

    # Avoid igraph:: usage: call exported function via namespace
    as_adj <- getExportedValue("igraph", "as_adjacency_matrix")
    A <- as_adj(graph, sparse = TRUE)

    # Ensure a Csparse numeric Matrix class if possible
    # Use base as() (no methods::) to avoid CRAN dependency warnings.
    A <- as(A, "dgCMatrix")

    A <- .Matrix_diag_set(A, 0)
    if (!directed) {
      A <- .sparse_pmax(A, M$t(A))
      A <- .Matrix_diag_set(A, 0)
    }
    if (!allow_weights) {
      s <- M$summary(A)
      if (nrow(s)) {
        A <- M$sparseMatrix(
          i = s$i, j = s$j,
          x = rep(1, nrow(s)),
          dims = dim(A),
          giveCsparse = TRUE
        )
        A <- .Matrix_diag_set(A, 0)
      }
      A <- M$drop0(A)
    }
    return(A)
  }

  stop("Unsupported `graph` type: provide a matrix, Matrix, adjacency list, or igraph object.")
}

# Sparse elementwise max(A, B) without densifying
.sparse_pmax <- function(A, B) {
  M <- .matrix_ns()
  if (nrow(A) != nrow(B) || ncol(A) != ncol(B)) stop("Dimension mismatch in .sparse_pmax().")

  sa <- M$summary(A)
  sb <- M$summary(B)

  # Merge keys (i,j) taking max(x)
  key_a <- paste(sa$i, sa$j, sep = ":")
  key_b <- paste(sb$i, sb$j, sep = ":")

  vals <- sa$x
  names(vals) <- key_a

  for (k in seq_along(sb$x)) {
    kk <- key_b[k]
    xv <- sb$x[k]
    if (!is.finite(xv)) stop("Non-finite value encountered in sparse matrix.")
    if (kk %in% names(vals)) {
      vals[kk] <- max(vals[kk], xv)
    } else {
      vals[kk] <- xv
    }
  }

  parts <- strsplit(names(vals), ":", fixed = TRUE)
  ii <- as.integer(vapply(parts, `[[`, character(1), 1L))
  jj <- as.integer(vapply(parts, `[[`, character(1), 2L))
  xx <- as.numeric(unname(vals))

  M$drop0(M$sparseMatrix(i = ii, j = jj, x = xx, dims = dim(A), giveCsparse = TRUE))
}

#' Validate adjacency matrix invariants
#'
#' @param A Numeric square adjacency matrix (base or Matrix).
#' @param directed Logical; if FALSE, requires symmetry within tolerance.
#' @param tol Numeric tolerance for symmetry check.
#' @return Invisibly returns TRUE on success.
#' @keywords internal
validate_adjacency <- function(A, directed = FALSE, tol = 0) {
  if (!is.matrix(A) && !inherits(A, "Matrix")) stop("Adjacency must be a matrix (base or Matrix).")
  if (nrow(A) != ncol(A)) stop("Adjacency must be square.")

  if (inherits(A, "Matrix")) {
    M <- .matrix_ns()
    s <- M$summary(A)
    if (nrow(s) && any(!is.finite(s$x))) stop("Adjacency matrix contains non-finite values.")
    if (nrow(s) && any(s$x < 0)) stop("Adjacency/weight matrix must be nonnegative.")
    if (any(M$diag(A) != 0)) stop("Adjacency matrix must have zero diagonal (no self-loops).")

    if (!directed) {
      D <- A - M$t(A)
      sd <- M$summary(D)
      if (nrow(sd) && max(abs(sd$x)) > tol) stop("Adjacency matrix must be symmetric for an undirected graph.")
    }
    invisible(TRUE)
  } else {
    if (!is.numeric(A)) stop("Adjacency matrix must be numeric.")
    if (any(!is.finite(A))) stop("Adjacency matrix contains non-finite values.")
    if (any(A < 0)) stop("Adjacency/weight matrix must be nonnegative.")
    if (any(diag(A) != 0)) stop("Adjacency matrix must have zero diagonal (no self-loops).")
    if (!directed) {
      if (max(abs(A - t(A))) > tol) stop("Adjacency matrix must be symmetric for an undirected graph.")
    }
    invisible(TRUE)
  }
}

#' Convert adjacency matrix to neighbor index lists
#'
#' Neighbors are returned as integer vectors. If \code{include_weights = TRUE},
#' returns a list with components \code{neighbors} and \code{weights}, where
#' \code{weights[[i]]} aligns with \code{neighbors[[i]]}.
#'
#' Works for base matrices and sparse \pkg{Matrix} matrices without densifying.
#'
#' @param A Numeric adjacency/weight matrix (base or Matrix).
#' @param include_weights Logical; include edge weights per node.
#' @return Either a list of integer neighbor vectors, or a list with
#'   \code{neighbors} and \code{weights}.
#' @keywords internal
adjacency_to_neighbors <- function(A, include_weights = TRUE) {
  validate_adjacency(A, directed = TRUE)

  n <- nrow(A)
  nbrs <- vector("list", n)
  wts  <- if (include_weights) vector("list", n) else NULL

  if (inherits(A, "Matrix")) {
    M <- .matrix_ns()
    s <- M$summary(A)
    keep <- s$i != s$j
    s <- s[keep, , drop = FALSE]

    if (nrow(s) == 0) {
      for (i in seq_len(n)) {
        nbrs[[i]] <- integer(0)
        if (include_weights) wts[[i]] <- numeric(0)
      }
      return(if (include_weights) list(neighbors = nbrs, weights = wts) else nbrs)
    }

    rows <- split(seq_len(nrow(s)), s$i)
    for (i in seq_len(n)) {
      idx <- rows[[as.character(i)]]
      if (is.null(idx)) {
        nbrs[[i]] <- integer(0)
        if (include_weights) wts[[i]] <- numeric(0)
      } else {
        nbrs[[i]] <- as.integer(s$j[idx])
        if (include_weights) wts[[i]] <- as.numeric(s$x[idx])
      }
    }
  } else {
    for (i in seq_len(n)) {
      j <- which(A[i, ] != 0)
      j <- j[j != i]
      nbrs[[i]] <- as.integer(j)
      if (include_weights) wts[[i]] <- as.numeric(A[i, j])
    }
  }

  if (include_weights) list(neighbors = nbrs, weights = wts) else nbrs
}

#' Standardize graph input to a ready-to-use structure
#'
#' Returns a list with adjacency \code{A}, neighbor lists \code{nbrs},
#' weight lists \code{wts}, and degrees (weighted and unweighted).
#'
#' @param graph Graph representation.
#' @param n_nodes Optional number of nodes (for adjacency lists).
#' @param directed Logical; if FALSE, symmetrize.
#' @param allow_weights Logical; if FALSE, treat any nonzero as 1.
#' @return List with components: \code{A}, \code{nbrs}, \code{wts},
#'   \code{deg}, \code{wdeg}, \code{n_nodes}.
#' @keywords internal
standardize_graph <- function(graph,
                              n_nodes = NULL,
                              directed = FALSE,
                              allow_weights = TRUE) {
  A <- as_adjacency_matrix(
    graph = graph,
    n_nodes = n_nodes,
    allow_weights = allow_weights,
    directed = directed
  )

  validate_adjacency(A, directed = directed, tol = 0)

  nn <- adjacency_to_neighbors(A, include_weights = TRUE)
  nbrs <- nn$neighbors
  wts  <- nn$weights

  deg  <- vapply(nbrs, length, integer(1))
  wdeg <- vapply(wts, function(x) if (length(x)) sum(x) else 0, numeric(1))

  list(
    A = A,
    nbrs = nbrs,
    wts = wts,
    deg = deg,
    wdeg = wdeg,
    n_nodes = nrow(A)
  )
}

