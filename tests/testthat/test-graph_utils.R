test_that("as_adjacency_matrix() rejects unsupported types", {
  expect_error(as_adjacency_matrix(1:3), "Unsupported `graph` type")
  expect_error(as_adjacency_matrix("x"), "Unsupported `graph` type")
  expect_error(as_adjacency_matrix(data.frame(a = 1)), "Unsupported `graph` type")
})


test_that("as_adjacency_matrix() validates square base matrices", {
  A <- matrix(0, 2, 3)
  expect_error(as_adjacency_matrix(A), "must be square")

  B <- matrix(c(0, 1, 0, 0), 2, 2)
  out <- as_adjacency_matrix(B, directed = TRUE)
  expect_true(is.matrix(out))
  expect_equal(diag(out), c(0, 0))
})


test_that("as_adjacency_matrix() removes self-loops and symmetrizes when undirected", {
  A <- matrix(c(1, 2,
                0, 3), 2, 2, byrow = TRUE)
  out <- as_adjacency_matrix(A, directed = FALSE, allow_weights = TRUE)
  expect_equal(diag(out), c(0, 0))
  expect_equal(out[1, 2], 2)
  expect_equal(out[2, 1], 2) # symmetrized via max
})


test_that("as_adjacency_matrix() with allow_weights=FALSE coerces to 0/1", {
  A <- matrix(c(0, 2,
                0, 0), 2, 2, byrow = TRUE)
  out <- as_adjacency_matrix(A, allow_weights = FALSE, directed = TRUE)
  expect_equal(out, matrix(c(0, 1, 0, 0), 2, 2, byrow = TRUE))
})


test_that("as_adjacency_matrix() replaces NA with 0 (base) and rejects Inf", {
  A <- matrix(c(0, NA, 0, 0), 2, 2)
  out <- as_adjacency_matrix(A, directed = TRUE)
  expect_equal(out[1, 2], 0)

  B <- matrix(c(0, Inf, 0, 0), 2, 2)
  expect_error(as_adjacency_matrix(B), "non-finite")
})


test_that("validate_adjacency() rejects negatives, self-loops, non-square", {
  A <- matrix(c(0, -1,
                0,  0), 2, 2, byrow = TRUE)
  expect_error(validate_adjacency(A, directed = TRUE), "nonnegative")

  B <- matrix(c(1, 0,
                0, 0), 2, 2)
  expect_error(validate_adjacency(B, directed = TRUE), "zero diagonal")

  C <- matrix(0, 2, 3)
  expect_error(validate_adjacency(C, directed = TRUE), "square")
})


test_that("validate_adjacency() symmetry enforcement works", {
  A <- matrix(c(0, 1,
                0, 0), 2, 2, byrow = TRUE)
  expect_error(validate_adjacency(A, directed = FALSE), "symmetric")
  expect_silent(validate_adjacency(A, directed = TRUE))
})


test_that("as_adjacency_matrix() supports adjacency lists with range checking", {
  g <- list(c(2L), integer(0))
  A <- as_adjacency_matrix(g, n_nodes = 2, directed = TRUE)
  expect_equal(A, matrix(c(0, 1, 0, 0), 2, 2, byrow = TRUE))

  g_bad <- list(c(3L), integer(0))
  expect_error(as_adjacency_matrix(g_bad, n_nodes = 2), "out-of-range")

  g_logical <- list(c(TRUE, FALSE), integer(0))
  expect_error(as_adjacency_matrix(g_logical, n_nodes = 2), "not logical")
})



test_that("adjacency_to_neighbors() works for base matrices", {
  A <- matrix(c(0, 2, 0,
                2, 0, 3,
                0, 0, 0), 3, 3, byrow = TRUE)
  nn <- adjacency_to_neighbors(A, include_weights = TRUE)
  expect_equal(nn$neighbors[[1]], c(2L))
  expect_equal(nn$weights[[1]], c(2))
  expect_equal(nn$neighbors[[2]], c(1L, 3L))
  expect_equal(nn$weights[[2]], c(2, 3))
  expect_equal(nn$neighbors[[3]], integer(0))
  expect_equal(nn$weights[[3]], numeric(0))
})


test_that("standardize_graph() returns consistent structure", {
  A <- matrix(c(0, 1, 0,
                0, 0, 2,
                0, 0, 0), 3, 3, byrow = TRUE)

  s <- standardize_graph(A, directed = TRUE, allow_weights = TRUE)
  expect_true(is.list(s))
  expect_equal(s$n_nodes, 3)
  expect_equal(s$deg, c(1L, 1L, 0L))
  expect_equal(s$wdeg, c(1, 2, 0))
})


test_that("standardize_graph() symmetrizes when directed=FALSE", {
  A <- matrix(c(0, 0, 0,
                5, 0, 0,
                0, 0, 0), 3, 3, byrow = TRUE)
  s <- standardize_graph(A, directed = FALSE, allow_weights = TRUE)
  expect_equal(s$A[1, 2], 5)
  expect_equal(s$A[2, 1], 5)
  expect_equal(s$deg[1], 1L)
  expect_equal(s$deg[2], 1L)
})


test_that("Matrix sparse adjacency is supported without densifying and checks are enforced", {
  skip_if_not_installed("Matrix")
  A <- Matrix::sparseMatrix(
    i = c(1, 2, 2),
    j = c(2, 1, 3),
    x = c(1, 1, 4),
    dims = c(3, 3),
    giveCsparse = TRUE
  )
  # add an explicit diagonal nonzero to ensure it gets removed
  A[1, 1] <- 9

  out <- as_adjacency_matrix(A, directed = TRUE, allow_weights = TRUE)
  expect_true(inherits(out, "Matrix"))
  expect_equal(as.numeric(Matrix::diag(out)), c(0, 0, 0))

  # validate rejects negative nonzeros
  B <- A
  B[2, 3] <- -1
  expect_error(validate_adjacency(B, directed = TRUE), "nonnegative")
})


test_that("sparse undirected symmetrization uses elementwise max", {
  skip_if_not_installed("Matrix")
  A <- Matrix::sparseMatrix(
    i = c(1, 2),
    j = c(2, 3),
    x = c(2, 7),
    dims = c(3, 3),
    giveCsparse = TRUE
  )
  # t(A) has different nonzeros; max should include both
  out <- as_adjacency_matrix(A, directed = FALSE, allow_weights = TRUE)
  expect_equal(out[2, 1], 2)
  expect_equal(out[1, 2], 2)
  expect_equal(out[2, 3], 7)
  expect_equal(out[3, 2], 7)
})


test_that("sparse allow_weights=FALSE coerces all nonzeros to 1", {
  skip_if_not_installed("Matrix")
  A <- Matrix::sparseMatrix(
    i = c(1, 1, 2),
    j = c(2, 3, 3),
    x = c(2, 5, 9),
    dims = c(3, 3),
    giveCsparse = TRUE
  )
  out <- as_adjacency_matrix(A, directed = TRUE, allow_weights = FALSE)
  expect_true(inherits(out, "Matrix"))
  expect_equal(sort(as.numeric(out@x)), c(1, 1, 1))
})


test_that("adjacency_to_neighbors() works for sparse Matrix", {
  skip_if_not_installed("Matrix")
  A <- Matrix::sparseMatrix(
    i = c(1, 2, 2),
    j = c(2, 1, 3),
    x = c(2, 2, 4),
    dims = c(3, 3),
    giveCsparse = TRUE
  )
  nn <- adjacency_to_neighbors(A, include_weights = TRUE)
  expect_equal(nn$neighbors[[1]], c(2L))
  expect_equal(nn$weights[[1]], c(2))
  expect_equal(nn$neighbors[[2]], c(1L, 3L))
  expect_equal(nn$weights[[2]], c(2, 4))
  expect_equal(nn$neighbors[[3]], integer(0))
})


test_that("igraph inputs are supported when igraph is installed", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("Matrix")

  g <- igraph::make_ring(5, directed = FALSE)
  A <- as_adjacency_matrix(g, directed = FALSE, allow_weights = FALSE)
  expect_true(inherits(A, "Matrix") || is.matrix(A))
  validate_adjacency(A, directed = FALSE)

  s <- standardize_graph(g, directed = FALSE, allow_weights = FALSE)
  expect_equal(s$n_nodes, 5)
  expect_true(all(s$deg == 2L))
})


test_that("test extreme edge cases: empty graph, isolated nodes, and duplicates in list", {
  # Empty edges, 4 isolated nodes
  A <- matrix(0, 4, 4)
  s <- standardize_graph(A, directed = FALSE)
  expect_equal(s$deg, rep(0L, 4))
  expect_equal(s$wdeg, rep(0, 4))

  # adjacency list with duplicate neighbors should still produce 1s
  g <- list(c(2L, 2L, 3L), integer(0), integer(0))
  A2 <- as_adjacency_matrix(g, n_nodes = 3, directed = TRUE, allow_weights = FALSE)
  expect_equal(A2[1, 2], 1)
  expect_equal(A2[1, 3], 1)
})


test_that("randomized adversarial base matrices do not stress algorithm", {
  set.seed(1)

  fails <- character()
  for (iter in 1:50) {
    n <- sample(2:10, 1)
    A <- matrix(sample(c(0, 0, 0, 1, 2, 5), n * n, replace = TRUE), n, n)
    diag(A) <- sample(c(0, 0, 0, 1), n, replace = TRUE)

    if (any(diag(A) != 0)) {
      ok <- tryCatch({
        validate_adjacency(A, directed = TRUE)
        FALSE
      }, error = function(e) grepl("zero diagonal", conditionMessage(e)))
      if (!ok) fails <- c(fails, sprintf("iter=%d: expected zero-diagonal error", iter))
      next
    }

    ok <- tryCatch({
      validate_adjacency(A, directed = TRUE)
      out1 <- as_adjacency_matrix(A, directed = TRUE, allow_weights = TRUE)
      out0 <- as_adjacency_matrix(A, directed = TRUE, allow_weights = FALSE)

      if (!all(out0 %in% c(0, 1))) stop("allow_weights=FALSE not binary")
      if (!all(diag(out1) == 0)) stop("diagonal not zeroed")

      TRUE
    }, error = function(e) {
      fails <<- c(fails, sprintf("iter=%d: %s", iter, conditionMessage(e)))
      FALSE
    })

    if (!ok) next
  }

  expect_true(length(fails) == 0, info = paste(fails, collapse = "\n"))
})


test_that("randomized adversarial sparse matrices do not stress algorithm", {
  skip_if_not_installed("Matrix")
  set.seed(2)

  fails <- character()
  for (iter in 1:50) {
    n <- sample(3:20, 1)
    nnz <- sample(0:(n * 2), 1)

    A <- if (nnz == 0) {
      Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(n, n))
    } else {
      i <- sample(1:n, nnz, replace = TRUE)
      j <- sample(1:n, nnz, replace = TRUE)
      x <- sample(c(0.5, 1, 2, 5), nnz, replace = TRUE)
      Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n, n), giveCsparse = TRUE)
    }

    # inject a diagonal entry; should be stripped
    A[sample.int(n, 1), sample.int(n, 1)] <- 3

    tryCatch({
      out <- as_adjacency_matrix(A, directed = FALSE, allow_weights = TRUE)
      if (!inherits(out, "Matrix")) stop("output not Matrix")
      if (!all(as.numeric(Matrix::diag(out)) == 0)) stop("diagonal not zeroed")

      nn <- adjacency_to_neighbors(out, include_weights = TRUE)
      if (length(nn$neighbors) != n) stop("neighbors length mismatch")
      if (length(nn$weights) != n) stop("weights length mismatch")
    }, error = function(e) {
      fails <<- c(fails, sprintf("iter=%d: %s", iter, conditionMessage(e)))
    })
  }

  expect_true(length(fails) == 0, info = paste(fails, collapse = "\n"))
})
