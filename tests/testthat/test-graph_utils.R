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
