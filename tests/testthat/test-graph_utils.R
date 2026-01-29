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

