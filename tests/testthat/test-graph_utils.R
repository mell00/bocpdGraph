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
