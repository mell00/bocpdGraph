test_that("as_adjacency_matrix() rejects unsupported types", {
  expect_error(as_adjacency_matrix(1:3), "Unsupported `graph` type")
  expect_error(as_adjacency_matrix("x"), "Unsupported `graph` type")
  expect_error(as_adjacency_matrix(data.frame(a = 1)), "Unsupported `graph` type")
})
