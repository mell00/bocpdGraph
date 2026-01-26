
test_that("hazard_constant / hazard_geometric output correct length and bounds", {
  h <- hazard_constant(0.3)
  out <- h(0:9, t = 1)
  expect_equal(length(out), 10)
  expect_true(all(out == 0.3))
  expect_true(all(out >= 0 & out <= 1))

  hg <- hazard_geometric(0.7)
  outg <- hg(0:3, t = 10)
  expect_equal(outg, rep(0.7, 4))

  expect_error(hazard_constant(NA_real_), "finite")
  expect_error(hazard_constant(-0.1), "in \\[0,1\\]")
  expect_error(hazard_constant(1.1), "in \\[0,1\\]")
})


test_that("hazard_piecewise_time behaves correctly at boundaries and validates input", {
  h <- hazard_piecewise_time(times = c(5, 10), p = c(0.1, 0.2, 0.3))

  expect_equal(h(0:2, t = 1), rep(0.1, 3))
  expect_equal(h(0:2, t = 5), rep(0.1, 3))
  expect_equal(h(0:2, t = 6), rep(0.2, 3))
  expect_equal(h(0:2, t = 10), rep(0.2, 3))
  expect_equal(h(0:2, t = 11), rep(0.3, 3))

  # empty times => always p[1]
  h2 <- hazard_piecewise_time(times = integer(0), p = 0.25)
  expect_equal(h2(0:4, t = 999), rep(0.25, 5))

  expect_error(hazard_piecewise_time(times = c(2, 1), p = c(0.1, 0.2, 0.3)), "strictly increasing")
  expect_error(hazard_piecewise_time(times = c(5, 10), p = c(0.1, 0.2)), "length")
  expect_error(hazard_piecewise_time(times = c(5, 10), p = c(-0.1, 0.2, 0.3)), "\\[0,1\\]")
  expect_error(h(0:2, t = NA_real_), "finite")
})
