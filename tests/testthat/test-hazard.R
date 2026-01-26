
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
