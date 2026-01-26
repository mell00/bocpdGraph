
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


test_that("hazard_logistic_run_length is monotone when b>0 and stable for extremes", {
  h <- hazard_logistic_run_length(a = -10, b = 0.5)
  r <- 0:20
  out <- h(r, t = 1)

  expect_equal(length(out), length(r))
  expect_true(all(out >= 0 & out <= 1))
  # monotone nondecreasing (allow tiny numeric noise)
  expect_true(all(diff(out) >= -1e-12))

  # extreme positive eta -> ~1
  h_hi <- hazard_logistic_run_length(a = 1000, b = 0)
  out_hi <- h_hi(0:3, t = 1)
  expect_true(all(out_hi > 1 - 1e-12))

  # extreme negative eta -> ~0
  h_lo <- hazard_logistic_run_length(a = -1000, b = 0)
  out_lo <- h_lo(0:3, t = 1)
  expect_true(all(out_lo < 1e-12))

  expect_error(hazard_logistic_run_length(a = NA_real_, b = 1), "finite")
  expect_error(hazard_logistic_run_length(a = 1, b = NA_real_), "finite")
})



test_that("hazard_from_pmf matches definition and handles tails", {
  # uniform on {1,2,3}
  h <- hazard_from_pmf(c(1, 1, 1), tail_hazard = 1)

  # survival: S1=1, S2=2/3, S3=1/3
  # hazard r=0 -> pmf1/S1 = 1/3
  # hazard r=1 -> pmf2/S2 = (1/3)/(2/3) = 1/2
  # hazard r=2 -> pmf3/S3 = (1/3)/(1/3) = 1
  out <- h(0:5, t = 1)
  expect_equal(out[1], 1/3, tolerance = 1e-12)
  expect_equal(out[2], 1/2, tolerance = 1e-12)
  expect_equal(out[3], 1, tolerance = 1e-12)
  # tail hazard forces cp for r>=3
  expect_equal(out[4:6], rep(1, 3))

  # tail_hazard not 1
  h2 <- hazard_from_pmf(c(0.2, 0.8), tail_hazard = 0.25)
  out2 <- h2(0:10, t = 1)
  expect_equal(out2[3], 0.25) # r=2 >= K
  expect_true(all(out2 >= 0 & out2 <= 1))

  # pmf normalization invariant
  h3 <- hazard_from_pmf(c(2, 2, 2), tail_hazard = 0.5)
  expect_equal(h3(0:2, t = 1), h(0:2, t = 1), tolerance = 1e-12)

  # validations
  expect_error(hazard_from_pmf(numeric(0)), "non-empty")
  expect_error(hazard_from_pmf(c(NA_real_)), "finite")
  expect_error(hazard_from_pmf(c(-1, 1)), "nonnegative")
  expect_error(hazard_from_pmf(c(0, 0, 0)), "positive sum")
  expect_error(hazard_from_pmf(c(1, 1), tail_hazard = 1.5), "\\[0,1\\]")
})



test_that("All hazards return length(run_length) and stay in [0,1] under random stress", {
  set.seed(123)

  hazards <- list(
    hazard_constant(0.01),
    hazard_piecewise_time(times = c(10, 20, 30), p = c(0.1, 0.2, 0.3, 0.4)),
    hazard_logistic_run_length(a = -5, b = 0.1),
    hazard_from_pmf(c(0.1, 0.2, 0.3, 0.4), tail_hazard = 0.05)
  )

  hazard_names <- c("constant", "piecewise_time", "logistic_run_length", "from_pmf")

  for (i in seq_along(hazards)) {
    h <- hazards[[i]]
    name <- hazard_names[[i]]

    fail_msg <- NULL

    for (iter in 1:200) {
      n <- sample.int(50, 1)
      r <- sample.int(200, n, replace = TRUE) - 1L
      t <- sample.int(1000, 1)

      out <- h(r, t = t)

      if (length(out) != length(r)) {
        fail_msg <- sprintf("iter=%d length(out)=%d != length(r)=%d", iter, length(out), length(r))
        break
      }
      if (!all(is.finite(out))) {
        fail_msg <- sprintf("iter=%d non-finite values present", iter)
        break
      }
      if (!all(out >= 0 & out <= 1)) {
        bad <- which(!(out >= 0 & out <= 1))[1]
        fail_msg <- sprintf("iter=%d out of bounds at index %d: %g", iter, bad, out[[bad]])
        break
      }
    }

    expect_true(
      is.null(fail_msg),
      info = paste0("Hazard '", name, "' failed: ", fail_msg)
    )
  }
})
