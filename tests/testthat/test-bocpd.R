test_that("bocpd validates hazard contract", {
  model <- gaussian_iid_model()

  expect_error(
    bocpd(as.list(1:3), model, hazard = 1),
    "hazard.*function"
  )

  bad_hazard1 <- function(r) rep(1.5, length(r))
  expect_error(
    bocpd(as.list(1:3), model, bad_hazard1),
    "\\[0,1\\]"
  )

  bad_hazard2 <- function(r) rep(NA_real_, length(r))
  expect_error(
    bocpd(as.list(1:3), model, bad_hazard2),
    "\\[0,1\\]"
  )
})


test_that("bocpd run-length posterior invariants hold (nonneg, sums to 1, finite)", {
  set.seed(123)

  model <- gaussian_iid_model(mu0 = 0, sigma = 1)
  h <- hazard_constant(1/50)

  x <- rnorm(200)
  fit <- bocpd(x, model, h, control = list(r_max = 250, prune_eps = 0))

  expect_length(fit$rl, length(x))

  for (t in c(1, 2, 5, 50, 200)) {
    rl <- fit$rl[[t]]
    expect_true(is.numeric(rl))
    expect_true(all(is.finite(rl)))
    expect_true(all(rl >= 0))
    expect_equal(sum(rl), 1, tolerance = 1e-12)
  }
})

