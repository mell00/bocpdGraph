test_that("bocpd validates hazard contract", {
  model <- gaussian_iid_model(mu0 = 0, sigma = 1)

  expect_error(
    bocpd(as.list(1:3), model, hazard = 1),
    "hazard.*function"
  )

  bad_hazard1 <- function(run_length, t, ...) rep(1.5, length(run_length))
  expect_error(
    bocpd(as.list(1:3), model, bad_hazard1),
    "hazard.*\\[0,1\\]"
  )

  bad_hazard2 <- function(run_length, t, ...) rep(NA_real_, length(run_length))
  expect_error(
    bocpd(as.list(1:3), model, bad_hazard2),
    "hazard.*\\[0,1\\]"
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



test_that("bocpd respects r_max truncation exactly", {
  set.seed(1)
  model <- gaussian_iid_model(mu0 = 0, sigma = 1)
  h <- hazard_constant(1/100)

  x <- rnorm(80)

  # r_max = 10  => length(rl_t) <= 11 for all t
  fit1 <- bocpd(x, model, h, control = list(r_max = 10, prune_eps = 0))
  lens1 <- vapply(fit1$rl, length, integer(1))
  expect_true(all(lens1 <= 11L))

  # r_max = 0 => only r=0 kept: length == 1 and value == 1
  fit2 <- bocpd(x, model, h, control = list(r_max = 0, prune_eps = 0))
  lens2 <- vapply(fit2$rl, length, integer(1))
  expect_true(all(lens2 == 1L))
  expect_true(all(vapply(fit2$rl, function(v) isTRUE(all.equal(v[1], 1, tolerance = 1e-12)), logical(1))))
})



test_that("bocpd pruning removes tiny mass but keeps r=0 and renormalizes", {
  set.seed(2)
  model <- gaussian_iid_model(mu0 = 0, sigma = 1)
  h <- hazard_constant(1/200)

  x <- rnorm(150)
  fit <- bocpd(x, model, h, control = list(r_max = 300, prune_eps = 1e-6))

  for (t in c(10, 50, 150)) {
    rl <- fit$rl[[t]]
    expect_true(length(rl) >= 1L)
    expect_equal(sum(rl), 1, tolerance = 1e-12)
    expect_true(rl[1] > 0)
    if (length(rl) > 1L) {
      expect_true(all(rl[-1] > 1e-6))
    }
  }
})


