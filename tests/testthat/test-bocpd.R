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


test_that("bocpd handles extreme observations without NaN/Inf posteriors", {
  model <- gaussian_iid_model(mu0 = 0, sigma = 1)
  h <- hazard_constant(1/100)

  x <- c(0, 1e6, -1e6, 1e12, -1e12, 0)
  fit <- bocpd(x, model, h, control = list(r_max = 50, prune_eps = 0))

  for (t in seq_along(x)) {
    rl <- fit$rl[[t]]
    expect_true(all(is.finite(rl)))
    expect_true(all(rl >= 0))
    expect_equal(sum(rl), 1, tolerance = 1e-12)
  }
})


test_that("near-zero and near-one hazards behave sensibly", {
  set.seed(3)
  model <- gaussian_iid_model(mu0 = 0, sigma = 1)

  x <- rnorm(30)

  h0 <- function(run_length, t, ...) rep(1e-12, length(run_length))
  fit0 <- bocpd(x, model, h0, control = list(r_max = 200, prune_eps = 0))
  expect_lt(fit0$rl[[30]][1], 1e-6)

  h1 <- function(run_length, t, ...) rep(1 - 1e-12, length(run_length))
  fit1 <- bocpd(x, model, h1, control = list(r_max = 200, prune_eps = 0))
  expect_gt(fit1$rl[[30]][1], 0.9)
})



test_that("bocpd supports list input and vector input equivalently", {
  set.seed(4)
  model <- gaussian_iid_model(mu0 = 0, sigma = 1)
  h <- hazard_constant(1/50)

  x <- rnorm(40)
  fit_vec <- bocpd(x, model, h, control = list(r_max = 80, prune_eps = 0))
  fit_lst <- bocpd(as.list(x), model, h, control = list(r_max = 80, prune_eps = 0))

  expect_equal(length(fit_vec$rl), length(fit_lst$rl))
  for (t in c(1, 10, 40)) {
    expect_equal(fit_vec$rl[[t]], fit_lst$rl[[t]], tolerance = 1e-12)
  }
})



test_that("bocpd_changepoints(cp_prob) returns valid indices and respects min_sep", {
  set.seed(5)
  x <- c(rnorm(60, 0), rnorm(60, 4))
  model <- gaussian_iid_model(mu0 = 0, sigma = 1)
  h <- hazard_constant(1/100)

  fit <- bocpd(x, model, h, control = list(r_max = 200, prune_eps = 1e-12))

  cps <- bocpd_changepoints(fit, method = "cp_prob", threshold = 0.2, min_sep = 10)
  expect_true(is.integer(cps) || is.numeric(cps))
  expect_true(all(cps >= 1 & cps <= length(x)))
  if (length(cps) > 1) {
    expect_true(all(diff(cps) >= 10))
  }
})

test_that("bocpd_changepoints(map_drop) returns valid indices and is robust", {
  set.seed(6)
  x <- c(rnorm(40, 0), rnorm(40, 2), rnorm(40, -2))
  model <- gaussian_iid_model(mu0 = 0, sigma = 1)
  h <- hazard_constant(1/80)

  fit <- bocpd(x, model, h, control = list(r_max = 250, prune_eps = 1e-14))

  cps <- bocpd_changepoints(fit, method = "map_drop", threshold = 8, min_sep = 5)
  expect_true(all(cps >= 1 & cps <= length(x)))
  if (length(cps) > 1) {
    expect_true(all(diff(cps) >= 5))
  }
})


test_that("bocpd returns expected structure (rl + state, optional rl_matrix)", {
  set.seed(7)
  model <- gaussian_iid_model(mu0 = 0, sigma = 1)
  h <- hazard_constant(1/30)
  x <- rnorm(25)

  fit <- bocpd(x, model, h, control = list(return_rl_matrix = FALSE))
  expect_true(is.list(fit))
  expect_true("rl" %in% names(fit))
  expect_true("state" %in% names(fit))
  expect_false("rl_matrix" %in% names(fit))

  fitm <- bocpd(x, model, h, control = list(r_max = 20, return_rl_matrix = TRUE))
  expect_true("rl_matrix" %in% names(fitm))
  expect_true(is.matrix(fitm$rl_matrix))
  expect_equal(nrow(fitm$rl_matrix), length(x))
  expect_equal(ncol(fitm$rl_matrix), 21)
  rs <- rowSums(fitm$rl_matrix)
  expect_true(all(abs(rs - 1) < 1e-10))
})

