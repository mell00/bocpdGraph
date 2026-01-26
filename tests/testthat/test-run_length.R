test_that("log_sum_exp is stable on extreme inputs", {
  # all -Inf
  expect_true(is.infinite(log_sum_exp(c(-Inf, -Inf))))
  expect_equal(log_sum_exp(c(-Inf, -Inf)), -Inf)

  # single finite
  expect_equal(log_sum_exp(0), 0)

  # large magnitude mix (no overflow)
  x <- c(1000, 999, 998)
  out <- log_sum_exp(x)
  expect_true(is.finite(out))
  expect_true(out >= max(x))
  expect_true(out <= max(x) + log(length(x)) + 1e-12)

  # mix with -Inf
  x2 <- c(-Inf, 0, -1000)
  out2 <- log_sum_exp(x2)
  expect_true(is.finite(out2))
  expect_equal(out2, log(sum(exp(x2))), tolerance = 1e-12)

  # very negative values (underflow-safe)
  x3 <- c(-10000, -10001, -9999)
  out3 <- log_sum_exp(x3)
  expect_true(is.finite(out3))

  # compare via normalization trick
  m <- max(x3)
  expect_equal(out3, m + log(sum(exp(x3 - m))), tolerance = 1e-12)
})


test_that("normalize_log_probs returns a proper log-prob vector", {
  logp <- c(0, -1, -2)
  out <- normalize_log_probs(logp)
  p <- exp(out)
  expect_equal(sum(p), 1, tolerance = 1e-12)
  expect_true(all(p >= 0))
  expect_true(all(is.finite(out)))

  # shift invariance
  out2 <- normalize_log_probs(logp + 1000)
  expect_equal(out, out2, tolerance = 1e-12)

  # degenerate: one mass at -Inf is okay if at least one finite
  logp3 <- c(-Inf, 0, -Inf)
  out3 <- normalize_log_probs(logp3)
  p3 <- exp(out3)
  expect_equal(p3, c(0, 1, 0), tolerance = 0)

  # all -Inf => undefined normalization; should return all NaN/-Inf? current impl returns -Inf
  # only assert it doesn't error and preserves shape
  out4 <- normalize_log_probs(c(-Inf, -Inf))
  expect_equal(length(out4), 2)
})


test_that("as_hazard_fun handles function hazards and constant hazards", {
  hf <- function(run_length, t, ...) rep(0.123, length(run_length))
  out <- as_hazard_fun(hf)
  expect_true(is.function(out))
  expect_equal(out(0:4, 1), rep(0.123, 5))

  out2 <- as_hazard_fun(0.2)
  expect_true(is.function(out2))
  expect_equal(out2(0:4, 99), rep(0.2, 5))

  expect_error(as_hazard_fun(1.2), "in \\[0,1\\]")
  expect_error(as_hazard_fun(-0.1), "in \\[0,1\\]")
  expect_error(as_hazard_fun(c(0.1, 0.2)))
  expect_error(as_hazard_fun(NA_real_))
})


test_that("truncate_run_length validates and truncates consistently", {
  log_R <- c(0, -1, -2, -3)
  stats <- list(a=1, b=2, c=3, d=4)

  tr <- truncate_run_length(log_R, stats, max_run_length = 2)
  expect_equal(length(tr$log_R), 3)   # 0..2
  expect_equal(length(tr$stats), 3)
  expect_equal(tr$log_R, log_R[1:3])
  expect_equal(tr$stats, stats[1:3])

  # no truncation
  tr2 <- truncate_run_length(log_R, stats, max_run_length = 10)
  expect_equal(tr2$log_R, log_R)
  expect_equal(tr2$stats, stats)

  # validations
  expect_error(truncate_run_length("x", stats, 2), "numeric")
  expect_error(truncate_run_length(log_R, "x", 2), "list")
  expect_error(truncate_run_length(log_R, stats[-1], 2), "same length")
  expect_error(truncate_run_length(log_R, stats, -1), "nonnegative")
})



test_that("update_run_length enforces invariants and is numerically stable", {
  set.seed(1)

  # basic well-formed
  log_R_prev <- normalize_log_probs(c(0, -1, -2))  # r=0..2
  loglik <- c(-0.5, -0.2, -1.0)
  loglik_prior <- -0.3

  res <- update_run_length(
    log_R_prev = log_R_prev,
    loglik = loglik,
    loglik_prior = loglik_prior,
    hazard = 0.1,
    t = 1,
    max_run_length = 5
  )
  expect_true(is.list(res))
  expect_equal(length(res$log_R), 6)
  expect_equal(length(res$R), 6)
  expect_equal(sum(res$R), 1, tolerance = 1e-12)
  expect_true(all(res$R >= 0))
  expect_true(all(is.finite(res$log_R)))

  # hazard = 0 => changepoint mass should be 0 (unless numerical), growth only
  res0 <- update_run_length(
    log_R_prev = log_R_prev,
    loglik = loglik,
    loglik_prior = loglik_prior,
    hazard = 0,
    t = 2,
    max_run_length = 5
  )
  expect_equal(res0$R[1], 0, tolerance = 1e-15)
  expect_equal(sum(res0$R), 1, tolerance = 1e-12)

  # hazard = 1 => all mass resets to r=0
  res1 <- update_run_length(
    log_R_prev = log_R_prev,
    loglik = loglik,
    loglik_prior = loglik_prior,
    hazard = 1,
    t = 3,
    max_run_length = 5
  )
  expect_equal(res1$R, c(1, rep(0, 5)), tolerance = 1e-15)

  # if max_run_length small, still normalized and stable
  res_tr <- update_run_length(
    log_R_prev = log_R_prev,
    loglik = loglik,
    loglik_prior = loglik_prior,
    hazard = 0.2,
    t = 4,
    max_run_length = 0
  )
  expect_equal(length(res_tr$R), 1)
  expect_equal(res_tr$R[1], 1, tolerance = 1e-12)

  # very large negative/positive values should not overflow
  log_R_prev2 <- normalize_log_probs(c(0, -1000, -2000, -3000))
  loglik2 <- c(1000, -1000, -1000, -1000)  # huge
  res_big <- update_run_length(
    log_R_prev = log_R_prev2,
    loglik = loglik2,
    loglik_prior = 1000,
    hazard = 0.5,
    t = 5,
    max_run_length = 10
  )
  expect_equal(sum(res_big$R), 1, tolerance = 1e-12)
  expect_true(all(is.finite(res_big$log_R)))

  # validations
  expect_error(update_run_length(numeric(0), numeric(0), 0, 0.1, 1, 5), "non-empty")
  expect_error(update_run_length(c(0, -1), c(-1), 0, 0.1, 1, 5), "same length")
  expect_error(update_run_length(c(0, -1), c(-1, -1), NA_real_, 0.1, 1, 5), "finite")
  expect_error(update_run_length(c(0, -1), c(-1, -1), 0, 0.1, 0, 5), ">= 1")
  expect_error(update_run_length(c(0, -1), c(-1, -1), 0, 0.1, 1, -1), "nonnegative")

  # hazard function wrong length
  bad_h <- function(run_length, t, ...) 0.2
  expect_error(
    update_run_length(log_R_prev, loglik, loglik_prior, bad_h, 1, 5),
    "wrong length"
  )

  # hazard function out of bounds
  bad_h2 <- function(run_length, t, ...) rep(2, length(run_length))
  expect_error(
    update_run_length(log_R_prev, loglik, loglik_prior, bad_h2, 1, 5),
    "in \\[0,1\\]"
  )
})
