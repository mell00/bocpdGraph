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
