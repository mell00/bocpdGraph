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
