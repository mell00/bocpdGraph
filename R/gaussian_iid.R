#' IID Gaussian mean model (reference)
#'
#' Minimal Gaussian model provided primarily for testing and examples.
#'
#' @param mu0 Prior mean.
#' @param sigma Known standard deviation.
#'
#' @return A model object usable with \code{\link{bocpd}}.
#'
#' @export
gaussian_iid_model <- function(mu0 = 0, sigma = 1) {

  prior <- list(nu = 1, chi = list(mu = mu0))

  model <- list(
    prior = prior,
    pred_loglik = function(nu, chi, x_t) {
      stats::dnorm(x_t, mean = chi$mu, sd = sigma, log = TRUE)
    }
  )

  class(model) <- "gaussian_iid"
  model
}

#' @export
model_init.gaussian_iid <- function(model) {
  list(nu = model$prior$nu, chi = list(model$prior$chi))
}

#' @export
model_update.gaussian_iid <- function(model, state, x_t) {
  chi_new <- lapply(state$chi, function(chi) {
    list(mu = (chi$mu * state$nu[1] + x_t) / (state$nu[1] + 1))
  })
  list(
    nu  = c(model$prior$nu, state$nu + 1),
    chi = c(list(model$prior$chi), chi_new)
  )
}
