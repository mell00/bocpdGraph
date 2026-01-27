#' Initialize model state
#'
#' @param model Model object.
#' @return Initial model state.
#'
#' @export
model_init <- function(model) {
  UseMethod("model_init")
}

#' Update model state after observing data
#'
#' @param model Model object.
#' @param state Current model state.
#' @param x_t New observation.
#' @return Updated model state.
#'
#' @export
model_update <- function(model, state, x_t) {
  UseMethod("model_update")
}
