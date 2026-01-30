#' Graph model interface for bocpdGraph
#'
#' These are S3 generics used by \code{\link{bocpdGraph}}. Implement methods for
#' your model class (e.g., \code{graph_model_init.myclass}).
#'
#' @param model A graph-aware model object.
#' @param graph Graph structure (standardized adjacency or related object).
#' @param data Node-by-time data matrix.
#' @param control Optional control list.
#' @param y_t Observation at time \code{t} for a given node.
#' @param node Integer node index.
#' @param run_length Integer vector of run lengths.
#' @param neighbors Integer vector of neighbor node indices.
#' @param t Integer time index.
#' @param ... Additional arguments.
#'
#' @name graph_model_interface
#' @rdname graph_model_interface
#' @export
graph_model_init <- function(model, graph, data, control = list(), ...) UseMethod("graph_model_init")

#' @rdname graph_model_interface
#' @export
graph_pred_loglik <- function(model, y_t, node, run_length, neighbors = integer(), t = NULL, ...) UseMethod("graph_pred_loglik")

#' @rdname graph_model_interface
#' @export
graph_model_update <- function(model, y_t, node, run_length, neighbors = integer(), t = NULL, ...) UseMethod("graph_model_update")


#' @export
graph_model_init <- function(model, graph, data, control = list(), ...) {
  UseMethod("graph_model_init")
}

#' @export
graph_pred_loglik <- function(model, y_t, node, run_length, neighbors = integer(), t = NULL, ...) {
  UseMethod("graph_pred_loglik")
}

#' @export
graph_model_update <- function(model, y_t, node, run_length, neighbors = integer(), t = NULL, ...) {
  UseMethod("graph_model_update")
}

#' @export
graph_model_init.default <- function(model, graph, data, control = list(), ...) {
  stop("No graph_model_init() method for class: ", paste(class(model), collapse = "/"))
}

#' @export
graph_pred_loglik.default <- function(model, y_t, node, run_length, neighbors = integer(), t = NULL, ...) {
  stop("No graph_pred_loglik() method for class: ", paste(class(model), collapse = "/"))
}

#' @export
graph_model_update.default <- function(model, y_t, node, run_length, neighbors = integer(), t = NULL, ...) {
  stop("No graph_model_update() method for class: ", paste(class(model), collapse = "/"))
}
