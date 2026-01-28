#' bocpdGraph: Bayesian Online Changepoint Detection
#'
#' This package implements Bayesian online changepoint detection (BOCPD)
#' following Adams and MacKay (2007), with additional implementation-level
#' safeguards and extensions for practical use.
#'
#' @details
#' \strong{Implementation notes and deviations from the idealized algorithm}
#'
#' The core BOCPD recursion implemented in \code{\link{bocpd}} follows the
#' two-branch run-length message passing scheme of Adams and MacKay (2007).
#' The following differences are intentional and arise from numerical
#' stability and software engineering considerations:
#'
#' \itemize{
#'   \item \emph{Hazard interface.} Hazard functions are called as
#'   \code{hazard(run_length, t, ...)} to allow time-varying hazards. The
#'   original formulation typically treats the hazard as a function of
#'   run length only.
#'
#'   \item \emph{Run-length truncation and pruning.} Optional truncation
#'   (\code{r_max}) and probability pruning (\code{prune_eps}) are provided
#'   to control computational complexity. These operations approximate the
#'   exact recursion but do not change the underlying probabilistic model.
#'
#'   \item \emph{Numerical safeguards.} When floating-point underflow or
#'   overflow prevents reliable normalization of the run-length posterior,
#'   the implementation falls back to a valid distribution concentrating
#'   all mass on \eqn{r_t = 0}. This preserves probabilistic invariants under
#'   extreme inputs but is not part of the mathematical algorithm.
#' }
#'
#' @references
#' Adams, R. P. and MacKay, D. J. C. (2007).
#' \emph{Bayesian Online Changepoint Detection}.
#' arXiv:0710.3742.
#'
#' Bretz, F. (2021).
#' \emph{Bayesian Online Changepoint Detection}.
#' Lecture notes / technical report.
#'
#' @keywords internal
"_PACKAGE"
