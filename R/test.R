#' Get test predictions for a 3D BaTFLED model.
#'
#' This is just a wrapper that calls test_CP or test_Tucker depending on the type of model provided.
#'
#' @param d data object created with input_data
#' @param m model object created with mk_model
#'
#' @export
#' @param d object of the class \code{input_data}
#' @param m object of the class \code{CP_model} or \code{Tucker_model}
#' @return An array of predicted responses the same size as \code{m$resp}.
#'
#' @examples value
#'

test <- function(d, m, ...) {
  if(length(m$core.mean)==0) {
    test_CP(d, m, ...)
  } else {
    test_Tucker(d, m, ...)
  }
}
