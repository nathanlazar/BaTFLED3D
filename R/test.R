#' Get test predictions for a 3D BaTFLED model.
#'
#' This is just a wrapper that calls test_CP or test_Tucker depending on the type of model provided.
#'
#' @export
#' @param d object of the class \code{input_data} created with input_data()
#' @param m object of the class \code{CP_model} or \code{Tucker_model} created with mk_model()
#' @param ... extra parameters passed to test_CP or test_Tucker
#' @return An array of predicted responses the same size as \code{m$resp}.

test <- function(d, m, ...) {
  if(length(m$core.mean)==0) {
    test_CP(d, m, ...)
  } else {
    test_Tucker(d, m, ...)
  }
}
