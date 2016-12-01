#' Make a new model object 
#' 
#' This function is a wrapper calling Tucker_model$new()  or CP_model$new() depending
#' on whether \code{params$decomp=='Tucker'} or \code{params$decomp=='CP'}
#' 
#' @export
#' @param d An \code{input_data} object. See \code{input_data}.
#' @param params A list of parameter values \code{get_model_params}.
#' @return \code{CP_model} or \code{Tucker_model} object
#' @seealso \link{Tucker_model}, \link{CP_model}
#'
mk_model <- function(d, params) {
  if(params$decomp=='CP') model <- CP_model$new(d, params)
  if(params$decomp=='Tucker') model <- Tucker_model$new(d, params)
  return(model)
}