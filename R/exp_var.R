#' Get the explained variance for a set of predictions
#' 
#' Calculates 1-var(obs-pred)/var(obs). If verbose == TRUE the result is printed.
#' 
#' @importFrom stats var
#' 
#' @export
#' @param obs data.frame, vector or matrix
#' @param pred data.frame, vector or matrix
#' @param verbose logical indicating whether to print result
#' 
#' @return numeric value of the explained variance
#' 
#' @examples
#' exp_var(rnorm(100) + seq(0,9.9,.1),  seq(0,9.9,.1))
#' 

exp_var <- function(obs, pred, verbose=F) {
  # Get the explained variance for a set of predictions
  if(is.data.frame(obs)) obs <- unlist(obs)
  if(!is.vector(obs)) obs <- as.vector(obs)
  if(is.data.frame(pred)) pred <- unlist(pred)
  if(!is.vector(pred)) pred <- as.vector(pred)
  exp.var <- 1 - var(obs - pred, na.rm=T)/var(obs, na.rm=T)
  
  if(verbose) print(sprintf("Explained variance: %.4f", exp.var))
  exp.var
}
