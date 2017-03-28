#' Computes the normalized root mean squared error
#' 
#' @export
#' @param obs observed vector, matrix or data.frame
#' @param pred predicted vector, matrix or data.frame
#' 
#' @return numeric value of the root mean squared error normalized to the standard deviation
#' of the observed data

nrmse <- function(obs, pred) {
  sqrt( mean((obs-pred)^2, na.rm=T) / mean((obs-mean(obs, na.rm=T))^2, na.rm=T) )
}
