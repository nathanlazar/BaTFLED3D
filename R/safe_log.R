#' Take logarithm avoiding underflow
#'
#' Returns the normal log if there is no underflow. If there is underflow, then
#' returns the minimum for which log can return (-744.4401)
#'
#' @export
#' @param x vector
#' @return vector log in base e of input or minimum possible log value of -744.4401
#'
#' @examples
#' log(c(1e-323, 1e-324))      # gives -Inf for the second value
#' safe_log(c(1e-323, 1e-324)) # gives the minimum value of -744.4401

safe_log <- function(x) {
  ret <- log(x)
#   if(sum(ret == -Inf) > 0) stop('Cutting log')
  ret[ret == -Inf] <- -744.4401
  return(ret)
}
