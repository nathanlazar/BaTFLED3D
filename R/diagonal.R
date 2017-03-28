#' Version of diag that has more consistent behavior
#' 
#' @export
#' @param x A vector, matrix or array with third mode length 1
#' @param len numeric dimensions of new diagonal matrix to me made. Recycles values in x.
#' @param ... parameters passed to diag
#' 
#' @return if x is a vector or integer, return a matrix with x on the diagonal. 
#' If x is a matrix, or degenerate array, return the diagonal of x.
#' 
#' @examples
#' diagonal(c(1,3))
#' diagonal(matrix(1:6, 2,3))
#' diagonal(5)
#' diagonal(c(5,2),3)
#' diagonal(array(1:12, dim=c(3,4,1)))

diagonal <- function(x, len=NA, ...) {
  if(is.matrix(x))
    return(diag(x, ...))
  if(is.array(x) && dim(x)[3]==1)
    return(diagonal(matrix(x[,,1], dim(x)[1], dim(x)[2])))
  if(length(x) == 1L && is.na(len))
    return(matrix(x, 1,1))
  else {
    if(is.na(len)) len <- length(x)
    return(diag(x, len, ...))
  }
}