#' Takes the product of two matrices adding a column of constants if necessary to the first matrix.
#'
#' @export
#' @param A matrix one
#' @param B matrix two
#' @return matrix product of A and B
#'

safe_prod <- function(A, B) {
  # take the product of two matrices adding a column of 1s if necessary
  if(ncol(A) == nrow(B)-1) {
    C <- cbind(1, A) %*% B
  } else C <- A %*% B
  return(C)
}