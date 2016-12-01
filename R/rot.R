#' Rotate a matrix for printing
#' 
#' Rotates a matrix so that when view is called the rows and columns appear in the 
#' same order as when looking at the matrix with print 
#' 
#' @export
#' @param m matrix
#' @return matrix that has been transposed and the columns reversed
#' 
#' @examples
#' 
#' # Normally image shows a matrix with the first entry in the bottom left
#' image((matrix(1:12, nrow=3, ncol=4))
#' 
#' # With rot the image is shown in the same order as print
#' image(rot(matrix(1:12, nrow=3, ncol=4)))

rot <- function(m) t(m)[,nrow(m):1,drop=F]
