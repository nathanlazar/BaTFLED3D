#' Multiply three matrices (or vectors) through a given core tensor to form a 
#' three dimensional tensor.
#' 
#' The package 'rTensor' is required and the number of columns of x, y and z must 
#' match the dimensions of core.
#' 
#' @export
#' @param core array
#' @param x matrix to multiply by the first mode of \code{core}
#' @param y matrix to multiply by the second mode of \code{core}
#' @param z matrix to multiply by the third mode of \code{core}
#' @param names 
#' 
#' @return Array with sizes given by the number of rows in x, y and z
#'   
#' @examples
#' mult_3d(array(1:24, dim=c(2,3,4)), matrix(1:4,2,2), matrix(1:6,2,3), matrix(1:8,2,4))

mult_3d <- function(core, x, y, z, names=T) {
  # Multiply three matrices through a core tensor
  # Requires rTensor library
<<<<<<< HEAD
  # resp <- rTensor::ttm(rTensor::ttm(rTensor::ttm(rTensor::as.tensor(core), x, 1), y, 2), z, 3)@data
  
  # If one of the dimensions of the core is degenerate, return a zero tensor of the appropriate size
  if(sum(dim(core)==0)) {
    resp <- array(0, dim=c(nrow(x), nrow(y), nrow(z)))
  } else {
    resp <- rTensor::ttl(rTensor::as.tensor(core), list(x, y, z), c(1,2,3))@data
  }
  
  if(names) dimnames(resp) <- list(rownames(x), rownames(y), rownames(z))
=======
  resp <-  rTensor::ttm(rTensor::ttm(rTensor::ttm(rTensor::as.tensor(core), x, 1), y, 2), z, 3)
>>>>>>> b452907b1fbd5d5edc1cccb97e72c2c1d3cebb36

  return(resp)
}
