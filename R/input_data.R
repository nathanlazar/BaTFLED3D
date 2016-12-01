#' Object storing input data for BaTFLED algorithm with 3-D response tensor.
#' 
#' @export
#' @slot mode1.X matrix of predictors for mode 1
#' @slot mode2.X matrix of predictors for mode 2
#' @slot mode3.X matrix of predictors for mode 3
#' @slot resp three dimensional array of responses with dimensions matching the number 
#' of rows in mode1.X, mode2.X and mode3.X
#' @examples
#' a <- input_data$new(mode1.X = matrix(rnorm(30), nrow=3, ncol=10),
#'                     mode2.X = matrix(rnorm(36), nrow=4, ncol=9), 
#'                     mode3.X = matrix(rnorm(40), nrow=5, ncol=8),
#'                     resp = array(rnorm(60), dim=c(3,4,5)))
#' im_mat(a$mode1.X)
#' im_mat(a$mode2.X)
#' im_mat(a$mode3.X)
#' im_mat(a$resp[,,1])

# Class definition for Tensor factorization run. The input data is not stored 
# here to avoid it being repeated unnecessarily.
input_data <- R6Class("input_data",
  portable=F,
  class=F,
  cloneable=T,
  public = list(
    mode1.X = matrix(0, 2, 2),
    mode2.X = matrix(0, 2, 2), 
    mode3.X = matrix(0, 2, 2), 
    resp = array(0, dim=c(2,2,2)),
    delta = array(0, dim=c(2,2,2)),
    initialize = function(mode1.X, mode2.X, mode3.X, resp=NA) {
      mode1.X <<- mode1.X
      mode2.X <<- mode2.X
      mode3.X <<- mode3.X
      resp <<- resp
      delta <<- (!is.na(resp)) * 1 
    }
  )
)
