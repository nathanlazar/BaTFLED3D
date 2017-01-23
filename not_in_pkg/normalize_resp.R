#' z-scale responses by params$norm.mode and store means and sds as attributes
#' 
#' @param resp.tens array
#' @param params params object
#' 
#' @return normalized array

normalize_resp <- function(resp.tens, params) {
  # If params$norm.mode is not null, normalize responses for training
  # TODO: move this to inside making a input.data object and 
  #       make an object for the input data before splitting?
  # OR: Make a normalize and split function?

  attr(resp.tens, 'm1.means') <- apply(resp.tens, c(2,3), mean, na.rm=T)
  attr(resp.tens, 'm1.sds')   <- apply(resp.tens, c(2,3), sd, na.rm=T)
  attr(resp.tens, 'm2.means') <- apply(resp.tens, c(1,3), mean, na.rm=T)
  attr(resp.tens, 'm2.sds')   <- apply(resp.tens, c(1,3), sd, na.rm=T)
  attr(resp.tens, 'm3.means') <- apply(resp.tens, c(1,2), mean, na.rm=T)
  attr(resp.tens, 'm3.sds')   <- apply(resp.tens, c(1,2), sd, na.rm=T)

  if(params$verbose) print(paste('** Normalizing responses to train for mode', 
                                 params$norm.mode, ' **'))
  # Normalize response values so they have mean 0 and sd 1
  if(params$norm.mode == 0) {
    resp.tens <- resp.tens - mean(resp.tens, na.rm=T)
    resp.tens <- resp.tens / sd(resp.tens, na.rm=T)
    attr(resp.tens, 'm1.means') <- matrix(0, dim(resp.tens)[2], dim(resp.tens)[3])
    attr(resp.tens, 'm1.sds') <- apply(resp.tens, c(2,3), sd, na.rm=T)
  } else if(params$norm.mode == 1) {
    resp.tens <- sweep(resp.tens, c(2,3), m1.means, FUN='-')
    resp.tens <- sweep(resp.tens, c(2,3), m1.sds, FUN='/')
    attr(resp.tens, 'm2.means') <- apply(resp.tens, c(2,3), mean, na.rm=T)
    attr(resp.tens, 'm2.sds') <- apply(resp.tens, c(2,3), sd, na.rm=T)
  } else if(params$norm.mode == 2) {
    resp.tens <- sweep(resp.tens, c(1,3), m2.means, FUN='-')
    resp.tens <- sweep(resp.tens, c(1,3), m2.sds, FUN='/')
    attr(resp.tens, 'm3.means') <- apply(resp.tens, c(1,3), mean, na.rm=T)
    attr(resp.tens, 'm3.sds') <- apply(resp.tens, c(1,3), sd, na.rm=T)
  } else if(params$norm.mode == 3) {
    resp.tens <- sweep(resp.tens, c(1,2), m3.means, FUN='-')
    resp.tens <- sweep(resp.tens, c(1,2), m3.sds, FUN='/')
    m3.means <- apply(resp.tens, c(1,2), mean, na.rm=T)
    m3.sds <- apply(resp.tens, c(1,2), sd, na.rm=T)
  }
  return(resp.tens)
}