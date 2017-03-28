#' Updates the root mean squared error for training data.  Predicting both from data
#' and from just the latent (H) matrices.
#' 
#' @export
#' @param m training object
#' @param d data object
#' @param verbose Logical indicating whether to print the results (TRUE)
#' 
#' @return numeric value of the explained variance
#' 

rmse <- function(m, d, verbose=T) {
  m$RMSE[m$iter] <- nrmse(d$resp, m$resp)
  
  if(length(m$core.mean)) { # Tucker model
    H.resp <- mult_3d(m$core.mean, m$mode1.H.mean, m$mode2.H.mean, m$mode3.H.mean)
    m$H.RMSE[m$iter] <- nrmse(d$resp, H.resp)
  } else {  # CP model
    R <- ncol(m$mode1.H.mean)
    core <- array(0, dim=c(R,R,R))
    for(r in 1:R) core[r,r,r] <- 1
    H.resp <- mult_3d(core, m$mode1.H.mean, m$mode2.H.mean, m$mode3.H.mean)
    m$H.RMSE[m$iter] <- nrmse(d$resp, H.resp)
  }
  if(verbose) {
    print(sprintf("A RMSE: %.4f, H RMSE: %.4f", m$RMSE[m$iter], m$H.RMSE[m$iter]))
  }
  m$RMSE[m$iter]
}
