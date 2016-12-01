#' Plot matrices from a model object with im_mat
#'  
#' @param m model object created with mk_model
#' @param d input data object created with get_input_data
#' @param show.mode vector of modes whose projection and latent matrices are 
#' to be displayed
#' @param scale Logical should the columns of matrices be scaled

# Show A and H matrices
show_mat <- function(m, d, show.mode, scale=F) {
  if(1 %in% show.mode) {
    if(nrow(m$mode1.A.mean)!=0) {
      im_mat(m$mode1.A.mean, main="Mode 1 A", scale=scale,
             ylab=paste(nrow(m$mode1.A.mean), 'predictors'))
    }
    im_mat(m$mode1.H.mean, main="Mode 1 H")
    #    im_2_mat(m$mode1.A.mean, toy$mode1.A, scale=scale)
    #    im_2_mat(m$mode1.H.mean, toy$mode1.H, scale=scale)
  }
  if(2 %in% show.mode) {
    if(nrow(m$mode2.A.mean)!=0) {
      im_mat(m$mode2.A.mean, main="Mode 2 A", scale=scale,
             ylab=paste(nrow(m$mode2.A.mean), 'predictors'))
    }
    im_mat(m$mode2.H.mean, main="Mode 2 H")
    #     im_2_mat(m$mode2.A.mean, toy$mode2.A, scale=scale)
    #     im_2_mat(m$mode2.H.mean, toy$mode2.H, scale=scale)
  }
  if(3 %in% show.mode) {
    if(nrow(m$mode3.A.mean)!=0) {
      im_mat(m$mode3.A.mean, main="Mode 3 A", scale=scale,
             ylab=paste(nrow(m$mode3.A.mean), 'predictors'))
    }
    im_mat(m$mode3.H.mean, main="Mode 3 H")
    #     im_2_mat(m$mode3.A.mean, toy$mode3.A, scale=scale)
    #     im_2_mat(m$mode3.H.mean, toy$mode3.H, scale=scale)
  }
}
