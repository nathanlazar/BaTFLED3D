#' Plot a heatmap of a matrix in red and blue
#' 
#' Displays a heatmap of a matrix using red and blue colors. Options to scale
#' and sort as well as any other graphical parameters with ...
#' 
#' @export
#' @param x matrix
#' @param high string of either 'red' or 'blue' used to show higher values
#' @param xaxt string indicating how to display the x axis. Suppress x axis with 'n'
#' @param yaxt string indicating how to display the y axis. Suppress y axis with 'n'
#' @param sort logical indicating whether the columns of the matrix should 
#' be sorted in decreasing order of their means
#' @param scale logical indicating whether the matrix should be z scaled to have 
#' columns with norm zero and standard deviation one.
#' @param ballance logical indicating whether to expand the range so it stays centered at zero
#' @param zlim numeric bounds on the max and min range for colors.
#' @param ... other graphical parameters passed to image
#' 
#' @return none
#' 
#' @examples
#' im_mat(matrix(1:12, nrow=3, ncol=4), sort=FALSE, scale=FALSE)
#' im_mat(matrix(1:12, nrow=3, ncol=4), sort=TRUE, scale=FALSE)
#' im_mat(matrix(1:12, nrow=3, ncol=4), sort=FALSE, scale=TRUE)
#' im_mat(matrix(1:12, nrow=3, ncol=4), sort=TRUE, scale=TRUE)

im_mat <- function(x, high='red', xaxt='n', yaxt='n', sort=FALSE, scale=FALSE, ballance=FALSE, zlim=NA, ...) {
  # display image of rotated matrix with the option to sort columns by 
  # absolute magnitude
  rwb <- colorRampPalette(c('red', 'white', 'blue'), space='rgb')
  
  if(scale & nrow(x)>1) {
    # Safe version of scale (if sd of columns = 0, set to minimum value)
    x.col.means <- apply(x, 2, mean)
    x.col.sds <- apply(x, 2, sd)
    x.col.sds[x.col.sds==0] <- 1e-300
    x <- sweep(x, MARGIN=2, STATS=x.col.means, FUN='-')
    x <- sweep(x, MARGIN=2, STATS=x.col.sds, FUN='/')
  }
  
  # Rearrange columns sorted from largest mean value to smallest
  if(sort==TRUE) {
    sorted <- sort(apply(x, 2, mean), decreasing=TRUE, index.return=TRUE)
    x <- x[,sorted$ix,drop=FALSE]
  }

  # Expand range if ballance = TRUE
  if(length(zlim) != 2) {
    if(ballance) {
      zlim=c(-max(abs(range(x, na.rm=TRUE))), max(abs(range(x, na.rm=TRUE))))
    } else zlim=range(x, na.rm=TRUE)
  }

  if(high=='blue') 
    image(rot(x), col=rwb(256), xaxt=xaxt, yaxt=xaxt, zlim=zlim, ...)
  else if(high=='red')
    image(rot(x), col=rev(rwb(256)), xaxt=xaxt, yaxt=xaxt, zlim=zlim, ...)
}  
