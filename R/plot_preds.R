#' Make a scatterplot of observed vs. predicted values
#' 
#' If there are more than 25,0000 points then they are subsampled down to 25,000.
#' Observed values are on the x axis predicted values on the y. A blue line shows the 
#' diagonal. Points are transparent to show dense clusters. Predictions for points where
#' the true value is not known are plotted at zero in blue.
#' 
#' @export
#' @param pred matrix or vector of predicted values
#' @param true matrix or vector of predicted values
#' @param show.na logical, display NA values as blue dots at the mean for the x or y axis (def: T)
#' @param ... other parameters passed to plot
#' 
#' @return none
#' 
#' @examples
#' x <- seq(-10,10, 0.01)+rnorm(2001)
#' y <- seq(-10,10, 0.01)+rnorm(2001)
#' x[sample(2001, 100)] <- NA
#' plot_preds(y, x)

plot_preds <- function(pred, true, show.na=T, ...) {
  # If there are over 25,000 points, sample them down
  if(prod(dim(pred)) > 25000) {
    subs <- sample(prod(dim(pred)), 25000)
    pred <- pred[subs]
    true <- true[subs]
  }
  
  plot(pred ~ true, pch=20, col=rgb(0,0,0,0.1),
       xlab='Observed responses', ylab='Predicted responses', ...)
  abline(a=0,b=1, lwd=2, col='blue')
  # Add horizontal lines for predictions of NA values
  if(show.na) {
    points(rep(mean(true, na.rm=T), sum(is.na(true))), pred[is.na(true)], 
           pch=20, col=rgb(0,0,1,0.1))
  }
}
