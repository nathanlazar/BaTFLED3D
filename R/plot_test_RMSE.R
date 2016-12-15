#' Plot RMSE results from test data
#'
#' @export
#' @param test.results An object created with \code{test_results}
#' @param ylim Limits of the y-axis
#' @param mean Logical to print horizontal lines for mean predictions
#' @param main Main title of the plot

plot_test_RMSE <- function(test.results, ylim=NA, mean=F, main='Test RMSEs') {
  # Sort by row names so colors are consistent for different response measures
  test.results <- test.results[,sort(names(test.results))]

  if(is.na(ylim)) 
    ylim=c(0,max(test.results[-1,grepl('RMSE', names(test.results))], na.rm=T))

  rmse.cols <- grep('RMSE$', names(test.results))
  types <- sub('.RMSE', '', names(test.results)[rmse.cols])
  types <- sub('^m', 'Mode ', types)
  types <- gsub('m', ' & ', types)
  types <- sub('war & ', 'Warm', types)
  types <- sub('1 & 2 & 3', '1,2 & 3', types)
  
  n <- length(types)
  colrs <- RColorBrewer::brewer.pal(n, "Set1")

  if(mean) {
    legend <- c(types, 'predicting the mean')
    colrs = c(colrs, 'black')
    pch=c(rep(20,4),45)
  } else {
    legend <- types
    pch=rep(20,4)
  }
  
  # if(is.na(main)) {
  #   main <- sprintf('Warm min at %.0f, cold mode 1 min at %.0f \n Cold mode 2 min at %.0f Both cold min at %.0f',
  #                   which.min(test.results$warm.RMSE), which.min(test.results$m1.RMSE),
  #                   which.min(test.results$m2.RMSE), which.min(test.results$m1m2.RMSE))
  # }
  
  plot(c(0, nrow(test.results)), ylim, xlab='Iteration', ylab='RMSE', type='n', main=main)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#ededed")
  for(i in 1:n) {
    points(test.results[,rmse.cols[i]], col=colrs[i], pch=20, type='b')
  }
  legend(x='topright', bty='n', legend=legend, col=colrs, pch=pch, pt.cex=2)
}

