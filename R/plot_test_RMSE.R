#' Plot RMSE results from test data
#'
#' @importFrom graphics plot rect points par
#'
#' @export
#' @param test.results An object created with \code{test_results}
#' @param baselines named vector of baseline values to draw as dotted horizontal 
#' lines e.g. c('warm'=0, 'm1'=0, 'm1m2'=0, 'm1m2m3'=0)
#' @param ylim Limits of the y-axis (default is (0, 1.5))
#' @param main Main title of the plot

plot_test_RMSE <- function(test.results, ylim='default', main='Test RMSEs',
  baselines=c('warm'=NA, 'm1'=NA, 'm2'=NA, 'm3'=NA, 'm1m2'=NA, 'm1m3'=NA, 
              'm2m3'=NA, 'm1m2m3'=NA)) {
  # Sort by row names so colors are consistent for different response measures
  sub.results <- test.results[,grepl('RMSE', names(test.results))]
  types <- paste0(c('warm', 'm1', 'm2', 'm3', 'm1m2', 'm1m3', 'm2m3', 'm1m2m3'), '.RMSE')
  sub.results <- sub.results[,types[types %in% names(test.results)]]

  bl2 <- rep(NA, length(sub.results))
  names(bl2) <- sub('.RMSE', '', names(sub.results))
  for(i in 1:length(baselines))
    bl2[names(baselines)[i]] <- baselines[i]

  if(typeof(ylim)=='character' && ylim=='default') {
    ylim <- c(0, 1.5)
  } else if(is.na(ylim[1])) ylim=c(0,max(sub.results, na.rm=T))

  types <- sub('.RMSE', '', names(sub.results))
  types <- sub('^m', 'Mode ', types)
  types <- gsub('m', ' & ', types)
  types <- sub('war & ', 'Warm', types)
  types <- sub('1 & 2 & 3', '1,2 & 3', types)
  
  n <- length(types)
    
  colrs <- RColorBrewer::brewer.pal(n, "Set1")

  plot(c(0, nrow(sub.results)), ylim, xlab='Iteration', ylab='RMSE', type='n', main=main)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#ededed")
  for(i in 1:n) {
    points(sub.results[,i], col=colrs[i], pch=20, type='b')
  }

  for(i in 1:n)
    abline(h=bl2[i], col=colrs[i], lty=2, lwd=2)

  keep <- apply(sub.results, 2, function(x) sum(!is.na(x)))!=0
  types <- types[keep]
  colrs <- colrs[keep]
  n <- sum(keep)
  
  if(sum(!is.na(bl2))) {
    legend <- c(types, 'predicting the mean')
    colrs = c(colrs, 'black')
    pch=c(rep(20,n),45)
  } else {
    legend <- types
    pch=rep(20,n)
  }
  
  legend(x='topright', bty='n', legend=legend, col=colrs, pch=pch, pt.cex=2)
  
  # if(is.na(main)) {
  #   main <- sprintf('Warm min at %.0f, cold mode 1 min at %.0f \n Cold mode 2 min at %.0f Both cold min at %.0f',
  #                   which.min(sub.results$warm.RMSE), which.min(sub.results$m1.RMSE),
  #                   which.min(sub.results$m2.RMSE), which.min(sub.results$m1m2.RMSE))
  # }

}

