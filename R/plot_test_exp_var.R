#' Plot explained variance results from test data
#'
#' @importFrom graphics plot rect points par
#'
#' @export
#' @param test.results an object generated with \code{test_results}
#' @param ylim Limits of the y-axis.
#' @param mean Logical to plot horizontal lines for mean prediction (FALSE)
#' @param main Main title of the plot

plot_test_exp_var <- function(test.results, ylim='default', mean=F, main=NA) {
  # Sort by row names so colors are consistent for different response measures
  test.results <- test.results[,sort(names(test.results))]
  
  exp.var.cols <- grep('exp.var$', names(test.results))

  if(ylim=='default')  ylim=c(-1,1)
  if(is.na(ylim)) 
    ylim=range(test.results[-1,exp.var.cols], na.rm=T)

  types <- sub('.exp.var', '', names(test.results)[exp.var.cols])
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
  #   main <- sprintf('Warm max at %.0f, cold mode 1 max at %d \n Cold mode 2 max at %d Both cold max at %d',
  #                   which.min(test.results$warm.exp.var), which.min(test.results$m1.exp.var),
  #                   which.min(test.results$m2.exp.var), which.min(test.results$m1m2.exp.var))
  # }
  
  plot(c(0, nrow(test.results)), ylim, xlab='Iteration', 
       ylab='Explained variance', type='n', main=main)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#ededed")
  for(i in 1:n) 
    points(test.results[,exp.var.cols[i]], col=colrs[i], pch=20, type='b')
  legend(x='bottomright', bty='n', legend=legend, col=colrs, pch=pch, pt.cex=2)
}

