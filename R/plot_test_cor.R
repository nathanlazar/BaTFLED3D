#' Plot correlation results from test data
#' 
#' @export
#' @param test.results results \code{generated with test_results}
#' @param ylim limits for the y-axis (NA)
#' @param mean Logical indicating whether to plot horizontal lines for mean predictions
#' @param main Main title of the plot
#' @param method Either 'pearson' or 'spearman' correlations

plot_test_cor <- function(test.results, ylim=NA, mean=F, main=NA, method='pearson') {
  # Sort by row names so colors are consistent for different response measures
  test.results <- test.results[,sort(names(test.results))]

  if(method=='pearson')
    cor.cols <- grep('p.cor$', names(test.results))
  if(method=='spearman')
    cor.cols <- grep('s.cor$', names(test.results))

  if(is.na(ylim))
    ylim=range(test.results[-1,cor.cols], na.rm=T)
  
  if(method=='pearson')
    types <- sub('.p.cor', '', names(test.results)[cor.cols])
  if(method=='spearman')
    types <- sub('.s.cor', '', names(test.results)[cor.cols])
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
  #   if(method=='pearson')
  #     main <- sprintf('Warm max at %.0f, cold mode 1 max at %d \n Cold mode 2 max at %d Both cold max at %d',
  #                     which.min(test.results$warm.p.cor), which.min(test.results$m1.p.cor),
  #                     which.min(test.results$m2.p.cor), which.min(test.results$m1m2.p.cor))
  #   if(method=='spearman')
  #     main <- sprintf('Warm max at %.0f, cold mode 1 max at %d \n Cold mode 2 max at %d Both cold max at %d',
  #                     which.min(test.results$warm.s.cor), which.min(test.results$m1.s.cor),
  #                     which.min(test.results$m2.s.cor), which.min(test.results$m1m2.s.cor))
  # }
  
  if(method=='pearson') ylab='Pearson correlation'
  if(method=='spearman') ylab='Spearman correlation'

  plot(c(0, nrow(test.results)), ylim, xlab='Iteration', 
       ylab=ylab, type='n', main=main)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#ededed")
  for(i in 1:n) 
    points(test.results[,cor.cols[i]], col=colrs[i], pch=20, type='b')
  legend(x='bottomright', bty='n', legend=legend, col=colrs, pch=pch, pt.cex=2)
}

