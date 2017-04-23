#' Plot correlation results from test data
#' 
#' @importFrom graphics plot rect points par
#' 
#' @export
#' @param test.results results \code{generated with test_results}
#' @param baselines named vector of baseline values to draw as dotted horizontal 
#' lines e.g. c('warm'=0, 'm1'=0, 'm1m2'=0, 'm1m2m3'=0)
#' @param ylim limits for the y-axis (NA)
#' @param main Main title of the plot
#' @param method Either 'pearson' or 'spearman' correlations

plot_test_cor <- function(test.results, ylim='default', main=NA, method='pearson',
  baselines=c('warm'=NA, 'm1'=NA, 'm2'=NA, 'm3'=NA, 'm1m2'=NA, 'm1m3'=NA, 'm2m3'=NA, 'm1m2m3'=NA)) {

  # Subset to just correlation data data
  if(method=='pearson') type='p.cor'
  if(method=='spearman') type='s.cor'
  sub.results <- test.results[,grepl(type, names(test.results))]
  types <- paste0(c('warm', 'm1', 'm2', 'm3', 'm1m2', 'm1m3', 'm2m3', 'm1m2m3'), '.' ,type)
  sub.results <- sub.results[,types[types %in% names(test.results)]]
  
  if(ylim=='default') { 
    ylim=c(-1, 1)
  } else if(is.na(ylim)) ylim=range(sub.results[-1,], na.rm=T)
  
  if(method=='pearson')
    types <- sub('.p.cor', '', names(sub.results))
  if(method=='spearman')
    types <- sub('.s.cor', '', names(sub.results))
  types <- sub('^m', 'Mode ', types)
  types <- gsub('m', ' & ', types)
  types <- sub('war & ', 'Warm', types)
  types <- sub('1 & 2 & 3', '1,2 & 3', types)
  
  n <- length(types)
  colrs <- RColorBrewer::brewer.pal(n, "Set1")

  if(method=='pearson') ylab='Pearson correlation'
  if(method=='spearman') ylab='Spearman correlation'

  plot(c(0, nrow(sub.results)), ylim, xlab='Iteration', 
       ylab=ylab, type='n', main=main)
  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "#ededed")
  for(i in 1:n) points(sub.results[,i], col=colrs[i], pch=20, type='b')
  for(i in 1:n) abline(h=baselines[i], col=colrs[i], lty=2, lwd=2)
  
  if(sum(!is.na(baselines))) {
    legend <- c(types, 'predicting the mean')
    colrs = c(colrs, 'black')
    pch=c(rep(20,n),45)
  } else {
    legend <- types
    pch=rep(20,n)
  }
  legend(x='bottomright', bty='n', legend=legend, col=colrs, pch=pch, pt.cex=2)
}

