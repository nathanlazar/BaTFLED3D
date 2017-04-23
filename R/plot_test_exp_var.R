#' Plot explained variance results from test data
#'
#' @importFrom graphics plot rect points par
#'
#' @export
#' @param test.results an object generated with \code{test_results}
#' @param baselines named vector of baseline values to draw as dotted horizontal 
#' lines e.g. c('warm'=0, 'm1'=0, 'm1m2'=0, 'm1m2m3'=0)
#' @param ylim Limits of the y-axis.
#' @param main Main title of the plot

plot_test_exp_var <- function(test.rersults, ylim='default', main=NA,
  baselines=c('warm'=NA, 'm1'=NA, 'm2'=NA, 'm3'=NA, 'm1m2'=NA, 'm1m3'=NA, 'm2m3'=NA, 'm1m2m3'=NA)) {
  # Subset to just explained variance data
  sub.results <- test.results[,grepl('exp.var', names(test.results))]
  types <- paste0(c('warm', 'm1', 'm2', 'm3', 'm1m2', 'm1m3', 'm2m3', 'm1m2m3'), '.exp.var')
  sub.results <- sub.results[,types[types %in% names(test.results)]]

  if(ylim=='default') {
    ylim <- c(-1,1)
  } else if(ylim=='auto') 
    ylim <- range(sub.results, na.rm=T)

  types <- sub('.exp.var', '', names(sub.results))
  types <- sub('^m', 'Mode ', types)
  types <- gsub('m', ' & ', types)
  types <- sub('war & ', 'Warm', types)
  types <- sub('1 & 2 & 3', '1,2 & 3', types)

  n <- length(types)
  colrs <- RColorBrewer::brewer.pal(n, "Set1")
  
  plot(c(0, nrow(sub.results)), ylim, xlab='Iteration', 
       ylab='Explained variance', type='n', main=main)
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

