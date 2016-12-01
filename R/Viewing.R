# some functions for viewing matrices
# TODO: Split these into individual files w/ help text


plot_gamma <- function(shape, scale, ...) {
  plot(1:3000/1000, dgamma(1:3000/1000, shape=shape, scale=scale), type='l', ...)
}

plot_left_out <- function(d, m, resp, ...) {
  if(sum(is.na(d@resp)) > 0) {
    actual <- resp[is.na(d@resp)]
    predicted <- m@resp[is.na(d@resp)]
    index <- sample(length(actual))
    plot(actual[index], pch=20, col='blue', ylim=range(actual, predicted, na.rm=T),
         xlab='', ylab='Response',...)
    points(predicted[index], pch=20, col='red')
    segments(1:length(actual), actual[index], 1:length(actual), predicted[index])
    abline(h=mean(actual), col='red')
    abline(h=mean(predicted), col='blue')
  } else {
    print("No missing data to predict")
  }
}

plot_left_out_curves <- function(d, m, resp, n=0, ...) {
  if(sum(is.na(d@resp)) > 0) {
    left.out <- which(is.na(d@resp[,,1]), arr.ind=T)
    if(n>0) {
      sampled <- sample(nrow(left.out), n)
    } else sampled <- 1:nrow(left.out)
    for(i in sampled) { 
      actual <- resp[left.out[i,1], left.out[i,2],]
      predicted <- m@resp[left.out[i,1], left.out[i,2],]
      plot(actual, col='blue', pch=20, ylim=range(actual, predicted), 
           xlab='Concentration', ylab='Mean response', 
           main=paste('Cell line:', dimnames(resp)[[1]][left.out[i,1]], 
                      '\nCompound:', dimnames(resp)[[2]][left.out[i,2]]),...)
      points(predicted, col='red', pch=20)
    }
  } else {
    print("No missing data to predict")
  }
}

plot_dr <- function(obs, pred, dose=0:9, cl="", dr="") {
  # Plots observed (blue) and predicted (red) responses 
  
  df <- data.frame(od=c(obs, pred), dose=dose, 
                   type=c(rep('Observed', 10), rep('Predicted', 10)),
                   cl=cl, dr=dr)
  
  p <- ggplot(df, aes(x=dose, y=od, colour=type)) + 
    geom_point(size=8) +
    ggtitle(paste("Cell line:", cl, "\nDrug:", dr)) +
    labs(x="Dose", y="Response") +
    theme_bw() +
    theme_classic() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.title = element_text(size=50, lineheight=.8, 
                                    face="bold", vjust=4),
          axis.text.x  = element_text(size=30),
          axis.text.y  = element_text(size=30),
          axis.title.x = element_text(size=35, vjust=-2.5),
          axis.title.y = element_text(size=35, vjust=3.5),
          plot.margin = unit(c(1, 1, 2, 2), "cm"),
          legend.title = element_blank(),
          legend.text = element_text(size = 25)) +
    scale_x_continuous(limits = c(0, 9), breaks=0:9) +
    scale_colour_manual(values=c("blue", "red"))
  p
}

## Function for arranging ggplots. use png(); arrange(p1, p2, ncol=1); dev.off() to save.
# From http://www.gettinggeneticsdone.com/2010/03/arrange-multiple-ggplot2-plots-in-same.html
require(grid)
vp.layout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)
arrange_ggplot2 <- function(..., nrow=NULL, ncol=NULL, as.table=FALSE) {
  dots <- list(...)
  n <- length(dots)
  if(is.null(nrow) & is.null(ncol)) { nrow = floor(n/2) ; ncol = ceiling(n/nrow)}
  if(is.null(nrow)) { nrow = ceiling(n/ncol)}
  if(is.null(ncol)) { ncol = ceiling(n/nrow)}
  ## NOTE see n2mfrow in grDevices for possible alternative
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(nrow,ncol) ) )
  ii.p <- 1
  for(ii.row in seq(1, nrow)){
    ii.table.row <- ii.row	
    if(as.table) {ii.table.row <- nrow - ii.table.row + 1}
    for(ii.col in seq(1, ncol)){
      ii.table <- ii.p
      if(ii.p > n) break
      print(dots[[ii.table]], vp=vp.layout(ii.table.row, ii.col))
      ii.p <- ii.p + 1
    }
  }
}

resp_density <- function(obs, pred, alpha=0.2) {
  df <- data.frame(Observed=as.vector(obs), Predicted=as.vector(pred))
  p <- ggplot(df,aes(x=Observed, y=Predicted)) + 
    geom_point(colour="blue", alpha=alpha) + 
    geom_density2d(colour="black") + 
    theme_bw() +
    theme_classic() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          plot.title = element_text(size=50, lineheight=.8, 
                                    face="bold", vjust=4),
          axis.text.x  = element_text(size=30),
          axis.text.y  = element_text(size=30),
          axis.title.x = element_text(size=35, vjust=-2.5),
          axis.title.y = element_text(size=35, vjust=3.5),
          plot.margin = unit(c(.5, .5, 1.5, 1.5), "cm"),
          legend.title = element_blank(),
          legend.text = element_text(size = 25)) +
    scale_x_continuous(limits = c(-50, 80), breaks=seq(-50,80,20)) +
    scale_y_continuous(limits = c(-40, 50), breaks=seq(-30,50,20)) +
    geom_segment(aes(x=-40, xend=50, y=-40, yend=50), size=1.5)
  p
}

plot_RMSE <- function(RMSE) {
  df <- data.frame(RMSE=RMSE, Iteration=1:length(RMSE))
  p1 <- ggplot(df,aes(x=Iteration, y=RMSE)) +
    geom_point(size=6) +
    theme_bw() +
    theme_classic() +
    theme(axis.text.x  = element_text(size=30),
          axis.text.y  = element_text(size=30),
          axis.title.x = element_text(size=35, vjust=-2.5),
          axis.title.y = element_text(size=35, vjust=3.5),
          plot.margin = unit(c(0.1, 0.1, 1.3, 1.3), "cm"),
          legend.title = element_blank())
  df2 <- df[10:nrow(df),]
  p2 <- ggplot(df2,aes(x=Iteration, y=RMSE)) +
    geom_point(size=6) +
    theme_bw() +
    theme_classic() +
    theme(axis.text.x  = element_text(size=30),
          axis.text.y  = element_text(size=30),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.margin = unit(c(0.1, 0.1, 0, 0), "cm"),
          legend.title = element_blank())
  subvp <- viewport(width = 0.76, height = 0.76, x = 1, y = 1, just=c(1,1))
  print(p1)
  print(p2, vp = subvp)
}

plot_times <- function(times) {
  df <- data.frame(Times=times, Iteration=1:length(times))
  p <- ggplot(df,aes(x=Iteration, y=Times)) +
    geom_point(size=6) +
    theme_bw() +
    theme_classic() +
    theme(axis.text.x  = element_text(size=30),
          axis.text.y  = element_text(size=30),
          axis.title.x = element_text(size=35, vjust=-2.5),
          axis.title.y = element_text(size=35, vjust=3.5),
          plot.margin = unit(c(0.1, 0.1, 1.3, 1.3), "cm"),
          legend.title = element_blank())
  print(p)
}

is_inc <- function(vec) vec[-1]>vec[-length(vec)]
