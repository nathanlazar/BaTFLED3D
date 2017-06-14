# Make plots of predictions comparing BaTFLED to mean for each sample

# Usage: head_to_head_mean.R <run_prefix>

# Example: head_to_head.R  run_2335512/run_2335512
# run_2324578/run_2324578

library(ggplot2)
library(BaTFLED3D)

args <- commandArgs(TRUE)

run.prefix <- args[1]
# run.prefix <- 'run_2335512/run_2335512'
# run.prefix <- 'run_2324578/run_2324578'

# see how many runs there are
n.files <- length(list.files(path = dirname(run.prefix),
                  pattern = paste0(basename(run.prefix), '.[0-9]+.out')))

# Make list to store results as long form data frames
sum.res <- list()
for(task in c('warm', 'm1', 'm2', 'm3', 'm1m2', 'm1m3', 'm2m3', 'm1m2m3'))
  if(exists(paste0(task, '.resp')))
    sum.res[[task]] <- data.frame(sample=character(0), model=character(0), 
                                  measure=character(0), value=numeric(0))

# Load in output from one run and get predictions
load_run <- function(run, sum.res) {
  load(run)
  
  if(exists('m1.resp')) {
    for(i in 1:dim(m1.resp)[[1]]) {
      samp <- dimnames(m1.resp)[[1]][i]
      sum.res[['m1']] <- rbind(sum.res[['m1']], 
        data.frame(sample=rep(samp, 4), model=c(rep('BaTFLED', 4), rep('mean', 4)), 
                   measure=rep(c('RMSE', 'exp.var', 'p.cor', 's.cor'), 2), 
                   value=c(nrmse(m1.resp[i,,], m1.pred.resp[i,,]),
                           exp_var(m1.resp[i,,], m1.pred.resp[i,,]),
                           p_cor(as.vector(m1.resp[i,,]), as.vector(m1.pred.resp[i,,])),
                           s_cor(as.vector(m1.resp[i,,]), as.vector(m1.pred.resp[i,,])),
                           nrmse(m1.resp[i,,], mean.tens.list[['m1']][i,,]),
                           exp_var(m1.resp[i,,], mean.tens.list[['m1']][i,,]),
                           p_cor(as.vector(m1.resp[i,,]), as.vector(mean.tens.list[['m1']][i,,])),
                           s_cor(as.vector(m1.resp[i,,]), as.vector(mean.tens.list[['m1']][i,,])))))
    }
  }

  if(exists('m2.resp')) {
    for(j in 1:dim(m2.resp)[[2]]) {
      samp <- dimnames(m2.resp)[[2]][j]
      sum.res[['m2']] <- rbind(sum.res[['m2']], 
        data.frame(sample=rep(samp, 4), model=c(rep('BaTFLED', 4), rep('mean', 4)), 
                   measure=rep(c('RMSE', 'exp.var', 'p.cor', 's.cor'), 2), 
                   value=c(nrmse(m2.resp[,j,], m2.pred.resp[,j,]),
                           exp_var(m2.resp[,j,], m2.pred.resp[,j,]),
                           p_cor(as.vector(m2.resp[,j,]), as.vector(m2.pred.resp[,j,])),
                           s_cor(as.vector(m2.resp[,j,]), as.vector(m2.pred.resp[,j,])),
                           nrmse(m2.resp[,j,], mean.tens.list[['m2']][,j,]),
                           exp_var(m2.resp[,j,], mean.tens.list[['m2']][,j,]),
                           p_cor(as.vector(m2.resp[,j,]), as.vector(mean.tens.list[['m2']][,j,])),
                           s_cor(as.vector(m2.resp[,j,]), as.vector(mean.tens.list[['m2']][,j,])))))
    }
  }

  if(exists('m3.resp')) {
    for(k in 1:dim(m3.resp)[[3]]) {
      samp <- dimnames(m3.resp)[[3]][k]
      sum.res[['m3']] <- rbind(sum.res[['m3']], 
         data.frame(sample=rep(samp, 4), model=c(rep('BaTFLED', 4), rep('mean', 4)), 
                    measure=rep(c('RMSE', 'exp.var', 'p.cor', 's.cor'), 2), 
                    value=c(nrmse(m3.resp[,,k], m3.pred.resp[,,k]),
                            exp_var(m3.resp[,,k], m3.pred.resp[,,k]),
                            p_cor(as.vector(m3.resp[,,k]), as.vector(m3.pred.resp[,,k])),
                            s_cor(as.vector(m3.resp[,,k]), as.vector(m3.pred.resp[,,k])),
                            nrmse(m3.resp[,,k], mean.tens.list[['m3']][,,k]),
                            exp_var(m3.resp[,,k], mean.tens.list[['m3']][,,k]),
                            p_cor(as.vector(m3.resp[,,k]), as.vector(mean.tens.list[['m3']][,,k])),
                            s_cor(as.vector(m3.resp[,,k]), as.vector(mean.tens.list[['m3']][,,k])))))
    }
  }
  return(sum.res)
}

# TODO: adjust themes ###########################

plot_boxes <- function(df) {
  p <- ggplot(df, aes(x=model, y=value, color=model)) +
    facet_grid(. ~ measure) +
    geom_boxplot() + 
    geom_jitter()
  return(p)
}

plot_scatters <- function(df) {
  # Reshape so models are in separate columns
  BaTFLED.df <- df[df$model=='BaTFLED',]
  BaTFLED.df <- BaTFLED.df[,c('sample', 'measure', 'value')]
  names(BaTFLED.df)[3] <- 'BaTFLED'
  mean.df <- df[df$model=='mean',]
  mean.df <- mean.df[,c('sample', 'measure', 'value')]
  names(mean.df)[3] <- 'mean'
  merge.df <- merge(BaTFLED.df, mean.df)

  for(measure in levels(merge.df$measure)) {
    p <- ggplot(merge.df[merge.df$measure==measure,], aes(x=mean, y=BaTFLED)) +
      geom_point(size=6) +
      geom_abline(slope=1, intercept=0, lwd=2, color='blue') +
      ggtitle(measure)
    # TODO: Make these square.
    # TODO: Make title bigger
    # TODO: Add percentages on plots
    print(p)
  }
}

for(i in 0:(n.files-1)) {
  run <- paste0(run.prefix, '.', i, '/image.Rdata')
  sum.res <- load_run(run, sum.res)
}

# TODO: make this generate pdfs.
for(task in names(sum.res)) {
  plot_boxes(sum.res[[task]])
  plot_scatters(sum.res[[task]])
}

