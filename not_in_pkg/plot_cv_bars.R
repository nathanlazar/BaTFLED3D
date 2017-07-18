# Examine predictors 

library(lattice)
library(reshape2)
library(ggplot2)
library(gtools)

# Usage: plot_CV_bars.R <excel_file.xlsx> <dir>  

# Assumes that the runs are divided into different folders depending on the algorithm
# i.e. Elastic net, Neural net, Random forest, BaTFLED
# These names are used as group labels for plots. 

args <- commandArgs(TRUE)

test <- T
if(test) args <- c('Results/CV/')
dir <- args[[1]]

classes <- dir(path=dir)
files <- list()
for(class in classes) 
  files[[class]] <- list.files(path=paste(dir, class, sep='/'))

# class <- c()
# for(i in 1:length(classes)) class <- append(class, rep(names(files)[i], length(files[[i]])))
# class <- sub('_', ' ', class)

# model <- c('LASSO a=1', 'E. net a=0.9', 'E. net a=0.5', 'Ridge a=0',
#            'R.F. 1Kx5', 'R.F. 5Kx5', 'R.F. 1Kx10',
#            'N.N. 1L', 'N.N. 2L', 'N.N. 3L',
#            'CP', 'CP retrain', 'Tucker', 'Tucker retrain')

loadRData <- function(x, type) {
  load(x)
  if(type=='Elastic net') {
    if(alpha==1) model <- 'LASSO a=1'
    if(alpha==0) model <- 'Ridge a=0'
    if(!(alpha %in% c(0,1))) model=paste0('E. net a=',alpha)
  }
  if(type=='Neural net') {
    
  }
  sum <- apply(results$summaries, c(1,2), mean)
  sum <- sum[apply(!is.na(sum), 1, sum)>0,]

  return(list(model=model, summary=sum))
}

loadMean <- function(x) {
  load(x)
  return(apply(results$mean.pred, c(1,2), mean))
}

for(type in names(files))
  for(i in 1:length(files[[type]])) {
    res.list <- loadRData(paste(dir, type, files[[type]][i], sep='/'), type)
    res <- melt(res.list$summary)
  
    names(res) <- c('mode', 'measure', 'value')
    res$model <- res.list$model
    res$class <- type
    
    if(!exists('results.df')) {
      results.df <- res
    } else {
      results.df <- smartbind(results.df, res)
    }
}

# head(results.df)

# Also get results for predicting mean responses 
# from the LASSO run
mean.list <- loadMean(paste0(dir, 'run_', runs[1], '_summary.Rdata'))
mean.melt <- melt(mean.list)
names(mean.melt)=c('mode', 'measure', 'value')
mean.melt <- mean.melt[mean.melt$mode!='train.m2',]
mean.melt <- mean.melt[mean.melt$mode!='warm.m2',]
mean.melt$mode <- as.character(mean.melt$mode)
mean.melt$mode[mean.melt$mode=='train.m1'] <- 'train'
mean.melt$mode[mean.melt$mode=='warm.m1'] <- 'warm'
mean.melt$model <- 'Mean'
mean.melt$class <- 'Mean'
# head(mean.melt)

results.df <- smartbind(results.df, mean.melt)
results.df <- filter(results.df, measure %in%
  c('RMSE', 'exp.var', 'p.cor', 's.cor'))
results.df$measure <- droplevels(results.df$measure)
results.df <- filter(results.df, mode!='H')
results.df$mode[results.df$mode == 'train'] <- 'Training' 
results.df$mode[results.df$mode == 'warm'] <- 'Warm start' 
results.df$mode[results.df$mode == 'm1'] <- 'Cell lines' 
results.df$mode[results.df$mode == 'm2'] <- 'Drugs' 
results.df$mode[results.df$mode == 'm1m2'] <- 'Cell line/Drug' 

results.df$measure <- as.character(results.df$measure)
results.df$measure[results.df$measure == 'RMSE'] <- 'NRMSE' 
results.df$measure[results.df$measure == 'exp.var'] <- 'Exp. var.' 
results.df$measure[results.df$measure == 'p.cor'] <- 'Pearson cor.' 
results.df$measure[results.df$measure == 's.cor'] <- 'Spearman cor.' 

results.df$class <- factor(results.df$class, 
  levels=c('Mean', 'E. net', 'R. forest', 'N. network', 'CP', 'Tucker'))
results.df$model <- factor(results.df$model, 
  levels=c('Mean', model))
results.df$mode <- factor(results.df$mode, 
  levels=c('Training', 'Warm start', 'Cell lines', 'Drugs', 'Cell line/Drug'))
results.df$measure <- factor(results.df$measure,
  levels=c('NRMSE', 'Exp. var.', 'Pearson cor.', 'Spearman cor.'))

p <- ggplot(results.df[results.df$measure!='Exp. var.',],
            aes(x=model, y=value, fill=class)) +
  geom_col() +
  facet_grid(measure~mode, scales='free') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_blank())
# p

pdf(paste0(dir, 'DREAM_cv_bars.pdf'), height=7, width=7*1.62)
p
dev.off()

