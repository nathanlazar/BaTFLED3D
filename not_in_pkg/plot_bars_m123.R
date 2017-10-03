# Examine predictors 

library(lattice)
library(dplyr)
library(reshape2)
library(ggplot2)
library(gtools)

# Usage: plot_bars.R <design_mat.txt> <dir> <plot_title> 

# Reads the .txt design matrix to get annotations (which method each run is, etc.)
# and then reads the files from dir

args <- commandArgs(TRUE)

test <- T
if(test) {
  args <- c('Results/design_matrix.txt', 'Results/',
            'MEMA data: Cross-validation results with kernelized inputs')
  if(exists('results.df')) rm(results.df)  
}
design.df <- data.frame(read.table(args[[1]], sep='\t', header=T, stringsAsFactors=F))
dir <- args[2]
title <- args[3]

design.df$file <- sapply(design.df$run, function(x)
  grep('summary', grep(x, list.files(path=dir, recursive=T), value=T), value=T))

loadRData <- function(x) {
  load(x)
  sum <- apply(results$summaries, c(1,2), mean)
  sum <- sum[apply(!is.na(sum), 1, sum)>0,]
  if(!sum(grepl('warm', rownames(sum))))
    sum <- sum[rownames(sum) != 'train',,drop=F]
  return(sum)
}

loadMean <- function(x) {
  load(x)
  if('mean' %in% names(results)) {
    sum <- apply(results$mean, c(1,2), mean)
    sum <- sum[,apply(!is.na(sum), 2, sum)>0]
  }
  if('mean.pred' %in% names(results)) {
    sum <- apply(results$mean.pred, c(1,2), mean)
    sum <- sum[apply(!is.na(sum), 1, sum)>0,]
  }
  if(!sum(grepl('warm', colnames(sum))))
    sum <- sum[,!grepl('train', colnames(sum)),drop=F]
  if(sum(grepl('m1m2m3', colnames(sum))))
    sum <- sum[,!(colnames(sum) %in% c('m1', 'm2', 'm3', 'm1m2', 'm1m3', 'm2m3')),drop=F]
  if(sum(grepl('m1m2', colnames(sum))))
    sum <- sum[,!(colnames(sum) %in% c('m1', 'm2')),drop=F]
  if(sum(grepl('m1m3', colnames(sum))))
    sum <- sum[,!(colnames(sum) %in% c('m1', 'm3')),drop=F]
  if(sum(grepl('m2m3', colnames(sum))))
    sum <- sum[,!(colnames(sum) %in% c('m2', 'm3')),drop=F]
  return(sum)
}

for(i in 1:nrow(design.df)) {
  if(design.df$class[i] == 'Mean') {
    res <- loadMean(paste0(dir, design.df$file[i]))
  } else {
    res <- loadRData(paste0(dir, design.df$file[i]))
  }
  res <- melt(res)
  res <- res[!is.na(res$value),]
  
  if(design.df$class[i] == 'Mean') {
    names(res) <- c('measure', 'mode', 'value')
    if(design.df$model[i] != 'Mean')
      res <- filter(res, !(grepl('m1m|m2m', res$mode)))
    if(design.df$model[i]=='Cl. mean')
      res <- filter(res, !(mode %in% c('train.m2', 'train.m3', 'warm.m2', 'warm.m3')))
    if(design.df$model[i]=='Lig. mean')
      res <- filter(res, !(mode %in% c('train.m1', 'train.m3', 'warm.m1', 'warm.m3')))
    if(design.df$model[i]=='ECMP mean') {
      res <- filter(res, !(mode %in% c('train.m1', 'train.m2', 'warm.m1', 'warm.m2')))
    }
  } else {
    names(res) <- c('mode', 'measure', 'value')
    if(('m1' %in% res$mode) & ('m2' %in% res$mode) & ('m3' %in% res$mode)) {
     res <- filter(res, mode %in% 'm1m2m3')
    } else if(('m1' %in% res$mode) & ('m2' %in% res$mode)) {
      res <- filter(res, mode %in% 'm1m2')
    } else if(('m1' %in% res$mode) & ('m3' %in% res$mode)) {
      res <- filter(res, mode %in% 'm1m3')
    } else if(('m2' %in% res$mode) & ('m3' %in% res$mode)) {
    res <- filter(res, mode %in% 'm2m3')
  }
}
  res$model <- design.df$model[i]
  res$class <- design.df$class[i]
  res$run <- design.df$run[i]
  
  if(!exists('results.df')) {
    results.df <- res
  } else {
    results.df <- smartbind(results.df, res)
  }
}

# head(results.df)

results.df$mode <- sub('.m1m2m3', '', results.df$mode, fixed=T)
results.df$mode <- sub('.m1m2', '', results.df$mode, fixed=T)
results.df$mode <- sub('.m1m3', '', results.df$mode, fixed=T)
results.df$mode <- sub('.m2m3', '', results.df$mode, fixed=T)
results.df$mode <- sub('.m1', '', results.df$mode, fixed=T)
results.df$mode <- sub('.m2', '', results.df$mode, fixed=T)
results.df$mode <- sub('.m3', '', results.df$mode, fixed=T)
results.df <- filter(results.df, mode!='H')

results.df <- filter(results.df, measure %in%
  c('RMSE', 'exp.var', 'p.cor', 's.cor'))
results.df$measure <- droplevels(results.df$measure)

# results.df <- results.df[!grepl('m3', results.df$mode),]

results.df$mode[results.df$mode == 'train'] <- 'Training' 
results.df$mode[results.df$mode == 'warm'] <- 'Warm start' 
results.df$mode[results.df$mode == 'm1'] <- 'Cell lines' 
results.df$mode[results.df$mode == 'm2'] <- 'Ligands' 
results.df$mode[results.df$mode == 'm3'] <- 'ECMPs'
results.df$mode[results.df$mode == 'm1m2'] <- 'Cell line/Ligand' 
results.df$mode[results.df$mode == 'm1m3'] <- 'Cell line/ECMP' 
results.df$mode[results.df$mode == 'm2m3'] <- 'Ligand/ECMP' 
results.df$mode[results.df$mode == 'm1m2m3'] <- 'All' 
results.df$mode <- factor(results.df$mode, 
  levels=c('Training', 'Warm start', 'Cell lines', 'Ligands', 'ECMPs',
           'Cell line/Ligand', 'Cell line/ECMP' , 'Ligand/ECMP', 'All'))

results.df$measure <- as.character(results.df$measure)
results.df$measure[results.df$measure == 'RMSE'] <- 'NRMSE' 
results.df$measure[results.df$measure == 'exp.var'] <- 'Exp. var.' 
results.df$measure[results.df$measure == 'p.cor'] <- 'Pearson cor.' 
results.df$measure[results.df$measure == 's.cor'] <- 'Spearman cor.' 
results.df$measure <- factor(results.df$measure,
  levels=c('NRMSE', 'Exp. var.', 'Pearson cor.', 'Spearman cor.'))

results.df$class <- factor(results.df$class, 
  levels=c('Mean', 'E. net', 'R. forest', 'N. network', 'BaTFLED CP', 'BaTFLED Tucker'))
results.df$model <- factor(results.df$model,
  levels=c('Mean', 'Cl. mean', 'Lig. mean', 'ECMP mean', 
           'LASSO a=1', 'E. net a=0.9', 'E. net a=0.5', 'Ridge a=0',
           'R.F. 1000x5', 'R.F. 5000x5', 'R.F. 1000x10', 
           'N.N. 1L', 'N.N. 2L', 'N.N. 3L',
           'CP', 'CP retrain', 'CP selected', 
           'Tucker', 'Tucker retrain', 'Tucker selected'))

# Remove explained variation results
results.df <- results.df[results.df$measure != 'Exp. var.',]

# Set values above 1.5 to 1.5
# results.df$value[results.df$value > 1.5] <- 1.5

p <- ggplot(results.df, aes(x=model, y=value, fill=class)) +
  geom_col() +
  facet_grid(measure~mode, scales='free') +
  # scale_fill_hue(breaks=unique(results.df$class), drop=F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0),
        axis.title = element_blank()) +
  ggtitle(title)
p

pearson.df <- results.df[results.df$measure=='Pearson cor.',]
p.pearson <- ggplot(pearson.df, aes(x=model, y=value, fill=class)) +
  geom_col() +
  facet_wrap(measure~mode, scales='free_x') +
  # scale_fill_hue(breaks=unique(results.df$class), drop=F) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0),
        axis.title = element_blank()) +
  ggtitle(title)
p.pearson

# if(min(results.df$value) < 0 |
#    max(results.df$value) > 1.2) {
#   p  <- p + coord_cartesian(ylim = c(0, 1.2))
# } else
#   p <- p + coord_cartesian(ylim=c(0,1)) 

pdf(paste0(dir, 'bars.pdf'), height=7, width=14)
p
dev.off()

png(paste0(dir, 'pearson_bars.png'), height=480*2, width=480*2)
p.pearson
dev.off()

save(results.df, file=paste0(dir, 'results_df.Rdata'))
