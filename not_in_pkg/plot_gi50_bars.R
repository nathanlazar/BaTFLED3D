# Examine predictors 

library(lattice)
library(dplyr)
library(reshape2)
library(ggplot2)
library(gtools)

# Usage: plot_gi50_bars.R <design_mat.txt> <dir> <mean.csv> <plot_title> 

# Reads the .txt design matrix to get annotations (which method each run is, etc.)
# and then reads the files from dir

args <- commandArgs(TRUE)

test <- T
if(test) {
  # args <- c('Results/CV_design_matrix.txt', 'Results/gi50_predictions/',
  #           'Results/gi50_predictions/all_cv_mean_gi50s.csv',
  #           'Heiser data: GR50 prediction, cross-validation results')
  args <- c('Results/Final_design_matrix.txt', 'Results/gi50_predictions/',
            'Results/gi50_predictions/all_final_mean_gi50s.csv',
            'Heiser data: GR50 prediction, final test results')
  if(exists('results.df')) rm(results.df)  
}
design.df <- data.frame(read.table(args[[1]], sep='\t', header=T, stringsAsFactors=F))
dir <- args[2]
mean.file <- args[3]
title <- args[4]

types <- c('train', 'train.sd', 'warm', 'warm.sd', 'm1', 'm1.sd',
           'm2', 'm2.sd', 'm1m2', 'm1m2.sd')

loadData <- function(x) {
  csv <- read.csv(x, header=F, stringsAsFactors = F)
  names(csv) <- c('run', paste0('rmse.', types), paste0('exp.var.', types),
                  paste0('p.cor.', types), paste0('s.cor.', types))
  csv <- transform(csv, run=as.character(run))
  return(csv)
}

# Load mean results
all.df <- loadData(mean.file)
all.df$class <- 'Mean'
all.df$model <- 'Mean'

design.df <- filter(design.df, class != 'Mean')

# Load E. net results
e.net <- loadData(paste0(dir, grep('LASSO', list.files(path=dir), value=T)))
e.net <- e.net[e.net$run %in% design.df$run,]
e.net$class <- design.df$class[match(e.net$run, design.df$run, nomatch=0)]
e.net$model <- design.df$model[match(e.net$run, design.df$run, nomatch=0)]
all.df <- smartbind(all.df, e.net)

# Load neural net results
n.net <- loadData(paste0(dir, grep('NN', list.files(path=dir), value=T)))
n.net <- n.net[n.net$run %in% design.df$run,]
n.net$class <- design.df$class[match(n.net$run, design.df$run, nomatch=0)]
n.net$model <- design.df$model[match(n.net$run, design.df$run, nomatch=0)]
all.df <- smartbind(all.df, n.net)

# Load random forest results
r.forest <- loadData(paste0(dir, grep('RF', list.files(path=dir), value=T)))
r.forest <- r.forest[r.forest$run %in% design.df$run,]
r.forest$class <- design.df$class[match(r.forest$run, design.df$run, nomatch=0)]
r.forest$model <- design.df$model[match(r.forest$run, design.df$run, nomatch=0)]
all.df <- smartbind(all.df, r.forest)

# Load BaTFLED results
batfled <- loadData(paste0(dir, grep('BaTFLED', list.files(path=dir), value=T)))
batfled <- batfled[batfled$run %in% design.df$run,]
batfled$class <- design.df$class[match(batfled$run, design.df$run, nomatch=0)]
batfled$model <- design.df$model[match(batfled$run, design.df$run, nomatch=0)]
all.df <- smartbind(all.df, batfled)

# Remove training results for all but the warm runs
all.df[is.na(all.df$rmse.warm), grepl('train', colnames(all.df))] <- NA

# remove m1, m2, results for the m1m2 runs
all.df[!is.na(all.df$rmse.m1m2), grepl('m1$', colnames(all.df))] <- NA
all.df[!is.na(all.df$rmse.m1m2), grepl('m1.sd$', colnames(all.df))] <- NA
all.df[!is.na(all.df$rmse.m1m2), grepl('.m2', colnames(all.df), fixed=T)] <- NA

results.df <- melt(all.df, id.vars=c('class', 'model', 'run'))
# Remove NA rows
results.df <- results.df[!is.na(results.df$value),]

results.df$stat <- 'mean'
results.df$stat[grepl('sd', results.df$variable)] <- 'sd'
results.df$variable <- sub('.mean|.sd', '', results.df$variable)
results.df$measure <- 'rmse'
results.df$measure[grepl('exp.var', results.df$variable)] <- 'exp.var'
results.df$measure[grepl('p.cor', results.df$variable)] <- 'p.cor'
results.df$measure[grepl('s.cor', results.df$variable)] <- 's.cor'
results.df$mode <- sub('rmse.|exp.var.|p.cor.|s.cor.', '', results.df$variable)
  
results.df$mode[results.df$mode == 'train'] <- 'Training' 
results.df$mode[results.df$mode == 'warm'] <- 'Warm start' 
results.df$mode[results.df$mode == 'm1'] <- 'Cell lines' 
results.df$mode[results.df$mode == 'm2'] <- 'Drugs' 
results.df$mode[results.df$mode == 'm1m2'] <- 'Cell line/Drug' 

results.df$measure <- as.character(results.df$measure)
results.df$measure[results.df$measure == 'rmse'] <- 'NRMSE' 
results.df$measure[results.df$measure == 'exp.var'] <- 'Exp. var.' 
results.df$measure[results.df$measure == 'p.cor'] <- 'Pearson cor.' 
results.df$measure[results.df$measure == 's.cor'] <- 'Spearman cor.' 

results.df$class <- factor(results.df$class, 
  levels=c('Mean', 'E. net', 'R. forest', 'N. network', 'BaTFLED CP', 'BaTFLED Tucker'))
results.df$model <- factor(results.df$model,
  levels=c('Mean', 'LASSO a=1', 'E. net a=0.9', 'E. net a=0.5', 'Ridge a=0',
           'R.F. 1000x5', 'R.F. 5000x5', 'R.F. 1000x10', 
           'N.N. 1L', 'N.N. 2L', 'N.N. 3L',
           'CP', 'CP retrain', 'CP selected', 
           'Tucker', 'Tucker retrain', 'Tucker selected'))
results.df$mode <- factor(results.df$mode, 
  levels=c('Training', 'Warm start', 'Cell lines', 'Drugs', 'Cell line/Drug'))
results.df$measure <- factor(results.df$measure,
  levels=c('NRMSE', 'Exp. var.', 'Pearson cor.', 'Spearman cor.'))

if(grepl('cross', title))
  save(results.df, file=paste0(dir, 'cv_results.Rdata'))
if(grepl('final', title))
  save(results.df, file=paste0(dir, 'final_results.Rdata'))

# Remove the redundant data for training etc. 
results.df <- filter(results.df, stat=='mean')

# Set high rmse values to 1.5 for CV runs.
if(grepl('cross', title))
  results.df$value[results.df$value > 1.5] <- 1.5

p <- ggplot(results.df[results.df$measure!='Exp. var.',],
            aes(x=model, y=value, fill=class)) +
  geom_col() +
  facet_grid(measure~mode, scales='free') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0),
        axis.title = element_blank()) +
  ggtitle(title)
p

p.pearson <- ggplot(results.df[results.df$measure=='Pearson cor.',],
                    aes(x=model, y=value, fill=class)) +
  geom_col() +
  facet_grid(.~mode, scales='free') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0),
        axis.title = element_blank()) +
  ggtitle(title)
p.pearson


if(grepl('cross', title))
  pdf(paste0(dir, 'cv_bars.pdf'), height=7, width=7*1.62)
if(grepl('final', title))
  pdf(paste0(dir, 'final_bars.pdf'), height=7, width=7*1.62)
p
dev.off()

if(grepl('cross', title))
  png(paste0(dir, 'cv_pearson_bars.png'), width=480*2)
if(grepl('final', title))
  png(paste0(dir, 'final_pearson_bars.png'), width=480*2)
p.pearson
dev.off()

save.image(paste0(dir, 'gi50_summary.Rdata'))

#############################################
# for(i in 1:nrow(design.df)) {
#   if(design.df$class[i] == 'Mean') {
#     res <- loadMean(paste0(dir, design.df$file[i]))
#   } else {
#     res <- loadRData(paste0(dir, design.df$file[i]))
#   }
#   res <- melt(res)
#   res <- res[!is.na(res$value),]
#   
#   if(design.df$class[i] == 'Mean') {
#     names(res) <- c('measure', 'mode', 'value')
#     if(design.df$model[i]=='Cl. mean')
#       res <- filter(res, !grepl('m2|m3', res$mode))
#     if(design.df$model[i]=='Dr. mean')
#       res <- filter(res, !grepl('m1|m3', res$mode))
#     if(design.df$model[i]=='Cl. Dr. Mean') {
#       res <- filter(res, grepl('m1m2', res$mode))
#       res <- filter(res, !grepl('m1m2m3', res$mode))
#     }
#   } else {
#     names(res) <- c('mode', 'measure', 'value')
#     if(('m1' %in% res$mode) & ('m2' %in% res$mode)) 
#       res <- filter(res, mode %in% 'm1m2')
#   }
#   res$model <- design.df$model[i]
#   res$class <- design.df$class[i]
#   res$run <- design.df$run[i]
#   
#   if(!exists('results.df')) {
#     results.df <- res
#   } else {
#     results.df <- smartbind(results.df, res)
#   }
# }
