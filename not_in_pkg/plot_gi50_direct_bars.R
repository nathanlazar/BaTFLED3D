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
  args <- c('Results/final_design_matrix.txt', 'Results/gi50_predictions/',
            'Results/gi50_predictions/Heiser_gr50_direct_results.csv',
            'Heiser data: GR50 direct prediction final test results')
  if(exists('results.df')) rm(results.df)  
}
design.df <- data.frame(read.table(args[[1]], sep='\t', header=T, stringsAsFactors=F))
dir <- args[2]
data.file <- args[3]
title <- args[4]

types <- c('train', 'train.sd', 'warm', 'warm.sd', 'm1', 'm1.sd',
           'm2', 'm2.sd', 'm1m2', 'm1m2.sd')

csv <- read.csv(data.file, header=T, stringsAsFactors = F)
names(csv) <- c('task', 'model', 'class', 'run', paste0('rmse.', types), paste0('exp.var.', types),
                paste0('p.cor.', types), paste0('s.cor.', types))
all.df <- transform(csv, run=as.character(run))

# Remove training results for all but the warm runs
all.df[is.na(all.df$rmse.warm), grepl('train', colnames(all.df))] <- NA

# remove m1, m2, results for the m1m2 runs
all.df[!is.na(all.df$rmse.m1m2), grepl('m1$', colnames(all.df))] <- NA
all.df[!is.na(all.df$rmse.m1m2), grepl('m1.sd$', colnames(all.df))] <- NA
all.df[!is.na(all.df$rmse.m1m2), grepl('.m2', colnames(all.df), fixed=T)] <- NA

results.df <- melt(all.df, id.vars=c('task', 'model', 'class', 'run'))
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

save(results.df, file=paste0(dir, 'gi_direct_results.Rdata'))

# Remove the redundant data for training etc. 
results.df <- filter(results.df, stat=='mean')

# Set high rmse values to 2 for CV runs.
results.df$value[results.df$value > 2] <- 2e

p <- ggplot(results.df[results.df$measure!='Exp. var.',],
            aes(x=model, y=value, fill=class)) +
  geom_col() +
  scale_fill_hue(breaks=levels(results.df$class), drop=F) +
  facet_grid(measure~mode, scales='free') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0),
        axis.title = element_blank()) +
  ggtitle(title)
p

pdf(paste0(dir, 'Heiser_final_gi50_direct_bars.pdf'), height=7, width=7*1.62)
p
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
