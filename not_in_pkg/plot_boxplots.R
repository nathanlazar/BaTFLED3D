# Examine predictors 

library(reshape2)
library(dplyr)
library(ggplot2)
library(gtools)
library(grid)
library(gridExtra)

# Usage: plot_boxplots.R <design_mat.txt> <dir> <plot_title> 

# Reads the .txt design matrix to get annotations (which method each run is, etc.)
# and then reads the files from dir

args <- commandArgs(TRUE)

test <- F
if(test) {
  args <- c('Results/CV_design_matrix.txt', 'Results/CV/',
            'Heiser data: cross-validation results')
  # args <- c('Results/CV_kern_design_matrix.txt', 'Results/CV_kernels/',
  #           'Heiser data: cross-validation results with kernelized inputs')
  # args <- c('Results/Final_design_matrix.txt', 'Results/Final/',
  #           'Heiser data: final test results')
  # args <- c('Results/Final_kern_design_matrix.txt', 'Results/Final_kernels/',
  #           'Heiser data: final test results with kernelized inputs')
  if(exists('results.df')) rm(results.df)  
}
design.df <- data.frame(read.table(args[[1]], sep='\t', header=T, stringsAsFactors=F))
design.df <- design.df[design.df$task!='Warm',]
dir <- args[2]
title <- args[3]

design.df$file <- sapply(design.df$run, function(x)
  grep(x, list.files(path=dir, recursive=T), value=T))

loadRData <- function(x) {
  load(x)
  return(list(m1=results$by.m1, m2=results$by.m2, m1m2=results$by.m1m2))
}

loadMean <- function(x) {
  load(x)
  return(list(m1=results$by.m1.mean, m2=results$by.m2.mean, m1m2=results$by.m1m2.mean))
}

# for(i in 1:12) {
for(i in 1:nrow(design.df)) {
  if(design.df$class[i]=='Mean') {
    res.list <- loadMean(paste0(dir, design.df$file[i]))
  } else 
    res.list <- loadRData(paste0(dir, design.df$file[i]))

  if(sum(!is.na(res.list$m1$rmse))) {
    res.m1 <- melt(res.list$m1)
    names(res.m1)=c('name', 'measure', 'value')
    res.m1$model <- design.df$model[i]
    res.m1$class <- design.df$class[i]
    res.m1$mode <- 'Cell line'
    if(!exists('results.df')) {
      results.df <- res.m1
    } else
      results.df <- smartbind(results.df, res.m1)
  }
  
  if(sum(!is.na(res.list$m2$rmse))) {
    res.m2 <- melt(res.list$m2)
    names(res.m2)=c('name', 'measure', 'value')
    res.m2$model <- design.df$model[i]
    res.m2$class <- design.df$class[i]
    res.m2$mode <- 'Drug'
    results.df <- smartbind(results.df, res.m2)
  }
  
  if(sum(!is.na(res.list$m1m2$rmse))) {
    res.m1m2 <- melt(res.list$m1m2)
    res.m1m2$name <- paste(res.m1m2$name1, res.m1m2$name2, sep='_+_')
    res.m1m2 <- res.m1m2[,c('name', 'variable', 'value')]
    names(res.m1m2)=c('name', 'measure', 'value')
    res.m1m2$model <- design.df$model[i]
    res.m1m2$class <- design.df$class[i]
    res.m1m2$mode <- 'Cell line/Drug'
    results.df <- smartbind(results.df, res.m1m2)
  }  
}

# Remove any rows without data
results.df <- results.df[!is.na(results.df$value),] 

# TODO: get results for training and warm cell lines, drugs?.

results.df$measure <- as.character(results.df$measure)
results.df$measure[results.df$measure == 'rmse'] <- 'NRMSE' 
results.df$measure[results.df$measure == 'exp.var'] <- 'Exp. var.' 
results.df$measure[results.df$measure == 'p.cor'] <- 'Pearson cor.' 
results.df$measure[results.df$measure == 's.cor'] <- 'Spearman cor.' 

results.df$class <- factor(results.df$class, 
  levels=c('Mean', 'E. net', 'R. forest', 'N. network', 
           'BaTFLED CP', 'BaTFLED Tucker'))
results.df$model <- factor(results.df$model,
                           levels=c('Cl. mean', 'Dr. mean', 'Cl. Dr. mean', 
                                    'LASSO a=1', 'E. net a=0.9', 'E. net a=0.5', 'Ridge a=0',
                                    'R.F. 1000x5', 'R.F. 5000x5', 'R.F. 1000x10', 
                                    'N.N. 1L', 'N.N. 2L', 'N.N. 3L',
                                    'CP', 'CP retrain', 'CP selected', 
                                    'Tucker', 'Tucker retrain', 'Tucker selected'))
results.df$mode <- factor(results.df$mode, 
  levels=c('Cell line', 'Drug', 'Cell line/Drug'))
results.df$measure <- factor(results.df$measure,
  levels=c('NRMSE', 'Exp. var.', 'Pearson cor.', 'Spearman cor.'))

# Remove the results for explained variance
results.df <- results.df[results.df$measure != 'Exp. var.',]

# Remove cl and dr results from cl/dr runs
results.df <- results.df[!(results.df$model=='Cl. Dr. mean' & 
                          results.df$mode %in% c('Cell line', 'Drug')),]

# Squash values above 3 to 3 for plotting
res.tmp <- results.df
res.tmp$value[res.tmp$value > 3] <- 3

# Remove values greater than 3 for plotting
# res.tmp <- results.df[results.df$value < 3,]

# Plot where all the axes are free
p.free <- ggplot(res.tmp, aes(x=model, y=value, color=class)) +
  geom_jitter(width=0.2, alpha=0.5) +
  geom_boxplot(color='black', fill= "transparent", outlier.alpha=0) +
  facet_wrap(measure~mode, scales='free') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0),
        axis.title = element_blank())
# p.free

pdf(paste0(dir, 'boxplot.pdf'), height=12, width=14)
p.free
dev.off()

p.cl <- ggplot(res.tmp[res.tmp$mode=='Cell line',], 
                 aes(x=model, y=value, color=class)) +
  geom_jitter(width=0.2, alpha=0.5) +
  geom_boxplot(color='black', fill= "transparent", outlier.alpha=0) +
  stat_summary(fun.y=mean, colour="red", geom="point") +
  facet_wrap(measure~mode, scales='free') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0)) 
# p.cl

# Plot head to head for cell lines and drugs
head_to_head <- function(results.df, model1, model2, mode1, measure1,
                         xlim=c(0,1), ylim=c(0,1)) {
  sub.df <- results.df %>% 
    dplyr::filter(mode==mode1) %>% 
    dplyr::filter(model==model1 | model==model2) %>%
    dplyr::filter(measure==measure1)

  sub.df1 <- sub.df[sub.df$model==model1,]
  names(sub.df1)[names(sub.df1)=='value'] <- 'model1'
    # as.character(unique(sub.df1$model))
  sub.df1 <- select(sub.df1, -model, -class)
  sub.df2 <- sub.df[sub.df$model==model2,]
  names(sub.df2)[names(sub.df2)=='value'] <- 'model2'
    # as.character(unique(sub.df2$model))
  sub.df2 <- select(sub.df2, -model, -class)
  
  sub.df <- merge(sub.df1, sub.df2)
  sub.df$col <- 3
  if(measure1=='NRMSE') { # Color the method better than the mean red.
    sub.df$col[sub.df$model1 < sub.df$model2] <- 2
    sub.df$col[sub.df$model1 > sub.df$model2] <- 1
  } else {
    sub.df$col[sub.df$model1 < sub.df$model2] <- 1
    sub.df$col[sub.df$model1 > sub.df$model2] <- 2
  }
  sub.df$col <- factor(sub.df$col, levels=1:3)
  
  p <- ggplot(sub.df, aes(x=model1, y=model2, color=col)) +
    geom_point(cex=2) +
    geom_abline(slope=1, intercept=0, color='blue') +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(x=model1, y=model2) +
    guides(color=FALSE) + 
    scale_color_manual(values=c('1'='red', '2'='black', '3'='blue')) +
    annotate('text', x=c(0.1,xlim[2]*0.9), y=c(ylim[2]*0.9,0.1),
             label=c(sprintf('%.1f%%', mean(sub.df$model1 < sub.df$model2)*100),
                     sprintf('%.1f%%', mean(sub.df$model1 > sub.df$model2)*100)))
  return(p)
}

rep.models <- c('Cl. mean', 'LASSO a=1', 'E. net a=0.9', 'E. net a=0.5', 'Ridge a=0',
                'N.N. 1L', 'N.N. 2L', 'N.N. 3L', 'R.F. 1000x5', 'R.F. 5000x5', 'R.F. 1000x10',
                'CP', 'CP retrain', 'CP selected', 'Tucker', 'Tucker retrain', 'Tucker selected')

# Cell line prediction
######################
# Start list with an empty plot to insert
rmse.h.to.h.m1 <- list(ggplot(data.frame()) + geom_blank() +
                         theme(panel.background = element_blank()))
pearson.h.to.h.m1 <- rmse.h.to.h.m1

i <- 1
for(m1 in rep.models) for(m2 in rep.models) {
  if(sum(results.df$model==m1) & sum(results.df$model==m2)) {
    i <- i+1
    rmse.h.to.h.m1[[i]] <- head_to_head(results.df=results.df, 
      model1=m1, model2=m2, mode1='Cell line', measure1='NRMSE')
    pearson.h.to.h.m1[[i]] <- head_to_head(results.df=results.df, 
      model1=m1, model2=m2, mode1='Cell line', measure1='Pearson cor.')
  }
}

# Display head to head plots comparing with mean prediction
if(('CP retrain' %in% results.df$model) | ('CP selected' %in% results.df$model)){
  selected <- c(3:9,1,10:12,1,13:16)
} else selected <- c(3:9,1,10:12,1,13,14)

pdf(paste0(dir, 'cl_rmse_head2head.pdf'), height=10, width=13)
grid.arrange(grobs=rmse.h.to.h.m1[selected], ncol=4,
             top=textGrob('Normalized RMSE performance on cold-test cell lines',
                          gp=gpar(fontsize=20)))
dev.off()

pdf(paste0(dir, 'cl_pearson_head2head.pdf'), height=10, width=13)
grid.arrange(grobs=pearson.h.to.h.m1[selected], ncol=4,
             top=textGrob('Pearson correlation performance on cold-test cell lines',
                          gp=gpar(fontsize=20)))
dev.off()

# Drug prediction
#############################
# Start list with an empty plot to insert
rmse.h.to.h.m2 <- list(ggplot(data.frame()) + geom_blank() +
                         theme(panel.background = element_blank()))
pearson.h.to.h.m2 <- rmse.h.to.h.m2

rep.models[1] <- 'Dr. mean'
  
i <- 1
for(m1 in rep.models) for(m2 in rep.models) {
  if(sum(results.df$model==m1) & sum(results.df$model==m2)) {
    i <- i+1
    rmse.h.to.h.m2[[i]] <- head_to_head(results.df, model1=m1, model2=m2, 
      mode1='Drug', measure1='NRMSE', xlim=c(0,1.5), ylim=c(0,2))
    pearson.h.to.h.m2[[i]] <- head_to_head(results.df, 
      model1=m1, model2=m2, mode1='Drug', measure1='Pearson cor.')
  }
}

pdf(paste0(dir, 'dr_rmse_head2head.pdf'), height=10, width=13)
grid.arrange(grobs=rmse.h.to.h.m2[selected], ncol=4,
             top=textGrob('Normalized RMSE performance on cold-test drugs',
                          gp=gpar(fontsize=20)))
dev.off()

pdf(paste0(dir, 'dr_pearson_head2head.pdf'), height=10, width=13)
grid.arrange(grobs=pearson.h.to.h.m2[selected], ncol=4,
             top=textGrob('Pearson correlation performance on cold-test drugs',
                          gp=gpar(fontsize=20)))
dev.off()

# Cell line/drug prediction
#############################
# Start list with an empty plot to insert
rmse.h.to.h.m1m2 <- list(ggplot(data.frame()) + geom_blank() +
                         theme(panel.background = element_blank()))
pearson.h.to.h.m1m2 <- rmse.h.to.h.m1m2

rep.models[1] <- 'Cl. Dr. mean'

i <- 1
for(m1 in rep.models) for(m2 in rep.models) {
  if(sum(results.df$model==m1) & sum(results.df$model==m2)) {
    i <- i+1
    rmse.h.to.h.m1m2[[i]] <- head_to_head(results.df=results.df, 
      model1=m1, model2=m2, mode1='Cell line/Drug', measure1='NRMSE',
      xlim=c(0,3), ylim=c(0,3))
    pearson.h.to.h.m1m2[[i]] <- head_to_head(results.df=results.df, 
      model1=m1, model2=m2, mode1='Cell line/Drug', measure1='Pearson cor.')
  }
}

# Display head to head plots comparing with mean prediction
if(('CP retrain' %in% results.df$model) | ('CP selected' %in% results.df$model)){
  selected <- c(3:9,1,10:12,1,13:16)
} else selected <- c(3:9,1,10:12,1,13,14)

pdf(paste0(dir, 'cldr_rmse_head2head.pdf'), height=10, width=13)
grid.arrange(grobs=rmse.h.to.h.m1m2[selected], ncol=4,
  top=textGrob('Normalized RMSE performance on cold-test cell line/drug combinations',
               gp=gpar(fontsize=20)))
dev.off()

pdf(paste0(dir, 'cldr_pearson_head2head.pdf'), height=10, width=13)
grid.arrange(grobs=pearson.h.to.h.m1m2[selected], ncol=4,
  top=textGrob('Pearson correlation performance on cold-test cell line/drug combinations',
               gp=gpar(fontsize=20)))
dev.off()
