#!/usr/bin/env Rscript

library(dplyr)
library(methods)
library(R6)
library(BaTFLED3D)
# devtools::document()

args <- commandArgs(TRUE)
run_prefix <- args[1]

# Cut offs used when determining sparsity
m1.rem.cut <- 1e-5
m2.rem.cut <- 1e-5
m3.rem.cut <- 1e-5
core.rem.cut <- 1e-3

# Functions
#############################################################

# Read in a file to get the number of iterations
get_iter_binary <- function(f1) {
  load(f1)
  binary <- F
  if(identical(sort(unique(train.data$resp[!is.na(train.data$resp)])), c(-1,1)))
    binary <- T
  ret <- list(iters=trained$iter, binary=binary)
  return(ret)
}

# Function to load in data from runs
loadData <- function(f, all.results){
  #loads an RData file, and returns a list with the 
  #trained model and warm & cold RMSE vectors.

  load(f)

  # Fix warm response measures
  warm.preds <- train.pred.resp[is.na(train.data$resp)]
  if(warm.per > 0) {
    results$trained['RMSE', 'warm'] <- nrmse(all.resp.train[is.na(train.data$resp)],
                                             warm.preds)
    results$trained['exp.var', 'warm'] <- exp_var(all.resp.train[is.na(train.data$resp)],
                                                  warm.preds)
    results$trained['p.cor', 'warm'] <- cor(all.resp.train[is.na(train.data$resp)],
      warm.preds, use='complete.obs')
    results$trained['s.cor', 'warm'] <- cor(all.resp.train[is.na(train.data$resp)],
      warm.preds, use='complete.obs', method='spearman')
  }

  for(mode in colnames(results$mean)) {
    for(type in rownames(results$mean)) {
      all.results$mean[type, mode, fold] <- results$mean[type, mode]
    }
  }

  all.results$training['lower.bnd', fold, 1:trained$iter] <-
    trained$lower.bnd
  all.results$lower.bnd['final', fold] <- trained$lower.bnd[trained$iter]
  all.results$lower.bnd['max', fold] <- max(trained$lower.bnd)
  all.results$lower.bnd['which.max', fold] <- which.max(trained$lower.bnd)
  all.results$lower.bnd['not.mono', fold] <- sum((trained$lower.bnd[-1] -
                              trained$lower.bnd[-trained$iter]) < 0)

  all.results$training['A.RMSE', fold, 1:trained$iter] <- trained$RMSE
  all.results$training['H.RMSE', fold, 1:trained$iter] <- trained$H.RMSE
  all.results$training['exp.var', fold, 1:trained$iter] <- trained$exp.var
  if(length(trained$p.cor))
    all.results$training['p.cor', fold, 1:trained$iter] <- trained$p.cor
  if(length(trained$s.cor))
    all.results$training['s.cor', fold, 1:trained$iter] <- trained$s.cor

  # Remove NA or NaN columns from test.results
  test.results <- test.results[,apply(test.results, 2, function(x) sum(!is.na(x))>0)]

  for(mode in c('warm', 'm1', 'm2', 'm3', 'm1m2', 'm1m3', 'm2m3', 'm1m2m3')) {
    if(paste0(mode, '.RMSE') %in% names(test.results)) {
      all.results$training[paste0(mode, '.RMSE'), fold, 1:trained$iter] <- 
        test.results[[paste0(mode, '.RMSE')]]
      all.results$training[paste0(mode, '.exp.var'), fold, 1:trained$iter] <- 
        test.results[[paste0(mode, '.exp.var')]]
      if(length(test.results[[paste0(mode, '.p.cor')]]))
        all.results$training[paste0(mode, '.p.cor'), fold, 1:trained$iter] <- 
          test.results[[paste0(mode, '.p.cor')]]
      if(length(test.results[[paste0(mode, '.s.cor')]]))
        all.results$training[paste0(mode, '.s.cor'), fold, 1:trained$iter] <- 
          test.results[[paste0(mode, '.p.cor')]]
    }
  }

  all.results$summaries['H', 'RMSE', fold] <- nrmse(train.data$resp, train.H.resp)
  all.results$summaries['H', 'min.RMSE.iter', fold] <- which.min(trained$H.RMSE)
  all.results$summaries['H', 'min.RMSE', fold] <- min(trained$H.RMSE)

  all.results$summaries['train', 'RMSE', fold] <- nrmse(train.data$resp, trained$resp)
  all.results$summaries['train', 'min.RMSE.iter', fold] <- which.min(trained$RMSE)
  all.results$summaries['train', 'min.RMSE', fold] <- min(trained$RMSE)
  all.results$summaries['train', 'exp.var', fold] <- exp_var(train.data$resp, trained$resp)
  all.results$summaries['train', 'max.exp.var.iter', fold] <- which.max(trained$exp.var)
  all.results$summaries['train', 'max.exp.var', fold] <- max(trained$exp.var)
  if(length(trained$p.cor)) {
    all.results$summaries['train', 'p.cor', fold] <- p_cor(train.data$resp, trained$resp)
    all.results$summaries['train', 'max.p.cor.iter', fold] <- which.max(trained$p.cor)
    all.results$summaries['train', 'max.p.cor', fold] <- max(trained$p.cor)
    all.results$summaries['train', 's.cor', fold] <- s_cor(train.data$resp, trained$resp)
    all.results$summaries['train', 'max.s.cor.iter', fold] <- which.max(trained$s.cor)
    all.results$summaries['train', 'max.s.cor', fold] <- max(trained$s.cor)
  }

  for(mode in c('warm', 'm1', 'm2', 'm3', 'm1m2', 'm1m3', 'm2m3', 'm1m2m3')) {
    if(paste0(mode, '.RMSE') %in% names(test.results)) {
      if(mode=='warm') {
        preds <- warm.preds
      } else {
        preds <- get(paste0(mode, '.pred.resp'))
      }
      resp <- get(paste0(mode, '.resp'))

      all.results$summaries[mode, 'RMSE', fold] <- nrmse(resp, preds)
      all.results$summaries[mode, 'min.RMSE.iter', fold] <- 
        which.min(test.results[[paste0(mode, '.RMSE')]])
      all.results$summaries[mode, 'min.RMSE', fold] <- 
        min(test.results[[paste0(mode, '.RMSE')]])
      all.results$summaries[mode, 'clip.RMSE', fold] <- 
        test.results[[paste0(mode, '.RMSE.clip')]][trained$iter]
      all.results$summaries[mode, 'min.clip.RMSE.iter', fold] <- 
        which.min(test.results[[paste0(mode, '.RMSE.clip')]])
      all.results$summaries[mode, 'min.clip.RMSE', fold] <- 
        min(test.results[[paste0(mode, '.RMSE.clip')]])
      all.results$summaries[mode, 'exp.var', fold] <- exp_var(resp, preds)
      all.results$summaries[mode, 'max.exp.var.iter', fold] <- 
        which.max(test.results[[paste0(mode, '.exp.var')]])
      all.results$summaries[mode, 'max.exp.var', fold] <- 
        max(test.results[[paste0(mode, '.exp.var')]])
      if(length(test.results[[paste0(mode, '.p.cor')]])) {
        all.results$summaries[mode, 'p.cor', fold] <- p_cor(resp, preds)
        all.results$summaries[mode, 'max.p.cor.iter', fold] <- 
          which.max(test.results[[paste0(mode, '.p.cor')]])
        all.results$summaries[mode, 'max.p.cor', fold] <- 
          max(test.results[[paste0(mode, '.p.cor')]])
        all.results$summaries[mode, 's.cor', fold] <- s_cor(resp, preds)
        all.results$summaries[mode, 'max.s.cor.iter', fold] <- 
          which.max(test.results[[paste0(mode, '.s.cor')]])
        all.results$summaries[mode, 'max.s.cor', fold] <- 
          max(test.results[[paste0(mode, '.s.cor')]])
      }
    }
  }

  # If responses are binary, then map to -1 or 1. and get accuracy measures
  if(binary) {
    trained.resp <- trained$resp
    trained.resp[trained.resp < 0] <- -1
    trained.resp[trained.resp > 0] <- 1
    all.results$summaries['train', 'acc', fold] <- mean(train.data$resp == trained.resp, na.rm=T)


    for(mode in c('warm', 'm1', 'm2', 'm3', 'm1m2', 'm1m3', 'm2m3', 'm1m2m3')) {
      if(paste0(mode, '.RMSE') %in% names(test.results)) {
        if(mode=='warm') {
          preds <- warm.preds
        } else {
          preds <- get(paste0(mode, '.pred.resp'))
        }
        resp <- get(paste0(mode, '.resp'))
        resp[resp < 0] <- -1
        resp[resp > 0] <- 1
        pred[pred < 0] <- -1
        pred[pred > 0] <- 1

        all.results$summaries[mode, 'acc', fold] <- mean(resp == pred, na.rm=T)
      }
    }
  }

  # Testing just comparing the number called correctly w/ BaTFLED vs. mean
  # mean(train.data$resp == trained$resp, na.rm=T)
  # mean(train.data$resp == mean.tens.list[[1]], na.rm=T)
  
  all.results$sparsity['m1', fold] <- sum(1/(trained$mode1.lambda.shape * 
    trained$mode1.lambda.scale) > m1.rem.cut)/dim(train.data$mode1.X)[2]
  all.results$sparsity['m2', fold] <- sum(1/(trained$mode2.lambda.shape * 
    trained$mode2.lambda.scale) > m2.rem.cut)/dim(train.data$mode2.X)[2]
  all.results$sparsity['m3', fold] <- sum(1/(trained$mode3.lambda.shape * 
    trained$mode3.lambda.scale) > m3.rem.cut)/dim(train.data$mode3.X)[3]
  all.results$sparsity['core', fold] <- sum(1/(trained$core.lambda.shape * 
    trained$core.lambda.scale) > core.rem.cut)/prod(dim(trained$core.mean))

  # print(paste('Read fold', fld))
  return(all.results)
}

########### MAIN ###############

# Determine the number of runs with this prefix
n.files <- length(list.files(path = dirname(run_prefix),
  pattern = paste0(basename(run_prefix), '.[0-9]+.out')))

f1 <- paste0(run_prefix, '.0/image.Rdata')
tmp <- get_iter_binary(f1)
iters <- tmp$iters
binary <- tmp$binary
####################################################### fix the binary stuff ######
# binary <- F
rm(tmp)

training <- array(NA, dim=c(38, n.files, iters),
  dimnames=list(c('lower.bnd', 'A.RMSE', 'H.RMSE', 'warm.RMSE', 
                  'm1.RMSE', 'm2.RMSE', 'm3.RMSE', 
                  'm1m2.RMSE', 'm1m3.RMSE', 'm2m3.RMSE', 'm1m2m3.RMSE',
                  'exp.var', 'warm.exp.var', 
                  'm1.exp.var', 'm2.exp.var', 'm3.exp.var',
                  'm1m2.exp.var', 'm1m3.exp.var', 'm2m3.exp.var', 'm1m2m3.exp.var',
                  'p.cor', 'warm.p.cor', 
                  'm1.p.cor', 'm2.p.cor', 'm3.p.cor',
                  'm1m2.p.cor', 'm1m3.p.cor', 'm2m3.p.cor', 'm1m2m3.p.cor',
                  's.cor', 'warm.s.cor', 
                  'm1.s.cor', 'm2.s.cor', 'm3.s.cor',
                  'm1m2.s.cor', 'm1m3.s.cor', 'm2m3.s.cor', 'm1m2m3.s.cor'),
                paste0('fold.', 1:n.files), 1:iters))

summaries <- array(NA, dim=c(11, 15, n.files),
  dimnames=list(c('A' ,'H', 'train', 'warm', 'm1', 'm2', 'm3', 
                  'm1m2', 'm1m3', 'm2m3', 'm1m2m3'),
                c('RMSE', 'min.RMSE', 'min.RMSE.iter',
                  'clip.RMSE', 'min.clip.RMSE', 'min.clip.RMSE.iter', 
                  'exp.var', 'max.exp.var.iter', 'max.exp.var',
                  'p.cor', 'max.p.cor.iter', 'max.p.cor',
                  's.cor', 'max.s.cor.iter', 'max.s.cor'),
                paste0('fold.', 1:n.files)))

if(binary) 
  summaries <- array(NA, dim=c(11, 18, n.files),
    dimnames=list(c('A' ,'H', 'train', 'warm', 'm1', 'm2', 'm3',
                    'm1m2', 'm1m3', 'm2m3', 'm1m2m3'),
                  c('RMSE', 'min.RMSE', 'min.RMSE.iter',
                    'clip.RMSE', 'min.clip.RMSE', 'min.clip.RMSE.iter',
                    'exp.var', 'max.exp.var.iter', 'max.exp.var',
                    'p.cor', 'max.p.cor.iter', 'max.p.cor',
                    's.cor', 'max.s.cor.iter', 'max.s.cor',
                    'acc', 'max.acc.iter', 'max.acc'),
                  paste0('fold.', 1:n.files)))

mean <- array(NA, dim=c(4, 21, n.files), 
  dimnames=list(c('RMSE', 'exp.var', 'p.cor', 's.cor'),
                c('train.m1', 'train.m2', 'train.m3',
                  'train.m1m2', 'train.m1m3', 'train.m2m3', 'train.m1m2m3',
                  'warm.m1', 'warm.m2', 'warm.m3', 
                  'warm.m1m2', 'warm.m1m3', 'warm.m2m3', 'warm.m1m2m3',
                  'm1', 'm2', 'm3', 'm1m2', 'm1m3', 'm2m3', 'm1m2m3'),
                paste0('fold.', 1:n.files)))

lower.bnd <- matrix(NA, 4, n.files, dimnames=
  list(c('final', 'max', 'which.max', 'not.mono'),
       paste0('fold.', 1:n.files)))

sparsity <- matrix(NA, 4, n.files,
  dimnames=list(c('m1', 'm2', 'm3', 'core'),
                paste0('fold.', 1:n.files)))

all.results <- list(training=training, summaries=summaries,
   mean=mean, lower.bnd=lower.bnd, sparsity=sparsity)

rm(training, summaries, mean, lower.bnd, sparsity)

for(fld in 1:n.files) {
  # Load in the run data
  f <- paste0(run_prefix, '.', (fld-1), '/image.Rdata')
  all.results <- loadData(f, all.results)
}

results <- all.results

# Remove NA results for tests that weren't performed
# results$training <- results$training[apply(results$training, 1, function(x) sum(!is.na(x))>0),,]
# results$summaries <- results$summaries[apply(results$summaries, 1, function(x) sum(!is.na(x))>0),,]
# results$mean <- results$mean[apply(results$mean, 1, function(x) sum(!is.na(x))>0),,]

# Make data frame counting how many folds peform better than the mean
#####################################################################

# Decide which mean to compare to depending on which test data we have results for
# default to m1
which.mean <- 'm1'
if(sum(!is.na(results$summaries['m1',,]))) which.mean <- 'm1'
if(sum(!is.na(results$summaries['m2',,]))) which.mean <- 'm2'
if(sum(!is.na(results$summaries['m3',,]))) which.mean <- 'm3'
if(sum(!is.na(results$summaries['m1m2',,]))) which.mean <- 'm1m2'
if(sum(!is.na(results$summaries['m1m3',,]))) which.mean <- 'm1m3'
if(sum(!is.na(results$summaries['m2m3',,]))) which.mean <- 'm2m3'
if(sum(!is.na(results$summaries['m1m2m3',,]))) which.mean <- 'm1m2m3'

better <- matrix(NA, dim(results$summaries)[1], 
  length(dimnames(results$summaries)[[2]][!grepl('iter', dimnames(results$summaries)[[2]])]), 
  dimnames=list(dimnames(results$summaries)[[1]],
                dimnames(results$summaries)[[2]][!grepl('iter', dimnames(results$summaries)[[2]])]))

for(type in c('A', 'H', 'train')) for(resp in c('RMSE', 'min.RMSE', 'clip.RMSE', 'min.clip.RMSE')) 
  if(type %in% dimnames(results$summaries)[[1]])
    better[type, resp] <- mean(results$summaries[type, resp,] < results$mean['RMSE', paste0('train.', which.mean),],na.rm=T)
for(resp in c('RMSE', 'min.RMSE', 'clip.RMSE', 'min.clip.RMSE')) 
  better['warm', resp] <- mean(results$summaries['warm', resp,] < results$mean['RMSE', paste0('warm.', which.mean),],na.rm=T)
for(type in c('m1', 'm2', 'm3', 'm1m2', 'm1m3', 'm2m3', 'm1m2m3')) 
  for (resp in c('RMSE', 'min.RMSE', 'clip.RMSE', 'min.clip.RMSE')) 
    better[type, resp] <- mean(results$summaries[type, resp,] < results$mean['RMSE', type,],na.rm=T)
for(type in c('A', 'H', 'train')) for(resp in c('exp.var', 'max.exp.var')) 
  if(type %in% dimnames(results$summaries)[[1]])
    better[type, resp] <- mean(results$summaries[type, resp,] > results$mean['exp.var', paste0('train.', which.mean),],na.rm=T)
for(resp in c('exp.var', 'max.exp.var')) 
  better['warm', resp] <- mean(results$summaries['warm', resp,] > results$mean['exp.var', paste0('warm.', which.mean),],na.rm=T)
for(type in c('m1', 'm2', 'm3', 'm1m2', 'm1m3', 'm2m3', 'm1m2m3')) for (resp in c('exp.var', 'max.exp.var')) 
  better[type, resp] <- mean(results$summaries[type, resp,] > results$mean['exp.var', type,],na.rm=T)
for(type in c('A', 'H', 'train')) for(resp in c('p.cor', 'max.p.cor')) 
  if(type %in% dimnames(results$summaries)[[1]])
    better[type, resp] <- mean(results$summaries[type, resp,] > results$mean['p.cor', paste0('train.', which.mean),],na.rm=T)
for(resp in c('p.cor', 'max.p.cor')) 
  better['warm', resp] <- mean(results$summaries['warm', resp,] > results$mean['p.cor', paste0('warm.', which.mean),],na.rm=T)
for(type in c('m1', 'm2', 'm3', 'm1m2', 'm1m3', 'm2m3', 'm1m2m3')) for (resp in c('p.cor', 'max.p.cor')) 
  better[type, resp] <- mean(results$summaries[type, resp,] > results$mean['p.cor', type,],na.rm=T)
for(type in c('A', 'H', 'train')) for(resp in c('s.cor', 'max.s.cor')) 
  if(type %in% dimnames(results$summaries)[[1]])
    better[type, resp] <- mean(results$summaries[type, resp,] > results$mean['s.cor', paste0('train.', which.mean),],na.rm=T)
for(resp in c('s.cor', 'max.s.cor')) 
  better['warm', resp] <- mean(results$summaries['warm', resp,] > results$mean['s.cor', paste0('warm.', which.mean),],na.rm=T)
for(type in c('m1', 'm2', 'm3', 'm1m2', 'm1m3', 'm2m3', 'm1m2m3')) for (resp in c('s.cor', 'max.s.cor')) 
  better[type, resp] <- mean(results$summaries[type, resp,] > results$mean['s.cor', type,],na.rm=T)

if(binary) {
  for(type in c('A', 'H', 'train')) for(resp in c('acc', 'max.acc'))
    if(type %in% dimnames(results$summaries)[[1]])
      better[type, resp] <- mean(results$summaries[type, resp,] > 0.5, na.rm=T)
  for(resp in c('acc', 'max.acc'))
    better['warm', resp] <- mean(results$summaries['warm', resp,] > .5, na.rm=T)
  for(type in c('m1', 'm2', 'm3', 'm1m2', 'm1m3', 'm2m3', 'm1m2m3')) for (resp in c('acc', 'max.acc'))
    better[type, resp] <- mean(results$summaries[type, resp,] > 0.5, na.rm=T)
}

# Read log files to get run time if the runs finished
if(length(system2('grep', c('"Job terminated"', paste0(run_prefix, '.*.log')), stdout=T)>0)) {
  months <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  month.days <- cumsum(months)
  log.starts <- system2('grep', c('"Job submitted"', paste0(run_prefix, '.*.log')), stdout=T)
  log.ends <- system2('grep', c('"Job terminated"', paste0(run_prefix, '.*.log')), stdout=T)
  month.start <- as.numeric(sapply(strsplit(sapply(strsplit(log.starts, split=' '), 
    '[', 3), split='/'), '[', 1))
  month.end <- as.numeric(sapply(strsplit(sapply(strsplit(log.ends, split=' '), 
    '[', 3), split='/'), '[', 1))
  day.start <- month.days[month.start] + 
    as.numeric(sapply(strsplit(sapply(strsplit(log.starts, split=' '), '[', 3), split='/'), '[', 2))
  day.end <- month.days[month.end] + 
    as.numeric(sapply(strsplit(sapply(strsplit(log.ends, split=' '), '[', 3), split='/'), '[', 2))
  days <- day.end - day.start
  hour.start <- as.numeric(sapply(strsplit(sapply(strsplit(log.starts, split=' '), 
    '[', 4), split=':'), '[', 1))
  hour.end <- as.numeric(sapply(strsplit(sapply(strsplit(log.ends, split=' '), 
    '[', 4), split=':'), '[', 1))
  min.start <- as.numeric(sapply(strsplit(sapply(strsplit(log.starts, split=' '), 
    '[', 4), split=':'), '[', 2))
  min.end <- as.numeric(sapply(strsplit(sapply(strsplit(log.ends, split=' '), 
    '[', 4), split=':'), '[', 2))
  
  hours <- days * 24 + (hour.end - hour.start) + (min.end - min.start)/60
  rm(log.starts, log.ends)
  print("Run time statistics (hours):")
  summary(hours)
}

# Save all the data
save.image(paste0(run_prefix, '_summary.Rdata'))

print("######################################################")
print('## Sparsity ##')
print(apply(results$sparsity, 1, summary))
print("######################################################")
print("Lower bounds")
print(apply(results$lower.bnd, 1, summary, na.rm=T))
print("Standard deviations")
print(apply(results$lower.bnd, 1, sd, na.rm=T))
print("######################################################")
print('## Means ##')
print(apply(results$summaries, c(1,2), mean, na.rm=T))
print('## Standard deviations ##')
print(apply(results$summaries, c(1,2), sd, na.rm=T))
print("######################################################")
print('## Better than mean ##')
print(better)

# Plotting
#####################################################################

cols <- rainbow(n.files)

pdf(file=paste0(run_prefix, '_lower_bounds.pdf'), height=14)
par(mfrow=c(2,1))
plot(results$training['lower.bnd',1,], type='l', lwd=2, col=cols[1], 
  ylim=range(results$training['lower.bnd',,]))
for(i in 2:n.files) points(results$training['lower.bnd',i,], type='l', lwd=2, col=cols[i])
plot(results$training['lower.bnd',1,], type='l', lwd=2, col=cols[1], 
     ylim=range(results$training['lower.bnd',,5:iters]))
for(i in 2:n.files) points(results$training['lower.bnd',i,], type='l', lwd=2, col=cols[i])
dev.off()

######### Barplot of accuracies if binary responses #################
if(binary) {
  pdf(file=paste0(run_prefix, '_acc_bars.pdf'))
  par(mfrow=c(2,1))
  barplot(results$summaries['train', 'acc',], main='training accuracy',
          ylim=c(0,1))
  abline(h=0.5, lty=2, col='red')

  type <- dimnames(results$summaries)[[1]][which(!is.na(results$summaries[,'acc', 1]))[2]]
  barplot(results$summaries[type, 'acc',], main=paste(type, 'accuracy'),
          ylim=c(0,1))
  abline(h=0.5, lty=2, col='red')
}

######### Plot RMSEs ##############

pdf(file=paste0(run_prefix, '_training_RMSE.pdf'))
plot(results$training['A.RMSE',1,] - results$mean['RMSE', 'train.m1',1], 
  type='l', lwd=2, col=cols[1], 
  ylim=range(results$training['A.RMSE',,] - results$mean['RMSE', 'train.m1',]),
  main="Training RMSEs relative \n to predicting mean response",
  xlab="Iteration", ylab="Relative RMSE")
for(i in 2:n.files)
  points(results$training['A.RMSE',i,] - results$mean['RMSE', 'train.m1',i],
     type='l', lwd=2, col=cols[i])
abline(h=0, lty=2, lwd=2)
dev.off()

if(sum(!is.na(results$training['warm.RMSE',1,]))) {
  pdf(file=paste0(run_prefix, '_warm_RMSEs.pdf'))
  plot(results$training['warm.RMSE',1,] - results$mean['RMSE', 'warm.m1',1], 
    type='l', lwd=2, col=cols[1], 
    ylim=range(results$training['warm.RMSE',,] - results$mean['RMSE', 'warm.m1',]),
    main="Warm RMSEs relative to predicting mean response",
    xlab="Iteration", ylab="Relative RMSE")
  for(i in 2:n.files)
    points(results$training['warm.RMSE',i,] - results$mean['RMSE', 'warm.m1',i],
       type='l', lwd=2, col=cols[i])
  abline(h=0, lty=2, lwd=2)
  dev.off()
}
  
if(sum(!is.na(results$training['m1.RMSE',1,]))) {
  pdf(file=paste0(run_prefix, '_cold_m1_RMSEs.pdf'))
  plot(results$training['m1.RMSE',1,] - results$mean['RMSE', 'm1',1], 
    type='l', lwd=2, col=cols[1], 
    ylim=range(results$training['m1.RMSE',,] - results$mean['RMSE', 'm1',]),
    main="Cold mode 1 RMSEs relative to predicting mean response",
    xlab="Iteration", ylab="Relative RMSE")
  for(i in 2:n.files)
    points(results$training['m1.RMSE',i,] - results$mean['RMSE', 'm1',i],
       type='l', lwd=2, col=cols[i])
  abline(h=0, lty=2, lwd=2)
  dev.off()
}

if(sum(!is.na(results$training['m2.RMSE',1,]))) {
  pdf(file=paste0(run_prefix, '_cold_m2_RMSEs.pdf'))
  plot(results$training['m2.RMSE',1,] - results$mean['RMSE', 'm2',1], 
    type='l', lwd=2, col=cols[1], 
    ylim=range(results$training['m2.RMSE',,] - results$mean['RMSE', 'm2',]),
    main="Cold mode 2 RMSEs relative to predicting mean response",
    xlab="Iteration", ylab="Relative RMSE")
  for(i in 2:n.files)
    points(results$training['m2.RMSE',i,] - results$mean['RMSE', 'm2',i],
       type='l', lwd=2, col=cols[i])
  abline(h=0, lty=2, lwd=2)
  dev.off()
}

if(sum(!is.na(results$training['m1m2.RMSE',1,]))) {
  pdf(file=paste0(run_prefix, '_cold_m1m2_RMSEs.pdf'))
  plot(results$training['m1m2.RMSE',1,] - results$mean['RMSE', 'm1m2',1], 
    type='l', lwd=2, col=cols[1], 
    ylim=range(results$training['m1m2.RMSE',,] - results$mean['RMSE', 'm1m2',], na.rm=T),
    main="Cold mode 1/mode 2 RMSEs relative to predicting mean response",
    xlab="Iteration", ylab="Relative RMSE")
  for(i in 2:n.files)
    points(results$training['m1m2.RMSE',i,] - results$mean['RMSE', 'm1m2',i],
       type='l', lwd=2, col=cols[i])
  abline(h=0, lty=2, lwd=2)
  dev.off()
}

########## Plot explained variances ################

pdf(file=paste0(run_prefix, '_training_exp_var.pdf'))
plot(results$training['exp.var',1,] - results$mean['exp.var', 'train.m1',1], 
     type='l', lwd=2, col=cols[1], 
     ylim=range(results$training['exp.var',,] - results$mean['exp.var', 'train.m1',]),
     main="Training explained variance \n relative to predicting mean response",
     xlab="Iteration", ylab="Relative exp. var.")
for(i in 2:n.files)
  points(results$training['exp.var',i,] - results$mean['exp.var', 'train.m1',i],
         type='l', lwd=2, col=cols[i])
abline(h=0, lty=2, lwd=2)
dev.off()

if(sum(!is.na(results$training['warm.exp.var',1,]))) {
  pdf(file=paste0(run_prefix, '_warm_exp_var.pdf'))
  plot(results$training['warm.exp.var',1,] - results$mean['exp.var', 'warm.m1',1], 
       type='l', lwd=2, col=cols[1], 
       ylim=range(results$training['warm.exp.var',,] - results$mean['exp.var', 'warm.m1',]),
       main="Warm explained variance \n relative to predicting mean response",
       xlab="Iteration", ylab="Relative exp. var.")
  for(i in 2:n.files)
    points(results$training['warm.exp.var',i,] - results$mean['exp.var', 'warm.m1',i],
           type='l', lwd=2, col=cols[i])
  abline(h=0, lty=2, lwd=2)
  dev.off()
}

if(sum(!is.na(results$training['m1.exp.var',1,]))) {
  pdf(file=paste0(run_prefix, '_cold_m1_exp_var.pdf'))
  plot(results$training['m1.exp.var',1,] - results$mean['exp.var', 'm1',1], 
       type='l', lwd=2, col=cols[1], 
       ylim=range(results$training['m1.exp.var',,] - results$mean['exp.var', 'm1',]),
       main="Cold mode 1 explained variance \n relative to predicting mean response",
       xlab="Iteration", ylab="Relative exp. var.")
  for(i in 2:n.files)
    points(results$training['m1.exp.var',i,] - results$mean['exp.var', 'm1',i],
           type='l', lwd=2, col=cols[i])
  abline(h=0, lty=2, lwd=2)
  dev.off()
}
  
if(sum(!is.na(results$training['m2.exp.var',1,]))) {
  pdf(file=paste0(run_prefix, '_cold_m2_exp_var.pdf'))
  plot(results$training['m2.exp.var',1,] - results$mean['exp.var', 'm2',1], 
       type='l', lwd=2, col=cols[1], 
       ylim=range(results$training['m2.exp.var',,] - results$mean['exp.var', 'm2',]),
       main="Cold mode 2 explained variance \n relative to predicting mean response",
       xlab="Iteration", ylab="Relative exp. var.")
  for(i in 2:n.files)
    points(results$training['m2.exp.var',i,] - results$mean['exp.var', 'm2',i],
           type='l', lwd=2, col=cols[i])
  abline(h=0, lty=2, lwd=2)
  dev.off()
}

if(sum(!is.na(results$training['m1m2.exp.var',1,]))) {
  pdf(file=paste0(run_prefix, '_cold_m1m2_exp_var.pdf'))
  plot(results$training['m1m2.exp.var',1,] - results$mean['exp.var', 'm1m2',1], 
       type='l', lwd=2, col=cols[1], 
       ylim=range(results$training['m1m2.exp.var',,] - results$mean['exp.var', 'm1m2',], na.rm=T),
       main="Cold mode 1/mode 2 explained variance \n relative to predicting mean response",
       xlab="Iteration", ylab="Relative exp. var.")
  for(i in 2:n.files)
    points(results$training['m1m2.exp.var',i,] - results$mean['exp.var', 'm1m2',i],
           type='l', lwd=2, col=cols[i])
  abline(h=0, lty=2, lwd=2)
  dev.off()
}

########## Plot Pearson correlations ################
if(sum(!is.na(results$training['p.cor',1,]))) {
  pdf(file=paste0(run_prefix, '_training_p_cor.pdf'))
  plot(results$training['p.cor',1,] - results$mean['p.cor', 'train.m1',1], 
       type='l', lwd=2, col=cols[1], 
       ylim=range(results$training['p.cor',,] - results$mean['p.cor', 'train.m1',], na.rm=T),
       main="Training Pearson correlation \n relative to predicting mean response",
       xlab="Iteration", ylab="Relative Pearson correlation")
  for(i in 2:n.files)
    points(results$training['p.cor',i,] - results$mean['p.cor', 'train.m1',i],
           type='l', lwd=2, col=cols[i])
  abline(h=0, lty=2, lwd=2)
  dev.off()
}

if(sum(!is.na(results$training['warm.p.cor',1,]))) {
  pdf(file=paste0(run_prefix, '_warm_p_cor.pdf'))
  plot(results$training['warm.p.cor',1,] - results$mean['p.cor', 'warm.m1',1], 
       type='l', lwd=2, col=cols[1], 
       ylim=range(results$training['warm.p.cor',,] - results$mean['p.cor', 'warm.m1',], na.rm=T),
       main="Warm Pearson correlation \n relative to predicting mean response",
       xlab="Iteration", ylab="Relative Pearson correlation")
  for(i in 2:n.files)
    points(results$training['warm.p.cor',i,] - results$mean['p.cor', 'warm.m1',i],
           type='l', lwd=2, col=cols[i])
  abline(h=0, lty=2, lwd=2)
  dev.off()
}

if(sum(!is.na(results$training['m1.p.cor',1,]))) {
  pdf(file=paste0(run_prefix, '_cold_m1_p_cor.pdf'))
  plot(results$training['m1.p.cor',1,] - results$mean['p.cor', 'm1',1], 
       type='l', lwd=2, col=cols[1], 
       ylim=range(results$training['m1.p.cor',,] - results$mean['p.cor', 'm1',], na.rm=T),
       main="Cold mode 1 Pearson correlation \n relative to predicting mean response",
       xlab="Iteration", ylab="Relative Pearson correlation")
  for(i in 2:n.files)
    points(results$training['m1.p.cor',i,] - results$mean['p.cor', 'm1',i],
           type='l', lwd=2, col=cols[i])
  abline(h=0, lty=2, lwd=2)
  dev.off()
}

if(sum(!is.na(results$training['m2.p.cor',1,]))) {
  pdf(file=paste0(run_prefix, '_cold_m2_p_cor.pdf'))
  plot(results$training['m2.p.cor',1,] - results$mean['p.cor', 'm2',1], 
       type='l', lwd=2, col=cols[1], 
       ylim=range(results$training['m2.p.cor',,] - results$mean['p.cor', 'm2',], na.rm=T),
       main="Cold mode 2 Pearson correlation \n relative to predicting mean response",
       xlab="Iteration", ylab="Relative Pearson correlation")
  for(i in 2:n.files)
    points(results$training['m2.p.cor',i,] - results$mean['p.cor', 'm2',i],
           type='l', lwd=2, col=cols[i])
  abline(h=0, lty=2, lwd=2)
  dev.off()
}

if(sum(!is.na(results$training['m1m2.p.cor',1,]))) {
  pdf(file=paste0(run_prefix, '_cold_m1m2_p_cor.pdf'))
  plot(results$training['m1m2.p.cor',1,] - results$mean['p.cor', 'm1m2',1], 
       type='l', lwd=2, col=cols[1], 
       ylim=range(results$training['m1m2.p.cor',,] - results$mean['p.cor', 'm1m2',], na.rm=T),
       main="Cold mode 1/mode 2 Pearson correlation \n relative to predicting mean response",
       xlab="Iteration", ylab="Relative Pearson correlation")
  for(i in 2:n.files)
    points(results$training['m1m2.p.cor',i,] - results$mean['p.cor', 'm1m2',i],
           type='l', lwd=2, col=cols[i])
  abline(h=0, lty=2, lwd=2)
  dev.off()
}

########## Plot Spearman correlations ################
if(sum(!is.na(results$training['s.cor',1,]))) {
  pdf(file=paste0(run_prefix, '_training_s_cor.pdf'))
  plot(results$training['s.cor',1,] - results$mean['s.cor', 'train.m1',1], 
       type='l', lwd=2, col=cols[1], 
       ylim=range(results$training['s.cor',,] - results$mean['s.cor', 'train.m1',], na.rm=T),
       main="Training Spearman correlation \n relative to predicting mean response",
       xlab="Iteration", ylab="Relative Spearman correlation")
  for(i in 2:n.files)
    points(results$training['s.cor',i,] - results$mean['s.cor', 'train.m1',i],
           type='l', lwd=2, col=cols[i])
  abline(h=0, lty=2, lwd=2)
  dev.off()
}

if(sum(!is.na(results$training['warm.s.cor',1,]))) {
  pdf(file=paste0(run_prefix, '_warm_s_cor.pdf'))
  plot(results$training['warm.s.cor',1,] - results$mean['s.cor', 'warm.m1',1], 
       type='l', lwd=2, col=cols[1], 
       ylim=range(results$training['warm.s.cor',,] - results$mean['s.cor', 'warm.m1',], na.rm=T),
       main="Warm Spearman correlation \n relative to predicting mean response",
       xlab="Iteration", ylab="Relative Spearman correlation")
  for(i in 2:n.files)
    points(results$training['warm.s.cor',i,] - results$mean['s.cor', 'warm.m1',i],
           type='l', lwd=2, col=cols[i])
  abline(h=0, lty=2, lwd=2)
  dev.off()
}

if(sum(!is.na(results$training['m1.s.cor',1,]))) {
  pdf(file=paste0(run_prefix, '_cold_m1_s_cor.pdf'))
  plot(results$training['m1.s.cor',1,] - results$mean['s.cor', 'm1',1], 
       type='l', lwd=2, col=cols[1], 
       ylim=range(results$training['m1.s.cor',,] - results$mean['s.cor', 'm1',], na.rm=T),
       main="Cold mode 1 Spearman correlation \n relative to predicting mean response",
       xlab="Iteration", ylab="Relative Spearman correlation")
  for(i in 2:n.files)
    points(results$training['m1.s.cor',i,] - results$mean['s.cor', 'm1',i],
           type='l', lwd=2, col=cols[i])
  abline(h=0, lty=2, lwd=2)
  dev.off()
}

if(sum(!is.na(results$training['m2.s.cor',1,]))) {
  pdf(file=paste0(run_prefix, '_cold_m2_s_cor.pdf'))
  plot(results$training['m2.s.cor',1,] - results$mean['s.cor', 'm2',1], 
       type='l', lwd=2, col=cols[1], 
       ylim=range(results$training['m2.s.cor',,] - results$mean['s.cor', 'm2',], na.rm=T),
       main="Cold mode 2 Spearman correlation \n relative to predicting mean response",
       xlab="Iteration", ylab="Relative Spearman correlation")
  for(i in 2:n.files)
    points(results$training['m2.s.cor',i,] - results$mean['s.cor', 'm2',i],
           type='l', lwd=2, col=cols[i])
  abline(h=0, lty=2, lwd=2)
  dev.off()
}

if(sum(!is.na(results$training['m1m2.s.cor',1,]))) {
  pdf(file=paste0(run_prefix, '_cold_m1m2_s_cor.pdf'))
  plot(results$training['m1m2.s.cor',1,] - results$mean['s.cor', 'm1m2',1], 
       type='l', lwd=2, col=cols[1], 
       ylim=range(results$training['m1m2.s.cor',,] - results$mean['s.cor', 'm1m2',], na.rm=T),
       main="Cold mode 1/mode 2 Spearman correlation \n relative to predicting mean response",
       xlab="Iteration", ylab="Relative Spearman correlation")
  for(i in 2:n.files)
    points(results$training['m1m2.s.cor',i,] - results$mean['s.cor', 'm1m2',i],
           type='l', lwd=2, col=cols[i])
  abline(h=0, lty=2, lwd=2)
  dev.off()
}

# pdf(file=paste0(run_prefix, '_warm_RMSE_barplots.pdf'), height=14/1.62)
# par(mfrow=c(2,1))
# barplot(matrix(c(warm.rmse.final.list, warm.rmse.for.mean.list),
#                2, n.files, byrow=T), beside = TRUE,
#         names.arg = test.m1.list, las=2,
#         legend.text = c('Pred', 'Mean'),
#         main="Warm RMSEs")
# barplot(matrix(c(warm.rmse.min.list, warm.rmse.for.mean.list),
#                2, n.files, byrow=T), beside = TRUE,
#         names.arg = test.m1.list, las=2,
#         legend.text = c('Pred', 'Mean'),
#         main="Minimum warm RMSEs")
# dev.off()
# 
# pdf(file=paste0(run_prefix, '_cold_RMSE_barplots.pdf'), height=21/1.62)
# par(mfrow=c(3,1))
# barplot(matrix(c(cold.rmse.final.list, cold.rmse.for.mean.list),
#                2, n.files, byrow=T), beside = TRUE,
#         names.arg = test.m1.list, las=2,
#         legend.text = c('Pred', 'Mean'),
#         main="Cold RMSEs")
# barplot(matrix(c(cold.clip.rmse.list, cold.rmse.for.mean.list),
#                2, n.files, byrow=T), beside = TRUE,
#         names.arg = test.m1.list, las=2,
#         legend.text = c('Pred', 'Mean'),
#         main="Cold clipped RMSEs")
# barplot(matrix(c(cold.rmse.min.list, cold.rmse.for.mean.list),
#                2, n.files, byrow=T), beside = TRUE,
#         names.arg = test.m1.list, las=2,
#         legend.text = c('Pred', 'Mean'),
#         main="Minimum cold RMSEs")
# dev.off()

# pdf(file=paste0(run_prefix, '_lb_vs_cold_rmse.pdf'))
# plot(lower.bnd.final.list, cold.rmse.final.list)
# dev.off()
