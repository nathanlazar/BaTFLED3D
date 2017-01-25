#!/usr/bin/env Rscript

# nathan dot lazar at gmail dot com

# Usage: mean_pred_normalize.R  <response data> 
#  <mode 1 inputs OR 'none'> <mode2 inputs OR 'none'> <mode3 inputs OR 'none'>
#  <mode1 folds OR 'none'> <mode2 folds OR 'none'> <mode3 folds OR 'none'>
#  OPTIONAL: normalize.for=<'mode1', 'mode2', 'mode12', etc.> warm.per=<num>

print('Opened file')

library(dplyr)          # General data frame manipulation
library(ggplot2)        # Used for plotting
library(BaTFLED3D)

print('Packages loaded')
args <- commandArgs(TRUE)

#****** Example running options (using expression and mutation data for COSMIC genes) ******* #
# args <- c(#'../CTRP2/Removed_NAs/log2resp.Rdata',
#           #'../CTRP2/Removed_NAs/param_resp.Rdata',
#           '../CTRP2/Removed_NAs/raw_train_resp.Rdata',
#           # '../CTRP2/Rem_NAs_COSMIC/CTRP2_cl_train_gr.Rdata',
#           '../CTRP2/Removed_NAs/cl_train_mat.Rdata',
#           '../CTRP2/Removed_NAs/dr_train_mat.Rdata',
#           'none',
#           '../CTRP2/Removed_NAs_10fold_cv/cv_folds_cl.Rdata',
#           '../CTRP2/Removed_NAs_10fold_cv/cv_folds_dr.Rdata',
#           'none',
#           'normalize.for=mode1',
#           'warm.per=0.01')

# DREAM7 gi50 data predicting for cold test cell lines
# args <- list('../DREAM7/DREAM7_cl_dr_gi50s_nonorm.Rdata',
#              '../DREAM7/DREAM7_cl.Rdata',
#              '../DREAM7/DREAM7_dr_nodt.Rdata',
#              'none',
#              '../DREAM7/cv_folds_cl_35.Rdata',
#              'none', 'none', 'normalize.for=mode1', 'warm.per=0')

# Report arguments and summary stats of the input objects
print(unlist(args))

# Read in options
resp.file <-      args[[1]]
m1.file <-        args[[2]] # Used to subset mode 1 responses
m2.file <-        args[[3]] # Used to subset mode 2 responses
m3.file <-        args[[4]]
m1.fold.file <-   args[[5]]
m2.fold.file <-   args[[6]]
m3.fold.file <-   args[[7]]
if(sum(grepl('^normalize.for=', args)))
  normalize.for <- sub('normalize.for=', '', args[grepl('^normalize.for=', args)])
if(sum(grepl('^warm.per=', args))) {
  warm.per <- as.numeric(sub('warm.per=', '', args[grepl('^warm.per=', args)]))
} else warm.per <- 0

# Functions
######################################################
rmse <- function(obs, pred) return(sqrt(mean((obs - pred)^2, na.rm=T)))

exp_var <- function(obs, pred, verbose=F) {
  # Get the explained variance for a set of predictions
  if(is.data.frame(obs)) obs <- unlist(obs)
  if(!is.vector(obs)) obs <- as.vector(obs)
  if(is.data.frame(pred)) pred <- unlist(pred)
  if(!is.vector(pred)) pred <- as.vector(pred)
  exp.var <- 1 - var(obs - pred, na.rm=T)/var(obs, na.rm=T)
  if(is.na(exp.var)) return(exp.var)

  if(verbose) print(sprintf("Explained variance: %.4f", exp.var))
  if(abs(exp.var) < 1e-10) return(0)
  else(return(exp.var))
}

corr <- function(obs, pred, ...) {
  if(sum(!is.na(obs + pred))==0) return(NA)
  return(cor(obs, pred, use='complete.obs', ...))
}

# Load data
############################################################
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Load response data
resp.tens <- loadRData(resp.file)

# Load m1, m2 & m3 matrices
if(m1.file != 'none') {
  m1.mat <- loadRData(m1.file)
} else m1.mat <- matrix(0, nrow=dim(resp.tens)[1], ncol=0,
                        dimnames=list(dimnames(resp.tens)[[1]],c()))
if(m2.file != 'none') {
  m2.mat <- loadRData(m2.file)
} else m2.mat <- matrix(0, nrow=dim(resp.tens)[2], ncol=0,
                        dimnames=list(dimnames(resp.tens)[[2]],c()))
if(m3.file != 'none') {
  m3.mat <- loadRData(m3.file)
} else m3.mat <- matrix(0, nrow=dim(resp.tens)[3], ncol=0,
                        dimnames=list(dimnames(resp.tens)[[3]],c()))

# Load the matrices with names for cross validation
if(m1.fold.file != 'none') {
  m1.cv.fold <- loadRData(m1.fold.file)
} else m1.cv.fold <- matrix('', 0, 0)
if(m2.fold.file != 'none') {
  m2.cv.fold <- loadRData(m2.fold.file)
} else m2.cv.fold <- matrix('', 0, 0)
if(m3.fold.file != 'none') {
  m3.cv.fold <- loadRData(m3.fold.file)
} else m3.cv.fold <- matrix('', 0, 0)

m1.folds <- nrow(m1.cv.fold)
m2.folds <- nrow(m2.cv.fold)
m3.folds <- nrow(m3.cv.fold)

# Normalize for the specific prediction task (mode1, 2 or 3)
if(normalize.for == 'mode1') norm.dims <- c(2,3)
if(normalize.for == 'mode2') norm.dims <- c(1,3)
if(normalize.for == 'mode3') norm.dims <- c(1,2)
if(normalize.for == 'mode12') norm.dims <- 3
if(normalize.for == 'mode13') norm.dims <- 2
if(normalize.for == 'mode23') norm.dims <- 1

# Subset responses using m1.mat and m2.mat
resp.tens <- resp.tens[match(rownames(m1.mat), dimnames(resp.tens)[[1]], nomatch=0),
                       match(rownames(m2.mat), dimnames(resp.tens)[[2]], nomatch=0),
                       match(rownames(m3.mat), dimnames(resp.tens)[[3]], nomatch=0),
                       drop=F]

# Preallocate result data frame
folds <- max(m1.folds, 1) * max(m2.folds, 1) * max(m3.folds, 1)
results <- data.frame(warm.m1.RMSE=rep(NA, folds),
                      warm.m2.RMSE=rep(NA, folds),
                      warm.m3.RMSE=rep(NA, folds),
                      m1.RMSE=rep(NA, folds),
                      m2.RMSE=rep(NA, folds),
                      m3.RMSE=rep(NA, folds),
                      m1m2.RMSE=rep(NA, folds),
                      m1m3.RMSE=rep(NA, folds),
                      m2m3.RMSE=rep(NA, folds),
                      m1m2m3.RMSE=rep(NA, folds),
                      warm.m1.exp.var=rep(NA, folds),
                      warm.m2.exp.var=rep(NA, folds),
                      warm.m3.exp.var=rep(NA, folds),
                      m1.exp.var=rep(NA, folds),
                      m2.exp.var=rep(NA, folds),
                      m3.exp.var=rep(NA, folds),
                      m1m2.exp.var=rep(NA, folds),
                      m1m3.exp.var=rep(NA, folds),
                      m2m3.exp.var=rep(NA, folds),
                      m1m2m3.exp.var=rep(NA, folds),
                      warm.m1.p.cor=rep(NA, folds),
                      warm.m2.p.cor=rep(NA, folds),
                      warm.m3.p.cor=rep(NA, folds),
                      m1.p.cor=rep(NA, folds),
                      m2.p.cor=rep(NA, folds),
                      m3.p.cor=rep(NA, folds),
                      m1m2.p.cor=rep(NA, folds),
                      m1m3.p.cor=rep(NA, folds),
                      m2m3.p.cor=rep(NA, folds),
                      m1m2m3.p.cor=rep(NA, folds),
                      warm.m1.s.cor=rep(NA, folds),
                      warm.m2.s.cor=rep(NA, folds),
                      warm.m3.s.cor=rep(NA, folds),
                      m1.s.cor=rep(NA, folds),
                      m2.s.cor=rep(NA, folds),
                      m3.s.cor=rep(NA, folds),
                      m1m2.s.cor=rep(NA, folds),
                      m1m3.s.cor=rep(NA, folds),
                      m2m3.s.cor=rep(NA, folds),
                      m1m2m3.s.cor=rep(NA, folds))

if(m1.folds) {m1.loop <- 1:m1.folds} else m1.loop <- 0
if(m2.folds) {m2.loop <- 1:m2.folds} else m2.loop <- 0
if(m3.folds) {m3.loop <- 1:m3.folds} else m3.loop <- 0

n <- 0
for(m1.fold in m1.loop) 
  for(m2.fold in m2.loop) 
    for(m3.fold in m3.loop) {
  n <- n + 1 
  # Remove the samples for this cross validation fold to be test
  # data for this run
  ##############################################################
  resp.train <- resp.tens[!(dimnames(resp.tens)[[1]] %in% m1.cv.fold[m1.fold,]),
                          !(dimnames(resp.tens)[[2]] %in% m2.cv.fold[m2.fold,]),
                          !(dimnames(resp.tens)[[3]] %in% m3.cv.fold[m3.fold,]),drop=F]

  # Normalize all responses using just the training data mean and sd
  if(exists('norm.dims')) {
    means <- apply(resp.train, norm.dims, mean, na.rm=T)
    sds <- apply(resp.train, norm.dims, sd, na.rm=T)
    # If the standard deviation is zero (e.g. drug PS-1145), don't divide by
    # anything. We don't remove it since something may be learned by other
    # similar samples.
    sds[sds==0] <- 1
    resp.tens <- sweep(resp.tens, norm.dims, means, FUN='-')
    resp.tens <- sweep(resp.tens, norm.dims, sds, FUN='/')
    resp.train <- sweep(resp.train, norm.dims, means, FUN='-')
    resp.train <- sweep(resp.train, norm.dims, sds, FUN='/')
  } else if(normalize.for == 'mode123') {
      means <- mean(resp.train, na.rm=T)
      sds <- sd(resp.train, na.rm=T)
      resp.tens <- resp.tens - means
      resp.tens <- resp.tens/sds
      resp.train <- resp.tens - means
      resp.train <- resp.tens/sds
  }

  resp.test.m1 <- resp.tens[dimnames(resp.tens)[[1]] %in% m1.cv.fold[m1.fold,],
                            !(dimnames(resp.tens)[[2]] %in% m2.cv.fold[m2.fold,]),
                            !(dimnames(resp.tens)[[3]] %in% m3.cv.fold[m3.fold,]),drop=F]
  resp.test.m2 <- resp.tens[!(dimnames(resp.tens)[[1]] %in% m1.cv.fold[m1.fold,]),
                            dimnames(resp.tens)[[2]] %in% m2.cv.fold[m2.fold,],
                            !(dimnames(resp.tens)[[3]] %in% m3.cv.fold[m3.fold,]),drop=F]
  resp.test.m3 <- resp.tens[!(dimnames(resp.tens)[[1]] %in% m1.cv.fold[m1.fold,]),
                            !(dimnames(resp.tens)[[2]] %in% m2.cv.fold[m2.fold,]),
                            dimnames(resp.tens)[[3]] %in% m3.cv.fold[m3.fold,],drop=F]
  resp.test.m1m2 <- resp.tens[dimnames(resp.tens)[[1]] %in% m1.cv.fold[m1.fold,],
                              dimnames(resp.tens)[[2]] %in% m2.cv.fold[m2.fold,],
                              !(dimnames(resp.tens)[[3]] %in% m3.cv.fold[m3.fold,]),drop=F]
  resp.test.m1m3 <- resp.tens[dimnames(resp.tens)[[1]] %in% m1.cv.fold[m1.fold,],
                              !(dimnames(resp.tens)[[2]] %in% m2.cv.fold[m2.fold,]),
                              dimnames(resp.tens)[[3]] %in% m3.cv.fold[m3.fold,],drop=F]
  resp.test.m2m3 <- resp.tens[!(dimnames(resp.tens)[[1]] %in% m1.cv.fold[m1.fold,]),
                              dimnames(resp.tens)[[2]] %in% m2.cv.fold[m2.fold,],
                              dimnames(resp.tens)[[3]] %in% m3.cv.fold[m3.fold,],drop=F]
  resp.test.m1m2m3 <- resp.tens[dimnames(resp.tens)[[1]] %in% m1.cv.fold[m1.fold,],
                                dimnames(resp.tens)[[2]] %in% m2.cv.fold[m2.fold,],
                                dimnames(resp.tens)[[3]] %in% m3.cv.fold[m3.fold,],drop=F]

  # Remove warm test data (per) percent
  all.resp.train <- resp.train
  if(warm.per > 0) {
    mask <- sample(prod(dim(resp.train)[1:2]),
                   round(warm.per*prod(dim(resp.train)[1:2])))
    for(i in 1:dim(resp.train)[3]) {
      resp.train[,,i][mask] <- NA
    }
  }

  # Make tensors of responses if predicting the mean across the training data
  m1.means <- apply(resp.train, c(2,3), mean, na.rm=T)
  m2.means <- apply(resp.train, c(1,3), mean, na.rm=T)
  m3.means <- apply(resp.train, c(1,2), mean, na.rm=T)
  m1m2.means <- apply(resp.train, 3, mean, na.rm=T)
  m1m3.means <- apply(resp.train, 2, mean, na.rm=T)
  m2m3.means <- apply(resp.train, 1, mean, na.rm=T)
  m1m2m3.means <- mean(resp.train, na.rm=T)

  if(warm.per) {
    mean.train.m1 <- array(0, dim=dim(resp.train), dimnames=dimnames(resp.train))
    mean.train.m1 <- sweep(mean.train.m1, c(2,3), m1.means, FUN='+')
    mean.train.m2 <- array(0, dim=dim(resp.train), dimnames=dimnames(resp.train))
    mean.train.m2 <- sweep(mean.train.m2, c(1,3), m2.means, FUN='+')
    mean.train.m3 <- array(0, dim=dim(resp.train), dimnames=dimnames(resp.train))
    mean.train.m3 <- sweep(mean.train.m3, c(1,2), m3.means, FUN='+')
  }
  
  mean.test.m1 <- array(0, dim=dim(resp.test.m1), dimnames=dimnames(resp.test.m1))
  mean.test.m1 <- sweep(mean.test.m1, c(2,3), m1.means, FUN='+')
  mean.test.m2 <- array(0, dim=dim(resp.test.m2), dimnames=dimnames(resp.test.m2))
  mean.test.m2 <- sweep(mean.test.m2, c(1,3), m2.means, FUN='+')
  mean.test.m3 <- array(0, dim=dim(resp.test.m3), dimnames=dimnames(resp.test.m3))
  mean.test.m3 <- sweep(mean.test.m3, c(1,2), m3.means, FUN='+')
  mean.test.m1m2 <- array(0, dim=dim(resp.test.m1m2), dimnames=dimnames(resp.test.m1m2))
  mean.test.m1m2 <- sweep(mean.test.m1m2, 3, m1m2.means, FUN='+')
  mean.test.m1m3 <- array(0, dim=dim(resp.test.m1m3), dimnames=dimnames(resp.test.m1m3))
  mean.test.m1m3 <- sweep(mean.test.m1m3, 2, m1m3.means, FUN='+')
  mean.test.m2m3 <- array(0, dim=dim(resp.test.m2m3), dimnames=dimnames(resp.test.m2m3))
  mean.test.m2m3 <- sweep(mean.test.m2m3, 1, m2m3.means, FUN='+')
  mean.test.m1m2m3 <- array(m1m2m3.means, dim=dim(resp.test.m1m2m3), dimnames=dimnames(resp.test.m1m2m3))

  # Calculate RMSE, exp.var., Pearson correlation & Spearman correlation
  if(warm.per) {
    results$warm.m1.RMSE[n] <- rmse(all.resp.train[is.na(resp.train)], mean.train.m1[is.na(resp.train)])/sd(resp.train, na.rm=T)
    results$warm.m2.RMSE[n] <- rmse(all.resp.train[is.na(resp.train)], mean.train.m2[is.na(resp.train)])/sd(resp.train, na.rm=T)
    results$warm.m3.RMSE[n] <- rmse(all.resp.train[is.na(resp.train)], mean.train.m3[is.na(resp.train)])/sd(resp.train, na.rm=T)
  }
  results$m1.RMSE[n] <- rmse(resp.test.m1, mean.test.m1)/sd(resp.train, na.rm=T)
  results$m2.RMSE[n] <- rmse(resp.test.m2, mean.test.m2)/sd(resp.train, na.rm=T)
  results$m3.RMSE[n] <- rmse(resp.test.m3, mean.test.m3)/sd(resp.train, na.rm=T)
  results$m1m2.RMSE[n] <- rmse(resp.test.m1m2, mean.test.m1m2)/sd(resp.train, na.rm=T)
  results$m1m3.RMSE[n] <- rmse(resp.test.m1m3, mean.test.m1m3)/sd(resp.train, na.rm=T)
  results$m2m3.RMSE[n] <- rmse(resp.test.m2m3, mean.test.m2m3)/sd(resp.train, na.rm=T)
  results$m1m2m3.RMSE[n] <- rmse(resp.test.m1m2m3, mean.test.m1m2m3)/sd(resp.train, na.rm=T)

  if(warm.per) {
    results$warm.m1.exp.var[n] <- exp_var(all.resp.train[is.na(resp.train)], mean.train.m1[is.na(resp.train)])
    results$warm.m2.exp.var[n] <- exp_var(all.resp.train[is.na(resp.train)], mean.train.m2[is.na(resp.train)])
    results$warm.m3.exp.var[n] <- exp_var(all.resp.train[is.na(resp.train)], mean.train.m3[is.na(resp.train)])
  }
  results$m1.exp.var[n] <- exp_var(resp.test.m1, mean.test.m1)
  results$m2.exp.var[n] <- exp_var(resp.test.m2, mean.test.m2)
  results$m3.exp.var[n] <- exp_var(resp.test.m3, mean.test.m3)
  results$m1m2.exp.var[n] <- exp_var(resp.test.m1m2, mean.test.m1m2)
  results$m1m3.exp.var[n] <- exp_var(resp.test.m1m3, mean.test.m1m3)
  results$m2m3.exp.var[n] <- exp_var(resp.test.m2m3, mean.test.m2m3)
  results$m1m2m3.exp.var[n] <- exp_var(resp.test.m1m2m3, mean.test.m1m2m3)

  if(warm.per) {
    results$warm.m1.p.cor[n] <- corr(all.resp.train[is.na(resp.train)], mean.train.m1[is.na(resp.train)])
    results$warm.m2.p.cor[n] <- corr(all.resp.train[is.na(resp.train)], mean.train.m2[is.na(resp.train)])
    results$warm.m3.p.cor[n] <- corr(all.resp.train[is.na(resp.train)], mean.train.m3[is.na(resp.train)])
  }
  results$m1.p.cor[n] <- corr(resp.test.m1, mean.test.m1)
  results$m2.p.cor[n] <- corr(resp.test.m2, mean.test.m2)
  results$m3.p.cor[n] <- corr(resp.test.m3, mean.test.m3)
  results$m1m2.p.cor[n] <- corr(resp.test.m1m2, mean.test.m1m2)
  results$m1m3.p.cor[n] <- corr(resp.test.m1m3, mean.test.m1m3)
  results$m2m3.p.cor[n] <- corr(resp.test.m2m3, mean.test.m2m3)
  results$m1m2m3.p.cor[n] <- corr(resp.test.m1m2m3, mean.test.m1m2m3)

  if(warm.per) {
    results$warm.m1.s.cor[n] <- corr(all.resp.train[is.na(resp.train)], mean.train.m1[is.na(resp.train)], method='spearman')
    results$warm.m2.s.cor[n] <- corr(all.resp.train[is.na(resp.train)], mean.train.m2[is.na(resp.train)], method='spearman')
    results$warm.m3.s.cor[n] <- corr(all.resp.train[is.na(resp.train)], mean.train.m3[is.na(resp.train)], method='spearman')
  }
  results$m1.s.cor[n] <- corr(resp.test.m1, mean.test.m1, method='spearman')
  results$m2.s.cor[n] <- corr(resp.test.m2, mean.test.m2, method='spearman')
  results$m3.s.cor[n] <- corr(resp.test.m3, mean.test.m3, method='spearman')
  results$m1m2.s.cor[n] <- corr(resp.test.m1m2, mean.test.m1m2, method='spearman')
  results$m1m3.s.cor[n] <- corr(resp.test.m1m3, mean.test.m1m3, method='spearman')
  results$m2m3.s.cor[n] <- corr(resp.test.m2m3, mean.test.m2m3, method='spearman')
  results$m1m2m3.s.cor[n] <- corr(resp.test.m1m2m3, mean.test.m1m2m3, method='spearman')
}

# Print table of mean results and sds.
print(paste('Means:'))
print(apply(results, 2, mean))
print(paste('Standard deviations:'))
print(apply(results, 2, sd))

#if(sum(is.nan(apply(results, 2, mean)))!=0) {
#  nan.counts <- apply(results, 2, function(x) sum(is.nan(x)))
#  nan.cols <- which(nan.counts > 0 & nan.counts != nrow(results))
#  print(sprintf('Warning!: %d folds have no responses', sum(is.nan(results[,nan.cols[1]]))))
#}

# Format printing for excel results sheet
result.mat <- matrix(apply(results, 2, mean), nrow=1)
result.mat <- rbind(result.mat, apply(results, 2, sd))
result.vec <- as.vector(result.mat)
tmp <- sapply(names(results), rep, 2)
tmp[1,] <- paste0(tmp[1,], '.mean')
tmp[2,] <- paste0(tmp[2,], '.sd')
names(result.vec) <- as.vector(tmp)
result.vec['warm.m2.RMSE.mean'] <- mean(result.vec['warm.m1.RMSE.mean'],
  result.vec['warm.m2.RMSE.mean'], result.vec['warm.m3.RMSE.mean'], na.rm=T)
result.vec['warm.m2.exp.var.mean'] <- mean(result.vec['warm.m1.exp.var.mean'],
  result.vec['warm.m2.exp.var.mean'], result.vec['warm.m3.exp.var.mean'], na.rm=T)
result.vec['warm.m2.p.cor.mean'] <- mean(result.vec['warm.m1.p.cor.mean'],
  result.vec['warm.m2.p.cor.mean'], result.vec['warm.m3.p.cor.mean'], na.rm=T)
result.vec['warm.m2.s.cor.mean'] <- mean(result.vec['warm.m1.s.cor.mean'],
  result.vec['warm.m2.s.cor.mean'], result.vec['warm.m3.s.cor.mean'], na.rm=T)

result.vec['warm.m2.RMSE.sd'] <- mean(result.vec['warm.m1.RMSE.sd'],
  result.vec['warm.m2.RMSE.sd'], result.vec['warm.m3.RMSE.sd'], na.rm=T)
result.vec['warm.m2.exp.var.sd'] <- mean(result.vec['warm.m1.exp.var.sd'],
  result.vec['warm.m2.exp.var.sd'], result.vec['warm.m3.exp.var.sd'], na.rm=T)
result.vec['warm.m2.p.cor.sd'] <- mean(result.vec['warm.m1.p.cor.sd'],
  result.vec['warm.m2.p.cor.sd'], result.vec['warm.m3.p.cor.sd'], na.rm=T)
result.vec['warm.m2.s.cor.sd'] <- mean(result.vec['warm.m1.s.cor.sd'],
  result.vec['warm.m2.s.cor.sd'], result.vec['warm.m3.s.cor.sd'], na.rm=T)

names(result.vec)[names(result.vec) == 'warm.m1.RMSE.mean'] <- 'train.RMSE.mean'
names(result.vec)[names(result.vec) == 'warm.m1.exp.var.mean'] <- 'train.exp.var.mean'
names(result.vec)[names(result.vec) == 'warm.m1.p.cor.mean'] <- 'train.p.cor.mean'
names(result.vec)[names(result.vec) == 'warm.m1.s.cor.mean'] <- 'train.s.cor.mean'

names(result.vec)[names(result.vec) == 'warm.m1.RMSE.sd'] <- 'train.RMSE.sd'
names(result.vec)[names(result.vec) == 'warm.m1.exp.var.sd'] <- 'train.exp.var.sd'
names(result.vec)[names(result.vec) == 'warm.m1.p.cor.sd'] <- 'train.p.cor.sd'
names(result.vec)[names(result.vec) == 'warm.m1.s.cor.sd'] <- 'train.s.cor.sd'

result.vec[grepl('train', names(result.vec))] <- NA

names(result.vec)[names(result.vec) == 'warm.m2.RMSE.mean'] <- 'warm.RMSE.mean'
names(result.vec)[names(result.vec) == 'warm.m2.exp.var.mean'] <- 'warm.exp.var.mean'
names(result.vec)[names(result.vec) == 'warm.m2.p.cor.mean'] <- 'warm.p.cor.mean'
names(result.vec)[names(result.vec) == 'warm.m2.s.cor.mean'] <- 'warm.s.cor.mean'

names(result.vec)[names(result.vec) == 'warm.m2.RMSE.sd'] <- 'warm.RMSE.sd'
names(result.vec)[names(result.vec) == 'warm.m2.exp.var.sd'] <- 'warm.exp.var.sd'
names(result.vec)[names(result.vec) == 'warm.m2.p.cor.sd'] <- 'warm.p.cor.sd'
names(result.vec)[names(result.vec) == 'warm.m2.s.cor.sd'] <- 'warm.s.cor.sd'

result.vec <- result.vec[!grepl('warm.m3', names(result.vec))]
result.vec[is.nan(result.vec)] <- NA

print(paste(result.vec, sep=',', collapse=','))