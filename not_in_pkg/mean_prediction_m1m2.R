#!/usr/bin/env Rscript

# nathan dot lazar at gmail dot com

# Usage: mean_prediction.R  <response data> <mode 1 inputs> <mode2 inputs> 
#                           <mode1 fold matrix> <mode2 fold matrix> <warm.per>

print('Opened file')

library(foreach)        # For parallel processing in loops
library(R6)             # For making memory efficent R6 objects
library(rTensor)        # Needed for multiplying matrices into tensors (could be removed)
library(dplyr)          # General data frame manipulation
library(ggplot2)        # Used for plotting
library(BaTFLED3D)
# devtools::document()    # Builds BaTFLED package
# library(BaTFLED_dev)

print('Packages loaded')
args <- commandArgs(TRUE)

#****** Example running options (using expression and mutation data for COSMIC genes) ******* #
# args <- c(#'../CTRP2/Removed_NAs/log2resp.Rdata',
#           #'../CTRP2/Removed_NAs/param_resp.Rdata',
#           '../CTRP2/Removed_NAs/raw_train_resp.Rdata',
#           # '../CTRP2/Rem_NAs_COSMIC/CTRP2_cl_train_gr.Rdata',
#           '../CTRP2/Removed_NAs/cl_train_mat.Rdata',
#           '../CTRP2/Removed_NAs/dr_train_mat.Rdata',
#           '../CTRP2/Removed_NAs_10fold_cv/cv_folds_cl.Rdata',
#           '../CTRP2/Removed_NAs_10fold_cv/cv_folds_dr.Rdata',
#           '0.01')

# DREAM7 dataset
# args <- list('../DREAM7/DREAM7_resp.Rdata',
#              '../DREAM7/DREAM7_cl.Rdata',
#              '../DREAM7/DREAM7_dr_nodt.Rdata',
#              '../DREAM7/cv_folds_cl_10.Rdata',
#              '../DREAM7/cv_folds_dr_10.Rdata',
#              '0.01')

# Report arguments and summary stats of the input objects
print(unlist(args))

# Read in options
resp.file <-      args[[1]]
m1.file <-        args[[2]] # Used to subset mode 1 responses
m2.file <-        args[[3]] # Used to subset mode 2 responses
m1.fold.file <-   args[[4]]
m2.fold.file <-   args[[5]]
warm.per <-       as.numeric(args[[6]]) # Percentage of warm data to remove

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

  if(verbose) print(sprintf("Explained variance: %.4f", exp.var))
  if(abs(exp.var) < 1e-10) return(0)
  else(return(exp.var))
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

# Load m1 and m2 matrices
m1.mat <- loadRData(m1.file)
m2.mat <- loadRData(m2.file)

# Load the matrices with names for cross validation
m1.cv.fold <- loadRData(m1.fold.file)
m2.cv.fold <- loadRData(m2.fold.file)

folds <- nrow(m1.cv.fold)

# Subset responses using m1.mat and m2.mat
resp.tens <- resp.tens[match(rownames(m1.mat), dimnames(resp.tens)[[1]], nomatch=0),
                       match(rownames(m2.mat), dimnames(resp.tens)[[2]], nomatch=0),
                       ,drop=F]

results <- data.frame(warm.RMSE=rep(0, folds),
                      m1.RMSE=rep(0, folds),
                      m2.RMSE=rep(0, folds),
                      m1m2.RMSE=rep(0, folds),
                      warm.exp.var=rep(0, folds),
                      m1.exp.var=rep(0, folds),
                      m2.exp.var=rep(0, folds),
                      m1m2.exp.var=rep(0, folds),
                      warm.p.cor=rep(0, folds),
                      m1.p.cor=rep(0, folds),
                      m2.p.cor=rep(0, folds),
                      m1m2.p.cor=rep(0, folds),
                      warm.s.cor=rep(0, folds),
                      m1.s.cor=rep(0, folds),
                      m2.s.cor=rep(0, folds),
                      m1m2.s.cor=rep(0, folds))

results.m3 <- array(NA, dim=c(folds, 16, dim(resp.tens)[[3]]),
                    dimnames=list(paste0('fold.', 1:folds),
                                  names(results),
                                  dimnames(resp.tens)[[3]]))

for(fold in 1:folds) {
  # Remove the mode 1 and mode 2 samples for this cross validation fold to be test
  # data for this run
  ###############################################################
  resp.train <- resp.tens[!(dimnames(resp.tens)[[1]] %in% m1.cv.fold[fold,]),
                          !(dimnames(resp.tens)[[2]] %in% m2.cv.fold[fold,]),,drop=F]
  resp.test.m1 <- resp.tens[dimnames(resp.tens)[[1]] %in% m1.cv.fold[fold,],
                            !(dimnames(resp.tens)[[2]] %in% m2.cv.fold[fold,]),,drop=F]
  resp.test.m2 <- resp.tens[!(dimnames(resp.tens)[[1]] %in% m1.cv.fold[fold,]),
                            dimnames(resp.tens)[[2]] %in% m2.cv.fold[fold,],,drop=F]
  resp.test.m1m2 <- resp.tens[dimnames(resp.tens)[[1]] %in% m1.cv.fold[fold,],
                              dimnames(resp.tens)[[2]] %in% m2.cv.fold[fold,],,drop=F]

  # Remove warm test data (per) percent
  all.resp.train <- resp.train
  if(warm.per > 0) {
    mask <- sample(prod(dim(resp.train)[1:2]),
                   round(warm.per*prod(dim(resp.train)[1:2])))
    for(i in 1:dim(resp.train)[3]) {
      resp.train[,,i][mask] <- NA
    }
  }

  # Make a tensor of the means for each value of mode2 & mode3
  m1.means <- apply(resp.train, c(2,3), mean, na.rm=T)
  m1.sds <- apply(resp.train, c(2,3), sd, na.rm=T)
  m2.means <- apply(resp.train, c(1,3), mean, na.rm=T)
  m2.sds <- apply(resp.train, c(1,3), sd, na.rm=T)
  m3.means <- apply(resp.train, c(1,2), mean, na.rm=T)
  m3.sds <- apply(resp.train, c(1,2), sd, na.rm=T)

  train.means <- array(NA, dim=dim(resp.train))
  for(i in 1:dim(resp.train)[1]) train.means[i,,] <- m1.means
  m1.test.means <- array(NA, dim=dim(resp.test.m1))
  for(i in 1:dim(resp.test.m1)[1]) m1.test.means[i,,] <- m1.means
  m2.test.means <- array(NA, dim=dim(resp.test.m2))
  for(j in 1:dim(resp.test.m2)[2]) m2.test.means[,j,] <- m2.means
  m1m2.test.means <- array(0, dim=dim(resp.test.m1m2))
  m1m2.means <- apply(resp.train, 3, mean, na.rm=T)
  for(i in 1:dim(resp.test.m1m2)[1]) 
    for(j in 1:dim(resp.test.m1m2)[2]) 
      m1m2.test.means[i,j,] <- m1m2.means

  # Calculate RMSE, exp.var., Pearson correlation & Spearman correlation for each mode3
  K <- dim(resp.tens)[[3]]
  for(k in 1:K) {
    results.m3[fold, 'warm.RMSE', k] <- rmse(all.resp.train[,,k][is.na(resp.train[,,k])], 
                                             train.means[,,k][is.na(resp.train[,,k])])/sd(resp.train[,,k], na.rm=T)
    results.m3[fold, 'm1.RMSE', k] <- rmse(resp.test.m1[,,k], m1.test.means[,,k])/sd(resp.train[,,k], na.rm=T)
    results.m3[fold, 'm2.RMSE', k] <- rmse(resp.test.m2[,,k] ,m2.test.means[,,k])/sd(resp.train[,,k], na.rm=T)
    results.m3[fold, 'm1m2.RMSE', k] <- rmse(resp.test.m1m2[,,k], m1m2.test.means[,,k])/sd(resp.train[,,k], na.rm=T)
    results.m3[fold, 'warm.exp.var', k] <- exp_var(all.resp.train[,,k][is.na(resp.train[,,k])], train.means[,,k][is.na(resp.train[,,k])])
    results.m3[fold, 'm1.exp.var', k] <- exp_var(resp.test.m1[,,k], m1.test.means[,,k])
    results.m3[fold, 'm2.exp.var', k] <- exp_var(resp.test.m2[,,k], m2.test.means[,,k])
    results.m3[fold, 'm1m2.exp.var', k] <- exp_var(resp.test.m1m2[,,k], m1m2.test.means[,,k])
    results.m3[fold, 'warm.p.cor', k] <- cor(all.resp.train[,,k][is.na(resp.train[,,k])], train.means[,,k][is.na(resp.train[,,k])], use='complete.obs')
    results.m3[fold, 'm1.p.cor', k] <- cor(as.vector(resp.test.m1[,,k]), as.vector(m1.test.means[,,k]), use='complete.obs')
    results.m3[fold, 'm2.p.cor', k] <- cor(as.vector(resp.test.m2[,,k]), as.vector(m2.test.means[,,k]), use='complete.obs')
    results.m3[fold, 'm1m2.p.cor', k] <- cor(as.vector(resp.test.m1m2[,,k]), as.vector(m1m2.test.means[,,k]), use='complete.obs')
    results.m3[fold, 'warm.s.cor', k] <- cor(all.resp.train[,,k][is.na(resp.train[,,k])], train.means[,,k][is.na(resp.train[,,k])], 
                                             use='complete.obs', method='spearman')
    results.m3[fold, 'm1.s.cor', k] <- cor(as.vector(resp.test.m1[,,k]), as.vector(m1.test.means[,,k]), use='complete.obs', method='spearman')
    results.m3[fold, 'm2.s.cor', k] <- cor(as.vector(resp.test.m2[,,k]), as.vector(m2.test.means[,,k]), use='complete.obs', method='spearman')
    results.m3[fold, 'm1m2.s.cor', k] <- cor(as.vector(resp.test.m1m2[,,k]), as.vector(m1m2.test.means[,,k]), use='complete.obs', method='spearman')
  }
  
  # Calculate RMSE, exp.var., Pearson correlation & Spearman correlation
  results$warm.RMSE[fold] <- rmse(all.resp.train[is.na(resp.train)], train.means[is.na(resp.train)])/sd(resp.train, na.rm=T)
  results$m1.RMSE[fold] <- rmse(resp.test.m1, m1.test.means)/sd(resp.train, na.rm=T)
  results$m2.RMSE[fold] <- rmse(resp.test.m2, m2.test.means)/sd(resp.train, na.rm=T)
  results$m1m2.RMSE[fold] <- rmse(resp.test.m1m2, m1m2.test.means)/sd(resp.train, na.rm=T)
  results$warm.exp.var[fold] <- exp_var(all.resp.train[is.na(resp.train)], train.means[is.na(resp.train)])
  results$m1.exp.var[fold] <- exp_var(resp.test.m1, m1.test.means)
  results$m2.exp.var[fold] <- exp_var(resp.test.m2, m2.test.means)
  results$m1m2.exp.var[fold] <- exp_var(resp.test.m1m2, m1m2.test.means)
  results$warm.p.cor[fold] <- cor(all.resp.train[is.na(resp.train)], train.means[is.na(resp.train)], use='complete.obs')
  results$m1.p.cor[fold] <- cor(resp.test.m1, m1.test.means, use='complete.obs')
  results$m2.p.cor[fold] <- cor(resp.test.m2, m2.test.means, use='complete.obs')
  results$m1m2.p.cor[fold] <- cor(resp.test.m1m2, m1m2.test.means, use='complete.obs')
  results$warm.s.cor[fold] <- cor(all.resp.train[is.na(resp.train)], train.means[is.na(resp.train)], use='complete.obs', method='spearman')
  results$m1.s.cor[fold] <- cor(resp.test.m1, m1.test.means, use='complete.obs', method='spearman')
  results$m2.s.cor[fold] <- cor(resp.test.m2, m2.test.means, use='complete.obs', method='spearman')
  results$m1m2.s.cor[fold] <- cor(resp.test.m1m2, m1m2.test.means, use='complete.obs', method='spearman')
}

# print(results)

# print('Results for each mode3 value:')
# print(apply(results.m3, c(2,3), mean))

print(paste('Means:'))
print(apply(results, 2, mean))
print(paste('Standard deviations:'))
print(apply(results, 2, sd))
