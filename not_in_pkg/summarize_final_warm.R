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
rmse <- function(obs, pred) sqrt(mean((obs-pred)^2,na.rm=T))

# Read in a file to get the number of iterations
get_iters <- function(f1) {
  load(f1)
  ret <- list()
  return(trained$iter)
}

# Function to load in data from runs
loadData <- function(f, results){
  #loads an RData file, and returns a list with the 
  #trained model and warm & cold RMSE vectors.

  load(f)

  # Get std. devs of responses
  stddev <- sd(train.data$resp, na.rm=T)

  # Make a tensor of the means for each value of mode2 & mode3 for comparisons
  # TODO: Add stuff for m3
  m1.mean.tens <- train.data$resp
  m2.mean.tens <- train.data$resp
  m1m2.mean.tens <- train.data$resp
  for(i in 1:dim(m1.mean.tens)[1]) m1.mean.tens[i,,] <- m1.means
  for(j in 1:dim(m2.mean.tens)[2]) m2.mean.tens[,j,] <- m2.means
  for(k in 1:dim(m1m2.mean.tens)[3]) m1m2.mean.tens[,,k] <- m3.means

  results$mean['RMSE', 'train.m1', fold] <- rmse(train.data$resp, m1.mean.tens)/stddev
  results$mean['RMSE', 'train.m2', fold] <- rmse(train.data$resp, m2.mean.tens)/stddev
  results$mean['RMSE', 'train.m1m2', fold] <- rmse(train.data$resp, m1m2.mean.tens)/stddev
  results$mean['RMSE', 'warm.m1', fold] <- rmse(all.resp.train[is.na(train.data$resp)],
    m1.mean.tens[is.na(train.data$resp)])/stddev
  results$mean['RMSE', 'warm.m2', fold] <- rmse(all.resp.train[is.na(train.data$resp)],
    m2.mean.tens[is.na(train.data$resp)])/stddev
  results$mean['RMSE', 'warm.m1m2', fold] <- rmse(all.resp.train[is.na(train.data$resp)], 
    m1m2.mean.tens[is.na(train.data$resp)])/stddev
  if(exists('test.data.m1'))
    results$mean['RMSE', 'm1', fold] <- rmse(test.data.m1$resp,
      m1.mean.tens[1:dim(test.data.m1$resp)[1],,])/stddev
  if(exists('test.data.m2'))
    results$mean['RMSE', 'm2', fold] <- rmse(test.data.m2$resp,
      m2.mean.tens[,1:dim(test.data.m2$resp)[2],])/stddev
  if(exists('test.data.m1m2'))
    results$mean['RMSE', 'm1m2', fold] <- rmse(test.data.m1m2$resp,
      m1m2.mean.tens[1:dim(test.data.m1m2$resp)[1], 1:dim(test.data.m1m2$resp)[2],])/stddev

  results$mean['exp.var', 'train.m1', fold] <- exp_var(train.data$resp, m1.mean.tens)
  results$mean['exp.var', 'train.m2', fold] <- exp_var(train.data$resp, m2.mean.tens)
  results$mean['exp.var', 'train.m1m2', fold] <- exp_var(train.data$resp, m1m2.mean.tens)
  results$mean['exp.var', 'warm.m1', fold] <- exp_var(all.resp.train[is.na(train.data$resp)],
    m1.mean.tens[is.na(train.data$resp)])
  results$mean['exp.var', 'warm.m2', fold] <- exp_var(all.resp.train[is.na(train.data$resp)],
    m2.mean.tens[is.na(train.data$resp)])
  results$mean['exp.var', 'warm.m1m2', fold] <- exp_var(all.resp.train[is.na(train.data$resp)], 
    m1m2.mean.tens[is.na(train.data$resp)])
  if(exists('test.data.m1'))
    results$mean['exp.var', 'm1', fold] <- exp_var(test.data.m1$resp,
      m1.mean.tens[1:dim(test.data.m1$resp)[1],,])
  if(exists('test.data.m2'))
    results$mean['exp.var', 'm2', fold] <- exp_var(test.data.m2$resp,
      m2.mean.tens[,1:dim(test.data.m2$resp)[2],])
  if(exists('test.data.m1m2'))
    results$mean['exp.var', 'm1m2', fold] <- exp_var(test.data.m1m2$resp,
      m1m2.mean.tens[1:dim(test.data.m1m2$resp)[1], 1:dim(test.data.m1m2$resp)[2],])

  results$mean['p.cor', 'train.m1', fold] <- cor(train.data$resp, m1.mean.tens, use='complete.obs')
  results$mean['p.cor', 'train.m2', fold] <- cor(train.data$resp, m2.mean.tens, use='complete.obs')
  results$mean['p.cor', 'train.m1m2', fold] <- cor(train.data$resp, m1m2.mean.tens, use='complete.obs')
  results$mean['p.cor', 'warm.m1', fold] <- cor(all.resp.train[is.na(train.data$resp)],
    m1.mean.tens[is.na(train.data$resp)], use='complete.obs')
  results$mean['p.cor', 'warm.m2', fold] <- cor(all.resp.train[is.na(train.data$resp)],
    m2.mean.tens[is.na(train.data$resp)], use='complete.obs')
  results$mean['p.cor', 'warm.m1m2', fold] <- tryCatch(cor(all.resp.train[is.na(train.data$resp)], 
    m1m2.mean.tens[is.na(train.data$resp)], use='complete.obs'), error=function(e) NA)
  if(exists('test.data.m1'))
    results$mean['p.cor', 'm1', fold] <- cor(test.data.m1$resp,
      m1.mean.tens[1:dim(test.data.m1$resp)[1],,], use='complete.obs')
  if(exists('test.data.m2'))
    results$mean['p.cor', 'm2', fold] <- cor(test.data.m2$resp,
      m2.mean.tens[,1:dim(test.data.m2$resp)[2],], use='complete.obs')
  if(exists('test.data.m1m2'))
    results$mean['p.cor', 'm1m2', fold] <- cor(test.data.m1m2$resp,
      m1m2.mean.tens[1:dim(test.data.m1m2$resp)[1], 1:dim(test.data.m1m2$resp)[2],], use='complete.obs')
  
  results$mean['s.cor', 'train.m1', fold] <- cor(train.data$resp, m1.mean.tens, use='complete.obs', method='spearman')
  results$mean['s.cor', 'train.m2', fold] <- cor(train.data$resp, m2.mean.tens, use='complete.obs', method='spearman')
  results$mean['s.cor', 'train.m1m2', fold] <- cor(train.data$resp, m1m2.mean.tens, use='complete.obs', method='spearman')
  results$mean['s.cor', 'warm.m1', fold] <- cor(all.resp.train[is.na(train.data$resp)],
    m1.mean.tens[is.na(train.data$resp)], use='complete.obs', method='spearman')
  results$mean['s.cor', 'warm.m2', fold] <- cor(all.resp.train[is.na(train.data$resp)],
    m2.mean.tens[is.na(train.data$resp)], use='complete.obs', method='spearman')
  results$mean['s.cor', 'warm.m1m2', fold] <- tryCatch(cor(all.resp.train[is.na(train.data$resp)], 
    m1m2.mean.tens[is.na(train.data$resp)], use='complete.obs', method='spearman'), error=function(e) NA)
  if(exists('test.data.m1'))
    results$mean['s.cor', 'm1', fold] <- cor(test.data.m1$resp,
      m1.mean.tens[1:dim(test.data.m1$resp)[1],,], use='complete.obs', method='spearman')
  if(exists('test.data.m2'))
    results$mean['s.cor', 'm2', fold] <- cor(test.data.m2$resp,
      m2.mean.tens[,1:dim(test.data.m2$resp)[2],], use='complete.obs', method='spearman')
  if(exists('test.data.m1m2'))
    results$mean['s.cor', 'm1m2', fold] <- cor(test.data.m1m2$resp,
      m1m2.mean.tens[1:dim(test.data.m1m2$resp)[1], 1:dim(test.data.m1m2$resp)[2],], use='complete.obs', method='spearman')

  results$training['lower.bnd', fold, 1:trained$iter] <-
    trained$lower.bnd

  results$lower.bnd['final', fold] <- trained$lower.bnd[trained$iter]
  results$lower.bnd['max', fold] <- max(trained$lower.bnd)
  results$lower.bnd['which.max', fold] <- which.max(trained$lower.bnd)
  results$lower.bnd['not.mono', fold] <- sum((trained$lower.bnd[-1] -
                              trained$lower.bnd[-trained$iter]) < 0)

  results$training['A.RMSE', fold, 1:trained$iter] <- trained$RMSE
  results$training['H.RMSE', fold, 1:trained$iter] <- trained$H.RMSE
  results$training['exp.var', fold, 1:trained$iter] <- trained$exp.var
  if(length(trained$p.cor))
    results$training['p.cor', fold, 1:trained$iter] <- trained$p.cor
  if(length(trained$s.cor))
    results$training['s.cor', fold, 1:trained$iter] <- trained$s.cor
  if('warm.RMSE' %in% names(test.results)) {
    results$training['warm.RMSE', fold, 1:trained$iter] <- test.results$warm.RMSE
    results$training['warm.exp.var', fold, 1:trained$iter] <- test.results$warm.exp.var
    if(length(test.results$warm.p.cor))
      results$training['warm.p.cor', fold, 1:trained$iter] <- test.results$warm.p.cor
    if(length(test.results$warm.s.cor))
      results$training['warm.s.cor', fold, 1:trained$iter] <- test.results$warm.s.cor
  }
  if('m1.RMSE' %in% names(test.results)) {
    results$training['m1.RMSE', fold, 1:trained$iter] <- test.results$m1.RMSE
    results$training['m1.exp.var', fold, 1:trained$iter] <- test.results$m1.exp.var
    if(length(test.results$m1.p.cor))
      results$training['m1.p.cor', fold, 1:trained$iter] <- test.results$m1.p.cor
    if(length(test.results$m1.s.cor))
      results$training['m1.s.cor', fold, 1:trained$iter] <- test.results$m1.s.cor
  }
  if('m2.RMSE' %in% names(test.results)) {
    results$training['m2.RMSE', fold, 1:trained$iter] <- test.results$m2.RMSE
    results$training['m2.exp.var', fold, 1:trained$iter] <- test.results$m2.exp.var
    if(length(test.results$m2.p.cor))
      results$training['m2.p.cor', fold, 1:trained$iter] <- test.results$m2.p.cor
    if(length(test.results$m2.s.cor))
      results$training['m2.s.cor', fold, 1:trained$iter] <- test.results$m2.s.cor
  }
  if('m3.RMSE' %in% names(test.results)) {
    results$training['m3.RMSE', fold, 1:trained$iter] <- test.results$m3.RMSE
    results$training['m3.exp.var', fold, 1:trained$iter] <- test.results$m3.exp.var
    if(length(test.results$m3.p.cor))
      results$training['m3.p.cor', fold, 1:trained$iter] <- test.results$m3.p.cor
    if(length(test.results$m3.s.cor))
      results$training['m3.s.cor', fold, 1:trained$iter] <- test.results$m3.s.cor
  }
  if('m1m2.RMSE' %in% names(test.results)) {
    results$training['m1m2.RMSE', fold, 1:trained$iter] <- test.results$m1m2.RMSE
    results$training['m1m2.exp.var', fold, 1:trained$iter] <- test.results$m1m2.exp.var
    if(length(test.results$m1m2.p.cor))
      results$training['m1m2.p.cor', fold, 1:trained$iter] <- test.results$m1m2.p.cor
    if(length(test.results$m1m2.s.cor))
      results$training['m1m2.s.cor', fold, 1:trained$iter] <- test.results$m1m2.s.cor
  }
  if('m1m3.RMSE' %in% names(test.results)) {
    results$training['m1m3.RMSE', fold, 1:trained$iter] <- test.results$m1m3.RMSE
    results$training['m1m3.exp.var', fold, 1:trained$iter] <- test.results$m1m3.exp.var
    if(length(test.results$m1m3.p.cor))
      results$training['m1m3.p.cor', fold, 1:trained$iter] <- test.results$m1m3.p.cor
    if(length(test.results$m1m3.s.cor))
      results$training['m1m3.s.cor', fold, 1:trained$iter] <- test.results$m1m3.s.cor
  }
  if('m2m3.RMSE' %in% names(test.results)) {
    results$training['m2m3.RMSE', fold, 1:trained$iter] <- test.results$m2m3.RMSE
    results$training['m2m3.exp.var', fold, 1:trained$iter] <- test.results$m2m3.exp.var
    if(length(test.results$m2m3.p.cor))
      results$training['m2m3.p.cor', fold, 1:trained$iter] <- test.results$m2m3.p.cor
    if(length(test.results$m2m3.s.cor))
      results$training['m2m3.s.cor', fold, 1:trained$iter] <- test.results$m2m3.s.cor
  }

  if('m1m2m3.RMSE' %in% names(test.results)) {
    results$training['m1m2m3.RMSE', fold, 1:trained$iter] <- test.results$m1m2m3.RMSE
    results$training['m1m2m3.exp.var', fold, 1:trained$iter] <- test.results$m1m2m3.exp.var
    if(length(test.results$m1m2m3.p.cor))
      results$training['m1m2m3.p.cor', fold, 1:trained$iter] <- test.results$m1m2m3.p.cor
    if(length(test.results$m1m2m3.s.cor))
      results$training['m1m2m3.s.cor', fold, 1:trained$iter] <- test.results$m1m2m3.s.cor
  }
  
  results$summaries['H', 'RMSE', fold] <- trained$H.RMSE[trained$iter]
  results$summaries['H', 'min.RMSE.iter', fold] <- which.min(trained$H.RMSE)
  results$summaries['H', 'min.RMSE', fold] <- min(trained$H.RMSE)

  results$summaries['train', 'RMSE', fold] <- trained$RMSE[trained$iter]
  results$summaries['train', 'min.RMSE.iter', fold] <- which.min(trained$RMSE)
  results$summaries['train', 'min.RMSE', fold] <- min(trained$RMSE)
  results$summaries['train', 'exp.var', fold] <- trained$exp.var[trained$iter]
  results$summaries['train', 'max.exp.var.iter', fold] <- which.max(trained$exp.var)
  results$summaries['train', 'max.exp.var', fold] <- max(trained$exp.var)
  if(length(trained$p.cor)) {
    results$summaries['train', 'p.cor', fold] <- trained$p.cor[trained$iter]
    results$summaries['train', 'max.p.cor.iter', fold] <- which.max(trained$p.cor)
    results$summaries['train', 'max.p.cor', fold] <- max(trained$p.cor)
    results$summaries['train', 's.cor', fold] <- trained$s.cor[trained$iter]
    results$summaries['train', 'max.s.cor.iter', fold] <- which.max(trained$s.cor)
    results$summaries['train', 'max.s.cor', fold] <- max(trained$s.cor)
  }

  if('warm.RMSE' %in% names(test.results)) {
    results$summaries['warm', 'RMSE', fold] <- test.results$warm.RMSE[trained$iter]
    results$summaries['warm', 'min.RMSE.iter', fold] <- which.min(test.results$warm.RMSE)
    results$summaries['warm', 'min.RMSE', fold] <- min(test.results$warm.RMSE)
    results$summaries['warm', 'clip.RMSE', fold] <- test.results$warm.RMSE.clip[trained$iter]
    results$summaries['warm', 'min.clip.RMSE.iter', fold] <- which.min(test.results$warm.RMSE.clip)
    results$summaries['warm', 'min.clip.RMSE', fold] <- min(test.results$warm.RMSE.clip)
    results$summaries['warm', 'exp.var', fold] <- test.results$warm.exp.var[trained$iter]
    results$summaries['warm', 'max.exp.var.iter', fold] <- which.max(test.results$warm.exp.var)
    results$summaries['warm', 'max.exp.var', fold] <- max(test.results$warm.exp.var)
    if(length(test.results$warm.p.cor)) {
      results$summaries['warm', 'p.cor', fold] <- test.results$warm.p.cor[trained$iter]
      results$summaries['warm', 'max.p.cor.iter', fold] <- which.max(test.results$warm.p.cor)
      results$summaries['warm', 'max.p.cor', fold] <- max(test.results$warm.p.cor)
      results$summaries['warm', 's.cor', fold] <- test.results$warm.s.cor[trained$iter]
      results$summaries['warm', 'max.s.cor.iter', fold] <- which.max(test.results$warm.s.cor)
      results$summaries['warm', 'max.s.cor', fold] <- max(test.results$warm.s.cor)
    }
  }

  if('m1.RMSE' %in% names(test.results)) {
    results$summaries['m1', 'RMSE', fold] <- test.results$m1.RMSE[trained$iter]
    results$summaries['m1', 'min.RMSE.iter', fold] <- which.min(test.results$m1.RMSE)
    results$summaries['m1', 'min.RMSE', fold] <- min(test.results$m1.RMSE)
    results$summaries['m1', 'clip.RMSE', fold] <- test.results$m1.RMSE.clip[trained$iter]
    results$summaries['m1', 'min.clip.RMSE.iter', fold] <- which.min(test.results$m1.RMSE.clip)
    results$summaries['m1', 'min.clip.RMSE', fold] <- min(test.results$m1.RMSE.clip)
    results$summaries['m1', 'exp.var', fold] <- test.results$m1.exp.var[trained$iter]
    results$summaries['m1', 'max.exp.var.iter', fold] <- which.max(test.results$m1.exp.var)
    results$summaries['m1', 'max.exp.var', fold] <- max(test.results$m1.exp.var)
    if(length(test.results$m1.p.cor)) {
      results$summaries['m1', 'p.cor', fold] <- test.results$m1.p.cor[trained$iter]
      results$summaries['m1', 'max.p.cor.iter', fold] <- which.max(test.results$m1.p.cor)
      results$summaries['m1', 'max.p.cor', fold] <- max(test.results$m1.p.cor)
      results$summaries['m1', 's.cor', fold] <- test.results$m1.s.cor[trained$iter]
      results$summaries['m1', 'max.s.cor.iter', fold] <- which.max(test.results$m1.s.cor)
      results$summaries['m1', 'max.s.cor', fold] <- max(test.results$m1.s.cor)
    }
  }

  if('m2.RMSE' %in% names(test.results)) {
    results$summaries['m2', 'RMSE', fold] <- test.results$m2.RMSE[trained$iter]
    results$summaries['m2', 'min.RMSE.iter', fold] <- which.min(test.results$m2.RMSE)
    results$summaries['m2', 'min.RMSE', fold] <- min(test.results$m2.RMSE)
    results$summaries['m2', 'clip.RMSE', fold] <- test.results$m2.RMSE.clip[trained$iter]
    results$summaries['m2', 'min.clip.RMSE.iter', fold] <- which.min(test.results$m2.RMSE.clip)
    results$summaries['m2', 'min.clip.RMSE', fold] <- min(test.results$m2.RMSE.clip)
    results$summaries['m2', 'exp.var', fold] <- test.results$m2.exp.var[trained$iter]
    results$summaries['m2', 'max.exp.var.iter', fold] <- which.max(test.results$m2.exp.var)
    results$summaries['m2', 'max.exp.var', fold] <- max(test.results$m2.exp.var)
    if(length(test.results$m2.p.cor)) {
      results$summaries['m2', 'p.cor', fold] <- test.results$m2.p.cor[trained$iter]
      results$summaries['m2', 'max.p.cor.iter', fold] <- which.max(test.results$m2.p.cor)
      results$summaries['m2', 'max.p.cor', fold] <- max(test.results$m2.p.cor)
      results$summaries['m2', 's.cor', fold] <- test.results$m2.s.cor[trained$iter]
      results$summaries['m2', 'max.s.cor.iter', fold] <- which.max(test.results$m2.s.cor)
      results$summaries['m2', 'max.s.cor', fold] <- max(test.results$m2.s.cor)
    }
  }

  if('m3.RMSE' %in% names(test.results)) {
    results$summaries['m3', 'RMSE', fold] <- test.results$m3.RMSE[trained$iter]
    results$summaries['m3', 'min.RMSE.iter', fold] <- which.min(test.results$m3.RMSE)
    results$summaries['m3', 'min.RMSE', fold] <- min(test.results$m3.RMSE)
    results$summaries['m3', 'clip.RMSE', fold] <- test.results$m3.RMSE.clip[trained$iter]
    results$summaries['m3', 'min.clip.RMSE.iter', fold] <- which.min(test.results$m3.RMSE.clip)
    results$summaries['m3', 'min.clip.RMSE', fold] <- min(test.results$m3.RMSE.clip)
    results$summaries['m3', 'exp.var', fold] <- test.results$m3.exp.var[trained$iter]
    results$summaries['m3', 'max.exp.var.iter', fold] <- which.max(test.results$m3.exp.var)
    results$summaries['m3', 'max.exp.var', fold] <- max(test.results$m3.exp.var)
    if(length(test.results$m3.p.cor)) {
      results$summaries['m3', 'p.cor', fold] <- test.results$m3.p.cor[trained$iter]
      results$summaries['m3', 'max.p.cor.iter', fold] <- which.max(test.results$m3.p.cor)
      results$summaries['m3', 'max.p.cor', fold] <- max(test.results$m3.p.cor)
      results$summaries['m3', 's.cor', fold] <- test.results$m3.s.cor[trained$iter]
      results$summaries['m3', 'max.s.cor.iter', fold] <- which.max(test.results$m3.s.cor)
      results$summaries['m3', 'max.s.cor', fold] <- max(test.results$m3.s.cor)
    }
  }

  if('m1m2.RMSE' %in% names(test.results)) {
    results$summaries['m1m2', 'RMSE', fold] <- test.results$m1m2.RMSE[trained$iter]
    results$summaries['m1m2', 'min.RMSE.iter', fold] <- which.min(test.results$m1m2.RMSE)
    results$summaries['m1m2', 'min.RMSE', fold] <- min(test.results$m1m2.RMSE)
    results$summaries['m1m2', 'clip.RMSE', fold] <- test.results$m1m2.RMSE.clip[trained$iter]
    results$summaries['m1m2', 'min.clip.RMSE.iter', fold] <- which.min(test.results$m1m2.RMSE.clip)
    results$summaries['m1m2', 'min.clip.RMSE', fold] <- min(test.results$m1m2.RMSE.clip)
    results$summaries['m1m2', 'exp.var', fold] <- test.results$m1m2.exp.var[trained$iter]
    results$summaries['m1m2', 'max.exp.var.iter', fold] <- which.max(test.results$m1m2.exp.var)
    results$summaries['m1m2', 'max.exp.var', fold] <- max(test.results$m1m2.exp.var)
    if(length(test.results$m1m2.p.cor)) {
      results$summaries['m1m2', 'p.cor', fold] <- test.results$m1m2.p.cor[trained$iter]
      results$summaries['m1m2', 'max.p.cor.iter', fold] <- which.max(test.results$m1m2.p.cor)
      results$summaries['m1m2', 'max.p.cor', fold] <- max(test.results$m1m2.p.cor)
      results$summaries['m1m2', 's.cor', fold] <- test.results$m1m2.s.cor[trained$iter]
      results$summaries['m1m2', 'max.s.cor.iter', fold] <- which.max(test.results$m1m2.s.cor)
      results$summaries['m1m2', 'max.s.cor', fold] <- max(test.results$m1m2.s.cor)
    }
  }

  if('m1m3.RMSE' %in% names(test.results)) {
    results$summaries['m1m3', 'RMSE', fold] <- test.results$m1m3.RMSE[trained$iter]
    results$summaries['m1m3', 'min.RMSE.iter', fold] <- which.min(test.results$m1m3.RMSE)
    results$summaries['m1m3', 'min.RMSE', fold] <- min(test.results$m1m3.RMSE)
    results$summaries['m1m3', 'clip.RMSE', fold] <- test.results$m1m3.RMSE.clip[trained$iter]
    results$summaries['m1m3', 'min.clip.RMSE.iter', fold] <- which.min(test.results$m1m3.RMSE.clip)
    results$summaries['m1m3', 'min.clip.RMSE', fold] <- min(test.results$m1m3.RMSE.clip)
    results$summaries['m1m3', 'exp.var', fold] <- test.results$m1m3.exp.var[trained$iter]
    results$summaries['m1m3', 'max.exp.var.iter', fold] <- which.max(test.results$m1m3.exp.var)
    results$summaries['m1m3', 'max.exp.var', fold] <- max(test.results$m1m3.exp.var)
    if(length(test.results$m1m3.p.cor)) {
      results$summaries['m1m3', 'p.cor', fold] <- test.results$m1m3.p.cor[trained$iter]
      results$summaries['m1m3', 'max.p.cor.iter', fold] <- which.max(test.results$m1m3.p.cor)
      results$summaries['m1m3', 'max.p.cor', fold] <- max(test.results$m1m3.p.cor)
      results$summaries['m1m3', 's.cor', fold] <- test.results$m1m3.s.cor[trained$iter]
      results$summaries['m1m3', 'max.s.cor.iter', fold] <- which.max(test.results$m1m3.s.cor)
      results$summaries['m1m3', 'max.s.cor', fold] <- max(test.results$m1m3.s.cor)
    }
  }

  if('m2m3.RMSE' %in% names(test.results)) {
    results$summaries['m2m3', 'RMSE', fold] <- test.results$m2m3.RMSE[trained$iter]
    results$summaries['m2m3', 'min.RMSE.iter', fold] <- which.min(test.results$m2m3.RMSE)
    results$summaries['m2m3', 'min.RMSE', fold] <- min(test.results$m2m3.RMSE)
    results$summaries['m2m3', 'clip.RMSE', fold] <- test.results$m2m3.RMSE.clip[trained$iter]
    results$summaries['m2m3', 'min.clip.RMSE.iter', fold] <- which.min(test.results$m2m3.RMSE.clip)
    results$summaries['m2m3', 'min.clip.RMSE', fold] <- min(test.results$m2m3.RMSE.clip)
    results$summaries['m2m3', 'exp.var', fold] <- test.results$m2m3.exp.var[trained$iter]
    results$summaries['m2m3', 'max.exp.var.iter', fold] <- which.max(test.results$m2m3.exp.var)
    results$summaries['m2m3', 'max.exp.var', fold] <- max(test.results$m2m3.exp.var)
    if(length(test.results$m2m3.p.cor)) {
      results$summaries['m2m3', 'p.cor', fold] <- test.results$m2m3.p.cor[trained$iter]
      results$summaries['m2m3', 'max.p.cor.iter', fold] <- which.max(test.results$m2m3.p.cor)
      results$summaries['m2m3', 'max.p.cor', fold] <- max(test.results$m2m3.p.cor)
      results$summaries['m2m3', 's.cor', fold] <- test.results$m2m3.s.cor[trained$iter]
      results$summaries['m2m3', 'max.s.cor.iter', fold] <- which.max(test.results$m2m3.s.cor)
      results$summaries['m2m3', 'max.s.cor', fold] <- max(test.results$m2m3.s.cor)
    }
  }

  if('m1m2m3.RMSE' %in% names(test.results)) {
    results$summaries['m1m2m3', 'RMSE', fold] <- test.results$m1m2m3.RMSE[trained$iter]
    results$summaries['m1m2m3', 'min.RMSE.iter', fold] <- which.min(test.results$m1m2m3.RMSE)
    results$summaries['m1m2m3', 'min.RMSE', fold] <- min(test.results$m1m2m3.RMSE)
    results$summaries['m1m2m3', 'clip.RMSE', fold] <- test.results$m1m2m3.RMSE.clip[trained$iter]
    results$summaries['m1m2m3', 'min.clip.RMSE.iter', fold] <- which.min(test.results$m1m2m3.RMSE.clip)
    results$summaries['m1m2m3', 'min.clip.RMSE', fold] <- min(test.results$m1m2m3.RMSE.clip)
    results$summaries['m1m2m3', 'exp.var', fold] <- test.results$m1m2m3.exp.var[trained$iter]
    results$summaries['m1m2m3', 'max.exp.var.iter', fold] <- which.max(test.results$m1m2m3.exp.var)
    results$summaries['m1m2m3', 'max.exp.var', fold] <- max(test.results$m1m2m3.exp.var)
    if(length(test.results$m1m2m3.p.cor)) {
      results$summaries['m1m2m3', 'p.cor', fold] <- test.results$m1m2m3.p.cor[trained$iter]
      results$summaries['m1m2m3', 'max.p.cor.iter', fold] <- which.max(test.results$m1m2m3.p.cor)
      results$summaries['m1m2m3', 'max.p.cor', fold] <- max(test.results$m1m2m3.p.cor)
      results$summaries['m1m2m3', 's.cor', fold] <- test.results$m1m2m3.s.cor[trained$iter]
      results$summaries['m1m2m3', 'max.s.cor.iter', fold] <- which.max(test.results$m1m2m3.s.cor)
      results$summaries['m1m2m3', 'max.s.cor', fold] <- max(test.results$m1m2m3.s.cor)
    }
  }
  
  results$sparsity['m1', fold] <- sum(1/(trained$mode1.lambda.shape * 
    trained$mode1.lambda.scale) > m1.rem.cut)/dim(train.data$mode1.X)[2]
  results$sparsity['m2', fold] <- sum(1/(trained$mode2.lambda.shape * 
    trained$mode2.lambda.scale) > m2.rem.cut)/dim(train.data$mode2.X)[2]
  results$sparsity['m3', fold] <- sum(1/(trained$mode3.lambda.shape * 
    trained$mode3.lambda.scale) > m3.rem.cut)/dim(train.data$mode2.X)[3]
  results$sparsity['core', fold] <- sum(1/(trained$core.lambda.shape * 
    trained$core.lambda.scale) > core.rem.cut)/prod(dim(trained$core.mean))

  return(results)
}

########### MAIN ###############

# Determine the number of runs with this prefix
n.files <- length(list.files(path = dirname(run_prefix),
  pattern = paste0(basename(run_prefix), '.[0-9]+.out')))

f1 <- paste0(run_prefix, '.0/image.Rdata')
iters <- get_iters(f1)

training <- array(NA, dim=c(22, n.files, iters),
  dimnames=list(c('lower.bnd', 'A.RMSE', 'H.RMSE', 'warm.RMSE', 
                  'm1.RMSE', 'm2.RMSE', 'm1m2.RMSE',
                  'exp.var', 'warm.exp.var', 'm1.exp.var', 
                  'm2.exp.var', 'm1m2.exp.var', 
                  'p.cor', 'warm.p.cor', 'm1.p.cor', 'm2.p.cor', 'm1m2.p.cor',
                  's.cor', 'warm.s.cor', 'm1.s.cor', 'm2.s.cor', 'm1m2.s.cor'),
                paste0('fold.', 1:n.files), 1:iters))

summaries <- array(NA, dim=c(7, 15, n.files),
  dimnames=list(c('A' ,'H', 'train', 'warm', 'm1', 'm2', 'm1m2'),
                c('RMSE', 'min.RMSE', 'min.RMSE.iter',
                  'clip.RMSE', 'min.clip.RMSE', 'min.clip.RMSE.iter', 
                  'exp.var', 'max.exp.var.iter', 'max.exp.var',
                  'p.cor', 'max.p.cor.iter', 'max.p.cor',
                  's.cor', 'max.s.cor.iter', 'max.s.cor'),
                paste0('fold.', 1:n.files)))

mean <- array(NA, dim=c(4, 9, n.files), 
  dimnames=list(c('RMSE', 'exp.var', 'p.cor', 's.cor'),
                c('train.m1', 'train.m2', 'train.m1m2', 
                  'warm.m1', 'warm.m2', 'warm.m1m2', 
                  'm1', 'm2', 'm1m2'),
                paste0('fold.', 1:n.files)))

lower.bnd <- matrix(NA, 4, n.files, dimnames=
  list(c('final', 'max', 'which.max', 'not.mono'),
       paste0('fold.', 1:n.files)))

sparsity <- matrix(NA, 4, n.files,
  dimnames=list(c('m1', 'm2', 'm3', 'core'),
                paste0('fold.', 1:n.files)))

results <- list(training=training, summaries=summaries,
   mean=mean, lower.bnd=lower.bnd, sparsity=sparsity)

rm(training, summaries, mean, lower.bnd, sparsity)

for(fld in 1:n.files) {
  # Load in the run data
  f <- paste0(run_prefix, '.', (fld-1), '/image.Rdata')
  results <- loadData(f, results)
}

# Make data frame counting how many folds peform better than the mean
better <- matrix(NA, 7, 10, dimnames=list(dimnames(results$summaries)[[1]],
  dimnames(results$summaries)[[2]][!grepl('iter', dimnames(results$summaries)[[2]])]))

for(type in c('A', 'H', 'train')) for(resp in c('RMSE', 'min.RMSE', 'clip.RMSE', 'min.clip.RMSE')) 
  better[type, resp] <- sum(results$summaries[type, resp,] < results$mean['RMSE', 'train.m1',])
for(resp in c('RMSE', 'min.RMSE', 'clip.RMSE', 'min.clip.RMSE')) 
  better['warm', resp] <- sum(results$summaries['warm', resp,] < results$mean['RMSE', 'warm.m1',])
for(type in c('m1', 'm2', 'm1m2')) for (resp in c('RMSE', 'min.RMSE', 'clip.RMSE', 'min.clip.RMSE')) 
  better[type, resp] <- sum(results$summaries[type, resp,] < results$mean['RMSE', type,])
for(type in c('A', 'H', 'train')) for(resp in c('exp.var', 'max.exp.var')) 
  better[type, resp] <- sum(results$summaries[type, resp,] > results$mean['exp.var', 'train.m1',])
for(resp in c('exp.var', 'max.exp.var')) 
  better['warm', resp] <- sum(results$summaries['warm', resp,] > results$mean['exp.var', 'warm.m1',])
for(type in c('m1', 'm2', 'm1m2')) for (resp in c('exp.var', 'max.exp.var')) 
  better[type, resp] <- sum(results$summaries[type, resp,] > results$mean['exp.var', type,])
for(type in c('A', 'H', 'train')) for(resp in c('p.cor', 'max.p.cor')) 
  better[type, resp] <- sum(results$summaries[type, resp,] > results$mean['p.cor', 'train.m1',])
for(resp in c('p.cor', 'max.p.cor')) 
  better['warm', resp] <- sum(results$summaries['warm', resp,] > results$mean['p.cor', 'warm.m1',])
for(type in c('m1', 'm2', 'm1m2')) for (resp in c('p.cor', 'max.p.cor')) 
  better[type, resp] <- sum(results$summaries[type, resp,] > results$mean['p.cor', type,])
for(type in c('A', 'H', 'train')) for(resp in c('s.cor', 'max.s.cor')) 
  better[type, resp] <- sum(results$summaries[type, resp,] > results$mean['s.cor', 'train.m1',])
for(resp in c('s.cor', 'max.s.cor')) 
  better['warm', resp] <- sum(results$summaries['warm', resp,] > results$mean['s.cor', 'warm.m1',])
for(type in c('m1', 'm2', 'm1m2')) for (resp in c('s.cor', 'max.s.cor')) 
  better[type, resp] <- sum(results$summaries[type, resp,] > results$mean['s.cor', type,])

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

########## Plot Pearson correlations ################

pdf(file=paste0(run_prefix, '_training_p_cor.pdf'))
plot(results$training['p.cor',1,] - results$mean['p.cor', 'train.m1',1], 
     type='l', lwd=2, col=cols[1], 
     ylim=range(results$training['p.cor',,] - results$mean['p.cor', 'train.m1',]),
     main="Training Pearson correlation \n relative to predicting mean response",
     xlab="Iteration", ylab="Relative Pearson correlation")
for(i in 2:n.files)
  points(results$training['p.cor',i,] - results$mean['p.cor', 'train.m1',i],
         type='l', lwd=2, col=cols[i])
abline(h=0, lty=2, lwd=2)
dev.off()

pdf(file=paste0(run_prefix, '_warm_p_cor.pdf'))
plot(results$training['warm.p.cor',1,] - results$mean['p.cor', 'warm.m1',1], 
     type='l', lwd=2, col=cols[1], 
     ylim=range(results$training['warm.p.cor',,] - results$mean['p.cor', 'warm.m1',]),
     main="Warm Pearson correlation \n relative to predicting mean response",
     xlab="Iteration", ylab="Relative Pearson correlation")
for(i in 2:n.files)
  points(results$training['warm.p.cor',i,] - results$mean['p.cor', 'warm.m1',i],
         type='l', lwd=2, col=cols[i])
abline(h=0, lty=2, lwd=2)
dev.off()

pdf(file=paste0(run_prefix, '_cold_m1_p_cor.pdf'))
plot(results$training['m1.p.cor',1,] - results$mean['p.cor', 'm1',1], 
     type='l', lwd=2, col=cols[1], 
     ylim=range(results$training['m1.p.cor',,] - results$mean['p.cor', 'm1',]),
     main="Cold mode 1 Pearson correlation \n relative to predicting mean response",
     xlab="Iteration", ylab="Relative Pearson correlation")
for(i in 2:n.files)
  points(results$training['m1.p.cor',i,] - results$mean['p.cor', 'm1',i],
         type='l', lwd=2, col=cols[i])
abline(h=0, lty=2, lwd=2)
dev.off()

pdf(file=paste0(run_prefix, '_cold_m2_p_cor.pdf'))
plot(results$training['m2.p.cor',1,] - results$mean['p.cor', 'm2',1], 
     type='l', lwd=2, col=cols[1], 
     ylim=range(results$training['m2.p.cor',,] - results$mean['p.cor', 'm2',]),
     main="Cold mode 2 Pearson correlation \n relative to predicting mean response",
     xlab="Iteration", ylab="Relative Pearson correlation")
for(i in 2:n.files)
  points(results$training['m2.p.cor',i,] - results$mean['p.cor', 'm2',i],
         type='l', lwd=2, col=cols[i])
abline(h=0, lty=2, lwd=2)
dev.off()

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

########## Plot Spearman correlations ################

pdf(file=paste0(run_prefix, '_training_s_cor.pdf'))
plot(results$training['s.cor',1,] - results$mean['s.cor', 'train.m1',1], 
     type='l', lwd=2, col=cols[1], 
     ylim=range(results$training['s.cor',,] - results$mean['s.cor', 'train.m1',]),
     main="Training Spearman correlation \n relative to predicting mean response",
     xlab="Iteration", ylab="Relative Spearman correlation")
for(i in 2:n.files)
  points(results$training['s.cor',i,] - results$mean['s.cor', 'train.m1',i],
         type='l', lwd=2, col=cols[i])
abline(h=0, lty=2, lwd=2)
dev.off()

pdf(file=paste0(run_prefix, '_warm_s_cor.pdf'))
plot(results$training['warm.s.cor',1,] - results$mean['s.cor', 'warm.m1',1], 
     type='l', lwd=2, col=cols[1], 
     ylim=range(results$training['warm.s.cor',,] - results$mean['s.cor', 'warm.m1',]),
     main="Warm Spearman correlation \n relative to predicting mean response",
     xlab="Iteration", ylab="Relative Spearman correlation")
for(i in 2:n.files)
  points(results$training['warm.s.cor',i,] - results$mean['s.cor', 'warm.m1',i],
         type='l', lwd=2, col=cols[i])
abline(h=0, lty=2, lwd=2)
dev.off()

pdf(file=paste0(run_prefix, '_cold_m1_s_cor.pdf'))
plot(results$training['m1.s.cor',1,] - results$mean['s.cor', 'm1',1], 
     type='l', lwd=2, col=cols[1], 
     ylim=range(results$training['m1.s.cor',,] - results$mean['s.cor', 'm1',]),
     main="Cold mode 1 Spearman correlation \n relative to predicting mean response",
     xlab="Iteration", ylab="Relative Spearman correlation")
for(i in 2:n.files)
  points(results$training['m1.s.cor',i,] - results$mean['s.cor', 'm1',i],
         type='l', lwd=2, col=cols[i])
abline(h=0, lty=2, lwd=2)
dev.off()

pdf(file=paste0(run_prefix, '_cold_m2_s_cor.pdf'))
plot(results$training['m2.s.cor',1,] - results$mean['s.cor', 'm2',1], 
     type='l', lwd=2, col=cols[1], 
     ylim=range(results$training['m2.s.cor',,] - results$mean['s.cor', 'm2',]),
     main="Cold mode 2 Spearman correlation \n relative to predicting mean response",
     xlab="Iteration", ylab="Relative Spearman correlation")
for(i in 2:n.files)
  points(results$training['m2.s.cor',i,] - results$mean['s.cor', 'm2',i],
         type='l', lwd=2, col=cols[i])
abline(h=0, lty=2, lwd=2)
dev.off()

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
