#/usr/bin/env Rscript

# Fits 4-parameter logistic curves to real data and predictions and 
# measures how well the predicted gi50s fit the true gi50s

# nathan dot lazar at gmail dot com

# Usage: pred_gi50.R <image.Rdata>

library(methods)
library(BaTFLED3D)
library(drc)

new.args <- commandArgs(TRUE)
test <- F
if(test)
  new.args <- c('HEISERresults/run_2409514.0/all.Rdata')

load(new.args[1])
plot <- F

# Functions
###########################################
gi50 <- function(fit) {
  co <- coef(fit)
  exp((1/co[1])*log((co[3]-co[2])/(.5-co[2])-1) + log(co[4]))
}

get_gi50s <- function(obs, pred, fit.obs=F, plot=F) {
  if(fit.obs) {
    obs.gi50s <- matrix(NA, nrow=dim(obs)[1], ncol=dim(obs)[2], dimnames=dimnames(obs)[1:2])
    obs.fit <- obs.gi50s 
  }
  pred.gi50s <- matrix(NA, nrow=dim(obs)[1], ncol=dim(obs)[2], dimnames=dimnames(obs)[1:2])
  pred.fit <- pred.gi50s
  
  for(i in 1:dim(obs)[1]) {
    for(j in 1:dim(obs)[2]) {
      if(sum(!is.na(obs[i,j,]))>1 && fit.obs) {
        df <- data.frame(resp=obs[i,j,], dose=exp(-(K:1)))
        log4.fit <- tryCatch(drm(resp~dose, data=df, fct = LL.4(), type="continuous"),
                             error = function(b) return(NA))
        if(!is.na(log4.fit) && log4.fit$fit$convergence) {
          obs.fit[i,j] <- log4.fit$fit$value
          obs.gi50s[i,j] <- K + log(gi50(log4.fit))+1
        }
      }
      
      pred.df <- data.frame(resp=pred[i,j,], dose=exp(-(K:1)))
      pred.log4.fit <- tryCatch(drm(resp~dose, data=pred.df, fct = LL.4(), type="continuous"),
                                error = function(b) return(NA))
      if(!is.na(pred.log4.fit) && pred.log4.fit$fit$convergence) {
        pred.fit[i,j] <- pred.log4.fit$fit$value
        pred.gi50s[i,j] <- K + log(gi50(pred.log4.fit)) + 1
      }
      
      if(plot) {
        plot(pred.df$resp, col='red', pch=20,
             main=paste(rownames(pred.gi50s)[i], colnames(pred.gi50s)[j]))
        if(!is.na(pred.log4.fit) && pred.log4.fit$fit$convergence) {
          points(seq(1,K,0.01), pred.log4.fit$curve[[1]](exp(seq(-K,-1,.01))), type='l', col='red')
          abline(v=pred.gi50s[i,j], col='red')
        }
        abline(h=.5, col='blue')
        if(fit.obs) {
          points(df$resp, ylim=c(-1,2), pch=20)
          if(!is.na(log4.fit) && log4.fit$fit$convergence) {
            points(seq(1,K,0.01), log4.fit$curve[[1]](exp(seq(-K,-1,.01))), type='l')
            abline(v=obs.gi50s[i,j])
          }
        }
      }
    }
  }
  if(fit.obs) {
    obs.gi50s[obs.gi50s < 1] <- NA
    obs.gi50s[obs.gi50s > K] <- K+1
    ret <- list(obs.gi50s=obs.gi50s, obs.fit=obs.fit)
  } else ret <- list()
  pred.gi50s[pred.gi50s < 1] <- NA
  pred.gi50s[pred.gi50s > K] <- K+1
  ret$pred.gi50s=pred.gi50s
  ret$pred.fit=pred.fit
  return(ret)
}

p_cor <- function(obs, pred) {
  if((sum(!is.na(obs))>0) & (sum(!is.na(pred))>0)) # Make sure there is data
    if(sum(!is.na(obs) & !is.na(pred)))            # Make sure there are pairs
      return(cor(as.vector(obs), as.vector(pred), use='complete.obs'))
  return(NA)
}

s_cor <- function(obs, pred) {
  if((sum(!is.na(obs))>0) & (sum(!is.na(pred))>0)) # Make sure there is data
    if(sum(!is.na(obs) & !is.na(pred)))            # Make sure there are pairs
      return(cor(as.vector(obs), as.vector(pred), method='spearman', use='complete.obs'))
  return(NA)
}

fns <- list(RMSE=nrmse, exp.var=exp_var, p.cor=p_cor, s.cor=s_cor)

# Fit logistic curves for the training and test data
####################################################
pred.train <- resp.train
pred.train[,,] <- NA
for(k in 1:K) { 
  # Need to rotate matrices to fill in correctly
  tmp <- t(pred.train[,,k])
  tmp[t(!is.na(resp.train[,,k]))] <- train.preds[[k]]
  tmp[t(is.na(resp.train[,,k]))] <- warm.test.preds[[k]]
  pred.train[,,k] <- t(tmp)
}

train.gi.list <- get_gi50s(resp.train, pred.train, fit.obs=T, plot=plot)

if(warm.per==0 && exists('all.resp')) {
  all.resp.train <- all.resp
  warm.per <- mean(!is.na(all.resp.train) & is.na(resp.train))
}

if(warm.per > 0) {
  warm.tens <- all.resp.train
  warm.tens[!is.na(resp.train)] <- NA
  warm.pred.tens <- pred.train
  warm.pred.tens[!is.na(resp.train)] <- NA
  warm.gi.list <- get_gi50s(warm.tens, warm.pred.tens, fit.obs=T, plot=plot)
}

if(exists('test.m1')) {
  m1.pred.resp <- resp.test.m1
  n <- 1
  for(i in 1:dim(resp.test.m1)[1]) for(j in 1:dim(resp.test.m1)[2]) {
    m1.pred.resp[i,j,] <- m1.test.preds[n,]
    n <- n+1
  }
  m1.gi.list <- get_gi50s(resp.test.m1, m1.pred.resp, fit.obs=T)
}
if(exists('test.m2')) {
  m2.pred.resp <- resp.test.m2
  n <- 1
  for(i in 1:dim(resp.test.m2)[1]) for(j in 1:dim(resp.test.m2)[2]) {
    m2.pred.resp[i,j,] <- m2.test.preds[n,]
    n <- n + 1
  }
  m2.gi.list <- get_gi50s(resp.test.m2, m2.pred.resp, fit.obs=T)
}
if(exists('test.m1') && exists('test.m2')) {
  m1m2.pred.resp <- resp.test.m1m2
  n <- 1
  for(i in 1:dim(resp.test.m1m2)[1]) for(j in 1:dim(resp.test.m1m2)[2]) {
    m1m2.pred.resp[i,j,] <- m1m2.test.preds[n,]
    n <- n+1
  }
  m1m2.gi.list <- get_gi50s(resp.test.m1m2, m1m2.pred.resp, fit.obs=T)
}

# Make scatter plots if plot==T
###############################
if(plot) {
  # TODO: make this output a pdf.
  plot_preds(train.gi.list$pred.gi50s, train.gi.list$obs.gi50s, main='Training gi50s')
  if(warm.per)
    plot_preds(warm.gi.list$pred.gi50s, warm.gi.list$obs.gi50s, main='Warm gi50s')
  if(exists('test.m1'))
    plot_preds(m1.gi.list$pred.gi50s, m1.gi.list$obs.gi50s, main='Mode 1 gi50s')
  if(exists('test.m2'))
    plot_preds(m2.gi.list$pred.gi50s, m2.gi.list$obs.gi50s, main='Mode 2 gi50s')
  if(exists('test.m1') && exists('test.m2'))
    plot_preds(m1m2.gi.list$pred.gi50s, m1m2.gi.list$obs.gi50s, main='Mode 1&2 gi50s')
}

# Get response measures for these predictions
#############################################
gi50.measures <- matrix(NA, 5, 4, 
  dimnames=list(c('train', 'warm', 'm1', 'm2', 'm1m2'),
                c('RMSE', 'exp.var', 'p.cor', 's.cor')))

for(i in 1:length(fns)) {
  nm <- names(fns)[[i]]
  fn <- fns[[i]]
  gi50.measures['train', nm] <- fn(train.gi.list$obs.gi50s, train.gi.list$pred.gi50s)
  if(warm.per)
    gi50.measures['warm', nm] <- fn(warm.gi.list$obs.gi50s, warm.gi.list$pred.gi50s)
  if(exists('test.m1'))
    gi50.measures['m1', nm] <- fn(m1.gi.list$obs.gi50s, m1.gi.list$pred.gi50s)
  if(exists('test.m2'))
    gi50.measures['m2', nm] <- fn(m2.gi.list$obs.gi50s, m2.gi.list$pred.gi50s)
  if(exists('test.m1') && exists('test.m2'))
    gi50.measures['m1m2', nm] <- fn(m1m2.gi.list$obs.gi50s, m1m2.gi.list$pred.gi50s)
}

# Get tensors for predicting the mean responses
###############################################
if(!exists('mean.tens.list'))
  mean.tens.list <- list()
if(exists('test.m1'))
  mean.tens.list$m1 <- mean.tens.list$train.m1[1:dim(resp.test.m1)[[1]],,]
if(exists('test.m2'))
  mean.tens.list$m2 <- mean.tens.list$train.m2[,1:dim(resp.test.m2)[[2]],]
if(exists('test.m1') && exists('test.m2')) {
  mean.tens.list$m1m2 <- resp.test.m1m2
  mean.tens.list$m1m2[,,] <- NA
  for(i in 1:I) for(j in 1:J) mean.tens.list$m1m2 <- apply(resp.train, c(3), mean, na.rm=T)
}

# Get gi50s for mean prediction
###############################
mean.tens.list <- list(train.m1=resp.train, train.m2=resp.train)
for(i in 1:I) mean.tens.list$train.m1[i,,] <- m1.means
for(j in 1:J) mean.tens.list$train.m2[,j,] <- m2.means

mean.m1.train.gi.list <- get_gi50s(resp.train, mean.tens.list$train.m1, fit.obs=F, plot=plot)
mean.m2.train.gi.list <- get_gi50s(resp.train, mean.tens.list$train.m2, fit.obs=F, plot=plot)

if(warm.per > 0) {
  mean.m1.warm.gi.list <- mean.m1.train.gi.list
  mean.m1.warm.gi.list$pred.gi50s[apply(!is.na(resp.train), c(1,2), sum)>0] <- NA
  mean.m1.warm.gi.list$pred.fit[apply(!is.na(resp.train), c(1,2), sum)>0] <- NA
  mean.m2.warm.gi.list <- mean.m2.train.gi.list
  mean.m2.warm.gi.list$pred.gi50s[apply(!is.na(resp.train), c(1,2), sum)>0] <- NA
  mean.m2.warm.gi.list$pred.fit[apply(!is.na(resp.train), c(1,2), sum)>0] <- NA
}

if(exists('test.m1'))
  mean.m1.gi.list <- get_gi50s(resp.test.m1, 
                               mean.tens.list$train.m1[1:dim(resp.test.m1)[1],,])
if(exists('test.m2'))
  mean.m2.gi.list <- get_gi50s(resp.test.m2, 
                               mean.tens.list$train.m2[,1:dim(resp.test.m2)[2],])
if(exists('test.m1') && exists('test.m2'))
  mean.m1m2.gi.list <- get_gi50s(resp.test.m1m2, 
    mean.tens.list$train.m1[1:dim(resp.test.m1)[1],1:dim(resp.test.m2)[2],])

# Make scatter plots for mean prediction if plot==T
###################################################
if(plot) {
  # TODO: make this output a pdf.
  plot_preds(mean.m1.train.gi.list$pred.gi50s, train.gi.list$obs.gi50s, main='Training mean m1 gi50s')
  plot_preds(mean.m2.train.gi.list$pred.gi50s, train.gi.list$obs.gi50s, main='Training mean m2 gi50s')
  if(warm.per) {
    plot_preds(mean.m1.warm.gi.list$pred.gi50s, warm.gi.list$obs.gi50s, main='Warm mean m1 gi50s')
    plot_preds(mean.m2.warm.gi.list$pred.gi50s, warm.gi.list$obs.gi50s, main='Warm mean m2 gi50s')
  }
  if(exists('test.m1'))
    plot_preds(mean.m1.gi.list$pred.gi50s, m1.gi.list$obs.gi50s, main='Mode 1 gi50s')
  if(exists('test.m2'))
    plot_preds(mean.m2.gi.list$pred.gi50s, m2.gi.list$obs.gi50s, main='Mode 2 gi50s')
  if(exists('test.m1') && exists('test.m2'))
    plot_preds(mean.m1m2.gi.list$pred.gi50s, m1m2.gi.list$obs.gi50s, main='Mode 1&2 gi50s')
}

# Get measures when predicting the mean
#######################################
gi50.mean.measures <- matrix(NA, 7, 4, 
  dimnames=list(c('train.m1', 'train.m2', 'warm.m1', 'warm.m2', 'm1', 'm2', 'm1m2'),
                c('RMSE', 'exp.var', 'p.cor', 's.cor')))

for(i in 1:length(fns)) {
  nm <- names(fns)[[i]]
  fn <- fns[[i]]
  gi50.mean.measures['train.m1', nm] <- fn(train.gi.list$obs.gi50s, mean.m1.train.gi.list$pred.gi50s)
  gi50.mean.measures['train.m2', nm] <- fn(train.gi.list$obs.gi50s, mean.m2.train.gi.list$pred.gi50s)
  if(warm.per) {
    gi50.mean.measures['warm.m1', nm] <- fn(warm.gi.list$obs.gi50s, mean.m1.warm.gi.list$pred.gi50s)
    gi50.mean.measures['warm.m2', nm] <- fn(warm.gi.list$obs.gi50s, mean.m2.warm.gi.list$pred.gi50s)
  }
  if(exists('test.m1'))
    gi50.mean.measures['m1', nm] <- fn(m1.gi.list$obs.gi50s, mean.m1.gi.list$pred.gi50s)
  if(exists('test.m2'))
    gi50.mean.measures['m2', nm] <- fn(m2.gi.list$obs.gi50s, mean.m2.gi.list$pred.gi50s)
  if(exists('test.m1') && exists('test.m2'))
    gi50.mean.measures['m1m2', nm] <- fn(m1m2.gi.list$obs.gi50s, mean.m1m2.gi.list$pred.gi50s)
}

# Print results and save
print('*** Results for predicting the mean ***')
print(gi50.mean.measures)
print('')
print('*** Results for model ***')
print(gi50.measures)

save.image(paste0(out.dir, 'image.Rdata'))