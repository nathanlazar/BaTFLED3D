#/usr/bin/env Rscript

# nathan dot lazar at gmail dot com

# Usage: BaTFLED3D_post.R <run.Rdata>

library(methods)
library(BaTFLED3D)      # Note: package must be installed from .tar file
# library(pryr)           # For monitoring memory usage
# library(microbenchmark) # For monitoring run time

print('Packages loaded')
args <- commandArgs(TRUE)

load(args[[1]])

# Transform everything back to the original space and then get performance 
# measures....
##########################################################################
train.resp <- train.data$resp
all.train.resp <- all.resp.train
train.pred.resp <- trained$resp
if(params$decomp=='Tucker') {
  train.H.resp <- mult_3d(trained$core.mean, trained$mode1.H.mean,
                          trained$mode2.H.mean, trained$mode3.H.mean)
} else 
  train.H.resp <- train.pred.resp

if(exists('norm.dims')) {
  resp.tens <- sweep(resp.tens, norm.dims, sds, FUN='*')
  resp.tens <- sweep(resp.tens, norm.dims, means, FUN='+')
  train.resp <- sweep(train.resp, norm.dims, sds, FUN='*')
  train.resp <- sweep(train.resp, norm.dims, means, FUN='+')
  all.train.resp <- sweep(all.train.resp, norm.dims, sds, FUN='*')
  all.train.resp <- sweep(all.train.resp, norm.dims, means, FUN='+')
  train.pred.resp <- sweep(train.pred.resp, norm.dims, sds, FUN='*')
  train.pred.resp <- sweep(train.pred.resp, norm.dims, means, FUN='+')
  train.H.resp <- sweep(train.H.resp, norm.dims, sds, FUN='*')
  train.H.resp <- sweep(train.H.resp, norm.dims, means, FUN='+')
}

for(m in c('m1', 'm2', 'm3', 'm1m2', 'm1m3', 'm2m3', 'm1m2m3')) {
  d <- get(paste0('test.data.', m))
  if(length(d)) {
    obs <- d$resp
    pred <- test(d, trained)
    if(exists('norm.dims')) {
      pred <- sweep(pred, norm.dims, sds, FUN='*')
      pred <- sweep(pred, norm.dims, means, FUN='+')
      obs <- sweep(obs, norm.dims, sds, FUN='*')
      obs <- sweep(obs, norm.dims, means, FUN='+')
    } else if(normalize.for == 'mode123') {
      pred <- pred * sds
      pred <- pred + means
      obs <- obs * sds
      obs <- obs + means
    }
    assign(paste0(m, '.pred.resp'), pred)
    assign(paste0(m, '.resp'), obs)
  }
}

## TODO ###################################################
# For parameters, need to transform back to original space, 
# since they've been log transformed etc.
###########################################################

# Get tensors of mean predictions for comparisons
#####################################################
m1.means <- apply(train.resp, c(2,3), mean, na.rm=T)
m1.sds <- apply(train.resp, c(2,3), sd, na.rm=T)
m2.means <- apply(train.resp, c(1,3), mean, na.rm=T)
m2.sds <- apply(train.resp, c(1,3), sd, na.rm=T)
m3.means <- apply(train.resp, c(1,2), mean, na.rm=T)
m3.sds <- apply(train.resp, c(1,2), sd, na.rm=T)

mean.tens.list <- list()
mean.tens.list[['train.m1']] <- train.resp
for(i in 1:dim(train.data$resp)[1]) mean.tens.list[['train.m1']][i,,] <- m1.means
mean.tens.list[['train.m2']] <- train.resp
for(j in 1:dim(train.data$resp)[2]) mean.tens.list[['train.m2']][,j,] <- m2.means
mean.tens.list[['train.m3']] <- train.resp
for(k in 1:dim(train.data$resp)[3]) mean.tens.list[['train.m3']][,,k] <- m3.means

if(length(test.data.m1)) {
  mean.tens.list[['m1']] <- test.data.m1$resp
  for(i in 1:dim(test.data.m1$resp)[1]) mean.tens.list[['m1']][i,,] <- m1.means
}
if(length(test.data.m2)) {
  mean.tens.list[['m2']] <- test.data.m2$resp
  for(j in 1:dim(test.data.m2$resp)[2]) mean.tens.list[['m2']][,j,] <- m2.means
}
if(length(test.data.m3)) {
  mean.tens.list[['m3']] <- test.data.m3$resp
  for(k in 1:dim(test.data.m3$resp)[3]) mean.tens.list[['m3']][,,k] <- m3.means
}
if(length(test.data.m1m2)) {
  mean.tens.list[['m1m2']] <- array(0, dim=dim(test.data.m1m2$resp))
  mean.tens.list[['m1m2']] <- apply(train.data$resp, 3, mean, na.rm=T)
}
if(length(test.data.m1m3)) {
  mean.tens.list[['m1m3']] <- array(0, dim=dim(test.data.m1m3$resp))
  mean.tens.list[['m1m3']] <- apply(train.data$resp, 2, mean, na.rm=T)
}
if(length(test.data.m2m3)) {
  mean.tens.list[['m2m3']] <- array(0, dim=dim(test.data.m2m3$resp))
  mean.tens.list[['m2m3']] <- apply(train.data$resp, 1, mean, na.rm=T)
}
if(length(test.data.m1m2m3)) {
  mean.tens.list[['m1m2m3']] <- array(0, dim=dim(test.data.m1m2m3$resp))
  mean.tens.list[['m1m2m3']][,,] <- mean(train.data$resp, na.rm=T)
}

results <- list()
results$mean <- matrix(NA, 4, 13,
  dimnames=list(c('RMSE', 'exp.var', 'p.cor', 's.cor'),
                c('train.m1', 'train.m2', 'train.m3',
                  'warm.m1', 'warm.m2', 'warm.m3',
                  'm1', 'm2', 'm3', 'm1m2', 'm1m3', 'm2m3', 'm1m2m3')))
results$trained <- matrix(NA, 4, 10,
  dimnames=list(c('RMSE', 'exp.var', 'p.cor', 's.cor'),
                c('train', 'train.H', 'warm',
                  'm1', 'm2', 'm3', 'm1m2', 'm1m3', 'm2m3', 'm1m2m3')))

p_cor <- function(obs, pred) {
  if(sum(!is.na(obs) & !is.na(pred))) {
    return(cor(obs, pred, use='complete.obs'))
  } else return(NA)
}
s_cor <- function(obs, pred) {
  if(sum(!is.na(obs) & !is.na(pred))) {
    return(cor(obs, pred, method='spearman', use='complete.obs'))
  } else return(NA)
}

fns <- list(RMSE=nrmse, exp.var=exp_var, p.cor=p_cor, s.cor=s_cor)

# Get results for predicting mean responses.
for(m in 1:length(fns)) {
  fn <- fns[[m]]
  name <- names(fns)[m]
  # Get results for training data
  for(mode in c('m1','m2','m3'))
    results$mean[name, paste0('train.', mode)] <-
      fn(train.resp, mean.tens.list[[paste0('train.', mode)]])

  # Get results for warm data
  for(mode in c('m1','m2','m3'))
    results$mean[name, paste0('warm.', mode)] <-
      fn(all.train.resp[is.na(train.data$resp)],
         mean.tens.list[[paste0('train.', mode)]][is.na(train.data$resp)])

  # Get results for modes
  for(mode in c('m1','m2','m3','m1m2','m1m3','m2m3','m1m2m3')) {
    if(exists(paste0(mode, '.pred.resp'))) {
      dat <- get(paste0(mode, '.resp'))
      if(length(dat))
        results$mean[name, mode] <- fn(dat, mean.tens.list[[mode]])
    }
  }
}

# Get results for BaTFLED predictions
for(m in 1:length(fns)) {
  fn <- fns[[m]]
  name <- names(fns)[m]

  results$trained[name, 'train'] <- fn(train.resp, train.pred.resp) 
  results$trained[name, 'train.H'] <- fn(train.resp, train.H.resp)
  results$trained[name, 'warm'] <- fn(train.resp[is.na(train.data$resp)],
                                      train.pred.resp[is.na(train.data$resp)])
  for(mode in c('m1', 'm2', 'm3', 'm1m2', 'm1m3', 'm2m3', 'm1m2m3')) {
    if(exists(paste0(mode, '.resp')))
      results$trained[name, mode] <- fn(get(paste0(mode, '.resp')), 
                                        get(paste0(mode, '.pred.resp')))
  }
}

print(results)

save(results, file=paste0(out.dir, 'results.Rdata'))

save.image(paste0(out.dir, 'image.Rdata'))

# Examine the fit.
####################################################
pdf(paste0(out.dir, 'training_plots.pdf'), height=7*2, width=7*2)
par(mfrow=c(2,2))
plot(trained$RMSE, main='Training RMSE')
plot(trained$lower.bnd, main=paste('Lower bound: max at', which.max(trained$lower.bnd)))
if(sum((trained$lower.bnd[-1]-trained$lower.bnd[-length(trained$lower.bnd)])<0)==0)
  print('Lower bounds are monotonically increasing')
plot(trained$times, main='Time for each iteration')
plot(trained$exp.var, main='Explained variance', ylim=c(0,1))
dev.off()

# RMSE 
pdf(paste0(out.dir, 'training_RMSEs.pdf'))
par(mfrow=c(1,1))
plot_test_RMSE(test.results, main="Test RMSEs", mean=T)
abline(h=warm.mean.RMSE, col='black', lty=2, lwd=2)
if(exists('m1.mean.RMSE'))
  abline(h=m1.mean.RMSE, col='red', lty=2, lwd=2)
if(exists('m2.mean.RMSE'))
  abline(h=m2.mean.RMSE, col='blue', lty=2, lwd=2)
if(exists('m3.mean.RMSE'))
  abline(h=m2.mean.RMSE, col='yellow', lty=2, lwd=2)
if(exists('m1m2.mean.RMSE'))
  abline(h=m1m2.mean.RMSE, col='purple', lty=3, lwd=2)
if(exists('m1m3.mean.RMSE'))
  abline(h=m1m3.mean.RMSE, col='orange', lty=3, lwd=2)
if(exists('m2m3.mean.RMSE'))
  abline(h=m2m3.mean.RMSE, col='green', lty=3, lwd=2)
if(exists('m1m2m3.mean.RMSE'))
  abline(h=m1m2m3.mean.RMSE, col='brown', lty=3, lwd=2)
dev.off()

save.image(paste0(out.dir, 'image.Rdata'))

## TODO: add in mode3 and combinations
pdf(paste0(out.dir, 'training_exp_var.pdf'))
plot_test_exp_var(test.results, mean=T,
     main=sprintf('Warm max at %.0f, cold mode 1 max at %d \n Cold mode 2 max at %d Both cold max at %d',
                  which.max(test.results$warm.exp.var), which.max(test.results$m1.exp.var),
                  which.max(test.results$m2.exp.var), which.max(test.results$m1m2.exp.var)))
abline(h=warm.mean.exp.var, col='red', lty=2, lwd=2)
if(exists('m1.mean.exp.var'))
  abline(h=m1.mean.exp.var, col='blue', lty=2, lwd=2)
if(exists('m2.mean.exp.var'))
  abline(h=m2.mean.exp.var, col='green', lty=2, lwd=2)
if(exists('m1m2.mean.exp.var'))
  abline(h=m1m2.mean.exp.var, col='purple', lty=2, lwd=2)
dev.off()

# Predictors
####################################################
# Plot the number of predictors at each iteration
# pdf(paste0(out.dir, 'num_preds.pdf'))
# par(mfrow=c(1,2))
# plot(sapply(m1.preds.list, length), pch=20)
# plot(sapply(m2.preds.list, length), pch=20)
# dev.off()

pdf(paste0(out.dir, 'mode1.matrices.pdf'), height=7, width=7)
par(mfrow=c(1,2))
if(ncol(trained$mode1.A.mean)) 
  im_mat(trained$mode1.A.mean, scale=T, 
    main=paste("Mode 1 A matrix\n", nrow(trained$mode1.A.mean), "predictors"))
im_mat(trained$mode1.H.mean, main="H matrix")
dev.off()

pdf(paste0(out.dir, 'mode2.matrices.pdf'), height=7, width=7)
par(mfrow=c(1,2))
if(ncol(trained$mode2.A.mean)) im_mat(trained$mode2.A.mean, scale=T,
                                      main=paste("Mode 2 A matrix\n", nrow(trained$mode2.A.mean), "predictors"))
im_mat(trained$mode2.H.mean, main="H matrix")
dev.off()

pdf(paste0(out.dir, 'mode3.matrices.pdf'), height=7, width=7)
par(mfrow=c(1,2))
if(nrow(trained$mode3.A.mean)) im_mat(trained$mode3.A.mean, scale=T,
                                      main=paste("Mode 3 A matrix\n", nrow(trained$mode3.A.mean), "predictors"))
im_mat(trained$mode3.H.mean, main="H matrix")
dev.off()

# Predictions
####################################################

# Plot warm start predictions ####################################################
pdf(paste0(out.dir, 'warm_preds.pdf'))
par(mfrow=c(1,1))
plot_preds(trained$resp, train.data$resp,
           main=sprintf('Warm Start: RMSE: preds: %.3f, mean:%.3f \n Explained variance: pred: %.3f mean: %.3f',
                        test.results$warm.RMSE[reps], warm.mean.RMSE, 
                        test.results$warm.exp.var[reps], warm.mean.exp.var))
points(warm.resp, warm.mean.preds, col='green')
points(warm.resp, warm.preds, col='red')
legend(x='topleft', legend=c('Preds', 'Mean'), text.col=c('red', 'green'), bty='n')
abline(a=0, b=1, lwd=2)
dev.off()

# Plot cold start predictions ####################################################
if(length(test.data.m1)) {
  pdf(paste0(out.dir, 'm1_cold_preds.pdf'))
  plot(test.data.m1$resp, m1.cold.preds, 
       main=sprintf('Mode 1 cold Start: \nRMSE: preds: %.3f, mean: %.3f\nExp. var.: preds: %.3f, mean %.3f',
                    test.results$m1.RMSE[trained$iter], m1.mean.RMSE, 
                    test.results$m1.exp.var[trained$iter], m1.mean.exp.var),
       col='red',  xlab='True responses', ylab='Predicted responses', pch=20,
       ylim=range(m1.cold.preds, m1.test.means, na.rm=T))
  points(test.data.m1$resp, m1.test.means, col='blue', pch=20)
  legend(x='topleft', legend=c('Preds', 'Mean'), text.col=c('red', 'blue'), bty='n')
  abline(a=0, b=1, lwd=2)
  dev.off()
}

if(length(test.data.m2)) {  
  pdf(paste0(out.dir, 'm2_cold_preds.pdf'))
  plot(test.data.m2$resp, m2.cold.preds, 
       main=sprintf('Mode 2 cold Start: \nRMSE: preds: %.3f, mean: %.3f\nExp. var.: preds: %.3f, mean %.3f',
                    test.results$m2.RMSE[trained$iter], m2.mean.RMSE, 
                    test.results$m2.exp.var[trained$iter], m2.mean.exp.var),
       col='red', xlab='True responses', ylab='Predicted responses', pch=20,
       ylim=range(m2.cold.preds, m2.test.means))
  points(test.data.m2$resp, m2.test.means, col='blue', pch=20)
  legend(x='topleft', legend=c('Preds', 'Mean'), text.col=c('red', 'blue'), bty='n')
  abline(a=0, b=1, lwd=2)
  dev.off()
}
 
if(length(test.data.m1m2)) { 
  pdf(paste0(out.dir, 'm1m2_cold_preds.pdf'))
  plot(test.data.m1m2$resp, m1m2.cold.preds, 
       main=sprintf('Mode 1&2 cold Start: \nRMSE: preds: %.3f, mean: %.3f\nExp. var.: preds: %.3f, mean %.3f',
                    test.results$m1m2.RMSE[trained$iter], m1m2.mean.RMSE, 
                    test.results$m1m2.exp.var[trained$iter], m1m2.mean.exp.var), 
       col='red', xlab='True responses', ylab='Predicted responses', pch=20,
       ylim=range(m1m2.cold.preds, m1m2.test.means))
  points(test.data.m1m2$resp, m1m2.test.means, col='blue', pch=20)
  legend(x='topleft', legend=c('Preds', 'Mean'), text.col=c('red', 'blue'), bty='n')
  abline(a=0, b=1, col='blue', lwd=2)
  dev.off()
}

# Which predictors are being used?
###################################################################################
# which(apply(trained$mode1.A.mean, 1, sd) > 1)
# sum(apply(trained$mode1.A.mean, 1, sd) > 1)
# which(apply(trained$mode2.A.mean, 1, sd) > 1)
# sum(apply(trained$mode2.A.mean, 1, sd) > 1)

# Look at predictions during training
##################################################################################
# Which iteration gives the smallest cold start RMSE?
print(sprintf('Warm RMSE min at iteration: %d of %.2f', 
              which.min(test.results$warm.RMSE), min(test.results$warm.RMSE, na.rm=T)))
print(sprintf('Cold mode 1 RMSE min at iteration: %d of %.2f', 
              which.min(test.results$m1.RMSE), min(test.results$m1.RMSE, na.rm=T)))
print(sprintf('Cold mode 2 RMSE min at iteration: %d of %.2f', 
              which.min(test.results$m2.RMSE), min(test.results$m2.RMSE, na.rm=T)))
print(sprintf('Cold mode 1/mode 2 combination RMSE min at iteration: %d of %.2f', 
              which.min(test.results$m1m2.RMSE), min(test.results$m1m2.RMSE, na.rm=T)))

# Which iteration gives the Largest explained variances?
print(sprintf('Warm explained variance max at iteration: %d of %.2f',
              which.max(test.results$warm.exp.var), max(test.results$warm.exp.var, na.rm=T)))
print(sprintf('Cold mode 1 explained variance max at iteration: %d of %.2f',
              which.max(test.results$m1.exp.var), max(test.results$m1.exp.var, na.rm=T)))
print(sprintf('Cold mode 2 explained variance max at iteration: %d of %.2f',
              which.max(test.results$m2.exp.var), max(test.results$m2.exp.var, na.rm=T)))
print(sprintf('Cold mode 1 & 2 combination explained variance max at iteration: %d of %.2f',
              which.max(test.results$m1m2.exp.var), max(test.results$m1m2.exp.var, na.rm=T)))

# Plot some response curves for left out mode 1
if(length(test.data.m1)) {
  pdf(paste0(out.dir, 'cl_resp_curves.pdf'), height=6, width=6*1.62)
  par(mfrow=c(2,3))
  for(n in 1:6) {
    i <- sample(1:dim(test.data.m1$resp)[[1]], 1)
    j <- sample(1:dim(test.data.m1$resp)[[2]], 1)
    # Unnormalize data
    if(multiplier != 0) {
      obs <- test.data.m1$resp[i,j,] / multiplier
      pred <- m1.cold.preds[i,j,] / multiplier
    } else {
      obs <- test.data.m1$resp[i,j,]
      pred <- m1.cold.preds[i,j,]
    }
    ylim=range(obs, pred, na.rm=T)
    plot(obs, pch=20, col='blue', ylim=ylim,
         main=sprintf("Mode 1: %s Drug: %s",
                      dimnames(test.data.m1$resp)[[1]][i],
                      dimnames(test.data.m1$resp)[[2]][j]))
    points(pred, pch=20, col='red')
  }
  dev.off()
}

if(length(test.data.m2)) { 
  # Plot some response curves for left out mode 2
  pdf(paste0(out.dir, 'dr_resp_curves.pdf'), height=6, width=6*1.62)
  par(mfrow=c(2,3))
  for(n in 1:6) {
    i <- sample(1:dim(test.data.m2$resp)[[1]], 1)
    j <- sample(1:dim(test.data.m2$resp)[[2]], 1)
    # Unnormalize data
    if(multiplier != 0) {
      obs <- test.data.m2$resp[i,j,] / multiplier
      pred <- m2.cold.preds[i,j,] / multiplier
    } else {
      obs <- test.data.m2$resp[i,j,]
      pred <- m2.cold.preds[i,j,]
    }
    ylim=range(obs, pred, na.rm=T)
    plot(obs, pch=20, col='blue', ylim=ylim,
         main=sprintf("Mode 1: %s Drug: %s",
                      dimnames(test.data.m2$resp)[[1]][i],
                      dimnames(test.data.m2$resp)[[2]][j]))
    points(pred, pch=20, col='red')
  }
  dev.off()
}

if(length(test.data.m1m2)) {
  # Plot some response curves for left out mode 1 and mode 2
  pdf(paste0(out.dir, 'm1m2_resp_curves.pdf'), height=6, width=6*1.62)
  par(mfrow=c(2,3))
  for(n in 1:6) {
    i <- sample(1:dim(test.data.m1m2$resp)[[1]], 1)
    j <- sample(1:dim(test.data.m1m2$resp)[[2]], 1)
    # Unnormalize data
    if(multiplier != 0) {
      obs <- test.data.m1m2$resp[i,j,] / multiplier
      pred <- m1m2.cold.preds[i,j,] / multiplier
    } else {
      obs <- test.data.m1m2$resp[i,j,]
      pred <- m1m2.cold.preds[i,j,]
    }
    ylim=range(obs, pred, na.rm=T)
    plot(obs, pch=20, col='blue', ylim=ylim,
         main=sprintf("Mode 1: %s Drug: %s",
                      dimnames(test.data.m1m2$resp)[[1]][i],
                      dimnames(test.data.m1m2$resp)[[2]][j]))
    points(pred, pch=20, col='red')
  }
  dev.off()
}

# Save everything again
save.image(paste0(out.dir, 'image.Rdata'))

#############################################################################################################