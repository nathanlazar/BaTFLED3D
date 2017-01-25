#!/usr/bin/env Rscript

# nathan dot lazar at gmail dot com

# Usage: BaTFLED_run_post.R <run_num/image.Rdata>

print('Opened file')

library(methods)
library(BaTFLED3D)      # Note: package must be installed from .tar file

args <- commandArgs(TRUE)
load(args[[1]])

# Make a tensor of the means for each value of mode2 & mode3 for comparisons
train.means <- train.data$resp
for(i in 1:dim(train.data$resp)[1]) train.means[i,,] <- m1.means
if(length(test.data.m1)) {
  m1.test.means <- test.data.m1$resp
  for(i in 1:dim(test.data.m1$resp)[1]) m1.test.means[i,,] <- m1.means
}
if(length(test.data.m2)) {
  m2.test.means <- test.data.m2$resp
  for(j in 1:dim(test.data.m2$resp)[2]) m2.test.means[,j,] <- m2.means
}
if(length(test.data.m3)) {
  m3.test.means <- test.data.m3$resp
  for(k in 1:dim(test.data.m3$resp)[3]) m3.test.means[,,k] <- m3.means
}
if(length(test.data.m1m2)) {
  m1m2.test.means <- array(0, dim=dim(test.data.m1m2$resp))
  m1m2.test.means <- apply(train.data$resp, 3, mean, na.rm=T)
}
if(length(test.data.m1m3)) {
  m1m3.test.means <- array(0, dim=dim(test.data.m1m3$resp))
  m1m3.test.means <- apply(train.data$resp, 2, mean, na.rm=T)
}
if(length(test.data.m2m3)) {
  m2m3.test.means <- array(0, dim=dim(test.data.m2m3$resp))
  m2m3.test.means <- apply(train.data$resp, 1, mean, na.rm=T)
}
if(length(test.data.m1m2m3)) {
  m1m2m3.test.means <- array(0, dim=dim(test.data.m1m2m3$resp))
  m1m2m3.test.means[,,] <- mean(train.data$resp, na.rm=T)
}

# Calculate predictions
warm.preds <- trained$resp[is.na(train.data$resp)]
warm.resp <- all.resp.train[is.na(train.data$resp)]
warm.mean.preds <- train.means[is.na(train.data$resp)]
warm.mean.exp.var <- exp_var(warm.resp, warm.mean.preds)
warm.mean.RMSE <- sqrt(mean((warm.resp-warm.mean.preds)^2, na.rm=T))/sd(train.data$resp, na.rm=T)

if(length(test.data.m1)) {
  m1.cold.preds <- test(d=test.data.m1, m=trained)
  m1.mean.RMSE <- sqrt(mean((test.data.m1$resp-m1.test.means)^2, na.rm=T))/sd(train.data$resp, na.rm=T)
  m1.mean.exp.var <- exp_var(test.data.m1$resp, m1.test.means)
}
if(length(test.data.m2)) {
  m2.cold.preds <- test(d=test.data.m2, m=trained)
  m2.mean.RMSE <- sqrt(mean((test.data.m2$resp-m2.test.means)^2, na.rm=T))/sd(train.data$resp, na.rm=T)
  m2.mean.exp.var <- exp_var(test.data.m2$resp, m2.test.means)
}
if(length(test.data.m3)) {
  m3.cold.preds <- test(d=test.data.m3, m=trained)
  m3.mean.RMSE <- sqrt(mean((test.data.m3$resp-m3.test.means)^2, na.rm=T))/sd(train.data$resp, na.rm=T)
  m3.mean.exp.var <- exp_var(test.data.m3$resp, m3.test.means)
}
if(length(test.data.m1m2)) {
  m1m2.cold.preds <- test(d=test.data.m1m2, m=trained)
  m1m2.mean.RMSE <- sqrt(mean((test.data.m1m2$resp-m1m2.test.means)^2, na.rm=T))/sd(train.data$resp, na.rm=T)
  m1m2.mean.exp.var <- exp_var(test.data.m1m2$resp, m1m2.test.means)
}
if(length(test.data.m1m3)) {
  m1m3.cold.preds <- test(d=test.data.m1m3, m=trained)
  m1m3.mean.RMSE <- sqrt(mean((test.data.m1m3$resp-m1m3.test.means)^2, na.rm=T))/sd(train.data$resp, na.rm=T)
  m1m3.mean.exp.var <- exp_var(test.data.m1m3$resp, m1m3.test.means)
}
if(length(test.data.m2m3)) {
  m2m3.cold.preds <- test(d=test.data.m2m3, m=trained)
  m2m3.mean.RMSE <- sqrt(mean((test.data.m2m3$resp-m2m3.test.means)^2, na.rm=T))/sd(train.data$resp, na.rm=T)
  m2m3.mean.exp.var <- exp_var(test.data.m2m3$resp, m2m3.test.means)
}
if(length(test.data.m1m2m3)) {
  m1m2m3.cold.preds <- test(d=test.data.m1m2m3, m=trained)
  m1m2m3.mean.RMSE <- sqrt(mean((test.data.m1m2m3$resp-m1m2m3.test.means)^2, na.rm=T))/sd(train.data$resp, na.rm=T)
  m1m2m3.mean.exp.var <- exp_var(test.data.m1m2m3$resp, m1m2m3.test.means)
}

warm.mean.exp.var <- exp_var(warm.resp, warm.mean.preds)

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
if(ncol(trained$mode1.A.mean)) im_mat(trained$mode1.A.mean, scale=T,
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

# Save everything again
save.image(paste0(out.dir, 'image.Rdata'))

##############################################################################################################

# Examine the fit.
####################################################
# Fit logistic curves to the predicted responses and get GI50 values
########################################################################################
# ndoses <- dim(train.data$resp)[[3]]
# 
# # Create empty matrices to store computed gi50s
# train.obs.gi50s <- matrix(NA, dim(train.data$resp)[[1]], dim(train.data$resp)[[2]],
#                           dimnames=dimnames(train.data$resp)[1:2])
# train.mean.gi50s <- train.obs.gi50s
# train.pred.gi50s <- train.obs.gi50s
# 
# test.m1.obs.gi50s <- matrix(NA, dim(test.data.m1$resp)[[1]], dim(test.data.m1$resp)[[2]],
#                             dimnames=dimnames(test.data.m1$resp)[1:2])
# test.m1.pred.gi50s <- test.m1.obs.gi50s
# test.m1.mean.gi50s <- test.m1.obs.gi50s
# 
# test.m2.obs.gi50s <- matrix(NA, dim(test.data.m2$resp)[[1]], dim(test.data.m2$resp)[[2]],
#                             dimnames=dimnames(test.data.m2$resp)[1:2])
# test.m2.pred.gi50s <- test.m2.obs.gi50s
# test.m2.mean.gi50s <- test.m2.obs.gi50s
# 
# test.m1m2.obs.gi50s <- matrix(NA, dim(test.data.m1m2$resp)[[1]], dim(test.data.m1m2$resp)[[2]],
#                               dimnames=dimnames(test.data.m1m2$resp)[1:2])
# test.m1m2.pred.gi50s <- test.m1m2.obs.gi50s
# test.m1m2.mean.gi50s <- test.m1m2.obs.gi50s
# 
# # Get observed and predicted GI50s for training data
# for(i in 1:nrow(train.obs.gi50s))
#   for(j in 1:ncol(train.obs.gi50s)) {
# 
#     # Just using 1...16 as doses now.
#     # TODO: Get actual doses
#     obs <- data.frame(dose=1:ndoses, resp=train.data$resp[i,j,])
#     pred <- data.frame(dose=1:ndoses, resp=trained$resp[i,j,])
#     me <- data.frame(dose=1:ndoses, resp=train.means[i,j,])
# 
#     # Fit a 4 parameter logistic curves
#     obs4.fit <- try(drm(resp~dose, data=obs, fct = LL.4(), type="continuous"), T)
#     pred4.fit <- try(drm(resp~dose, data=pred, fct = LL.4(), type="continuous"), T)
#     mean4.fit <- try(drm(resp~dose, data=me, fct = LL.4(), type="continuous"), T)
# 
#     # Add gi50s to matrices and plot curves if plot=T
#     if(plot) {
#       cl <- dimnames(train.obs.gi50s)[[1]][i]
#       dr <- dimnames(train.obs.gi50s)[[2]][j]
#       par(mfrow=c(1,3))
#     }
#     if(length(obs4.fit) > 1) if(obs4.fit$fit$convergence) {
#         train.obs.gi50s[i,j] <- ifelse(obs4.fit$coefficients[4] < max(obs$dose),
#                                  obs4.fit$coefficients[4], max(obs$dose))
#         if(plot) {
#           plot(obs4.fit, pch=20, col='blue',
#                       main=sprintf("Observed C.L.: %s \n Drug: %s Fit: %.2f",
#                                              cl, dr, obs4.fit$fit$value))
#           abline(v=train.obs.gi50s[i,j], lty=2, lwd=2)
#         }
#     }
#     if(length(pred4.fit) > 1) if(pred4.fit$fit$convergence) {
#       train.pred.gi50s[i,j] <- ifelse(pred4.fit$coefficients[4] < max(obs$dose),
#                                   pred4.fit$coefficients[4], max(obs$dose))
#       if(plot) {
#         plot(pred4.fit, pch=20, col='red',
#              main=sprintf("Predicted C.L.: %s \n Drug: %s Fit: %.2f",
#                                    cl, dr, pred4.fit$fit$value))
#         abline(v=train.pred.gi50s[i,j], lty=2, lwd=2)
#       }
#     }
#     if(length(mean4.fit) > 1) if(mean4.fit$fit$convergence) {
#       train.mean.gi50s[i,j] <- ifelse(mean4.fit$coefficients[4] < max(obs$dose),
#                                       mean4.fit$coefficients[4], max(obs$dose))
#       if(plot) {
#         plot(mean4.fit, pch=20,
#              main=sprintf("Mean C.L.: %s \n Drug: %s Fit: %.2f",
#                           cl, dr, mean4.fit$fit$value))
#         abline(v=train.mean.gi50s[i,j], lty=2, lwd=2)
#       }
#     }
# }
# 
# # Get observed and predicted GI50s for mode 1 test data
# for(i in 1:nrow(test.m1.obs.gi50s))
#   for(j in 1:ncol(test.m1.obs.gi50s)) {
# 
#     # Just using 1...16 as doses now.
#     # TODO: Get actual doses
#     obs <- data.frame(dose=1:ndoses, resp=test.data.m1$resp[i,j,])
#     pred <- data.frame(dose=1:ndoses, resp=m1.cold.preds[i,j,])
#     me <- data.frame(dose=1:ndoses, resp=m1.test.means[i,j,])
# 
#     # Fit a 4 parameter logistic curves
#     obs4.fit <- try(drm(resp~dose, data=obs, fct = LL.4(), type="continuous"), T)
#     pred4.fit <- try(drm(resp~dose, data=pred, fct = LL.4(), type="continuous"), T)
#     mean4.fit <- try(drm(resp~dose, data=me, fct = LL.4(), type="continuous"), T)
# 
#     # Add gi50s to matrices and plot curves if plot=T
#     if(plot) {
#       par(mfrow=c(1,3))
#       cl <- dimnames(test.m1.obs.gi50s)[[1]][i]
#       dr <- dimnames(test.m1.obs.gi50s)[[2]][j]
#     }
#     if(length(obs4.fit) > 1) if(obs4.fit$fit$convergence) {
#       test.m1.obs.gi50s[i,j] <- ifelse(obs4.fit$coefficients[4] < max(obs$dose),
#                                      obs4.fit$coefficients[4], max(obs$dose))
#       if(plot) {
#         plot(obs4.fit, pch=20, col='blue',
#              main=sprintf("Observed C.L.: %s \n Drug: %s Fit: %.2f",
#                           cl, dr, obs4.fit$fit$value))
#         abline(v=test.m1.obs.gi50s[i,j], lty=2, lwd=2)
#       }
#     }
#     if(length(pred4.fit) > 1) if(pred4.fit$fit$convergence) {
#       test.m1.pred.gi50s[i,j] <- ifelse(pred4.fit$coefficients[4] < max(obs$dose),
#                                       pred4.fit$coefficients[4], max(obs$dose))
#       if(plot) {
#         plot(pred4.fit, pch=20, col='red',
#              main=sprintf("Predicted C.L.: %s \n Drug: %s Fit: %.2f",
#                           cl, dr, pred4.fit$fit$value))
#         abline(v=test.m1.pred.gi50s[i,j], lty=2, lwd=2)
#       }
#     }
#     if(length(mean4.fit) > 1) if(mean4.fit$fit$convergence) {
#       test.m1.mean.gi50s[i,j] <- ifelse(mean4.fit$coefficients[4] < max(obs$dose),
#                                       mean4.fit$coefficients[4], max(obs$dose))
#       if(plot) {
#         plot(mean4.fit, pch=20,
#              main=sprintf("Mean C.L.: %s \n Drug: %s Fit: %.2f",
#                           cl, dr, mean4.fit$fit$value))
#         abline(v=test.m1.mean.gi50s[i,j], lty=2, lwd=2)
#       }
#     }
# }
# 
# # Get observed and predicted GI50s for mode 2 test data
# for(i in 1:nrow(test.m2.obs.gi50s))
#   for(j in 1:ncol(test.m2.obs.gi50s)) {
# 
#     # Just using 1...16 as doses now.
#     # TODO: Get actual doses
#     obs <- data.frame(dose=1:ndoses, resp=test.data.m2$resp[i,j,])
#     pred <- data.frame(dose=1:ndoses, resp=m2.cold.preds[i,j,])
#     me <- data.frame(dose=1:ndoses, resp=m2.test.means[i,j,])
# 
#     # Fit a 4 parameter logistic curves
#     obs4.fit <- try(drm(resp~dose, data=obs, fct = LL.4(), type="continuous"), T)
#     pred4.fit <- try(drm(resp~dose, data=pred, fct = LL.4(), type="continuous"), T)
#     mean4.fit <- try(drm(resp~dose, data=me, fct = LL.4(), type="continuous"), T)
# 
#     # Add gi50s to matrices and plot curves if plot=T
#     if(plot) {
#       par(mfrow=c(1,3))
#       cl <- dimnames(test.m2.obs.gi50s)[[1]][i]
#       dr <- dimnames(test.m2.obs.gi50s)[[2]][j]
#     }
#     if(length(obs4.fit) > 1) if(obs4.fit$fit$convergence) {
#       test.m2.obs.gi50s[i,j] <- ifelse(obs4.fit$coefficients[4] < max(obs$dose),
#                                        obs4.fit$coefficients[4], max(obs$dose))
#       if(plot) {
#         plot(obs4.fit, pch=20, col='blue',
#              main=sprintf("Observed C.L.: %s \n Drug: %s Fit: %.2f",
#                           cl, dr, obs4.fit$fit$value))
#         abline(v=test.m2.obs.gi50s[i,j], lty=2, lwd=2)
#       }
#     }
#     if(length(pred4.fit) > 1) if(pred4.fit$fit$convergence) {
#       test.m2.pred.gi50s[i,j] <- ifelse(pred4.fit$coefficients[4] < max(obs$dose),
#                                         pred4.fit$coefficients[4], max(obs$dose))
#       if(plot) {
#         plot(pred4.fit, pch=20, col='red',
#              main=sprintf("Predicted C.L.: %s \n Drug: %s Fit: %.2f",
#                           cl, dr, pred4.fit$fit$value))
#         abline(v=test.m2.pred.gi50s[i,j], lty=2, lwd=2)
#       }
#     }
#     if(length(mean4.fit) > 1) if(mean4.fit$fit$convergence) {
#       test.m2.mean.gi50s[i,j] <- ifelse(mean4.fit$coefficients[4] < max(obs$dose),
#                                         mean4.fit$coefficients[4], max(obs$dose))
#       if(plot) {
#         plot(mean4.fit, pch=20,
#              main=sprintf("Mean C.L.: %s \n Drug: %s Fit: %.2f",
#                           cl, dr, mean4.fit$fit$value))
#         abline(v=test.m2.mean.gi50s[i,j], lty=2, lwd=2)
#       }
#     }
# }
# 
# # Get observed and predicted GI50s for simultaneous mode 1 and mode 2 test data
# for(i in 1:nrow(test.m1m2.obs.gi50s))
#   for(j in 1:ncol(test.m1m2.obs.gi50s)) {
# 
#     # Just using 1...16 as doses now.
#     # TODO: Get actual doses
#     obs <- data.frame(dose=1:ndoses, resp=test.data.m1m2$resp[i,j,])
#     pred <- data.frame(dose=1:ndoses, resp=m1m2.cold.preds[i,j,])
#     me <- data.frame(dose=1:ndoses, resp=m1m2.test.means[i,j,])
# 
#     # Fit a 4 parameter logistic curves
#     obs4.fit <- try(drm(resp~dose, data=obs, fct = LL.4(), type="continuous"), T)
#     pred4.fit <- try(drm(resp~dose, data=pred, fct = LL.4(), type="continuous"), T)
#     mean4.fit <- try(drm(resp~dose, data=me, fct = LL.4(), type="continuous"), T)
# 
#     # Add gi50s to matrices and plot curves if plot=T
#     if(plot) {
#       par(mfrow=c(1,3))
#       cl <- dimnames(test.m1m2.obs.gi50s)[[1]][i]
#       dr <- dimnames(test.m1m2.obs.gi50s)[[2]][j]
#     }
#     if(length(obs4.fit) > 1) if(obs4.fit$fit$convergence) {
#       test.m1m2.obs.gi50s[i,j] <- ifelse(obs4.fit$coefficients[4] < max(obs$dose),
#                                          obs4.fit$coefficients[4], max(obs$dose))
#       if(plot) {
#         plot(obs4.fit, pch=20, col='blue',
#              main=sprintf("Observed C.L.: %s \n Drug: %s Fit: %.2f",
#                           cl, dr, obs4.fit$fit$value))
#         abline(v=test.m1m2.obs.gi50s[i,j], lty=2, lwd=2)
#       }
#     }
#     if(length(pred4.fit) > 1) if(pred4.fit$fit$convergence) {
#       test.m1m2.pred.gi50s[i,j] <- ifelse(pred4.fit$coefficients[4] < max(obs$dose),
#                                           pred4.fit$coefficients[4], max(obs$dose))
#       if(plot) {
#         plot(pred4.fit, pch=20, col='red',
#              main=sprintf("Predicted C.L.: %s \n Drug: %s Fit: %.2f",
#                           cl, dr, pred4.fit$fit$value))
#         abline(v=test.m1m2.pred.gi50s[i,j], lty=2, lwd=2)
#       }
#     }
#     if(length(mean4.fit) > 1) if(mean4.fit$fit$convergence) {
#       test.m1m2.mean.gi50s[i,j] <- ifelse(mean4.fit$coefficients[4] < max(obs$dose),
#                                           mean4.fit$coefficients[4], max(obs$dose))
#       if(plot) {
#         plot(mean4.fit, pch=20,
#              main=sprintf("Mean C.L.: %s \n Drug: %s Fit: %.2f",
#                           cl, dr, mean4.fit$fit$value))
#         abline(v=test.m1m2.mean.gi50s[i,j], lty=2, lwd=2)
#       }
#     }
# }
# 
# # Get RMSE and correlation between observed & predicted responses
# train.gi50.rmse <- sqrt(mean((train.obs.gi50s - train.pred.gi50s)^2, na.rm=T))
# train.gi50.cor <- cor(as.vector(train.obs.gi50s[!is.na(train.obs.gi50s)]),
#                       as.vector(train.pred.gi50s[!is.na(train.obs.gi50s)]))
# test.m1.gi50.rmse <- sqrt(mean((test.m1.obs.gi50s - test.m1.pred.gi50s)^2, na.rm=T))
# test.m1.gi50.cor <- cor(as.vector(test.m1.obs.gi50s[!is.na(test.m1.obs.gi50s)]),
#                       as.vector(test.m1.pred.gi50s[!is.na(test.m1.obs.gi50s)]))
# test.m2.gi50.rmse <- sqrt(mean((test.m2.obs.gi50s - test.m2.pred.gi50s)^2, na.rm=T))
# test.m2.gi50.cor <- cor(as.vector(test.m2.obs.gi50s[!is.na(test.m2.obs.gi50s)]),
#                         as.vector(test.m2.pred.gi50s[!is.na(test.m2.obs.gi50s)]))
# test.m1m2.gi50.rmse <- sqrt(mean((test.m1m2.obs.gi50s - test.m1m2.pred.gi50s)^2, na.rm=T))
# test.m1m2.gi50.cor <- cor(as.vector(test.m1m2.obs.gi50s[!is.na(test.m1m2.obs.gi50s)]),
#                         as.vector(test.m1m2.pred.gi50s[!is.na(test.m1m2.obs.gi50s)]))
# 
# # Compare to predicting the mean response at each mode 3
# train.mean.gi50.rmse <- sqrt(mean((train.obs.gi50s - train.mean.gi50s)^2, na.rm=T))
# train.mean.gi50.cor <- cor(as.vector(train.obs.gi50s[!is.na(train.obs.gi50s)]),
#                       as.vector(train.mean.gi50s[!is.na(train.obs.gi50s)]))
# test.m1.mean.gi50.rmse <- sqrt(mean((test.m1.obs.gi50s - test.m1.mean.gi50s)^2, na.rm=T))
# test.m1.mean.gi50.cor <- cor(as.vector(test.m1.obs.gi50s[!is.na(test.m1.obs.gi50s)]),
#                       as.vector(test.m1.mean.gi50s[!is.na(test.m1.obs.gi50s)]))
# test.m2.mean.gi50.rmse <- sqrt(mean((test.m2.obs.gi50s - test.m2.mean.gi50s)^2, na.rm=T))
# test.m2.mean.gi50.cor <- cor(as.vector(test.m2.obs.gi50s[!is.na(test.m2.obs.gi50s)]),
#                         as.vector(test.m2.mean.gi50s[!is.na(test.m2.obs.gi50s)]))
# test.m1m2.mean.gi50.rmse <- sqrt(mean((test.m1m2.obs.gi50s - test.m1m2.mean.gi50s)^2, na.rm=T))
# test.m1m2.mean.gi50.cor <- cor(as.vector(test.m1m2.obs.gi50s[!is.na(test.m1m2.obs.gi50s)]),
#                         as.vector(test.m1m2.mean.gi50s[!is.na(test.m1m2.obs.gi50s)]))
# 
# rmse.results <- data.frame(predicted_RMSE=c(train.gi50.rmse, test.m1.gi50.rmse,
#                                        test.m2.gi50.rmse, test.m1m2.gi50.rmse),
#                            mean_RMSE=c(train.mean.gi50.rmse, test.m1.mean.gi50.rmse,
#                                        test.m2.mean.gi50.rmse, test.m1m2.mean.gi50.rmse),
#                            predicted_cor=c(train.gi50.cor, test.m1.gi50.cor,
#                                            test.m2.gi50.cor, test.m1m2.gi50.cor),
#                            mean_cor=c(train.mean.gi50.cor, test.m1.mean.gi50.cor,
#                                       test.m2.mean.gi50.cor, test.m1m2.mean.gi50.cor))
# rownames(rmse.results) <- c('train', 'test.m1', 'test.m2', 'test.m1m2')
# print(rmse.results)
# 
# ## TODO: Add warm test results? ##
# 
# pdf('train_GI50_predictions.pdf')
# plot(train.obs.gi50s, train.pred.gi50s, pch=20, xlab="Observed", ylab="Predicted",
#      main=sprintf("Training GI50s \nPredictions RMSE : %.2f corr: %.3f \n RMSE for mean: %.2f, corr: %.3f",
#                   train.gi50.rmse, train.gi50.cor, train.mean.gi50.rmse, train.mean.gi50.cor))
# abline(a=0, b=1, col='blue', lwd=2)
# points(train.obs.gi50s, train.mean.gi50s, pch=20, xlab="Observed", ylab="Mean",
#        col='green')
# legend(x='bottomright', legend = c('predictions', 'mean'),
#        col = c('black', 'green'),  bty='n', pch=20)
# dev.off()
# 
# pdf('test_cl_GI50_predictions.pdf')
# plot(test.m1.obs.gi50s, test.m1.pred.gi50s, pch=20, xlab="Observed", ylab="Predicted",
#      main=sprintf("Mode 1 test GI50s \nPredictions RMSE : %.2f corr: %.3f \n RMSE for mean: %.2f, corr: %.3f",
#                   test.m1.gi50.rmse, test.m1.gi50.cor, test.m1.mean.gi50.rmse, test.m1.mean.gi50.cor))
# abline(a=0, b=1, col='blue', lwd=2)
# points(test.m1.obs.gi50s, test.m1.mean.gi50s, pch=20, xlab="Observed", ylab="Mean",
#        col='green')
# legend(x='bottomright', legend = c('predictions', 'mean'),
#        col = c('black', 'green'),  bty='n', pch=20)
# dev.off()
# 
# pdf('test_dr_GI50_predictions.pdf')
# plot(test.m2.obs.gi50s, test.m2.pred.gi50s, pch=20, xlab="Observed", ylab="Predicted",
#      main=sprintf("Drug test GI50s \nPredictions RMSE : %.2f corr: %.3f \n RMSE for mean: %.2f, corr: %.3f",
#                   test.m2.gi50.rmse, test.m2.gi50.cor, test.m2.mean.gi50.rmse, test.m2.mean.gi50.cor))
# abline(a=0, b=1, col='blue', lwd=2)
# points(test.m2.obs.gi50s, test.m2.mean.gi50s, pch=20, xlab="Observed", ylab="Mean",
#        col='green')
# legend(x='bottomright', legend = c('predictions', 'mean'),
#        col = c('black', 'green'),  bty='n', pch=20)
# dev.off()
# 
# pdf('test_m1m2_GI50_predictions.pdf')
# plot(test.m1m2.obs.gi50s, test.m1m2.pred.gi50s, pch=20, xlab="Observed", ylab="Predicted",
#      main=sprintf("Mode 1 and mode 2 test GI50s \nPredictions RMSE : %.2f corr: %.3f \n RMSE for mean: %.2f, corr: %.3f",
#                   test.m1m2.gi50.rmse, test.m1m2.gi50.cor, test.m1m2.mean.gi50.rmse, test.m1m2.mean.gi50.cor))
# abline(a=0, b=1, col='blue', lwd=2)
# points(test.m1m2.obs.gi50s, test.m1m2.mean.gi50s, pch=20, xlab="Observed", ylab="Mean",
#        col='green')
# legend(x='bottomright', legend = c('predictions', 'mean'),
#        col = c('black', 'green'),  bty='n', pch=20)
# dev.off()
# 
# pdf('train_GI50_matrices.pdf')
# par(mfrow=c(2,2))
# im_mat(train.obs.gi50s, main='Training observed')
# im_mat(train.pred.gi50s, main='Training predicted')
# im_mat(train.mean.gi50s, main='Predicting mean response')
# dev.off()
# 
# pdf('test_cl_GI50_matrices.pdf')
# par(mfrow=c(2,2))
# im_mat(test.m1.obs.gi50s, main='Mode 1 Test observed')
# im_mat(test.m1.pred.gi50s, main='Mode 1 Test predicted')
# im_mat(test.m1.mean.gi50s, main='Predicting mean response')
# dev.off()
# 
# pdf('test_dr_GI50_matrices.pdf')
# par(mfrow=c(2,2))
# im_mat(test.m2.obs.gi50s, main='Drug test observed')
# im_mat(test.m2.pred.gi50s, main='Drug test predicted')
# im_mat(test.m2.mean.gi50s, main='Predict mean response')
# dev.off()
# 
# pdf('test_m1m2_GI50_matrices.pdf')
# par(mfrow=c(2,2))
# im_mat(test.m1m2.obs.gi50s, main='C.L. & mode 2 test observed')
# im_mat(test.m1m2.pred.gi50s, main='C.L. & mode 2 test predicted')
# im_mat(test.m1m2.mean.gi50s, main='Predicting mean response')
# dev.off()
# 
# # Save everything again
# ###################################################################################
# save.image(paste0(out.dir, 'image.Rdata'))
