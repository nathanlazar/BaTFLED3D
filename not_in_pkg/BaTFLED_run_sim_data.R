#!/usr/bin/env Rscript

# nathan dot lazar at gmail dot com

# Usage: BaTFLED_run_sim_data.R <output dir> <options...>

# Example: BaTFLED_run_sim_data.R test decomp=Tucker \
#  row.share=T reps=30 warm.per=0.01 plot=F \
#  m1.rows=10 m1.cols=100 \
#  m2.rows=15 m2.cols=150 \
#  m3.rows=8 m3.cols=0 \
#  R1=8 R2=8 R3=8 \
#  m1.true=10 m2.true=15 m3.true=0 \
#  A1.intercept=T A2.intercept=T A3.intercept=F \
#  H1.intercept=T H2.intercept=T H3.intercept=T \
#  m1.test=2 m2.test=2 m3.test=0 \
#  m1.alpha=1e-10 m2.alpha=1e-10 m3.alpha=1 \
#  m1.beta=1e10 m2.beta=1e10 m3.beta=1 \
#  core.0D.alpha=1 core.0D.beta=1 \
#  core.1D.alpha=1e-4 core.1D.beta=1e4 \
#  core.2D.alpha=1e-8 core.2D.beta=1e8 \
#  core.3D.alpha=1e-12 core.3D.beta=1e12 \
#  m1.sigma2=0.01 m2.sigma2=0.01 m3.sigma2=1 \
#  sigma2=1 parallel=T cores=16 \
#  core.spar=.5 noise.sd=0.1

pkgTest <- function(x) { # load packages and install if needed
  if (!require(x,character.only = TRUE)) {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

pkgTest('foreach')      # For parallel processing in loops
pkgTest('R6')           # For making memory efficent R6 objects
pkgTest('iterators')    # Making iterator objects for efficient parallelization
pkgTest('rTensor')      # Needed for multiplying matrices into tensors (could be removed)
pkgTest('dplyr')        # General data frame manipulation
pkgTest('ggplot2')      # Pretty plotting
# Packages for registering parallel backend (depends on platform)
ifelse(.Platform$OS.type == "windows", library(doParallel), library(doMC))
# install.packages('BaTFLED3D_0.0.1.tar.gz', repos = NULL, type="source")  # Note: must be run from the directory containing the tar.gz file
library(BaTFLED3D)

## Main
#################################################
print('Packages loaded')
args <- commandArgs(TRUE)

# Report arguments and summary stats of the input objects
print(unlist(args))

# Read in options
out.dir <- paste0(args[[1]], '/')

# If the output directory doesn't exist, create it
print(paste('Output directory:', out.dir))
if (!file.exists(out.dir)){
  dir.create(file.path(out.dir))
}

# Set some defaults
reps <- 10; warm.per <- 0.05; cores <- 1
m1.test <- 0; m2.test <- 0; m3.test <- 0

# Override defaults if provided
if(sum(grepl('^reps=', args)))
  reps <- as.numeric(sub('reps=', '', args[grepl('^reps=', args)]))
if(sum(grepl('^warm.per=', args)))
  warm.per <- as.numeric(sub('warm.per=', '', args[grepl('^warm.per=', args)]))
if(sum(grepl('^cores=', args)))
  cores <- as.numeric(sub('cores=', '', args[grepl('^cores=', args)]))
if(sum(grepl('^m1.test=', args)))
  m1.test <- as.numeric(sub('m1.test=', '', args[grepl('^m1.test=', args)]))
if(sum(grepl('^m2.test=', args)))
  m2.test <- as.numeric(sub('m2.test=', '', args[grepl('^m2.test=', args)]))
if(sum(grepl('^m3.test=', args)))
  m3.test <- as.numeric(sub('m3.test=', '', args[grepl('^m3.test=', args)]))

# Set up backend for parallel execution
cores <- cores
if(.Platform$OS.type == "windows") {
  clust <- makeCluster(cores)
  registerDoParallel(clust)
} else registerDoMC(cores)

## Generating a simulated dataset
data.params <- get_data_params(args)
toy <- mk_toy(data.params)

## Make training data object
m1.rem <- sample(row.names(toy$mode1.X), m1.test)
m2.rem <- sample(row.names(toy$mode2.X), m2.test)
m3.rem <- sample(row.names(toy$mode3.X), m3.test)
train.data <- input_data$new(
  mode1.X=toy$mode1.X[!(row.names(toy$mode1.X) %in% m1.rem), 
                      !grepl('const', colnames(toy$mode1.X))],
  mode2.X=toy$mode2.X[!(row.names(toy$mode2.X) %in% m2.rem), 
                      !grepl('const', colnames(toy$mode2.X))],
  mode3.X=toy$mode3.X[!(row.names(toy$mode3.X) %in% m3.rem), 
                      !grepl('const', colnames(toy$mode3.X))],
  resp=toy$resp[!(row.names(toy$mode1.X) %in% m1.rem),
                !(row.names(toy$mode2.X) %in% m2.rem),
                !(row.names(toy$mode3.X) %in% m3.rem)]
)

all.resp.train <- train.data$resp

# Remove 'warm.per' percent of responses to be used as 'warm' test data
warm.idx <- sample(1:prod(dim(train.data$resp)), prod(dim(train.data$resp))*warm.per)
warm.resp <- train.data$resp[sort(warm.idx)]
train.data$resp[warm.idx] <- NA
train.data$delta[warm.idx] <- 0

# Make testing datasets
#######################
if(length(m1.rem)) {
  test.data.m1 <- input_data$new(
    mode1.X=toy$mode1.X[row.names(toy$mode1.X) %in% m1.rem, 
                        !grepl('const', colnames(toy$mode1.X))],
    mode2.X=toy$mode2.X[!(row.names(toy$mode2.X) %in% m2.rem), 
                        !grepl('const', colnames(toy$mode2.X))],
    mode3.X=toy$mode3.X[!(row.names(toy$mode3.X) %in% m3.rem), 
                        !grepl('const', colnames(toy$mode3.X))],
    resp=toy$resp[row.names(toy$mode1.X) %in% m1.rem,
                  !(row.names(toy$mode2.X) %in% m2.rem),
                  !(row.names(toy$mode3.X) %in% m3.rem)]
  )
} else test.data.m1 <- NULL

if(length(m2.rem)) {
  test.data.m2 <- input_data$new(
    mode1.X=toy$mode1.X[!(row.names(toy$mode1.X) %in% m1.rem), 
                        !grepl('const', colnames(toy$mode1.X))],
    mode2.X=toy$mode2.X[row.names(toy$mode2.X) %in% m2.rem, 
                        !grepl('const', colnames(toy$mode2.X))],
    mode3.X=toy$mode3.X[!(row.names(toy$mode3.X) %in% m3.rem), 
                        !grepl('const', colnames(toy$mode3.X))],
    resp=toy$resp[!(row.names(toy$mode1.X) %in% m1.rem),
                  row.names(toy$mode2.X) %in% m2.rem,
                  !(row.names(toy$mode3.X) %in% m3.rem)]
  )
} else test.data.m2 <- NULL

if(length(m3.rem)) {
  test.data.m3 <- input_data$new(
    mode1.X=toy$mode1.X[!(row.names(toy$mode1.X) %in% m1.rem), 
                        !grepl('const', colnames(toy$mode1.X))],
    mode2.X=toy$mode2.X[!(row.names(toy$mode2.X) %in% m2.rem), 
                        !grepl('const', colnames(toy$mode2.X))],
    mode3.X=toy$mode3.X[row.names(toy$mode3.X) %in% m3.rem, 
                        !grepl('const', colnames(toy$mode3.X))],
    resp=toy$resp[!(row.names(toy$mode1.X) %in% m1.rem),
                  !(row.names(toy$mode2.X) %in% m2.rem),
                  row.names(toy$mode3.X) %in% m3.rem]
  )
} else test.data.m3 <- NULL

if(length(m1.rem) & length(m2.rem)) {
  test.data.m1m2 <- input_data$new(
    mode1.X=toy$mode1.X[row.names(toy$mode1.X) %in% m1.rem, 
                        !grepl('const', colnames(toy$mode1.X))],
    mode2.X=toy$mode2.X[row.names(toy$mode2.X) %in% m2.rem, 
                        !grepl('const', colnames(toy$mode2.X))],
    mode3.X=toy$mode3.X[!(row.names(toy$mode3.X) %in% m3.rem), 
                        !grepl('const', colnames(toy$mode3.X))],
    resp=toy$resp[row.names(toy$mode1.X) %in% m1.rem,
                  row.names(toy$mode2.X) %in% m2.rem,
                  !(row.names(toy$mode3.X) %in% m3.rem)]
  )
} else test.data.m1m2 <- NULL

if(length(m1.rem) & length(m3.rem)) {
  test.data.m1m3 <- input_data$new(
    mode1.X=toy$mode1.X[row.names(toy$mode1.X) %in% m1.rem, 
                        !grepl('const', colnames(toy$mode1.X))],
    mode2.X=toy$mode2.X[!(row.names(toy$mode2.X) %in% m2.rem), 
                        !grepl('const', colnames(toy$mode2.X))],
    mode3.X=toy$mode3.X[row.names(toy$mode3.X) %in% m3.rem, 
                        !grepl('const', colnames(toy$mode3.X))],
    resp=toy$resp[row.names(toy$mode1.X) %in% m1.rem,
                  !(row.names(toy$mode2.X) %in% m2.rem),
                  row.names(toy$mode3.X) %in% m3.rem]
  )
} else test.data.m1m3 <- NULL

if(length(m2.rem) & length(m3.rem)) {
  test.data.m2m3 <- input_data$new(
    mode1.X=toy$mode1.X[!(row.names(toy$mode1.X) %in% m1.rem), 
                        !grepl('const', colnames(toy$mode1.X))],
    mode2.X=toy$mode2.X[row.names(toy$mode2.X) %in% m2.rem, 
                        !grepl('const', colnames(toy$mode2.X))],
    mode3.X=toy$mode3.X[row.names(toy$mode3.X) %in% m3.rem, 
                        !grepl('const', colnames(toy$mode3.X))],
    resp=toy$resp[!(row.names(toy$mode1.X) %in% m1.rem),
                  row.names(toy$mode2.X) %in% m2.rem,
                  row.names(toy$mode3.X) %in% m3.rem]
  )
} else test.data.m2m3 <- NULL

if(length(m1.rem) & length(m2.rem) & length(m3.rem)) {
  test.data.m1m2m3 <- input_data$new(
    mode1.X=toy$mode1.X[row.names(toy$mode1.X) %in% m1.rem, 
                        !grepl('const', colnames(toy$mode1.X))],
    mode2.X=toy$mode2.X[row.names(toy$mode2.X) %in% m2.rem, 
                        !grepl('const', colnames(toy$mode2.X))],
    mode3.X=toy$mode3.X[row.names(toy$mode3.X) %in% m3.rem, 
                        !grepl('const', colnames(toy$mode3.X))],
    resp=toy$resp[row.names(toy$mode1.X) %in% m1.rem,
                  row.names(toy$mode2.X) %in% m2.rem,
                  row.names(toy$mode3.X) %in% m3.rem]
  )
} else test.data.m1m2m3 <- NULL

# Get model parameters
model.params <- get_model_params(args[2:length(args)])
model <- mk_model(train.data, model.params)
model$rand_init(model.params)

# Data frame to store RMSEs & explained variances for test data while training 
test.results <- numeric(0)

# Copy of model object to train
trained <- model$clone()

# Save before training for debugging
save.image(file=paste0(out.dir, 'image.Rdata'))

for(i in 1:reps) {
  train(d=train.data, m=trained, params=model.params)

  # Get cold results
  test.results <- test_results(m=trained, d=train.data, test.results=test.results,
                               warm.resp=warm.resp, test.m1=test.data.m1, 
                               test.m2=test.data.m2, test.m3=test.data.m3,
                               test.m1m2=test.data.m1m2, test.m1m3=test.data.m1m3,
                               test.m2m3=test.data.m2m3, test.m1m2m3=test.data.m1m2m3)

  # Save every 10 iterations
  if((i %% 10) ==0) save.image(paste0(out.dir, 'image.Rdata'))
}

# Stop cluster (if using parallel package)
if(.Platform$OS.type == "windows") stopCluster(clust)

save.image(file=paste0(out.dir, 'image.Rdata'))

# Make a tensor of the means for each value of mode2 & mode3 for comparisons
m1.means <- apply(train.data$resp, c(2,3), mean, na.rm=T)
m1.sds <- apply(train.data$resp, c(2,3), sd, na.rm=T)
m2.means <- apply(train.data$resp, c(1,3), mean, na.rm=T)
m2.sds <- apply(train.data$resp, c(1,3), sd, na.rm=T)
m3.means <- apply(train.data$resp, c(1,2), mean, na.rm=T)
m3.sds <- apply(train.data$resp, c(1,2), sd, na.rm=T)

train.means <- train.data$resp
for(i in 1:dim(train.data$resp)[1]) train.means[i,,] <- m1.means
if(length(m1.rem)) {
  m1.test.means <- test.data.m1$resp
  for(i in 1:dim(test.data.m1$resp)[1]) m1.test.means[i,,] <- m1.means
  m1.cold.preds <- test(d=test.data.m1, m=trained)
  m1.mean.RMSE <- sqrt(mean((test.data.m1$resp-m1.test.means)^2, na.rm=T))/sd(train.data$resp, na.rm=T)
  m1.mean.exp.var <- exp_var(test.data.m1$resp, m1.test.means)
}

if(length(m2.rem)) {
  m2.test.means <- test.data.m2$resp
  for(j in 1:dim(test.data.m2$resp)[2]) m2.test.means[,j,] <- m2.means
  m2.cold.preds <- test(d=test.data.m2, m=trained)
  m2.mean.RMSE <- sqrt(mean((test.data.m2$resp-m2.test.means)^2, na.rm=T))/sd(train.data$resp, na.rm=T)
  m2.mean.exp.var <- exp_var(test.data.m2$resp, m2.test.means)
}

if(length(m3.rem)) {
  m3.test.means <- test.data.m3$resp
  for(k in 1:dim(test.data.m3$resp)[3]) m3.test.means[,,k] <- m3.means
  m3.cold.preds <- test(d=test.data.m3, m=trained)
  m3.mean.RMSE <- sqrt(mean((test.data.m3$resp-m3.test.means)^2, na.rm=T))/sd(train.data$resp, na.rm=T)
  m3.mean.exp.var <- exp_var(test.data.m3$resp, m3.test.means)
}

if(length(m1.rem) & length(m2.rem)) {
  m1m2.test.means <- array(0, dim=dim(test.data.m1m2$resp))
  m1m2.means <- apply(train.data$resp, 3, mean, na.rm=T)
  for(i in 1:dim(test.data.m1m2$resp)[1])
    for(j in 1:dim(test.data.m1m2$resp)[2])
      m1m2.test.means[i,j,] <- m1m2.means
  m1m2.cold.preds <- test(d=test.data.m1m2, m=trained)
  m1m2.mean.RMSE <- sqrt(mean((test.data.m1m2$resp-m1m2.test.means)^2, na.rm=T))/sd(train.data$resp, na.rm=T)
  m1m2.mean.exp.var <- exp_var(test.data.m1m2$resp, m1m2.test.means)
}

if(length(m1.rem) & length(m3.rem)) {
  m1m3.test.means <- array(0, dim=dim(test.data.m1m3$resp))
  m1m3.means <- apply(train.data$resp, 2, mean, na.rm=T)
  for(i in 1:dim(test.data.m1m3$resp)[1])
    for(k in 1:dim(test.data.m1m3$resp)[3])
      m1m3.test.means[i,,k] <- m1m3.means
  m1m3.cold.preds <- test(d=test.data.m1m3, m=trained)
  m1m3.mean.RMSE <- sqrt(mean((test.data.m1m3$resp-m1m3.test.means)^2, na.rm=T))/sd(train.data$resp, na.rm=T)
  m1m3.mean.exp.var <- exp_var(test.data.m1m3$resp, m1m3.test.means)
}

if(length(m2.rem) & length(m3.rem)) {
  m2m3.test.means <- array(0, dim=dim(test.data.m2m3$resp))
  m2m3.means <- apply(train.data$resp, 1, mean, na.rm=T)
  for(j in 1:dim(test.data.m2m3$resp)[2])
    for(k in 1:dim(test.data.m2m3$resp)[3])
      m2m3.test.means[,j,k] <- m2m3.means
  m2m3.cold.preds <- test(d=test.data.m2m3, m=trained)
  m2m3.mean.RMSE <- sqrt(mean((test.data.m2m3$resp-m2m3.test.means)^2, na.rm=T))/sd(train.data$resp, na.rm=T)
  m2m3.mean.exp.var <- exp_var(test.data.m2m3$resp, m2m3.test.means)
}

if(length(m1.rem) & length(m2.rem) & length(m3.rem)) {
  m1m2m3.test.means <- array(0, dim=dim(test.data.m1m2m3$resp))
  m1m2m3.mean <- mean(train.data$resp, na.rm=T)
  m1m2m3.test.means[,,] <- m1m2m3.mean
  m1m2m3.cold.preds <- test(d=test.data.m1m2m3, m=trained)
  m1m2m3.mean.RMSE <- sqrt(mean((test.data.m1m2m3$resp-m1m2m3.test.means)^2, na.rm=T))/sd(train.data$resp, na.rm=T)
  m1m2m3.mean.exp.var <- exp_var(test.data.m1m2m3$resp, m1m2m3.test.means)
}

# Calculate predictions
warm.preds <- trained$resp[is.na(train.data$resp)]
warm.mean.preds <- train.means[is.na(train.data$resp)]
warm.mean.exp.var <- exp_var(warm.resp, warm.mean.preds)
warm.mean.RMSE <- sqrt(mean((warm.resp-warm.mean.preds)^2, na.rm=T))/sd(train.data$resp, na.rm=T)
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

## TODO: add lines for mode 3 and combinations thereof

# RMSE
pdf(paste0(out.dir, 'training_RMSEs.pdf'))
par(mfrow=c(1,1))
plot_test_RMSE(test.results, main="Test RMSEs", mean=T)
abline(h=warm.mean.RMSE, col='red', lty=2, lwd=2)
abline(h=m1.mean.RMSE, col='blue', lty=2, lwd=2)
abline(h=m2.mean.RMSE, col='green', lty=2, lwd=2)
abline(h=m1m2.mean.RMSE, col='purple', lty=2, lwd=2)
dev.off()

pdf(paste0(out.dir, 'training_exp_var.pdf'))
plot_test_exp_var(test.results, mean=T,
     main=sprintf('Warm max at %.0f, cold mode 1 max at %d \n Cold mode 2 max at %d Both cold max at %d',
                  which.max(test.results$warm.exp.var), which.max(test.results$m1.exp.var),
                  which.max(test.results$m2.exp.var), which.max(test.results$m1m2.exp.var)))
abline(h=warm.mean.exp.var, col='red', lty=2, lwd=2)
abline(h=m1.mean.exp.var, col='blue', lty=2, lwd=2)
abline(h=m2.mean.exp.var, col='green', lty=2, lwd=2)
abline(h=m1m2.mean.exp.var, col='purple', lty=2, lwd=2)
dev.off()

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
  obs <- test.data.m1$resp[i,j,]
  pred <- m1.cold.preds[i,j,]
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
  obs <- test.data.m2$resp[i,j,]
  pred <- m2.cold.preds[i,j,]
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

  obs <- test.data.m1m2$resp[i,j,]
  pred <- m1m2.cold.preds[i,j,]
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
