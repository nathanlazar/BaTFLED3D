# ---
# title: "BaTFLED (Bayesian Tensor Factorization Linked to External Data) toy example"
# author: "Nathan Lazar"
# date: "May 14, 2016"
# output: 
#   html_document:
#     fig_caption: FALSE
# ---

args <- commandArgs(TRUE)
out.dir <- args[[1]]

# If the output directory doesn't exist, create it
print(paste('Output directory:', out.dir))
if (!file.exists(out.dir)){
  dir.create(file.path(out.dir))
}

## Example usage

# Setting global options
.libPaths("/home/users/lazar/R/x86_64-redhat-linux-gnu-library/3.2")

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
pkgTest('RColorBrewer') # Colors for plotting
# Packages for registering parallel backend (depends on platform)
ifelse(.Platform$OS.type == "windows", library(doParallel), library(doMC))
# Build the package: must be run from the directory containing the tar.gz file
# install.packages('BaTFLED3D_0.0.1.tar.gz', repos = NULL, type="source")  
library(BaTFLED3D)

# Set up backend for parallel execution
cores <- 16
if(.Platform$OS.type == "windows") {
  clust <- makeCluster(cores)
  registerDoParallel(clust)
} else registerDoMC(cores)

## Generating a simulated dataset

args <- list('decomp=Tucker',             # Should be either CP or Tucker factorization
             'row.share=T',               # Share variance across rows of projection matrices
             'm1.rows=30', 'm1.cols=100', # Dimensions of input matrices
             'm2.rows=30', 'm2.cols=100',
             'm3.rows=30', 'm3.cols=100',
             'A1.intercept=T',            # Add intercepts to projection matrices
             'A2.intercept=T', 
             'A3.intercept=T',
             'H1.intercept=T',            # Add intercepts to latent matrices
             'H2.intercept=T', 
             'H3.intercept=T',
             'm1.true=10',                 # Number of predictors affecting output 
             'm2.true=10', 
             'm3.true=10',
             'R1=4', 'R2=4', 'R3=4',      # Core dimension
             'core.spar=0.5',             # Amount of sparsity in the core [0-1], 1=no sparsity
             'noise.sd=0.1'               # Amount of Gaussian noise added to responses
)
data.params <- get_data_params(args)
toy <- mk_toy(data.params)

png(paste0(out.dir, '/model_mats.png'), height=800, width=800)
par(mfrow=c(2,2), mar=c(2,2,2,2), mgp=c(1,1,1))
im_mat(toy$mode1.X, main='Mode 1 X', ylab='Samples', xlab='Predictors')
im_mat(toy$mode1.A, main='Mode 1 A', ylab='Predictors', xlab='Latent factors')
im_mat(toy$mode1.H, main='Mode 1 H', ylab='Samples', xlab='Constant + latent factors')
im_mat(toy$resp[,,1], main='Slice of response tensor',
       xlab='Mode 2 samples', ylab='Mode 1 samples')
dev.off()

## Training a BaTFLED model

train.data <- input_data$new(
  mode1.X=toy$mode1.X[-c(1,2), !grepl('const', colnames(toy$mode1.X))],
  mode2.X=toy$mode2.X[-c(1,2), !grepl('const', colnames(toy$mode2.X))],
  mode3.X=toy$mode3.X[-c(1,2), !grepl('const', colnames(toy$mode3.X))],
  resp=toy$resp[-c(1,2),-c(1,2),-c(1,2)])

# Remove 'warm.per' percent of responses to be used as 'warm' test data
all.resp.train <- train.data$resp
warm.per <- 0.01
warm.idx <- sort(sample(1:prod(dim(train.data$resp)), prod(dim(train.data$resp))*warm.per))
warm.resp <- train.data$resp[warm.idx]
train.data$resp[warm.idx] <- NA
train.data$delta[warm.idx] <- 0

test.data.m1 <- input_data$new(
  mode1.X=toy$mode1.X[c(1,2), !grepl('const', colnames(toy$mode1.X))],
  mode2.X=toy$mode2.X[-c(1,2), !grepl('const', colnames(toy$mode2.X))],
  mode3.X=toy$mode3.X[-c(1,2), !grepl('const', colnames(toy$mode3.X))],
  resp=toy$resp[c(1,2),-c(1,2),-c(1,2)])
test.data.m2 <- input_data$new(
  mode1.X=toy$mode1.X[-c(1,2), !grepl('const', colnames(toy$mode1.X))],
  mode2.X=toy$mode2.X[c(1,2), !grepl('const', colnames(toy$mode2.X))],
  mode3.X=toy$mode3.X[-c(1,2), !grepl('const', colnames(toy$mode3.X))],
 resp=toy$resp[-c(1,2),c(1,2),-c(1,2)])
test.data.m3 <- input_data$new(
  mode1.X=toy$mode1.X[-c(1,2), !grepl('const', colnames(toy$mode1.X))],
  mode2.X=toy$mode2.X[-c(1,2), !grepl('const', colnames(toy$mode2.X))],
  mode3.X=toy$mode3.X[c(1,2), !grepl('const', colnames(toy$mode3.X))],
 resp=toy$resp[-c(1,2),-c(1,2),c(1,2)])

test.data.m1m2 <- input_data$new(
  mode1.X=toy$mode1.X[c(1,2), !grepl('const', colnames(toy$mode1.X))],
  mode2.X=toy$mode2.X[c(1,2), !grepl('const', colnames(toy$mode2.X))],
  mode3.X=toy$mode3.X[-c(1,2), !grepl('const', colnames(toy$mode3.X))],
  resp=toy$resp[c(1,2),c(1,2),-c(1,2)])
test.data.m1m3 <- input_data$new(
  mode1.X=toy$mode1.X[c(1,2), !grepl('const', colnames(toy$mode1.X))],
  mode2.X=toy$mode2.X[-c(1,2), !grepl('const', colnames(toy$mode2.X))],
  mode3.X=toy$mode3.X[c(1,2), !grepl('const', colnames(toy$mode3.X))],
  resp=toy$resp[c(1,2),-c(1,2),c(1,2)])
test.data.m2m3 <- input_data$new(
  mode1.X=toy$mode1.X[-c(1,2), !grepl('const', colnames(toy$mode1.X))],
  mode2.X=toy$mode2.X[c(1,2), !grepl('const', colnames(toy$mode2.X))],
  mode3.X=toy$mode3.X[c(1,2), !grepl('const', colnames(toy$mode3.X))],
  resp=toy$resp[-c(1,2),c(1,2),c(1,2)])

test.data.m1m2m3 <- input_data$new(
  mode1.X=toy$mode1.X[c(1,2), !grepl('const', colnames(toy$mode1.X))],
  mode2.X=toy$mode2.X[c(1,2), !grepl('const', colnames(toy$mode2.X))],
  mode3.X=toy$mode3.X[c(1,2), !grepl('const', colnames(toy$mode3.X))],
  resp=toy$resp[c(1,2),c(1,2),c(1,2)])

## Make the model

args <- list('decomp=Tucker', 'row.share=T',
             'A1.intercept=T', 'A2.intercept=T', 'A3.intercept=T',
             'H1.intercept=T', 'H2.intercept=T', 'H3.intercept=T',
             'plot=T', 'verbose=F',
             'R1=4', 'R2=4', 'R3=4',
             paste0('sigma2=', var(toy$resp)),
             'm1.alpha=1e-5', 'm2.alpha=1e-5', 'm3.alpha=1e-5',
             'm1.beta=1e5', 'm2.beta=1e5', 'm3.beta=1e5',
             'core.3D.alpha=1e-5', 'core.3D.beta=1e5',
             'parallel=T', 'cores=4', 'lower.bnd=T', 
             'update.order=c(3,1,2)', 'show.mode=c(1,2,3)')
model.params <- get_model_params(args)
model <- mk_model(train.data, model.params)
model$rand_init(model.params)

## Train model

reps <- 100

# Data frame to store RMSEs & explained variances for test data while training 
test.results <- numeric(0)

trained <- model$clone()  # Copy the model so it won't be updated in place

for(i in 1:reps) {
  train(d=train.data, m=trained, params=model.params)

  if(model.params$decomp=='Tucker') {
    rng <- range(trained$core.mean)
    im_mat(trained$core.mean[1,,], zlim=rng, main="core[1,,]")
  }

  # Get cold results
  test.results <- test_results(m=trained, d=train.data, test.results=test.results,
                     warm.resp=warm.resp, test.m1=test.data.m1, test.m2=test.data.m2,
                     test.m3=test.data.m3, test.m1m2=test.data.m1m2, 
                     test.m1m3=test.data.m1m3, test.m2m3=test.data.m2m3,
                     test.m1m2m3=test.data.m1m2m3)
}

# Stop cluster (if using parallel package)
if(.Platform$OS.type == "windows") stopCluster(clust)

## Exploring the resulting model

png(paste0(out.dir, '/compare_mats.png'), height=800, width=800)
par(mfrow=c(3,4), mar=c(2,2,2,2))
im_2_mat(toy$mode1.A, trained$mode1.A.mean, main1='Mode 1 A true', main2='Mode 1 A learned')
m1.H.order <- im_2_mat(toy$mode1.H[-(1:2),], trained$mode1.H.mean, 
                       main1='Mode 1 H true', main2='Mode 1 H learned')
im_2_mat(toy$mode2.A, trained$mode2.A.mean, main1='Mode 2 A true', main2='Mode 2 A learned')
m2.H.order <- im_2_mat(toy$mode2.H[-(1:2),], trained$mode2.H.mean, 
                       main1='Mode 2 H true', main2='Mode 2 H learned')
im_2_mat(toy$mode3.A, trained$mode3.A.mean, main1='Mode 3 A true', main2='Mode 3 A learned')
m3.H.order <- im_2_mat(toy$mode3.H, trained$mode3.H.mean, 
                       main1='Mode 3 H true', main2='Mode 3 H learned')
dev.off()

png(paste0(out.dir, '/compare_core.png'), height=800, width=800)
par(mfrow=c(2,4))
if(model.params$decomp=='Tucker') {
  core.reorder <- trained$core.mean[abs(m1.H.order), abs(m2.H.order), abs(m3.H.order)]
  core.reorder <- core.reorder * outer(sign(m1.H.order), outer(sign(m2.H.order), sign(m3.H.order)))
  im_2_mat(toy$core[,,1], core.reorder[,,1], sort=F, center=T, scale='all',
           main1='Core slice 1 true', main2='Core slice 1 learned')
  im_2_mat(toy$core[,,2], core.reorder[,,2], sort=F, center=T, scale='all',
           main1='Core slice 2 true', main2='Core slice 2 learned')
  im_2_mat(toy$core[,,3], core.reorder[,,3], sort=F, center=T, scale='all',
           main1='Core slice 3 true', main2='Core slice 3 learned')
  im_2_mat(toy$core[,,4], core.reorder[,,4], sort=F, center=T, scale='all',
           main1='Core slice 4 true', main2='Core slice 4 learned')
}
dev.off()

png(paste0(out.dir, '/roc.png'), height=800, width=800)
par(mfrow=c(1,3), mar=c(3,3,2,2), mgp=c(2,1,0))
plot_roc(toy$mode1.A, trained$mode1.A.mean, main='Mode 1')
plot_roc(toy$mode2.A, trained$mode2.A.mean, main='Mode 2')
plot_roc(toy$mode3.A, trained$mode3.A.mean, main='Mode 3')
dev.off()

png(paste0(out.dir, '/training.png'), height=800, width=800)
par(mfrow=c(2,2), mar=c(3,3,2,2), mgp=c(2,1,0))
if(model.params$lower.bnd) plot(trained$lower.bnd, main='Lower bound')
if(model.params$RMSE)      plot(trained$RMSE, main='Training RMSE')
if(model.params$exp.var)   plot(trained$exp.var, main='Training explained variation')
if(model.params$cor)       plot(trained$s.cor, main='Spearman correlation')
dev.off()

png(paste0(out.dir, '/test.png'), height=800, width=800)
par(mfrow=c(2,2))
plot_test_RMSE(test.results, main="Test RMSEs")
plot_test_exp_var(test.results, main="Test explained variances")
plot_test_cor(test.results, main="Pearson correlation")
plot_test_cor(test.results, main="Spearman correlation", method='spearman')
dev.off()

save.image(paste0(out.dir, '/image.Rdata'))