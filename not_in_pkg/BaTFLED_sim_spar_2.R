# ---
# title: "BaTFLED (Bayesian Tensor Factorization Linked to External Data) toy example"
# author: "Nathan Lazar"
# date: "Oct 12, 2016"
# ---

library(foreach)      # For parallel processing in loops
library(R6)           # For making memory efficent R6 objects
library(iterators)    # Making iterator objects for efficient parallelization
library(rTensor)      # Needed for multiplying matrices into tensors (could be removed)
library(dplyr)        # General data frame manipulation
library(RColorBrewer) # Colors for plotting
# Packages for registering parallel backend (depends on platform)
ifelse(.Platform$OS.type == "windows", library(doParallel), library(doMC))
library(BaTFLED3D)

args <- commandArgs(TRUE)
in.data <- args[[1]]
out.dir <- args[[2]]

if(sum(grepl('^warm.per=', args)))
  warm.per.new <- as.numeric(sub('warm.per=', '', args[grepl('^warm.per=', args)]))
if(sum(grepl('^cores=', args)))
  cores <- as.numeric(sub('cores=', '', args[grepl('^cores=', args)]))

# If the output directory doesn't exist, create it
print(paste('Output directory:', out.dir))
if (!file.exists(out.dir)){
  dir.create(file.path(out.dir))
}

# Set up backend for parallel execution
# cores <- 16
if(.Platform$OS.type == "windows") {
  clust <- makeCluster(cores)
  registerDoParallel(clust)
} else registerDoMC(cores)

# Function to load the input data and make the training data objects
# available globally
loadRData <- function(in.data, warm.per.new) {
  load(in.data)
  # If the warm percentage is the same as the run that generated the data,
  # then keep don't change train.data and assign warm responses to global env.
  if(warm.per == warm.per.new) {
    assign('warm.per', warm.per, envir=.GlobalEnv)
    assign('warm.resp', warm.resp, envir=.GlobalEnv)
  }
  objs <- mget(c('all.resp.train', 'train.data',
       'test.data.m1', 'test.data.m2', 'test.data.m3',
       'test.data.m1m2', 'test.data.m1m3', 'test.data.m2m3', 'test.data.m1m2m3',
       'toy'))
  for(i in 1:length(objs)) 
    assign(names(objs)[i], objs[[i]], envir=.GlobalEnv)
}

# Load the data
loadRData(in.data, warm.per.new)

# Remove 'warm.per' percent of responses to be used as 'warm' test data
# if the new warm percentage differs from that used in the run generating data
if(!exists('warm.per')) {
  warm.per <- warm.per.new
  n.train.resp <- prod(dim(train.data$resp))
  warm.idx <- sort(sample(1:n.train.resp, n.train.resp*warm.per))
  warm.resp <- all.resp.train[warm.idx]
  train.data$resp <- all.resp.train
  train.data$resp[warm.idx] <- NA
  train.data$delta[,,] <- 1
  train.data$delta[warm.idx] <- 0
}

# Get parameters for the model
##############################
# Set some defaults
reps <- 100
# Override defaults if provided
if(sum(grepl('^reps=', args)))
  reps <- as.numeric(sub('reps=', '', args[grepl('^reps=', args)]))

model.params <- get_model_params(args[3:length(args)])

# Make the model
################

model <- mk_model(train.data, model.params)
model$rand_init(model.params)
model$sigma2 <- var(toy$resp, na.rm=T)

# Train model
######################
# Data frame to store RMSEs & explained variances for test data while training 
test.results <- numeric(0)

trained <- model$clone()  # Copy the model so it won't be updated in place

save.image(paste0(out.dir, '/image.Rdata'))

for(i in 1:reps) {
  train(d=train.data, m=trained, params=model.params)

  if(model.params$decomp=='Tucker') 
    if(model.params$plot) {
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

save.image(paste0(out.dir, '/image.Rdata'))

# Explore the resulting model
###############################

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