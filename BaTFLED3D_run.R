#!/usr/bin/env Rscript

# nathan dot lazar at gmail dot com

# Usage: BaTFLED_run.R <mode1 data> <mode2 data> <mode3 data>
#                      <response> 
#                      <mode1 fold matrix> <mode2 fold matrix> <mode3 fold matrix>
#                      <output dir>
#                      <options...>

# The first 8 arguments must be in order
# The order of the rest of the arguments doesn't mattter and they can be absent

print('Opened file')

library(methods)
library(BaTFLED3D)      # Note: package must be installed from .tar file
# library(pryr)           # For monitoring memory usage
# library(microbenchmark) # For monitoring run time

print('Packages loaded')
args <- commandArgs(TRUE)

# DREAM7 dataset
# args <- list('../DREAM7/DREAM7_cl.Rdata',
#              '../DREAM7/DREAM7_dr.Rdata',
#              'none',
#              '../DREAM7/DREAM7_resp.Rdata',
#              '../DREAM7/cv_folds_cl_10.Rdata',
#              '../DREAM7/cv_folds_dr_10.Rdata',
#              'none',
#              'DREAM7results/test/',
#              'decomp=Tucker',
#              'row.share=F',
#              'reps=30',
#              'warm.per=0.05',
#              'seed=NA',
#              'multiplier=0',
#              'normalize.for=mode1',
#              'early.stop=NA', 'plot=T', 'verbose=F',
#              'remove.per=0', 'remove.start=Inf', 
#              'm1.remove.lmt=0', 'm2.remove.lmt',
#              'R1=6', 'R2=3', 'R3=3',
#              'A1.intercept=T', 'A2.intercept=T', 'A3.intercept=F',
#              'H1.intercept=T', 'H2.intercept=T', 'H3.intercept=T',
#              'm1.alpha=1e-10', 'm2.alpha=1e-10', 'm3.alpha=1',
#              'm1.beta=1e10', 'm2.beta=1e10', 'm3.beta=1',
#              'core.0D.alpha=1', 'core.0D.beta=1',
#              'core.1D.alpha=1e-4', 'core.1D.beta=1e4',
#              'core.2D.alpha=1e-8', 'core.2D.beta=1e8',
#              'core.3D.alpha=1e-12', 'core.3D.beta=1e12',
#              'm1.sigma2=0.01', 'm2.sigma2=0.01', 'm3.sigma2=1',
#              'sigma2=auto', 'parallel=T', 'cores=4',
#              'exp.var=T', 'lower.bnd=T', 'history=F', 'fold=0')

#'../cell_line_drug_datasets/BC_CellLines/DREAM7/Costello2014/norm_mean_tensor.Rdata',

# Report arguments and summary stats of the input objects
print(unlist(args))

# Read in options
m1.file <-        args[[1]]
m2.file <-        args[[2]]
m3.file <-        args[[3]]
resp.file <-      args[[4]]
m1.fold.file <-   args[[5]]
m2.fold.file <-   args[[6]]
m3.fold.file <-   args[[7]]
out.dir <-        paste0(args[[8]], '/')

# If the output directory doesn't exist, create it
print(paste('Output directory:', out.dir))
if (!file.exists(out.dir)){
  dir.create(file.path(out.dir))
}
save.image(paste0(out.dir, 'image.Rdata'))

# Get parameters
params <- get_model_params(args[9:length(args)])

# Set some defaults 
reps <- 10; warm.per <- 0.05
multiplier <- 0
normalize.for='none'
history <- F
fold <- 1
kernelize <- F

# Override defaults if provided
if(sum(grepl('^reps=', args)))
  reps <- as.numeric(sub('reps=', '', args[grepl('^reps=', args)]))
if(sum(grepl('^warm.per=', args)))
  warm.per <- as.numeric(sub('warm.per=', '', args[grepl('^warm.per=', args)]))
if(sum(grepl('^multiplier=', args)))
  multiplier <- as.numeric(sub('multiplier=', '', args[grepl('^multiplier=', args)]))
if(sum(grepl('^history=', args)))
  history <- as.logical(sub('history=', '', args[grepl('^history=', args)]))
if(sum(grepl('^normalize.for=', args)))
  normalize.for <- sub('normalize.for=', '', args[grepl('^normalize.for=', args)])
if(sum(grepl('^fold=', args)))
  fold <- as.numeric(sub('fold=', '', args[grepl('^fold=', args)]))+1
if(sum(grepl('^kernelize=', args)))
  kernelize <- as.logical(sub('kernelize=', '', args[grepl('^kernelize=', args)]))

if(params$parallel) {
  # Packages for registering parallel backend (depends on platform)
  if(.Platform$OS.type == "windows") {
    library(doParallel)
  } else {
    library(doMC)
  }
}

# rm(args)

# Load input data
############################################################
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Load and rename input data if supplied
if(m1.file != 'none') m1.mat <- loadRData(m1.file)
if(m2.file != 'none') m2.mat <- loadRData(m2.file)
if(m3.file != 'none') m3.mat <- loadRData(m3.file)

# Load the response data
resp.tens <- loadRData(resp.file)

# Load the matrices with names for cross validation
if(m1.fold.file != 'none') m1.cv.fold <- loadRData(m1.fold.file)
if(m2.fold.file != 'none') m2.cv.fold <- loadRData(m2.fold.file)
if(m3.fold.file != 'none') m3.cv.fold <- loadRData(m3.fold.file)

# Make empty matrices if no data is present
if(m1.file == 'none') {
  m1.mat <- matrix(nrow=dim(resp.tens)[1], ncol=0)
  dimnames(m1.mat) <- list(dimnames(resp.tens)[[1]])
}

if(m2.file == 'none') {
  m2.mat <- matrix(nrow=dim(resp.tens)[2], ncol=0)
  dimnames(m2.mat) <- list(dimnames(resp.tens)[[2]])
}

if(m3.file == 'none') {
  m3.mat <- matrix(nrow=dim(resp.tens)[3], ncol=0)
  dimnames(m3.mat) <- list(dimnames(resp.tens)[[3]])
}

# Subset predictors to those that have responses
m1.mat <- m1.mat[rownames(m1.mat) %in% dimnames(resp.tens)[[1]],]
m2.mat <- m2.mat[rownames(m2.mat) %in% dimnames(resp.tens)[[2]],]
m3.mat <- m3.mat[rownames(m3.mat) %in% dimnames(resp.tens)[[3]],]

# Transform predictors to kernel versions if 'kernerlize=T'
if(kernelize) {
  m1.mat <- kernelize(m1.mat)
  m2.mat <- kernelize(m2.mat)
  m3.mat <- kernelize(m3.mat)
}

# Scale responses if multiplier non-zero
if(multiplier != 0) {
  resp.tens <- resp.tens * multiplier
}

# Remove the samples for this cross validation fold to be test
# data for this run and make test.data object.
###############################################################
if(exists('m1.cv.fold')) {
  test.m1 <- m1.cv.fold[fold,]
  test.m1.mat <- m1.mat[row.names(m1.mat) %in% test.m1,,drop=F]
  train.m1.mat <- m1.mat[!(row.names(m1.mat) %in% test.m1),,drop=F]
} else {
  test.m1.mat <- m1.mat[0,]
  train.m1.mat <- m1.mat
}

if(exists('m2.cv.fold')) {
  test.m2 <- m2.cv.fold[fold,]
  test.m2.mat <- m2.mat[row.names(m2.mat) %in% test.m2,,drop=F]
  train.m2.mat <- m2.mat[!(row.names(m2.mat) %in% test.m2),,drop=F]
} else {
  test.m2.mat <- m2.mat[0,]
  train.m2.mat <- m2.mat
}

if(exists('m3.cv.fold')) {
  test.m3 <- m3.cv.fold[fold,]
  test.m3.mat <- m3.mat[row.names(m3.mat) %in% test.m3,,drop=F]
  train.m3.mat <- m3.mat[!(row.names(m3.mat) %in% test.m3),,drop=F]
} else {
  test.m3.mat <- m3.mat[0,]
  train.m3.mat <- m3.mat
}

# Remove predictors that have no variance in the training data
# or have only one non-zero value
if(ncol(train.m1.mat)) {
  train.m1.mat <- train.m1.mat[,apply(train.m1.mat, 2, function(x) sd(x) > 0)]
  m1.one <- apply(train.m1.mat, 2, function(x) sum(x==1)==1) &
    apply(train.m1.mat, 2, function(x) sum(x==0)==(nrow(train.m1.mat)-1))
  train.m1.mat <- train.m1.mat[,!m1.one]
  # Also remove these predictors from testing data
  test.m1.mat <- test.m1.mat[,match(dimnames(train.m1.mat)[[2]], 
                             dimnames(test.m1.mat)[[2]]),drop=F]
}

# Same for mode 2
if(ncol(train.m2.mat)) {
  train.m2.mat <- train.m2.mat[,apply(train.m2.mat, 2, function(x) sd(x) > 0)]
  m2.one <- apply(train.m2.mat, 2, function(x) sum(x==1)==1) &
    apply(train.m2.mat, 2, function(x) sum(x==0)==(nrow(train.m2.mat)-1))
  train.m2.mat <- train.m2.mat[,!m2.one]
  # Also remove these predictors from testing data
  test.m2.mat <- test.m2.mat[,match(dimnames(train.m2.mat)[[2]], 
                             dimnames(test.m2.mat)[[2]]), drop=F]
}

# Same for mode 3
if(ncol(train.m3.mat)) {
  train.m3.mat <- train.m3.mat[,apply(train.m3.mat, 2, function(x) sd(x) > 0)]
  m3.one <- apply(train.m3.mat, 2, function(x) sum(x==1)==1) &
    apply(train.m3.mat, 2, function(x) sum(x==0)==(nrow(train.m3.mat)-1))
  train.m3.mat <- train.m3.mat[,!m3.one]
  # Also remove these predictors from testing data
  test.m3.mat <- test.m3.mat[,match(dimnames(train.m3.mat)[[2]], 
                             dimnames(test.m3.mat)[[2]]),drop=F]
}

# Reorder responses to match the inputs
resp.tens <- resp.tens[match(rownames(m1.mat), dimnames(resp.tens)[[1]], nomatch=0),
                       match(rownames(m2.mat), dimnames(resp.tens)[[2]], nomatch=0),
                       match(rownames(m3.mat), dimnames(resp.tens)[[3]], nomatch=0), drop=F]

resp.train <- resp.tens[match(rownames(train.m1.mat), dimnames(resp.tens)[[1]], nomatch=0),
                        match(rownames(train.m2.mat), dimnames(resp.tens)[[2]], nomatch=0),
                        match(rownames(train.m3.mat), dimnames(resp.tens)[[3]], nomatch=0), drop=F]

# Normalize for the specific prediction task (mode1, 2 or 3)
if(normalize.for == 'mode1') norm.dims <- c(2,3)
if(normalize.for == 'mode2') norm.dims <- c(1,3)
if(normalize.for == 'mode3') norm.dims <- c(1,2)
if(normalize.for == 'mode12') norm.dims <- 3
if(normalize.for == 'mode13') norm.dims <- 2
if(normalize.for == 'mode23') norm.dims <- 1

if(exists('norm.dims')) {
  means <- apply(resp.train, norm.dims, mean, na.rm=T)
  sds <- apply(resp.train, norm.dims, sd, na.rm=T)
  # If the standard deviation is zero (e.g. drug PS-1145), don't divide by
  # anything. We don't remove it since something may be learned by other
  # similar samples.
  sds[sds==0] <- 1
  resp.tens <- sweep(resp.tens, norm.dims, means, FUN='-')
  resp.tens <- sweep(resp.tens, norm.dims, sds, FUN='/')
} else if(normalize.for == 'mode123') {
    means <- mean(resp.train, na.rm=T)
    sds <- sd(resp.train, na.rm=T)
    resp.tens <- resp.tens - means
    resp.tens <- resp.tens/sds
}

# If params$sigma2 is 'auto' set this value to the variance of the response data
if(params$sigma2=='auto')
  params$sigma2 <- var(resp.tens, na.rm=T)

# Clean up
# rm(m1.one, m2.one)
# rm(m1.cv.fold, m2.cv.fold, m3.cv.fold)

# Make the response objects
resp.train <- resp.tens

if(nrow(test.m1.mat)) {
  resp.train <- resp.train[dimnames(resp.train)[[1]] %in% row.names(train.m1.mat),,,drop=F]
  resp.test.m1 <- resp.tens
  resp.test.m1 <- resp.tens[dimnames(resp.tens)[[1]] %in% test.m1,,,drop=F]
}

if(nrow(test.m2.mat)) {
  resp.train <- resp.train[,dimnames(resp.train)[[2]] %in% row.names(train.m2.mat),,drop=F]
  resp.test.m2 <- resp.tens
  resp.test.m2 <- resp.tens[,dimnames(resp.tens)[[2]] %in% test.m2,,drop=F]
}

if(nrow(test.m3.mat)) {
  resp.train <- resp.train[,,dimnames(resp.train)[[3]] %in% row.names(train.m3.mat),drop=F]
  resp.test.m3 <- resp.tens
  resp.test.m3 <- resp.tens[,,dimnames(resp.tens)[[3]] %in% test.m3,drop=F]
}

if(nrow(test.m1.mat) & nrow(test.m2.mat)) {
  resp.test.m1 <- resp.test.m1[,dimnames(resp.test.m1)[[2]] %in% row.names(train.m2.mat),,drop=F]
  resp.test.m2 <- resp.test.m2[dimnames(resp.test.m2)[[1]] %in% row.names(train.m1.mat),,,drop=F]
  resp.test.m1m2 <- resp.tens[dimnames(resp.tens)[[1]] %in% test.m1,
                              dimnames(resp.tens)[[2]] %in% test.m2,,drop=F]
}

if(nrow(test.m1.mat) & nrow(test.m3.mat)) {
  resp.test.m1 <- resp.test.m1[,,dimnames(resp.test.m1)[[3]] %in% row.names(train.m3.mat),drop=F]
  resp.test.m3 <- resp.test.m3[dimnames(resp.test.m3)[[1]] %in% row.names(train.m1.mat),,,drop=F]
  resp.test.m1m3 <- resp.tens[dimnames(resp.tens)[[1]] %in% test.m1,,
                              dimnames(resp.tens)[[3]] %in% test.m3,drop=F]
}

if(nrow(test.m2.mat) & nrow(test.m3.mat)) {
  resp.test.m2 <- resp.test.m2[,,dimnames(resp.test.m2)[[3]] %in% row.names(train.m3.mat),drop=F]
  resp.test.m3 <- resp.test.m3[,dimnames(resp.test.m3)[[2]] %in% row.names(train.m2.mat),,drop=F]
  resp.test.m2m3 <- resp.tens[,dimnames(resp.tens)[[2]] %in% test.m2,
                              dimnames(resp.tens)[[3]] %in% test.m3,drop=F]
}

if(nrow(test.m1.mat) & nrow(test.m2.mat) & nrow(test.m3.mat)) {
  resp.test.m1m2 <- resp.test.m1m2[,,dimnames(resp.test.m1m2)[[3]] %in% row.names(train.m3.mat),drop=F]
  resp.test.m1m3 <- resp.test.m1m3[,dimnames(resp.test.m1m3)[[2]] %in% row.names(train.m2.mat),,drop=F]
  resp.test.m2m3 <- resp.test.m2m3[dimnames(resp.test.m2m3)[[1]] %in% row.names(train.m1.mat),,,drop=F]
  resp.test.m1m2m3 <- resp.tens[dimnames(resp.tens)[[1]] %in% row.names(test.m1.mat),
                                dimnames(resp.tens)[[2]] %in% row.names(test.m2.mat),
                                dimnames(resp.tens)[[3]] %in% row.names(test.m3.mat),drop=F]
}

# Remove warm.per percent of training responses for warm prediction
all.resp.train <- resp.train   # Stores responses before warms removed

if(warm.per > 0) {
  I <- dim(resp.train)[1]
  J <- dim(resp.train)[2]
  K <- dim(resp.train)[3]
  # If all modes are different, assume the third mode is dose and remove 
  # warm data across all doses
  if((m1.file != m2.file) & (m1.file != m3.file) & (m2.file != m3.file)) {0
    mask <- sample(I*J, round(warm.per*I*J))
    for(k in 1:K) 
      resp.train[,,k][mask] <- NA
  } 
  # otherwise remove randomly but the same for matching modes
  if(m1.file == m2.file) {
    mask <- sample(I*K, (round(warm.per*I*K)))
    for(m in mask) {
      i <- m %% I + 1
      k <- (m-1) %/% I +1
      resp.train[i,i,k] <- NA
    }
  }
  if(m2.file == m3.file) {
    mask <- sample(I*J, (round(warm.per*I*J)))
    for(m in mask) {
      i <- m %% I + 1
      j <- (m-1) %/% I +1
    resp.train[i,j,j] <- NA
    }
  }
}

m1.means <- apply(resp.train, c(2,3), mean, na.rm=T)
m1.sds <- apply(resp.train, c(2,3), sd, na.rm=T)
m2.means <- apply(resp.train, c(1,3), mean, na.rm=T)
m2.sds <- apply(resp.train, c(1,3), sd, na.rm=T)
m3.means <- apply(resp.train, c(1,2), mean, na.rm=T)
m3.sds <- apply(resp.train, c(1,2), sd, na.rm=T)

print(paste0(round(sum(is.na(resp.train))/prod(dim(resp.train))*100, 2),
            "% of responses are missing, ",
            round(sum(is.na(resp.train) & !is.na(all.resp.train))/prod(dim(resp.train))*100,2),
            "% are used for warm predictions"))

# Make the training data object (input_data_3d):
if(!exists('train.m1.mat')) train.m1.mat <- m1.mat
if(!exists('train.m2.mat')) train.m2.mat <- m2.mat
if(!exists('train.m3.mat')) train.m3.mat <- m3.mat

train.data <- input_data$new(mode1.X=train.m1.mat, 
                             mode2.X=train.m2.mat,
                             mode3.X=train.m3.mat, 
                             resp=resp.train)

# Make the testing data objects (input_data_3d):
if(nrow(test.m1.mat)) {
  test.data.m1 <- input_data$new(mode1.X=test.m1.mat,
                                 mode2.X=train.m2.mat,
                                 mode3.X=train.m3.mat,
                                 resp=resp.test.m1)
} else test.data.m1 <- numeric(0)

if(nrow(test.m2.mat)) {
  test.data.m2 <- input_data$new(mode1.X=train.m1.mat,
                                 mode2.X=test.m2.mat,
                                 mode3.X=train.m3.mat,
                                 resp=resp.test.m2)
} else test.data.m2 <- numeric(0)

if(nrow(test.m3.mat)) {
  test.data.m3 <- input_data$new(mode1.X=train.m1.mat,
                                 mode2.X=train.m2.mat,
                                 mode3.X=test.m3.mat,
                                 resp=resp.test.m3)
} else test.data.m3 <- numeric(0)

if(nrow(test.m1.mat) & nrow(test.m2.mat)) {
  test.data.m1m2 <- input_data$new(mode1.X=test.m1.mat,
                                   mode2.X=test.m2.mat,
                                   mode3.X=train.m3.mat,
                                   resp=resp.test.m1m2)
} else test.data.m1m2 <- numeric(0)

if(nrow(test.m1.mat) & nrow(test.m3.mat)) {
  test.data.m1m3 <- input_data$new(mode1.X=test.m1.mat,
                                   mode2.X=train.m2.mat,
                                   mode3.X=test.m3.mat,
                                   resp=resp.test.m1m3)
} else test.data.m1m3 <- numeric(0)

if(nrow(test.m2.mat) & nrow(test.m3.mat)) {
  test.data.m2m3 <- input_data$new(mode1.X=train.m1.mat,
                                   mode2.X=test.m2.mat,
                                   mode3.X=test.m3.mat,
                                   resp=resp.test.m2m3)
} else test.data.m2m3 <- numeric(0)

if(nrow(test.m1.mat) & nrow(test.m2.mat) & nrow(test.m3.mat)) {
  test.data.m1m2m3 <- input_data$new(mode1.X=test.m1.mat,
                                     mode2.X=test.m2.mat,
                                     mode3.X=test.m3.mat,
                                     resp=resp.test.m1m2m3)
} else test.data.m1m2m3 <- numeric(0)

print(sprintf('BaTFLED is being trained with %d mode1, %d mode2 and %d mode3 values.',
              dim(train.data$resp)[1], dim(train.data$resp)[2], dim(train.data$resp)[3]))
print(sprintf('Cold start predictions will be made for %d mode1 %d mode2 and %d mode3 values,',
              dim(resp.tens)[[1]] - dim(resp.train)[[1]],
              dim(resp.tens)[[2]] - dim(resp.train)[[2]],
              dim(resp.tens)[[3]] - dim(resp.train)[[3]]))
print(sprintf('using %d mode 1 predictors, %d mode 2 and %d mode3 predictors.',
              ncol(m1.mat), ncol(m2.mat), ncol(m3.mat)))

# Clean up:
# rm(train.m1.mat, train.m2.mat, train.m3.mat, test.m1.mat, test.m2.mat, test.m3.mat, resp.train, 
#    resp.test.m1, resp.test.m2, resp.test.m3, resp.test.m1m2, resp.test.m1m3, resp.test.m2m3,
#    resp.test.m1m2m3, m1.mat, m2.mat, m3.mat, resp.tens)

# Make the factorization model object (either CP or Tucker)
if(params$decomp=='CP') model <- CP_model$new(d=train.data, params=params)
if(params$decomp=='Tucker') model <- Tucker_model$new(d=train.data, params=params)

# Initialize the model with random values in the A and H matrices
model$rand_init(params)

########################################################################
################## Train model and explore results #####################
########################################################################

# Open a connection to store results while training if history == T
if(history) out <- file(paste0(out.dir, 'training_models.serialize'), "w")

# Set up backend for parallel execution
if(params$parallel) {
  if(.Platform$OS.type == "windows") {
    clust <- parallel::makeCluster(params$cores)
    doParallel::registerDoParallel(clust)
  } else {
    doMC::registerDoMC(params$cores)
  }
}

# Data frame to store RMSEs & explained variances while training 
# (filled in by cold_results function)
test.results <- numeric(0)

# Calculate warm resposes
warm.resp <- all.resp.train[is.na(train.data$resp)]

trained <- model$clone()
save.image(paste0(out.dir, 'image.Rdata'))

for(i in (trained$iter+1):reps) {
  train(d=train.data, m=trained, new.iter=1, params=params)
  
  # write out a serialized version of the trained model
  if(history) {print('writing out model'); serialize(trained, out)}
  
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
if(params$parallel) {
  if(.Platform$OS.type == "windows") stopCluster(clust)
}

# close the connection to the serialized data file
if(history) close(out)

save.image(paste0(out.dir, 'image.Rdata'))

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
