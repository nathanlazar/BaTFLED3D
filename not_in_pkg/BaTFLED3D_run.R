#/usr/bin/env Rscript

# nathan dot lazar at gmail dot com

# Usage: BaTFLED3D_run.R <mode1 data> <mode2 data> <mode3 data>
#                        <response>
#                        <mode1 fold matrix>
#                        <mode2 fold matrix>
#                        <mode3 fold matrix>
#                        <output dir>
#                        <options...>

# The first 8 arguments must be in order
# The order of the rest of the arguments doesn't mattter and they can be absent

print('Opened file')

library(methods)
library(BaTFLED3D)      # Note: package must be installed from .tar file
sessionInfo()
# library(pryr)           # For monitoring memory usage
# library(microbenchmark) # For monitoring run time

source('kern_combine.R')  # Script to kernelize input features

print('Packages loaded')
args <- commandArgs(TRUE)

save.a.lot <- 1  # Save after this many iterations

# DREAM data predicting for cell lines
testing <- F
if(testing) 
  args <- list('../HeiserData/Cell_Line_Data/Kerns/train_harmonic_kern.Rdata',
     # '../HeiserData/Cell_Line_Data/anno_mat.Rdata',
     # '../HeiserData/Cell_Line_Data/exp_mat.Rdata',
     # '../HeiserData/Drug_Data/dr_all_mat.Rdata',
     'none',
     'none',
     '../HeiserData/Responses/norm_tens1000.Rdata',
     '../HeiserData/Responses/Split/cl_10cv_mat.Rdata',
    'none', 'none',
    'HEISERresults/test_5x5x5_kern_nodr', 
    'decomp=Tucker', 'row.share=T',
    'reps=30', 'warm.per=0.01', 'plot=F', 'cores=12',
    'A1.intercept=T', 'A2.intercept=F', 'A3.intercept=F', 
    'H1.intercept=T', 'H2.intercept=T', 'H3.intercept=T',  
    'm1.sigma2=0.01', 'm2.sigma2=1','m3.sigma2=1',
    'm1.alpha=1', 'm2.alpha=1', 'm3.alpha=1',
    'm1.beta=1',   'm2.beta=1',   'm3.beta=1',
    'core.alpha=1',    'core.beta=1',
    'core.1D.alpha=1', 'core.1D.beta=1',
    'core.2D.alpha=1', 'core.2D.beta=1',
    'core.3D.alpha=1', 'core.3D.beta=1',
    'sigma2=auto',
    'R1=5', 'R2=5', 'R3=5', 'normalize.for=none', 
    'kern=F', 
    # 'kern.combine=linear', 'kern.scale=T',
    # 'exp.s=1', 'meth.s=1', 'rppa.s=1', 'cn.s=1', 'padel.s=1',
    'scale=F', 'fold=0', 'center=F')

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
parallel <- T
normalize.for='none'
fold <- 1
scale <- F
center <- F
kern <- F
kern.combine <- 'mean'
kern.scale <- F

# Override defaults if provided
if(sum(grepl('^reps=', args)))
  reps <- as.numeric(sub('reps=', '', args[grepl('^reps=', args)]))
if(sum(grepl('^warm.per=', args)))
  warm.per <- as.numeric(sub('warm.per=', '', args[grepl('^warm.per=', args)]))
if(sum(grepl('^multiplier=', args)))
  multiplier <- as.numeric(sub('multiplier=', '', args[grepl('^multiplier=', args)]))
if(sum(grepl('^normalize.for=', args)))
  normalize.for <- sub('normalize.for=', '', args[grepl('^normalize.for=', args)])
if(sum(grepl('^fold=', args)))
  fold <- as.numeric(sub('fold=', '', args[grepl('^fold=', args)]))+1
if(sum(grepl('^scale=', args)))
  scale <- as.logical(sub('scale=', '', args[grepl('^scale=', args)]))
if(sum(grepl('^center=', args)))
  center <- as.logical(sub('center=', '', args[grepl('^center=', args)]))
if(sum(grepl('^kern=', args)))
  kern <- as.logical(sub('kern=', '', args[grepl('^kern=', args)]))
if(sum(grepl('^kern.combine=', args)))
  kern.combine <- sub('kern.combine=', '', args[grepl('^kern.combine=', args)])
if(sum(grepl('^kern.scale=', args)))
  kern.scale <- as.logical(sub('kern.scale=', '', args[grepl('^kern.scale=', args)]))

if(params$parallel) {
  # Packages for registering parallel backend (depends on platform)
  if(.Platform$OS.type == "windows") {
    library(doParallel)
  } else {
    library(doMC)
  }
}

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
resp.tens[resp.tens == Inf] <- NA

# Load the matrices with names for cross validation
if(m1.fold.file != 'none') {
  m1.cv.fold <- loadRData(m1.fold.file)
  all.m1 <- unique(as.vector(m1.cv.fold))
  m1.mat <- m1.mat[rownames(m1.mat) %in% all.m1,]
}
if(m2.fold.file != 'none') {
  m2.cv.fold <- loadRData(m2.fold.file)
  all.m2 <- unique(as.vector(m2.cv.fold))
  m2.mat <- m2.mat[rownames(m2.mat) %in% all.m2,]
}
if(m3.fold.file != 'none') {
  m3.cv.fold <- loadRData(m3.fold.file)
  all.m3 <- unique(as.vector(m3.cv.fold))
  m3.mat <- m3.mat[rownames(m3.mat) %in% all.m3,]
}

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

# Subset samples to those that have responses
m1.mat <- m1.mat[rownames(m1.mat) %in% dimnames(resp.tens)[[1]],,drop=F]
m2.mat <- m2.mat[rownames(m2.mat) %in% dimnames(resp.tens)[[2]],,drop=F]
m3.mat <- m3.mat[rownames(m3.mat) %in% dimnames(resp.tens)[[3]],,drop=F]

# Set NA values to 0 in feature matrices
# m1.mat[is.na(m1.mat)] <- 0
# m2.mat[is.na(m2.mat)] <- 0
# m3.mat[is.na(m3.mat)] <- 0

# Scale responses if multiplier non-zero
if(multiplier != 0) {
  resp.tens <- resp.tens * multiplier
}

# Remove the samples for this cross validation fold to be test
# data for this run.
###############################################################
if(exists('m1.cv.fold')) {
  test.m1 <- m1.cv.fold[fold,]
  test.m1.mat <- m1.mat[row.names(m1.mat) %in% test.m1,
    !(colnames(m1.mat) %in% test.m1),drop=F]
  train.m1.mat <- m1.mat[!(row.names(m1.mat) %in% test.m1),
    !(colnames(m1.mat) %in% test.m1),drop=F]
} else {
  test.m1.mat <- m1.mat[0,,drop=F]
  train.m1.mat <- m1.mat
}

if(exists('m2.cv.fold')) {
  test.m2 <- m2.cv.fold[fold,]
  test.m2.mat <- m2.mat[row.names(m2.mat) %in% test.m2,
    !(colnames(m2.mat) %in% m2.mat),drop=F]
  train.m2.mat <- m2.mat[!(row.names(m2.mat) %in% test.m2),
    !(colnames(m2.mat) %in% m2.mat),drop=F]
} else {
  test.m2.mat <- m2.mat[0,,drop=F]
  train.m2.mat <- m2.mat
}

if(exists('m3.cv.fold')) {
  test.m3 <- m3.cv.fold[fold,]
  test.m3.mat <- m3.mat[row.names(m3.mat) %in% test.m3,
    !(colnames(m3.mat) %in% test.m3),drop=F]
  train.m3.mat <- m3.mat[!(row.names(m3.mat) %in% test.m3),
    !(colnames(m3.mat) %in% test.m3),drop=F]
} else {
  test.m3.mat <- m3.mat[0,,drop=F]
  train.m3.mat <- m3.mat
}

# Scale columns of input matrices and remove columns with no variance
# in the training data
if(scale & center) {
  train.m1.mat <- scale(train.m1.mat)
  test.m1.mat <- scale(test.m1.mat, scale=attr(train.m1.mat, 'scaled:scale'),
                                    center=attr(train.m1.mat, 'scaled:center'))
  test.mat.mat <- test.m1.mat[,apply(is.na(train.m1.mat), 2, sum)==0]
  train.m1.mat <- train.m1.mat[,apply(is.na(train.m1.mat), 2, sum)==0]
  if(ncol(train.m1.mat)) {
    test.m1.mat <- test.m1.mat[,apply(train.m1.mat==Inf, 2, sum)==0]
    train.m1.mat <- train.m1.mat[,apply(train.m1.mat==Inf, 2, sum)==0]
  }

  train.m2.mat <- scale(train.m2.mat)
  test.m2.mat <- scale(test.m2.mat, scale=attr(train.m2.mat, 'scaled:scale'),
                                    center=attr(train.m2.mat, 'scaled:center'))
  test.mat.mat <- test.m2.mat[,apply(is.na(train.m2.mat), 2, sum)==0]
  train.m2.mat <- train.m2.mat[,apply(is.na(train.m2.mat), 2, sum)==0]
  if(ncol(train.m2.mat)) {
    test.m2.mat <- test.m2.mat[,apply(train.m2.mat==Inf, 2, sum)==0]
    train.m2.mat <- train.m2.mat[,apply(train.m2.mat==Inf, 2, sum)==0]
  }

  train.m3.mat <- scale(train.m3.mat)
  test.m3.mat <- scale(test.m3.mat, scale=attr(train.m3.mat, 'scaled:scale'),
                                    center=attr(train.m3.mat, 'scaled:center'))
  test.mat.mat <- test.m3.mat[,apply(is.na(train.m3.mat), 2, sum)==0]
  train.m3.mat <- train.m3.mat[,apply(is.na(train.m3.mat), 2, sum)==0]
  if(ncol(train.m3.mat)) {
    test.m3.mat <- test.m3.mat[,apply(train.m3.mat==Inf, 2, sum)==0]
    train.m3.mat <- train.m3.mat[,apply(train.m3.mat==Inf, 2, sum)==0]
  }

} else if(scale & !center) {
  train.m1.mat <- scale(train.m1.mat, center=F, scale=apply(train.m1.mat, 2, sd, na.rm=T))
  test.m1.mat <- scale(test.m1.mat, scale=attr(train.m1.mat, 'scaled:scale'), center=F)
  test.m1.mat <- test.m1.mat[,apply(is.na(train.m1.mat), 2, sum)==0]
  train.m1.mat <- train.m1.mat[,apply(is.na(train.m1.mat), 2, sum)==0]
  if(ncol(train.m1.mat)) {
    test.m1.mat <- test.m1.mat[,apply(train.m1.mat==Inf, 2, sum)==0]
    train.m1.mat <- train.m1.mat[,apply(train.m1.mat==Inf, 2, sum)==0]
  }

  train.m2.mat <- scale(train.m2.mat, center=F, scale=apply(train.m2.mat, 2, sd, na.rm=T))
  test.m2.mat <- scale(test.m2.mat, scale=attr(train.m2.mat, 'scaled:scale'), center=F)
  test.m2.mat <- test.m2.mat[,apply(is.na(train.m2.mat), 2, sum)==0]
  train.m2.mat <- train.m2.mat[,apply(is.na(train.m2.mat), 2, sum)==0]
  if(ncol(train.m2.mat)) {
    test.m2.mat <- test.m2.mat[,apply(train.m2.mat==Inf, 2, sum)==0]
    train.m2.mat <- train.m2.mat[,apply(train.m2.mat==Inf, 2, sum)==0]
  }

  train.m3.mat <- scale(train.m3.mat, center=F, scale=apply(train.m3.mat, 2, sd, na.rm=T))
  test.m3.mat <- scale(test.m3.mat, scale=attr(train.m3.mat, 'scaled:scale'), center=F)
  test.m3.mat <- test.m3.mat[,apply(is.na(train.m3.mat), 2, sum)==0]
  train.m3.mat <- train.m3.mat[,apply(is.na(train.m3.mat), 2, sum)==0]
  if(ncol(train.m3.mat)) {
    test.m3.mat <- test.m3.mat[,apply(train.m3.mat==Inf, 2, sum)==0]
    train.m3.mat <- train.m3.mat[,apply(train.m3.mat==Inf, 2, sum)==0]
  }
} else if(center & !scale) {
  train.m1.mat <- scale(train.m1.mat, scale=F)
  test.m1.mat <- scale(test.m1.mat, center=attr(train.m1.mat, 'scaled:center'),
    scale=F)

  train.m2.mat <- scale(train.m2.mat, scale=F)
  test.m2.mat <- scale(test.m2.mat, center=attr(train.m2.mat, 'scaled:center'), 
    scale=F)

  train.m3.mat <-  scale(train.m3.mat, scale=F)
  test.m3.mat <- scale(test.m3.mat, center=attr(train.m3.mat, 'scaled:center'), 
    scale=F)
}

# Transform predictors to kernel versions if 'kern=T'
if(kern) {
  kern.list <- kern_combine(train.m1.mat, train.m2.mat, train.m3.mat, 
    test.m1.mat, test.m2.mat, test.m3.mat, kern.combine,
    args)
  train.m1.mat <- kern.list[['train.m1.mat']]
  test.m1.mat  <- kern.list[['test.m1.mat']]
  train.m2.mat <- kern.list[['train.m2.mat']]
  test.m2.mat  <- kern.list[['test.m2.mat']]
  train.m3.mat <- kern.list[['train.m3.mat']]
  test.m3.mat  <- kern.list[['test.m3.mat']]
  if(kern.combine=='linear') {
    if(ncol(train.m1.mat)) {
      orig.train.m1.mat <- kern.list[['orig.train.m1.mat']]
      orig.test.m1.mat <- kern.list[['orig.test.m1.mat']]
    }
    if(ncol(train.m2.mat)) {
      orig.train.m2.mat <- kern.list[['orig.train.m2.mat']]
      orig.test.m2.mat <- kern.list[['orig.test.m2.mat']]
    }
    if(ncol(train.m3.mat)) {
      orig.train.m3.mat <- kern.list[['orig.train.m3.mat']]
      orig.test.m3.mat <- kern.list[['orig.test.m3.mat']]
    }
  }
  if(kern.scale) {
    if(ncol(train.m1.mat)) {
      m1.m <- mean(train.m1.mat)
      m1.s <- sd(train.m1.mat)
      train.m1.mat <- (train.m1.mat - m1.m)/m1.s
      test.m1.mat <- (test.m1.mat - m1.m)/m1.s
    }
    if(ncol(train.m2.mat)) {
      m2.m <- mean(train.m2.mat)
      m2.s <- sd(train.m2.mat)
      train.m2.mat <- (train.m2.mat - m2.m)/m2.s
      test.m2.mat <- (test.m2.mat - m2.m)/m2.s
    }
    if(ncol(train.m3.mat)) {
      m3.m <- mean(train.m3.mat)
      m3.s <- sd(train.m3.mat)
      train.m3.mat <- (train.m3.mat - m3.m)/m3.s
      test.m3.mat <- (test.m3.mat - m3.m)/m3.s
    }
  }
}

# Remove predictors that have no variance in the training data
# or have only one non-zero value
if(ncol(train.m1.mat)) {
  train.m1.mat <- train.m1.mat[,apply(train.m1.mat, 2, function(x) sd(x) > 0),drop=F]
  m1.one <- apply(train.m1.mat, 2, function(x) sum(x==1)==1) &
    apply(train.m1.mat, 2, function(x) sum(x==0)==(nrow(train.m1.mat)-1))
  train.m1.mat <- train.m1.mat[,!m1.one,drop=F]
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
  train.m3.mat <- train.m3.mat[,apply(train.m3.mat, 2, function(x) sd(x) > 0),drop=F]
  m3.one <- apply(train.m3.mat, 2, function(x) sum(x==1)==1) &
    apply(train.m3.mat, 2, function(x) sum(x==0)==(nrow(train.m3.mat)-1))
  train.m3.mat <- train.m3.mat[,!m3.one,drop=F]
  # Also remove these predictors from testing data
  test.m3.mat <- test.m3.mat[,match(dimnames(train.m3.mat)[[2]], 
                             dimnames(test.m3.mat)[[2]]),drop=F]
}

# If there are NA values in feature matrices throw an error
if(sum(is.na(train.m1.mat))) stop('NAs in mode 1 feature matrix!')
if(sum(is.na(train.m2.mat))) stop('NAs in mode 2 feature matrix!')
if(sum(is.na(train.m3.mat))) stop('NAs in mode 3 feature matrix!')

# Reorder responses to match the inputs
resp.tens <- resp.tens[match(rownames(m1.mat), dimnames(resp.tens)[[1]], nomatch=0),
                       match(rownames(m2.mat), dimnames(resp.tens)[[2]], nomatch=0),
                       match(rownames(m3.mat), dimnames(resp.tens)[[3]], nomatch=0), drop=F]

resp.train <- resp.tens[match(rownames(train.m1.mat), dimnames(resp.tens)[[1]], nomatch=0),
                        match(rownames(train.m2.mat), dimnames(resp.tens)[[2]], nomatch=0),
                        match(rownames(train.m3.mat), dimnames(resp.tens)[[3]], nomatch=0), drop=F]

# Normalize responses for the specific prediction task given 'normalize.for' argument
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
  resp.train <- sweep(resp.train, norm.dims, means, FUN='-')
  resp.train <- sweep(resp.train, norm.dims, sds, FUN='/')
} else if(normalize.for == 'mode123') {
    means <- mean(resp.train, na.rm=T)
    sds <- sd(resp.train, na.rm=T)
    resp.tens <- resp.tens - means
    resp.tens <- resp.tens / sds
    resp.train <- resp.train - means
    resp.train <- resp.train / sds
}

# Make the test response objects
################################
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
  # If there is no predictor data for the third mode, 
  # assume the third mode is dose and remove warm data across all doses
  if(m3.file == 'none') {
    non.missing <- which(apply(is.na(resp.train), c(1,2), sum)!=K)
    mask <- sample(non.missing, round(warm.per*length(non.missing)))
    for(k in 1:K) 
      resp.train[,,k][mask] <- NA
  } else {
    # otherwise remove randomly but the same for matching modes
    if(m1.file == m2.file & m1.file != 'none') {
      mask <- sample(I*K, (round(warm.per*I*K)))
      for(m in mask) {
        i <- m %% I + 1
        k <- (m-1) %/% I +1
        resp.train[i,i,k] <- NA
      }
    }
    if(m2.file == m3.file & m2.file != 'none') {
      mask <- sample(I*J, (round(warm.per*I*J)))
      for(m in mask) {
        i <- m %% I + 1
        j <- (m-1) %/% I +1
      resp.train[i,j,j] <- NA
      }
    }
  }
}

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

# If params$sigma2 is 'auto' set this value to the standard deviation of the training response data
if(params$sigma2=='auto')
  params$sigma2 <- sd(resp.train, na.rm=T)

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
  
  # Get cold results
  if(warm.per | ncol(m1.mat) | ncol(m2.mat) | ncol(m3.mat)) 
    test.results <- test_results(m=trained, d=train.data, test.results=test.results,
                                 warm.resp=warm.resp, test.m1=test.data.m1, 
                                 test.m2=test.data.m2, test.m3=test.data.m3,
                                 test.m1m2=test.data.m1m2, test.m1m3=test.data.m1m3,
                                 test.m2m3=test.data.m2m3, test.m1m2m3=test.data.m1m2m3)
  
  # Save every save.a.lot iterations
  if((i %% save.a.lot) ==0) save.image(paste0(out.dir, 'image.Rdata'))

  # Collect garbage (may help with memory)
  gc()
}

# Stop cluster (if using parallel package)
if(params$parallel) {
  if(.Platform$OS.type == "windows") stopCluster(clust)
}

save.image(paste0(out.dir, 'image.Rdata'))

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
# if they've been log transformed etc.
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
  for(i in 1:dim(test.data.m1m2$resp)[[1]]) for(j in 1:dim(test.data.m1m2$resp)[[2]])
    mean.tens.list[['m1m2']][i,j,] <- apply(train.data$resp, 3, mean, na.rm=T)
}
if(length(test.data.m1m3)) {
  mean.tens.list[['m1m3']] <- array(0, dim=dim(test.data.m1m3$resp))
  for(i in 1:dim(test.data.m1m3$resp)[[1]]) for(k in 1:dim(test.data.m1m3$resp)[[3]])
    mean.tens.list[['m1m3']][i,,k] <- apply(train.data$resp, 2, mean, na.rm=T)
}
if(length(test.data.m2m3)) {
  mean.tens.list[['m2m3']] <- array(0, dim=dim(test.data.m2m3$resp))
  for(j in 1:dim(test.data.m2m3$resp)[[2]]) for(k in 1:dim(test.data.m2m3$resp)[[3]])
    mean.tens.list[['m2m3']][,j,k] <- apply(train.data$resp, 1, mean, na.rm=T)
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
  abline(h=m3.mean.RMSE, col='yellow', lty=2, lwd=2)
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