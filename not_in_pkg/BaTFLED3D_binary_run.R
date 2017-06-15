#/usr/bin/env Rscript

# nathan dot lazar at gmail dot com

# Usage: BaTFLED3D_binary_run.R <mode1 data> <mode2 data> <mode3 data>
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
# library(pryr)           # For monitoring memory usage
# library(microbenchmark) # For monitoring run time

print('Packages loaded')
args <- commandArgs(TRUE)

# HEISER data predicting for cell lines
test <- F
if(test) 
  args <- list('../HeiserData/Cell_Line_Data/exp_cna_meth_mat.Rdata',
    '../HeiserData/Drug_Data/dr_all_mat.Rdata', 
    'none', 
    '../HeiserData/Responses/Split/cl_train_param_trans_tens.Rdata',
    '../HeiserData/Responses/Split/cl_loocv_mat.Rdata', 
    'none', 'none',
    'HEISERresults/test/', 'decomp=Tucker', 'row.share=T',
    'reps=5', 'warm.per=0', 'plot=F', 'cores=1',
    'A1.intercept=T', 'A2.intercept=T', 'A3.intercept=F', 
    'H1.intercept=T', 'H2.intercept=T', 'H3.intercept=T',  
    'm1.sigma2=0.01','m2.sigma2=0.01','m3.sigma2=1',
    'R1=4', 'R2=4', 'R3=4', 'normalize.for=mode12', 
    'kern=T', 'kern.combine=prod', 
    'anno.s=2', 'exp.s=2', 'cna.s=2', 'meth.s=2',
    'scale=F', 'fold=1')

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
fold <- 1
scale <- T
kern <- F
kern.combine <- 'prod'

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
if(sum(grepl('^kern=', args)))
  kern <- as.logical(sub('kern=', '', args[grepl('^kern=', args)]))
if(sum(grepl('^kern.combine=', args)))
  kern.combine <- sub('kern.combine=', '', args[grepl('^kern.combine=', args)])

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

# Subset samples to those that have responses
m1.mat <- m1.mat[rownames(m1.mat) %in% dimnames(resp.tens)[[1]],,drop=F]
m2.mat <- m2.mat[rownames(m2.mat) %in% dimnames(resp.tens)[[2]],,drop=F]
m3.mat <- m3.mat[rownames(m3.mat) %in% dimnames(resp.tens)[[3]],,drop=F]

# Set NA values to 0 in feature matrices
m1.mat[is.na(m1.mat)] <- 0
m2.mat[is.na(m2.mat)] <- 0
m3.mat[is.na(m3.mat)] <- 0

# Scale responses if multiplier non-zero
if(multiplier != 0) {
  resp.tens <- resp.tens * multiplier
}

# Remove the samples for this cross validation fold to be test
# data for this run.
###############################################################
if(exists('m1.cv.fold')) {
  test.m1 <- m1.cv.fold[fold,]
  test.m1.mat <- m1.mat[row.names(m1.mat) %in% test.m1,,drop=F]
  train.m1.mat <- m1.mat[!(row.names(m1.mat) %in% test.m1),,drop=F]
} else {
  test.m1.mat <- m1.mat[0,,drop=F]
  train.m1.mat <- m1.mat
}

if(exists('m2.cv.fold')) {
  test.m2 <- m2.cv.fold[fold,]
  test.m2.mat <- m2.mat[row.names(m2.mat) %in% test.m2,,drop=F]
  train.m2.mat <- m2.mat[!(row.names(m2.mat) %in% test.m2),,drop=F]
} else {
  test.m2.mat <- m2.mat[0,,drop=F]
  train.m2.mat <- m2.mat
}

if(exists('m3.cv.fold')) {
  test.m3 <- m3.cv.fold[fold,]
  test.m3.mat <- m3.mat[row.names(m3.mat) %in% test.m3,,drop=F]
  train.m3.mat <- m3.mat[!(row.names(m3.mat) %in% test.m3),,drop=F]
} else {
  test.m3.mat <- m3.mat[0,,drop=F]
  train.m3.mat <- m3.mat
}

# Scale non-binary columns of input matrices
if(scale) {
  m1.bin <- apply(m1.mat, 2, min)==0 & apply(m1.mat, 2, max)==1 &
            apply(m1.mat, 2, function(x) length(unique(x)))==2
  scld <- scale(train.m1.mat[,!m1.bin,drop=F])
  train.m1.mat[,!m1.bin] <- scld
  test.m1.mat[,!m1.bin] <- scale(test.m1.mat[,!m1.bin,drop=F], 
    scale=attr(scld, 'scaled:scale'), center=attr(scld, 'scaled:center')) 

  m2.bin <- apply(m2.mat, 2, min)==0 & apply(m2.mat, 2, max)==1 &
            apply(m2.mat, 2, function(x) length(unique(x)))==2
  scld <- scale(train.m2.mat[,!m2.bin,drop=F])
  train.m2.mat[,!m2.bin] <- scld
  test.m2.mat[,!m2.bin] <- scale(test.m2.mat[,!m2.bin,drop=F],
    scale=attr(scld, 'scaled:scale'), center=attr(scld, 'scaled:center'))

  m3.bin <- apply(m3.mat, 2, min)==0 & apply(m3.mat, 2, max)==1 &
            apply(m3.mat, 2, function(x) length(unique(x)))==2
  scld <- scale(train.m3.mat[,!m3.bin,drop=F])
  train.m3.mat[,!m3.bin] <- scld
  test.m3.mat[,!m3.bin] <- scale(test.m3.mat[,!m3.bin,drop=F],
    scale=attr(scld, 'scaled:scale'), center=attr(scld, 'scaled:center'))
}

# Transform predictors to kernel versions if 'kern=T'
if(kern) {
  if(ncol(m1.mat)) {
    m1.types <- unique(sapply(strsplit(colnames(m1.mat), split='.', fixed=T), '[', 1))
    # Get the kernel width for each type from args
    m1.s <- rep(1, length(m1.types))
    names(m1.s) <- m1.types
    for(type in m1.types)
      if(sum(grepl(type, args[9:length(args)])))
        m1.s[[type]] <- as.numeric(
          strsplit(grep(type, args[9:length(args)], value=T), split='=')[[1]][2])

    m1.train.list <- list()
    if(exists('test.m1'))
      m1.test.list <- list()
    for(i in 1:length(m1.s)) {
      ind <- grep(m1.types[[i]], colnames(train.m1.mat))
      m1.train.list[[i]] <- kernelize(train.m1.mat[,ind,drop=F], 
        train.m1.mat[,ind,drop=F], m1.s[[i]])
      colnames(m1.train.list[[i]]) <- paste(m1.types[[i]], 
        colnames(m1.train.list[[i]]), sep='.')
      if(exists('test.m1')) {
        m1.test.list[[i]] <- kernelize(test.m1.mat[,ind,drop=F], 
          train.m1.mat[,ind,drop=F], m1.s[[i]])
        colnames(m1.test.list[[i]]) <- paste(m1.types[[i]], 
          colnames(m1.test.list[[i]]), sep='.')
      }
    }

    if(kern.combine == 'cat') {
    ########### CHANGE THIS? #############
    # Concatenate the unkernelized annotation data
    # and the other two kernels 
    ######################################
      train.m1.mat <- cbind(cbind(
        train.m1.mat[,grep('anno.', colnames(train.m1.mat)), drop=F],
                           m1.train.list[[2]]), m1.train.list[[3]])
      if(nrow(test.m1.mat))
        test.m1.mat <- cbind(cbind(
          test.m1.mat[,grep('anno.', colnames(test.m1.mat)), drop=F],
          m1.test.list[[2]]), m1.test.list[[3]])
    }
    if(kern.combine == 'mean') {
    # Take the average of the kernels
      train.m1.mat <- m1.train.list[[1]]
      if(nrow(test.m1.mat)) test.m1.mat <- m1.test.list[[1]]
      if(length(m1.train.list) > 1) {
        for(n in 2:length(m1.train.list)) {
          train.m1.mat <- train.m1.mat + m1.train.list[[n]]
          if(nrow(test.m1.mat))
            test.m1.mat <- test.m1.mat + m1.test.list[[n]]
        }
        train.m1.mat <- train.m1.mat/length(m1.train.list)
        if(nrow(test.m1.mat))  test.m1.mat <- test.m1.mat/length(m1.test.list)
      }
    }
    if(kern.combine == 'prod') {
    # Take the product of the kernels
      train.m1.mat <- m1.train.list[[1]]
      if(nrow(test.m1.mat)) test.m1.mat <- m1.test.list[[1]]
      if(length(m1.train.list) > 1) {
        for(n in 2:length(m1.train.list)) {
          train.m1.mat <- train.m1.mat * m1.train.list[[n]]
          if(nrow(test.m1.mat)) test.m1.mat <- test.m1.mat * m1.test.list[[n]]
        }
      }
    }
    if(!nrow(test.m1.mat)) test.m1.mat <- train.m1.mat[0,]
  }

  if(ncol(m2.mat)) {
    m2.types <- unique(sapply(strsplit(colnames(m2.mat), split='.', fixed=T), '[', 1))
    m2.s <- rep(1, length(m2.types))
    names(m2.s) <- m2.types
    for(type in m2.types)
      if(sum(grepl(type, args[9:length(args)])))
        m2.s[[type]] <- as.numeric(
          strsplit(grep(type, args[9:length(args)], value=T), split='=')[[1]][2])
    m2.train.list <- list()
    if(nrow(test.m2.mat)) m2.test.list <- list()
    for(i in 1:length(m2.s)) {
      ind <- grep(m2.types[[i]], colnames(train.m2.mat))
      m2.train.list[[i]] <- kernelize(train.m2.mat[,ind,drop=F],
        train.m2.mat[,ind,drop=F], m2.s[[i]])
      colnames(m2.train.list[[i]]) <- paste(m2.types[[i]],
        colnames(m2.train.list[[i]]), sep='.')
      if(nrow(test.m2.mat)) {
        m2.test.list[[i]] <- kernelize(test.m2.mat[,ind,drop=F],
          train.m2.mat[,ind,drop=F], m2.s[[i]])
        colnames(m2.test.list[[i]]) <- paste(m2.types[[i]],
          colnames(m2.test.list[[i]]), sep='.')
      }
    }

    if(kern.combine == 'cat') {
    ########### CHANGE THIS? #############
    # Concatenate the unkernelized target and class data
    # and the other two kernels
    #######################################
      train.m2.mat <- cbind(cbind(
        train.m2.mat[,grep('targ.|class.', colnames(train.m2.mat)), drop=F],
                           m2.train.list[[3]]), m2.train.list[[4]])
      if(nrow(test.m2.mat))
        test.m2.mat <- cbind(cbind(
          test.m2.mat[,grep('targ.|class.', colnames(test.m2.mat)), drop=F],
          m2.test.list[[3]]), m2.test.list[[4]])
    }
    if(kern.combine == 'mean') {
    # Take the average of the kernels
      train.m2.mat <- m2.train.list[[1]]
      if(nrow(test.m2.mat)) test.m2.mat <- m2.test.list[[1]]
      if(length(m2.train.list) > 1) {
        for(n in 2:length(m2.train.list)) {
          train.m2.mat <- train.m2.mat + m2.train.list[[n]]
          if(nrow(test.m2.mat)) test.m2.mat <- test.m2.mat + m2.test.list[[n]]
        }
        train.m2.mat <- train.m2.mat/length(m2.train.list)
        if(nrow(test.m2.mat)) test.m2.mat <- test.m2.mat/length(m2.test.list)
      }
    }
    if(kern.combine == 'prod') {
    # Take the product of the kernels
      train.m2.mat <- m2.train.list[[1]]
      if(nrow(test.m2.mat)) test.m2.mat <- m2.test.list[[1]]
      if(length(m2.train.list) > 1) {
        for(n in 2:length(m2.train.list)) {
          train.m2.mat <- train.m2.mat * m2.train.list[[n]]
          if(nrow(test.m2.mat)) test.m2.mat <- test.m2.mat * m2.test.list[[n]]
        }
      }
    }
    if(!nrow(test.m2.mat)) test.m2.mat <- train.m2.mat[0,]
  }

  if(ncol(m3.mat)) {
    m3.types <- unique(sapply(strsplit(colnames(m3.mat), split='.', fixed=T), '[', 1))
    m3.s <- rep(1, length(m3.types))
    names(m3.s) <- m3.types
    for(type in m3.types)
      if(sum(grepl(type, args[9:length(args)])))
        m3.s[[type]] <- as.numeric(strsplit(
          grep(type, args[9:length(args)], value=T), split='=')[[1]][2])

    m3.train.list <- list()
    if(nrow(test.m3.mat)) m3.test.list <- list()
    for(i in 1:length(m3.s)) {
      ind <- grep(m3.types[[i]], colnames(train.m3.mat))
      m3.train.list[[i]] <- kernelize(train.m3.mat[,ind,drop=F],
        train.m3.mat[,ind,drop=F], m3.s[[i]])
      colnames(m3.train.list[[i]]) <- paste(m3.types[[i]],
        colnames(m3.train.list[[i]]), sep='.')
      if(nrow(test.m3.mat)) {
        m3.test.list[[i]] <- kernelize(test.m3.mat[,ind,drop=F],
          train.m3.mat[,ind,drop=F], m3.s[[i]])
        colnames(m3.test.list[[i]]) <- paste(m3.types[[i]],
          colnames(m3.test.list[[i]]), sep='.')
      }
    }

    if(kern.combine == 'cat') {
    # Concatenate the kernels
      train.m3.mat <- m3.train.list[[1]]
      if(nrow(test.m3.mat)) test.m3.mat <- m3.test.list[[1]]
      if(length(m3.train.list) > 1) {
        for(n in 2:length(m3.train.list)) {
          train.m3.mat <- cbind(train.m3.mat, m3.train.list[[n]])
          if(nrow(test.m3.mat)) test.m3.mat <- cbind(test.m3.mat, m3.test.list[[n]])
        }
      }
    }

    if(kern.combine == 'mean') {
    # Take the average of the kernels
      train.m3.mat <- m3.train.list[[1]]
      if(nrow(test.m3.mat)) test.m3.mat <- m3.test.list[[1]]
      if(length(m3.train.list) > 1) {
        for(n in 2:length(m3.train.list)) {
          train.m3.mat <- train.m3.mat + m3.train.list[[n]]
          if(nrow(test.m3.mat)) test.m3.mat <- test.m3.mat + m3.test.list[[n]]
        }
        if(nrow(test.m3.mat)) train.m3.mat <- train.m3.mat/length(m3.train.list)
        test.m3.mat <- test.m3.mat/length(m3.test.list)
      }
    }
    if(kern.combine == 'prod') {
    # Take the product of the kernels
      train.m3.mat <- m3.train.list[[1]]
      if(nrow(test.m3.mat)) test.m3.mat <- m3.test.list[[1]]
      if(length(m3.train.list) > 1) {
        for(n in 2:length(m3.train.list)) {
          train.m3.mat <- train.m3.mat * m3.train.list[[n]]
          if(nrow(test.m3.mat)) test.m3.mat <- test.m3.mat * m3.test.list[[n]]
        }
      }
    }
    if(!nrow(test.m3.mat)) test.m3.mat <- train.m3.mat[0,]
  }
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
  train.m3.mat <- train.m3.mat[,apply(train.m3.mat, 2, function(x) sd(x) > 0),drop=F]
  m3.one <- apply(train.m3.mat, 2, function(x) sum(x==1)==1) &
    apply(train.m3.mat, 2, function(x) sum(x==0)==(nrow(train.m3.mat)-1))
  train.m3.mat <- train.m3.mat[,!m3.one,drop=F]
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

# If params$sigma2 is 'auto' set this value to the variance of the training response data
if(params$sigma2=='auto')
  params$sigma2 <- var(resp.train, na.rm=T)

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
  if(ncol(m1.mat) | ncol(m2.mat) | ncol(m3.mat)) 
    test.results <- test_results(m=trained, d=train.data, test.results=test.results,
                                 warm.resp=warm.resp, test.m1=test.data.m1, 
                                 test.m2=test.data.m2, test.m3=test.data.m3,
                                 test.m1m2=test.data.m1m2, test.m1m3=test.data.m1m3,
                                 test.m2m3=test.data.m2m3, test.m1m2m3=test.data.m1m2m3)
  
  # Save every 10 iterations
  # if((i %% 10) ==0) save.image(paste0(out.dir, 'image.Rdata'))
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