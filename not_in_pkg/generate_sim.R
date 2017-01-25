#!/usr/bin/env Rscript

# nathan dot lazar at gmail dot com

# Generates a simulated dataset and makes a file for separating training 
# examples into 10 cross-validation folds

# Usage: generate_sim.R <output dir> <options...>

# Example: generate_sim.R test decomp=Tucker \
#  row.share=T warm.per=0.01 plot=F \
#  m1.rows=10 m1.cols=100 \
#  m2.rows=15 m2.cols=150 \
#  m3.rows=8 m3.cols=0 \
#  R1=3 R2=3 R3=3 \
#  m1.true=10 m2.true=15 m3.true=0 \
#  A1.intercept=T A2.intercept=T A3.intercept=F \
#  H1.intercept=T H2.intercept=T H3.intercept=T \
#  core.spar=.5 noise.sd=0.1

# args <- list('test' , 'decomp=Tucker',
#              'row.share=T', 'warm.per=0.01', 'plot=F',
#              'm1.rows=10', 'm1.cols=100',
#              'm2.rows=15', 'm2.cols=150',
#              'm3.rows=8', 'm3.cols=0',
#              'R1=3', 'R2=3', 'R3=3',
#              'm1.true=10', 'm2.true=15', 'm3.true=0',
#              'A1.intercept=T', 'A2.intercept=T', 'A3.intercept=F',
#              'H1.intercept=T', 'H2.intercept=T', 'H3.intercept=T',
#              'core.spar=.5', 'noise.sd=0.1')

pkgTest <- function(x) { # load packages and install if needed
  if (!require(x,character.only = TRUE)) {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

pkgTest('foreach')      # For parallel processing in loops
pkgTest('R6')           # For making memory efficent R6 objects
pkgTest('rTensor')      # Needed for multiplying matrices into tensors (could be removed)
pkgTest('dplyr')        # General data frame manipulation
pkgTest('ggplot2')      # Pretty plotting
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
warm.per <- 0.05
plot <- F

# Override defaults if provided
if(sum(grepl('^warm.per=', args)))
  warm.per <- as.numeric(sub('warm.per=', '', args[grepl('^warm.per=', args)]))
if(sum(grepl('^plot=', args)))
  plot <- as.logical(sub('plot=', '', args[grepl('^plot=', args)]))

## Generating a simulated dataset
params <- get_data_params(args)
toy <- mk_toy(params)

# Make matrices for 10-fold cross validation
m1.per <- ceiling(params$m1.rows/10)
m2.per <- ceiling(params$m2.rows/10)
m3.per <- ceiling(params$m3.rows/10)
  
m1.cv.fold <- matrix(row.names(toy$mode1.X), 10, m1.per)
m2.cv.fold <- matrix(row.names(toy$mode2.X), 10, m2.per)
m3.cv.fold <- matrix(row.names(toy$mode3.X), 10, m3.per)

# Save everything
save.image(paste0(out.dir, 'image.Rdata'))

# Save input matrices, response tensor and cv fold matrices
###########################################################
if(params$m1.cols) {
  m1.X <- toy$mode1.X
  save(m1.X, file=paste0(out.dir, 'mode1.X.Rdata'))
  save(m1.cv.fold, file=paste0(out.dir, 'cv_folds_m1.Rdata'))
}
if(params$m2.cols) {
  m2.X <- toy$mode2.X
  save(m2.X, file=paste0(out.dir, 'mode2.X.Rdata'))
  save(m2.cv.fold, file=paste0(out.dir, 'cv_folds_m2.Rdata'))
}
if(params$m3.cols) {
  m3.X <- toy$mode3.X
  save(m3.X, file=paste0(out.dir, 'mode3.X.Rdata'))
  save(m3.cv.fold, file=paste0(out.dir, 'cv_folds_m3.Rdata'))  
}

train.resp <- toy$resp
save(train.resp, file=paste0(out.dir, 'resp.Rdata'))

