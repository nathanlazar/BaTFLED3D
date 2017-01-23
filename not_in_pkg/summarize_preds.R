# Reads in data from multiple runs builds a matrix of the weights on predictors

# Usage: summarize_preds.R <run_prefix>
# Output: <run_prefix>_preds.Rdata

# Resulting Rdata object stored in <run_prefix>_preds.Rdata

library(dplyr)
library(methods)
library(R6)
library(BaTFLED3D)
#devtools::document()

args <- commandArgs(TRUE)
run_prefix <- args[1]

############################################################################
# Functions

get_pred_weights <- function(file, pred.dfs=list()) {
  # Gets the mean of the absolute weights across rows of projection matrices
  # and adds a column to the data frames in pred.dfs
  # pred.dfs should be a list of data.frames m1.preds, m2.preds, m3.preds.

  load(file)

  if(length(pred.dfs)==0) {
    pred.dfs <- list(m1.preds=data.frame(row.names=rownames(trained$mode1.A.mean)),
                     m2.preds=data.frame(row.names=rownames(trained$mode2.A.mean)),
                     m3.preds=data.frame(row.names=rownames(trained$mode3.A.mean)))
  }

  m1.preds <- data.frame(apply(abs(trained$mode1.A.mean), 1, mean, na.rm=T))
  m2.preds <- data.frame(apply(abs(trained$mode2.A.mean), 1, mean, na.rm=T))
  m3.preds <- data.frame(apply(abs(trained$mode3.A.mean), 1, mean, na.rm=T))

  pred.dfs$m1.preds <- merge(pred.dfs$m1.preds, m1.preds, by=0, all=T)
  if(ncol(pred.dfs$m1.preds)) {
    rownames(pred.dfs$m1.preds) <- pred.dfs$m1.preds[,1]
    pred.dfs$m1.preds <- pred.dfs$m1.preds[,-1,drop=F]
    names(pred.dfs$m1.preds) <- paste0('run', 1:ncol(pred.dfs$m1.preds))
  }

  pred.dfs$m2.preds <- merge(pred.dfs$m2.preds, m2.preds, by=0, all=T)
  if(ncol(pred.dfs$m2.preds)) {
    rownames(pred.dfs$m2.preds) <- pred.dfs$m2.preds[,1]
    pred.dfs$m2.preds <- pred.dfs$m2.preds[,-1,drop=F]
    names(pred.dfs$m2.preds) <- paste0('run', 1:ncol(pred.dfs$m2.preds))
  }

  pred.dfs$m3.preds <- merge(pred.dfs$m3.preds, m3.preds, by=0, all=T)
  if(ncol(pred.dfs$m3.preds)) {
    rownames(pred.dfs$m3.preds) <- pred.dfs$m3.preds[,1]
    pred.dfs$m3.preds <- pred.dfs$m3.preds[,-1,drop=F]
    names(pred.dfs$m3.preds) <- paste0('run', 1:ncol(pred.dfs$m3.preds))
  }

  return(pred.dfs)
}

##############################################################################
# Main

n.files <- length(list.files(path = dirname(run_prefix),
  pattern = paste0(basename(run_prefix), '.[0-9]+.out')))

files <- paste0(run_prefix, '.', 0:(n.files-1), '/image.Rdata')

pred.dfs <- list()

for(file in files) {
  pred.dfs <- get_pred_weights(file, pred.dfs)
}

save(pred.dfs, file=paste0(run_prefix, '_preds.Rdata'))
