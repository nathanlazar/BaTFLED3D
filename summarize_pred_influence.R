# Reads in data from multiple runs and builds a matrices of the 
# influence of each predictor on the responses.

# Usage: Rscript summarize_pred_influence.R <run prefix>

# Example: Rscript summarize_pred_influence.R CTRP2results/run_1188453

.libPaths("/home/users/lazar/R/x86_64-redhat-linux-gnu-library/3.2")

library(dplyr)
library(methods)
library(R6)
library(rTensor)
library(BaTFLED3D)

args <- commandArgs(TRUE)
run_prefix <- args[1]

##############################################################################
# functions
load_res <- function(file) {
  load(file)
  return(list(m=trained, d=train.data))
}

##############################################################################
# Main

n.files <- length(list.files(path = dirname(run_prefix),
  pattern = paste0(basename(run_prefix), '.[0-9]+.out')))

files <- paste0(run_prefix, '.', 0:(n.files-1), '/image.Rdata')

pred.mats <- list()

for(i in 1:length(files)) {
  file <- files[i]
  res <- load_res(file)
  
  new.mats <- list(m1.preds=matrix(NA, nrow(res$m$mode1.A.mean), 1),
                   m2.preds=matrix(NA, nrow(res$m$mode2.A.mean), 1),
                   m3.preds=matrix(NA, nrow(res$m$mode3.A.mean), 1))
  dimnames(new.mats$m1.preds) <- list(rownames(res$m$mode1.A.mean), paste0('run.', i))
  dimnames(new.mats$m2.preds) <- list(rownames(res$m$mode2.A.mean), paste0('run.', i))
  dimnames(new.mats$m3.preds) <- list(rownames(res$m$mode3.A.mean), paste0('run.', i))

  infl <- get_influence(res$m, res$d, method='sub', interactions=F)
  new.mats$m1.preds[] <- infl$m1.inf
  new.mats$m2.preds[] <- infl$m2.inf
  new.mats$m3.preds[] <- infl$m3.inf
  
  if(length(pred.mats)) {
    pred.mats$m1.preds <- merge(pred.mats$m1.preds, new.mats$m1.preds, by=0, all=T)
    rownames(pred.mats$m1.preds) <- pred.mats$m1.preds$Row.names
    pred.mats$m1.preds <- pred.mats$m1.preds[,-1]
    pred.mats$m2.preds <- merge(pred.mats$m2.preds, new.mats$m2.preds, by=0, all=T)  
    rownames(pred.mats$m2.preds) <- pred.mats$m2.preds$Row.names
    pred.mats$m2.preds <- pred.mats$m2.preds[,-1]
    pred.mats$m3.preds <- merge(pred.mats$m3.preds, new.mats$m3.preds, by=0, all=T)  
    rownames(pred.mats$m3.preds) <- pred.mats$m3.preds$Row.names
    pred.mats$m3.preds <- pred.mats$m3.preds[,-1]
  } else pred.mats <- new.mats
}

save(pred.mats, file=paste0(run_prefix, '_pred_infl.Rdata'))
