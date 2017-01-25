#!/usr/bin/env Rscript

library(dplyr)
library(methods)
library(R6)
library(BaTFLED3D)

args <- commandArgs(TRUE)
run_prefix <- args[1]

# Read in a file to get ROCs
get_rocs <- function(f1) {
  load(f1)
  m1.roc <- plot_roc(toy$mode1.A, trained$mode1.A.mean)
  m2.roc <- plot_roc(toy$mode2.A, trained$mode2.A.mean)
  m3.roc <- plot_roc(toy$mode3.A, trained$mode3.A.mean)
  return(list(m1.roc, m2.roc, m3.roc))
}

n.files <- length(list.files(path = dirname(run_prefix),
  pattern = paste0(basename(run_prefix), '.[0-9]+.out')))

m1.rocs <- rep(0, n.files)
m2.rocs <- rep(0, n.files)
m3.rocs <- rep(0, n.files)

for(i in 1:n.files) {
  f <-  paste0(run_prefix, '.', i-1, '/image.Rdata')
  rocs <- get_rocs(f)
  m1.rocs[i] <- rocs[[1]]
  m2.rocs[i] <- rocs[[2]]
  m3.rocs[i] <- rocs[[3]]
}

print(sprintf('Mean ROC for mode 1 (sd): %.4f (%.4f)', mean(m1.rocs), sd(m1.rocs)))
print(sprintf('Mean ROC for mode 2 (sd): %.4f (%.4f)', mean(m2.rocs), sd(m2.rocs)))
print(sprintf('Mean ROC for mode 3 (sd): %.4f (%.4f)', mean(m3.rocs), sd(m3.rocs)))

print(paste(mean(m1.rocs), sd(m1.rocs),
            mean(m2.rocs), sd(m2.rocs),
            mean(m3.rocs), sd(m3.rocs), sep=','))