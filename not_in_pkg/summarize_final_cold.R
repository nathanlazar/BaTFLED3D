#!/usr/bin/env Rscript

library(dplyr)
library(methods)
library(R6)
library(BaTFLED3D)
# devtools::document()

args <- commandArgs(TRUE)
run_prefix <- args[1]

# Functions
############################################################

# Function to load in data from runs
loadData <- function(f, results){
  #loads an RData file, and returns a list with the 
  #trained model and warm & cold RMSE vectors.

  load(f)

  summaries['final', 'RMSE', fold] <- final.RMSE
  summaries['final', 'clip.RMSE', fold] <- final.RMSE.clip
  summaries['final', 'exp.var', fold] <- final.exp.var
  summaries['final', 'clip.exp.var', fold] <- final.exp.var.clip
  summaries['final', 'p.cor', fold] <- final.p.cor
  summaries['final', 'clip.p.cor', fold] <- final.p.cor.clip
  summaries['final', 's.cor', fold] <- final.s.cor
  summaries['final', 'clip.s.cor', fold] <- final.s.cor.clip

  summaries['mean', 'RMSE', fold] <- final.mean.RMSE
  summaries['mean', 'exp.var', fold] <- final.mean.exp.var
  summaries['mean', 'p.cor', fold] <- final.mean.p.cor
  summaries['mean', 's.cor', fold] <- final.mean.s.cor

  return(summaries)
}

########### MAIN ###############

# Determine the number of runs with this prefix
n.files <- length(list.files(path = dirname(run_prefix),
  pattern = paste0(basename(run_prefix), '.[0-9]+.out')))

summaries <- array(NA, dim=c(2, 8, n.files),
  dimnames=list(c('final', 'mean'),
                c('RMSE', 'clip.RMSE',
                  'exp.var', 'clip.exp.var',
                  'p.cor', 'clip.p.cor',
                  's.cor', 'clip.s.cor'),
                paste0('fold.', 1:n.files)))

for(fld in 1:n.files) {
  # Load in the run data
  f <- paste0(run_prefix, '.', (fld-1), '/final_test.Rdata')
  summaries <- loadData(f, results)
}

# Save all the data
# save.image(paste0(run_prefix, '_summary.Rdata'))

print("######################################################")
print('## Means ##')
print(apply(summaries, c(1,2), mean, na.rm=T))
print('## Standard deviations ##')
print(apply(summaries, c(1,2), sd, na.rm=T))

# print("######################################################")
# print('## Better than mean ##')
# print(better)


# # Make data frame counting how many folds peform better than the mean
# better <- matrix(NA, 7, 10, dimnames=list(dimnames(results$summaries)[[1]],
#   dimnames(results$summaries)[[2]][!grepl('iter', dimnames(results$summaries)[[2]])]))
# 
# for(type in c('A', 'H', 'train')) for(resp in c('RMSE', 'min.RMSE', 'clip.RMSE', 'min.clip.RMSE')) 
#   better[type, resp] <- sum(results$summaries[type, resp,] < results$mean['RMSE', 'train.m1',])
# for(resp in c('RMSE', 'min.RMSE', 'clip.RMSE', 'min.clip.RMSE')) 
#   better['warm', resp] <- sum(results$summaries['warm', resp,] < results$mean['RMSE', 'warm.m1',])
# for(type in c('m1', 'm2', 'm1m2')) for (resp in c('RMSE', 'min.RMSE', 'clip.RMSE', 'min.clip.RMSE')) 
#   better[type, resp] <- sum(results$summaries[type, resp,] < results$mean['RMSE', type,])
# for(type in c('A', 'H', 'train')) for(resp in c('exp.var', 'max.exp.var')) 
#   better[type, resp] <- sum(results$summaries[type, resp,] > results$mean['exp.var', 'train.m1',])
# for(resp in c('exp.var', 'max.exp.var')) 
#   better['warm', resp] <- sum(results$summaries['warm', resp,] > results$mean['exp.var', 'warm.m1',])
# for(type in c('m1', 'm2', 'm1m2')) for (resp in c('exp.var', 'max.exp.var')) 
#   better[type, resp] <- sum(results$summaries[type, resp,] > results$mean['exp.var', type,])
# for(type in c('A', 'H', 'train')) for(resp in c('p.cor', 'max.p.cor')) 
#   better[type, resp] <- sum(results$summaries[type, resp,] > results$mean['p.cor', 'train.m1',])
# for(resp in c('p.cor', 'max.p.cor')) 
#   better['warm', resp] <- sum(results$summaries['warm', resp,] > results$mean['p.cor', 'warm.m1',])
# for(type in c('m1', 'm2', 'm1m2')) for (resp in c('p.cor', 'max.p.cor')) 
#   better[type, resp] <- sum(results$summaries[type, resp,] > results$mean['p.cor', type,])
# for(type in c('A', 'H', 'train')) for(resp in c('s.cor', 'max.s.cor')) 
#   better[type, resp] <- sum(results$summaries[type, resp,] > results$mean['s.cor', 'train.m1',])
# for(resp in c('s.cor', 'max.s.cor')) 
#   better['warm', resp] <- sum(results$summaries['warm', resp,] > results$mean['s.cor', 'warm.m1',])
# for(type in c('m1', 'm2', 'm1m2')) for (resp in c('s.cor', 'max.s.cor')) 
#   better[type, resp] <- sum(results$summaries[type, resp,] > results$mean['s.cor', type,])
# 
# # Read log files to get run time if the runs finished
# if(length(system2('grep', c('"Job terminated"', paste0(run_prefix, '.*.log')), stdout=T)>0)) {
#   months <- c(31,28,31,30,31,30,31,31,30,31,30,31)
#   month.days <- cumsum(months)
#   log.starts <- system2('grep', c('"Job submitted"', paste0(run_prefix, '.*.log')), stdout=T)
#   log.ends <- system2('grep', c('"Job terminated"', paste0(run_prefix, '.*.log')), stdout=T)
#   month.start <- as.numeric(sapply(strsplit(sapply(strsplit(log.starts, split=' '), 
#     '[', 3), split='/'), '[', 1))
#   month.end <- as.numeric(sapply(strsplit(sapply(strsplit(log.ends, split=' '), 
#     '[', 3), split='/'), '[', 1))
#   day.start <- month.days[month.start] + 
#     as.numeric(sapply(strsplit(sapply(strsplit(log.starts, split=' '), '[', 3), split='/'), '[', 2))
#   day.end <- month.days[month.end] + 
#     as.numeric(sapply(strsplit(sapply(strsplit(log.ends, split=' '), '[', 3), split='/'), '[', 2))
#   days <- day.end - day.start
#   hour.start <- as.numeric(sapply(strsplit(sapply(strsplit(log.starts, split=' '), 
#     '[', 4), split=':'), '[', 1))
#   hour.end <- as.numeric(sapply(strsplit(sapply(strsplit(log.ends, split=' '), 
#     '[', 4), split=':'), '[', 1))
#   min.start <- as.numeric(sapply(strsplit(sapply(strsplit(log.starts, split=' '), 
#     '[', 4), split=':'), '[', 2))
#   min.end <- as.numeric(sapply(strsplit(sapply(strsplit(log.ends, split=' '), 
#     '[', 4), split=':'), '[', 2))
#   
#   hours <- days * 24 + (hour.end - hour.start) + (min.end - min.start)/60
#   rm(log.starts, log.ends)
#   print("Run time statistics (hours):")
#   summary(hours)
# }
# 
