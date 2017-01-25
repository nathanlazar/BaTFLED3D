# nathan dot lazar at gmail dot com

# Read in results from runs and collate the results from predicting the 
# mean response

# Usage: summarize_mean.R <run_summary.Rdata>

# Example: summarize_mean.R DAEMENresults/run_1731453_summary.Rdata

# Outputs results in a text file <run_prefix>_mean.txt

args <- commandArgs(TRUE)
summary_data <- args[1]

load(args[1])

sum.mean <- apply(results$mean, c(1,2), mean)
sum.sd <- apply(results$mean, c(1,2), sd)

print('Mean across replictes')
print('########################')
print(sum.mean)
print('Standard deviation across replicates')
print('########################')
print(sum.sd)

print('Predicting means for mode 1')
m1.vec <- c(NA, 18*4)
i <- 0
for(type in c('RMSE', 'exp.var', 'p.cor', 's.cor')) {
  for(resp in c('train.m1', 'warm.m1', 'm1', 'm2', 'm3', 
                'm1m2', 'm1m3', 'm2m3', 'm1m2m3')) {
    i <- i + 1
    if(resp %in% dimnames(sum.mean)[[2]])
      m1.vec[i] <- sum.mean[type, resp]
    i <- i + 1
    if(resp %in% dimnames(sum.sd)[[2]])
      m1.vec[i] <- sum.sd[type, resp]
  }
}
m1.vec[is.nan(m1.vec)] <- NA
print(paste(m1.vec, collapse=','))

print('Predicting means for mode 2')
m2.vec <- c(NA, 18*4)
i <- 0
for(type in c('RMSE', 'exp.var', 'p.cor', 's.cor')) {
  for(resp in c('train.m2', 'warm.m2', 'm1', 'm2', 'm3', 
                'm1m2', 'm1m3', 'm2m3', 'm1m2m3')) {
    i <- i + 1
    if(resp %in% dimnames(sum.mean)[[2]])
      m2.vec[i] <- sum.mean[type, resp]
    i <- i + 1
    if(resp %in% dimnames(sum.sd)[[2]])
      m2.vec[i] <- sum.sd[type, resp]
  }
}
m2.vec[is.nan(m2.vec)] <- NA
print(paste(m2.vec, collapse=','))

print('Predicting means for mode 3')
m3.vec <- c(NA, 18*4)
i <- 0
for(type in c('RMSE', 'exp.var', 'p.cor', 's.cor')) {
  for(resp in c('train.m3', 'warm.m3', 'm1', 'm2', 'm3', 
                'm1m2', 'm1m3', 'm2m3', 'm1m2m3')) {
    i <- i + 1
    if(resp %in% dimnames(sum.mean)[[2]])
      m3.vec[i] <- sum.mean[type, resp]
    i <- i + 1
    if(resp %in% dimnames(sum.sd)[[2]])
      m3.vec[i] <- sum.sd[type, resp]
  }
}
m3.vec[is.nan(m3.vec)] <- NA
print(paste(m3.vec, collapse=','))

print('Predicting means for modes 1&2')
m12.vec <- c(NA, 18*4)
i <- 0
for(type in c('RMSE', 'exp.var', 'p.cor', 's.cor')) {
  for(resp in c('train.m1m2', 'warm.m1m2', 'm1', 'm2', 'm3', 
                'm1m2', 'm1m3', 'm2m3', 'm1m2m3')) {
    i <- i + 1
    if(resp %in% dimnames(sum.mean)[[2]])
      m12.vec[i] <- sum.mean[type, resp]
    i <- i + 1
    if(resp %in% dimnames(sum.sd)[[2]])
      m12.vec[i] <- sum.sd[type, resp]
  }
}
m12.vec[is.nan(m12.vec)] <- NA
print(paste(m12.vec, collapse=','))

print('Predicting means for modes 1&3')
m13.vec <- c(NA, 18*4)
i <- 0
for(type in c('RMSE', 'exp.var', 'p.cor', 's.cor')) {
  for(resp in c('train.m1m3', 'warm.m1m3', 'm1', 'm2', 'm3', 
                'm1m2', 'm1m3', 'm2m3', 'm1m2m3')) {
    i <- i + 1
    if(resp %in% dimnames(sum.mean)[[2]])
      m13.vec[i] <- sum.mean[type, resp]
    i <- i + 1
    if(resp %in% dimnames(sum.sd)[[2]])
      m13.vec[i] <- sum.sd[type, resp]
  }
}
m13.vec[is.nan(m13.vec)] <- NA
print(paste(m13.vec, collapse=','))

print('Predicting means for modes 2&3')
m23.vec <- c(NA, 18*4)
i <- 0
for(type in c('RMSE', 'exp.var', 'p.cor', 's.cor')) {
  for(resp in c('train.m2m3', 'warm.m2m3', 'm1', 'm2', 'm3', 
                'm1m2', 'm1m3', 'm2m3', 'm1m2m3')) {
    i <- i + 1
    if(resp %in% dimnames(sum.mean)[[2]])
      m23.vec[i] <- sum.mean[type, resp]
    i <- i + 1
    if(resp %in% dimnames(sum.sd)[[2]])
      m23.vec[i] <- sum.sd[type, resp]
  }
}
m23.vec[is.nan(m23.vec)] <- NA
print(paste(m23.vec, collapse=','))

print('Predicting means for modes 1,2&3')
m123.vec <- c(NA, 18*4)
i <- 0
for(type in c('RMSE', 'exp.var', 'p.cor', 's.cor')) {
  for(resp in c('train.m1m2m3', 'warm.m1m2m3', 'm1', 'm2', 'm3', 
                'm1m2', 'm1m3', 'm2m3', 'm1m2m3')) {
    i <- i + 1
    if(resp %in% dimnames(sum.mean)[[2]])
      m123.vec[i] <- sum.mean[type, resp]
    i <- i + 1
    if(resp %in% dimnames(sum.sd)[[2]])
      m123.vec[i] <- sum.sd[type, resp]
  }
}
m123.vec[is.nan(m123.vec)] <- NA
print(paste(m123.vec, collapse=','))

