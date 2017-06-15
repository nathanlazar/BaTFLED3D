# nathan dot lazar at gmail dot com

# Subset predictors by fitting a mixture of two Gaussians to the density of log 
# predictor weights and keeping the Gaussian with the higher mean

# First argument is Rdata object created with summarize_preds.R 
# or summarize_pred_influence.R 

# Usage: subset_preds_fit_mixture.R <*_preds.Rdata> <m1_matrix.Rdata> <m2_matrix.Rdata>
#                       <m3_matrix.Rdata> <number of folds selecting to keep predictor> 
#                       <out directory>

# The output directory is filled with the input data matrices subsetted
# using both the intersection method with the given ovelap value 
# and the union method for the top 5%, 10%, 15%, 20% and 25% of predictors.

# args <- list('CTRP2results/run_225280_pred_infl.Rdata',
#              '../CTRP2/CTRP2_cl_train.Rdata',
#              '../CTRP2/CTRP2_dr_train.Rdata',
#              'none', '3', '../CTRP2/test')

# args <- list('DREAM7results/run_1047559_pred_infl.Rdata',
#              '../DREAM7/DREAM7_cl.Rdata',
#              '../DREAM7/DREAM7_dr_nodt.Rdata',
#              'none', '3', 'DREAM7results/run_1047559_selected_preds')

library(dplyr)
library(methods)
library(R6)
library(mixtools)
# library(BaTFLED3D)
# devtools::document()

########################################################
# Functions
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#######################################################
# Main
args <- commandArgs(TRUE)
print(args)

plot <- F

pred.mats <- loadRData(args[[1]])
folds <- ncol(pred.mats[[1]])

cuts <- list()
if(args[[2]] != 'none') {
  m1.mat <- loadRData(args[[2]])
  cuts[['mode1']] <- rep(0, folds)
} else m1.mat <- matrix(NA, 0, 0)
if(args[[3]] != 'none') {
  m2.mat <- loadRData(args[[3]])
  cuts[['mode2']] <- rep(0, folds)
} else m2.mat <- matrix(NA, 0, 0)
if(args[[4]] != 'none') {
  m3.mat <- loadRData(args[[4]])
  cuts[['mode3']] <- rep(0, folds)
} else m3.mat <- matrix(NA, 0, 0)

n_lap <- as.numeric(args[[5]])
out.dir <- args[[6]]

# If the output directory doesn't exist, create it
print(paste('Output directory:', out.dir))
if (!file.exists(out.dir)){
  dir.create(file.path(out.dir))
}

for(fold in 1:folds) {
  for(mode in 1:3) if(paste0('mode', mode) %in% names(cuts)) {
    n_gau <- 0
    iter <- 0
    while(n_gau !=2) { # Repeat if the number of Gaussians found is not 2
      iter <- iter + 1
      pred.dens <- log(pred.mats[[mode]][,fold])
      pred.dens <- pred.dens[!is.na(pred.dens)]
      # Check convergence?
      mixmdl = normalmixEM(pred.dens, maxit=10000)
      if(plot) {
        plot(mixmdl, which=2, breaks=50)
        lines(density(pred.dens), lty=2, lwd=2)
      }
    
      n_gau <- length(mixmdl$mu)
      print(sprintf("Mixture of %d Gaussians", n_gau))
      
      gau1 <- function(x) 
        mixmdl$lambda[1] * 1/sqrt(2*mixmdl$sigma[1]^2*pi) * exp(-(x-mixmdl$mu[1])^2/(2*mixmdl$sigma[1]^2))

      gau2 <- function(x) 
        mixmdl$lambda[2] * 1/sqrt(2*mixmdl$sigma[2]^2*pi) * exp(-(x-mixmdl$mu[2])^2/(2*mixmdl$sigma[2]^2))

      # split at the minimum of the mixture of Gaussians if there is one between
      # the two means. If there is no such point, split at the mean of the left Gaussian
      # gau_sum <- function(x) gau1(x) + gau2(x)
      # 
      # opt <- optimize(gau_sum, interval=range(mixmdl$mu[1], mixmdl$mu[2]))
      # cut <- opt$minimum
      
      # If the minum between the two gaussians is at the mean of the right 
      # Gaussian, set it to the peak of the left Gaussian
      # if(round(cut, 2) == round(max(c(mixmdl$mu[1], mixmdl$mu[2])), 2))
      #   cut <- min(c(mixmdl$mu[1], mixmdl$mu[2]))

      # Split where the two Gaussian functions cross (if there is such
      # a point) between the two means
      gau_diff <- function(x) gau1(x) - gau2(x)
      
      cross <- tryCatch(uniroot(gau_diff, range(mixmdl$mu[1], mixmdl$mu[2])), 
                        error=function(x) NULL)
      
      # If the Gaussians don't cross between their means, set the cut to the
      # smaller mean
      if(!length(cross)) {
        cut <- min(c(mixmdl$mu[1], mixmdl$mu[2]))
      } else cut <- cross$root

      if(plot) abline(v=cut, col='blue', lwd=2)
    
      cuts[[mode]][fold] <- cut
      if(iter >= 10)
        stop("One or more fold is not fit by two Gaussians")
    }
  }
}

if('m1.preds' %in% names(pred.mats)) {
  m1.preds <- pred.mats$m1.preds
  m1.mat <- m1.mat[,colnames(m1.mat) %in% row.names(m1.preds)]
} else m1.preds <- matrix(0, 0, folds)
if('m2.preds' %in% names(pred.mats)) {
  m2.preds <- pred.mats$m2.preds
  m2.mat <- m2.mat[,colnames(m2.mat) %in% row.names(m2.preds)]
} else m2.preds <- matrix(0, 0, folds)
if('m3.preds' %in% names(pred.mats)) {
  m3.preds <- pred.mats$m3.preds
  m3.mat <- m3.mat[,colnames(m3.mat) %in% row.names(m3.preds)]
} else m3.preds <- matrix(0, 0, folds)

# Make a matrix of predictors log(abs()) weights by fold. 
# Set non-selected predictor weights to -Inf
m1.select <- m1.preds; m2.select <- m2.preds; m3.select <- m3.preds
if(nrow(m1.preds)) m1.select[,] <- -Inf
if(nrow(m2.preds)) m2.select[,] <- -Inf
if(nrow(m3.preds)) m3.select[,] <- -Inf

for(fold in 1:folds) {
  if(nrow(m1.preds) > 0) {
    sel <- which(log(m1.preds[,fold]) > cuts[['mode1']][fold])
    m1.select[,fold][sel] <- log(m1.preds[sel, fold])
  }
  if(nrow(m2.preds) > 0) {
    sel <- which(log(m2.preds[,fold]) > cuts[['mode2']][fold])
    m2.select[,fold][sel] <- log(m2.preds[sel, fold])
  }
  if(nrow(m3.preds) > 0) { 
    sel <- which(log(m3.preds[,fold]) > cuts[['mode3']][fold])
    m3.select[,fold][sel] <- log(m3.preds[sel, fold])
  }
}

# Matrices of predictors for each fold
######################################

# Save matrices of predictors for each fold
m1.mat.upGau.list <- list()
m2.mat.upGau.list <- list()
m3.mat.upGau.list <- list()

m1.mat.top15.list <- list()
m2.mat.top15.list <- list()
m3.mat.top15.list <- list()

for(fold in 1:folds) {
  m1.mat.upGau.list[[fold]] <- m1.mat[,colnames(m1.mat) %in% rownames(m1.preds)[m1.select[,fold] > -Inf]]
  m2.mat.upGau.list[[fold]] <- m2.mat[,colnames(m2.mat) %in% rownames(m2.preds)[m2.select[,fold] > -Inf]]
  m3.mat.upGau.list[[fold]] <- m3.mat[,colnames(m3.mat) %in% rownames(m3.preds)[m3.select[,fold] > -Inf]]

  m1.mat.top15.list[[fold]] <- m1.mat[,colnames(m1.mat) %in% rownames(m1.preds[m1.preds[,fold] >= quantile(m1.preds[,fold], .85, na.rm=T),])]
  m2.mat.top15.list[[fold]] <- m2.mat[,colnames(m2.mat) %in% rownames(m2.preds[m2.preds[,fold] >= quantile(m2.preds[,fold], .85, na.rm=T),])]
  m3.mat.top15.list[[fold]] <- m3.mat[,colnames(m3.mat) %in% rownames(m3.preds[m3.preds[,fold] >= quantile(m3.preds[,fold], .85, na.rm=T),])]

  m1.to.save <- m1.mat.upGau.list[[fold]]
  m2.to.save <- m2.mat.upGau.list[[fold]]
  m3.to.save <- m3.mat.upGau.list[[fold]]
  if(ncol(m1.to.save))
    save(m1.to.save, file=paste0(out.dir, '/upGau_fold_', fold-1, '_m1_mat.Rdata'))
  if(ncol(m2.to.save))
    save(m2.to.save, file=paste0(out.dir, '/upGau_fold_', fold-1, '_m2_mat.Rdata'))
  if(ncol(m3.to.save))
    save(m3.to.save, file=paste0(out.dir, '/upGau_fold_', fold-1, '_m3_mat.Rdata'))

  m1.to.save <- m1.mat.top15.list[[fold]]
  m2.to.save <- m2.mat.top15.list[[fold]]
  m3.to.save <- m3.mat.top15.list[[fold]]
  if(ncol(m1.to.save))
    save(m1.to.save, file=paste0(out.dir, '/top15_fold_', fold-1, '_m1_mat.Rdata'))
  if(ncol(m2.to.save))
    save(m2.to.save, file=paste0(out.dir, '/top15_fold_', fold-1, '_m2_mat.Rdata'))
  if(ncol(m3.to.save))
    save(m3.to.save, file=paste0(out.dir, '/top15_fold_', fold-1, '_m3_mat.Rdata'))
}

# Print the number of predictors in the upper groups for each fold
print(sprintf('Top 15%% of predictors for mode 1 ranges from %d to %d',
  min(sapply(m1.mat.top15.list, ncol)), max(sapply(m1.mat.top15.list, ncol))))
print(sprintf('Top 15%% of predictors for mode 2 ranges from %d to %d',
  min(sapply(m2.mat.top15.list, ncol)), max(sapply(m2.mat.top15.list, ncol))))
print(sprintf('Top 15%% of predictors for mode 3 ranges from %d to %d',
  min(sapply(m3.mat.top15.list, ncol)), max(sapply(m3.mat.top15.list, ncol))))

for(fold in 1:folds)
  print(sprintf('For fold %d, keep %d predictors for mode 1', 
                fold-1, ncol(m1.mat.upGau.list[[fold]])))
print(sprintf('Range: %d - %d', min(sapply(m1.mat.upGau.list, ncol)),
  max(sapply(m1.mat.upGau.list, ncol))))
print('##############')
for(fold in 1:folds)
  print(sprintf('For fold %d, keep %d predictors for mode 2', 
                fold-1, ncol(m2.mat.upGau.list[[fold]])))
print(sprintf('Range: %d - %d', min(sapply(m2.mat.upGau.list, ncol)),
  max(sapply(m2.mat.upGau.list, ncol))))
print('##############')
for(fold in 1:folds)
  print(sprintf('For fold %d, keep %d predictors for mode 3', 
                fold-1, ncol(m3.mat.upGau.list[[fold]])))
print(sprintf('Range: %d - %d', min(sapply(m3.mat.upGau.list, ncol)),
  max(sapply(m3.mat.upGau.list, ncol))))

# Overlap of predictors chosen by at least n_lap folds
######################################################

# Get a list of predictor names to keep
if(nrow(m1.preds)) {
  m1.keep <- rownames(m1.select)[rowSums(m1.select > -Inf) >= n_lap]
} else m1.keep <- c()
if(nrow(m2.preds)) {
  m2.keep <- rownames(m2.select)[rowSums(m2.select > -Inf) >= n_lap]
} else m2.keep <- c()
if(nrow(m3.preds)) {
  m3.keep <- rownames(m3.select)[rowSums(m3.select > -Inf) >= n_lap]
} else m3.keep <- c()

print(sprintf('For predictors selected by at least %d folds:', n_lap))
print(sprintf("Keeping %d features of %d (%.2f%%) for mode 1", length(m1.keep), nrow(m1.preds), 
  length(m1.keep)/nrow(m1.preds)*100))
print(sprintf("Keeping %d features of %d (%.2f%%) for mode 2", length(m2.keep), nrow(m2.preds), 
  length(m2.keep)/nrow(m2.preds)*100))
print(sprintf("Keeping %d features of %d (%.2f%%) for mode 3", length(m3.keep), nrow(m3.preds), 
  length(m3.keep)/nrow(m3.preds)*100))

# Subset data matrices to just these predictors
m1.mat.sub <- m1.mat[,colnames(m1.mat) %in% m1.keep]
m2.mat.sub <- m2.mat[,colnames(m2.mat) %in% m2.keep]
m3.mat.sub <- m3.mat[,colnames(m3.mat) %in% m3.keep]

# Save matrices of intersecting predictors
if(ncol(m1.mat.sub))
  save(m1.mat.sub, file=paste0(out.dir, '/keep', n_lap, '_m1_mat.Rdata'))
if(ncol(m2.mat.sub))
  save(m2.mat.sub, file=paste0(out.dir, '/keep', n_lap, '_m2_mat.Rdata'))
if(ncol(m3.mat.sub))
  save(m3.mat.sub, file=paste0(out.dir, '/keep', n_lap, '_m3_mat.Rdata'))

# Union of the top per% of predictors for all folds (per=5%, 10%, 15%, 20%, 25%)
################################################################################
for(per in c(.05, .1, .15, .2, .25)) {
  per_str <- sub('0.', '_', sprintf('%.2f',per), fixed=T)

  if(nrow(m1.mat)) {
    m1.quantiles <- apply(m1.select, 2, function(x) quantile(x[x > -Inf], 1-per, na.rm=T))
    m1.keep.mat <- sweep(m1.select, 2, m1.quantiles, FUN='>')
    print(sprintf("A union of the top %d%% of predictors for mode 1 gives %d (%.1f%%)",
                  (per*100), sum(rowSums(m1.keep.mat)>0), 
                  sum(rowSums(m1.keep.mat)>0)/nrow(m1.keep.mat)*100))
  
    # Make the mode 1 matrix with the predictors chosen
    m1.keep <- rownames(m1.keep.mat)[rowSums(m1.keep.mat) > 0]
    m1.mat.sub <- m1.mat[,colnames(m1.mat) %in% m1.keep]

    save(m1.mat.sub, file=paste0(out.dir, '/union', per_str, '_m1_mat.Rdata'))
  }

  if(nrow(m2.mat)) {
    m2.quantiles <- apply(m2.select, 2, function(x) quantile(x[x > -Inf], 1-per, na.rm=T))
    m2.keep.mat <- sweep(m2.select, 2, m2.quantiles, FUN='>')
    print(sprintf("A union of the top %d%% of predictors for mode 2 gives %d (%.1f%%)",
                  (per*100), sum(rowSums(m2.keep.mat)>0),
                  sum(rowSums(m2.keep.mat)>0)/nrow(m2.keep.mat)*100))

    # Make the mode 2 matrix with the predictors chosen
    m2.keep <- rownames(m2.keep.mat)[rowSums(m2.keep.mat) > 0]
    m2.mat.sub <- m2.mat[,colnames(m2.mat) %in% m2.keep]

    save(m2.mat.sub, file=paste0(out.dir, '/union', per_str, '_m2_mat.Rdata'))
  }

  if(nrow(m3.mat)) {
    m3.quantiles <- apply(m3.select, 2, function(x) quantile(x[x > -Inf], 1-per, na.rm=T))
    m3.keep.mat <- sweep(m3.select, 2, m3.quantiles, FUN='>')
    print(sprintf("A union of the top %d%% of predictors for mode 3 gives %d (%.1f%%)",
                  (per*100), sum(rowSums(m3.keep.mat)>0),
                  sum(rowSums(m3.keep.mat)>0)/nrow(m3.keep.mat)*100))
  
    # Make the mode 3 matrix with the predictors chosen
    m3.keep <- rownames(m3.keep.mat)[rowSums(m3.keep.mat) > 0]
    m3.mat.sub <- m3.mat[,colnames(m3.mat) %in% m3.keep]

    save(m3.mat.sub, file=paste0(out.dir, '/union', per_str, '_m3_mat.Rdata'))
  }
}
