#!/usr/bin/env Rscript

# nathan dot lazar at gmail dot com

# Usage: BaTFLED_final_test.R <image from model.Rdata>
#                             <matrix of test input features .Rdata>
#                             <tensor of test responses>

# Example: BaTFLED_final_test.R DREAM7results/run_1223578.0/image.Rdata \
#                               ../DREAM7/DREAM7_cl_test.Rdata
#                               ../DREAM7/raw_median_tensor.Rdata \

library(R6)             # For making memory efficent R6 objects
library(rTensor)        # Needed for multiplying matrices into tensors (could be removed)
library(dplyr)          # General data frame manipulation
# library(ggplot2)        # Used for plotting
# library(drc)            # Used for fitting logistic curves
library(BaTFLED3D)      # Note: package must be installed from .tar file

print('Packages loaded')
args <- commandArgs(TRUE)

# Report arguments and summary stats of the input objects
print(unlist(args))

# Read in options
in.data <-        args[[1]]
test.features <-  args[[2]]
test.responses <- args[[3]]

# Load input data
############################################################
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

load(in.data)
final.m1.mat <- loadRData(test.features)
final.resp.tens <- loadRData(test.responses)

# Redefine RMSE function
rmse <- function(obs, pred, stdev=sd(train.data$resp, na.rm=T)) 
  sqrt(mean((obs-pred)^2, na.rm=T))/stdev

# Fix some drug names
dimnames(final.resp.tens)[[2]][dimnames(final.resp.tens)[[2]] == 'IKK 16'] <- 'IKK'
dimnames(final.resp.tens)[[2]][dimnames(final.resp.tens)[[2]] == 'Nilontinib'] <- 'Nilotinib'
dimnames(final.resp.tens)[[2]][dimnames(final.resp.tens)[[2]] == 'Methylglyoxol'] <- 'Methylglyoxal'

# Reorder/select predictors to match the training data
final.m1.mat <- final.m1.mat[,match(colnames(m1.mat), colnames(final.m1.mat), nomatch=0)]

# Reorder responses to match the inputs
final.resp.tens <- final.resp.tens[
  match(rownames(final.m1.mat), dimnames(final.resp.tens)[[1]], nomatch=0),
  match(rownames(m2.mat), dimnames(final.resp.tens)[[2]], nomatch=0),
  match(rownames(m3.mat), dimnames(final.resp.tens)[[3]], nomatch=0)]

# Make the testing data objects (input_data_3d):
final.test.data <- input_data$new(mode1.X=final.m1.mat,
                                  mode2.X=m2.mat,
                                  mode3.X=m3.mat,
                                  resp=final.resp.tens)

save.image(paste0(out.dir, 'final_test.Rdata'))

final.cold.preds <- test(final.test.data, trained)
final.RMSE <- rmse(final.test.data$resp, final.cold.preds)

final.cold.preds.clip <- final.cold.preds
final.cold.preds.clip[final.cold.preds.clip < min(train.data$resp, na.rm=T)] <- min(train.data$resp, na.rm=T)
final.cold.preds.clip[final.cold.preds.clip > max(train.data$resp, na.rm=T)] <- max(train.data$resp, na.rm=T)

final.RMSE.clip <- rmse(final.test.data$resp, final.cold.preds.clip)
final.exp.var <- exp_var(final.test.data$resp, final.cold.preds)
final.exp.var.clip <- exp_var(final.test.data$resp, final.cold.preds.clip)
final.p.cor <- cor(final.test.data$resp, final.cold.preds, use='complete.obs')
final.p.cor.clip <- cor(final.test.data$resp, final.cold.preds.clip, use='complete.obs')
final.s.cor <- cor(final.test.data$resp, final.cold.preds, use='complete.obs', method='spearman')
final.s.cor.clip <- cor(final.test.data$resp, final.cold.preds.clip, use='complete.obs', method='spearman')

print(sprintf('Final cold RMSE: %.2f, clipped: %.2f', final.RMSE, final.RMSE.clip))
print(sprintf("Final cold explained variance: %.2f, clipped %.2f", final.exp.var, final.exp.var.clip))
print(sprintf("Final cold Pearson correlation: %.2f, clipped %.2f", final.p.cor, final.p.cor.clip))
print(sprintf("Final cold Spearman correlation: %.2f, clipped %.2f", final.s.cor, final.s.cor.clip))

# Get results for predicting the mean response
m1.mean.tens <- final.resp.tens
for(i in 1:dim(final.resp.tens)[1])
  m1.mean.tens[i,,] <- m1.means

final.mean.RMSE <- rmse(final.test.data$resp, m1.mean.tens)
final.mean.exp.var <- exp_var(final.test.data$resp, m1.mean.tens)
final.mean.p.cor <- cor(final.test.data$resp, m1.mean.tens, use='complete.obs')
final.mean.s.cor <- cor(final.test.data$resp, m1.mean.tens, use='complete.obs', method='spearman')

save.image(paste0(out.dir, 'final_test.Rdata'))