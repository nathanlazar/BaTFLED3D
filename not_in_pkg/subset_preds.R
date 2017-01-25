# nathan dot lazar at gmail dot com

# Subset predictors by quantile (q) provided. Predictors above the q quantile in
# any run (if union is specified) will be used. If intersection is specified, 
# predictors above the q quantile in all runs will be used. 

# First argument is Rdata object created with summarize_preds.R 

# Usage: subset_preds.R <*_preds.Rdata> <m1_matrix.Rdata> <m2_matrix.Rdata>
#                       <m3_matrix.Rdata> <quantile> <folds> <union or intersection>
#                       <out directory>

# args <- list('CTRP2results/run_225280/run_225280_preds.Rdata',
#              '../CTRP2/CTRP2_cl_train.Rdata',
#              '../CTRP2/CTRP2_dr_train.Rdata',
#              'none', '0.95', '10', 'union', '../CTRP2/top5')

library(dplyr)
library(methods)
library(R6)
library(BaTFLED3D)
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

pred.mats <- loadRData(args[[1]])
if(args[[2]] != 'none') {
  m1.mat <- loadRData(args[[2]])
} else m1.mat <- matrix(NA, 0, 0)
if(args[[3]] != 'none') {
  m2.mat <- loadRData(args[[3]])
} else m2.mat <- matrix(NA, 0, 0)
if(args[[4]] != 'none') {
  m3.mat <- loadRData(args[[4]])
} else m3.mat <- matrix(NA, 0, 0)
q <- as.numeric(args[[5]])
folds <- as.numeric(args[[6]])
type <- args[[7]]
out.dir <- args[[8]]

if(!(type %in% c('union', 'intersection'))) 
  print("Type (argument 7) must be either 'union' or 'intersection'")

# If the output directory doesn't exist, create it
print(paste('Output directory:', out.dir))
if (!file.exists(out.dir)){
  dir.create(file.path(out.dir))
}

if('m1.preds' %in% names(pred.mats)) {
  m1.preds <- pred.mats$m1.preds
} else m1.preds <- matrix(0, 0, folds)
if('m2.preds' %in% names(pred.mats)) {
  m2.preds <- pred.mats$m2.preds
} else m2.preds <- matrix(0, 0, folds)
if('m3.preds' %in% names(pred.mats)) {
  m3.preds <- pred.mats$m3.preds
} else m3.preds <- matrix(0, 0, folds)
  
# Keep the top <q>% of predictors from any (or all) of the <folds> cross-validation runs
m1.keep <- c(); m2.keep <- c(); m3.keep <- c()

if(type == 'union') {
  for(fold in 1:folds) {
    if(nrow(m1.preds) > 0) 
      m1.keep <- c(m1.keep, row.names(m1.preds)[m1.preds[,fold] > quantile(m1.preds[,fold], q, na.rm=T)])
    if(nrow(m2.preds) > 0) 
      m2.keep <- c(m2.keep, row.names(m2.preds)[m2.preds[,fold] > quantile(m2.preds[,fold], q, na.rm=T)])
    if(nrow(m3.preds) > 0) 
      m3.keep <- c(m3.keep, row.names(m3.preds)[m3.preds[,fold] > quantile(m3.preds[,fold], q, na.rm=T)])
  }
  
  # End up keeping ~30% of predictors w/ q=.95.
  m1.keep <- unique(m1.keep)
  m2.keep <- unique(m2.keep)
  m3.keep <- unique(m3.keep)
}

if(type == 'intersection') {
  if(nrow(m1.preds) > 0)
    m1.keep <- row.names(m1.preds)[m1.preds[,1] > quantile(m1.preds[,1], q, na.rm=T)]
  if(nrow(m2.preds) > 0)
    m2.keep <- row.names(m2.preds)[m2.preds[,1] > quantile(m2.preds[,1], q, na.rm=T)]
  if(nrow(m3.preds) > 0)
    m3.keep <- row.names(m3.preds)[m3.preds[,1] > quantile(m3.preds[,1], q, na.rm=T)]

  for(fold in 2:folds) {
    if(nrow(m1.preds) > 0)
      m1.keep <- m1.keep[m1.keep %in% row.names(m1.preds)[m1.preds[,fold] > quantile(m1.preds[,fold], q, na.rm=T)]]
    if(nrow(m2.preds) > 0)
      m2.keep <- m2.keep[m2.keep %in% row.names(m2.preds)[m2.preds[,fold] > quantile(m2.preds[,fold], q, na.rm=T)]]
    if(nrow(m3.preds) > 0)
      m3.keep <- m3.keep[m3.keep %in% row.names(m3.preds)[m3.preds[,fold] > quantile(m3.preds[,fold], q, na.rm=T)]]
  }

  m1.keep <- m1.keep[!is.na(m1.keep)]
  m2.keep <- m2.keep[!is.na(m2.keep)]
  m3.keep <- m3.keep[!is.na(m3.keep)]
}

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

# Save data matrices
save(m1.mat.sub, file=paste0(out.dir, '/m1_mat.Rdata'))
save(m2.mat.sub, file=paste0(out.dir, '/m2_mat.Rdata'))
save(m3.mat.sub, file=paste0(out.dir, '/m3_mat.Rdata'))

##################################################################################

# Look at correlations between folds
# m1.cors <- cor(m1.preds)
# m2.cors <- cor(m2.preds, use="pairwise.complete.obs")

# im_mat(m1.cors)
# im_mat(m2.cors)

# range(m1.cors - diag(NA, 20), na.rm=T)
# range(m2.cors - diag(NA, 20), na.rm=T)

# barplot(apply(pred.mats$m1.preds, 1, mean), las=2)
# barplot(apply(pred.mats$m2.preds, 1, mean), las=2)

# show that by removing these predictors the influence of the remaining predictors
# is closer to a normal distribution.
# m1.preds <- as.matrix(m1.preds)
# m2.preds <- as.matrix(m2.preds)
# 
# par(mfrow=c(1,2))
# hist(m1.preds, breaks=100)
# hist(m1.preds[!(rownames(m1.preds) %in% m1.keep),], breaks=100)
# 
# # Alternatively use K-S test and maybe Student's T distribution.
# shapiro.test(sample(m1.preds, 1000))
# shapiro.test(sample(m1.preds[!(rownames(m1.preds) %in% m1.keep),], 1000))
# 
# qqnorm(m1.preds)
# qqnorm(m1.preds[!(rownames(m1.preds) %in% m1.keep),])

# Alternatively subset by means across all folds
# pred.mats <- lapply(pred.mats, scale)
# 
# m1.means <- apply(pred.mats$m1.preds, 1, mean, na.rm=T)
# m1.top <- m1.means[m1.means > quantile(m1.means, q)]
# 
# m2.means <- apply(pred.mats$m2.preds, 1, mean, na.rm=T)
# m2.top <- m2.means[m2.means > quantile(m2.means, q)]
# 
# m3.means <- apply(pred.mats$m3.preds, 1, mean, na.rm=T)
# m3.top <- m3.means[m3.means > quantile(m3.means, q)]
# 

