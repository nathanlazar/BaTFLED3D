# nathan dot lazar at gmail dot com

# Runs a neural net model on cell line/drug response curves using 
# all predictors for cell lines and drugs

# Usage: run_H2O.R <m1.Rdata> <m2.Rdata> <m3.Rdata> <resp.Rdata>
#                  <m1.cv.folds.Rdata> <m2.cv.folds.Rdata>
#                  <out.dir> <options>

library(methods)        # Need to include this for basic functioning on some systems
library(h2o)            # For training artificial neural nets and others
library(reshape2)       # For reshaping data frames
library(preprocessCore) # Used to quantile normalize columns of input data
library(BaTFLED3D)      # For evaluation and plotting functions
# h2o requires a Java run time environment

print('Packages loaded')
args <- commandArgs(TRUE)

test <- F
if(test)
  args <- list('../HeiserData/Cell_Line_Data/exp_mat.Rdata',
               '../HeiserData/Drug_Data/all_mat.Rdata',
               'none', 
               '../HeiserData/Responses/norm_tens.Rdata',
               '../HeiserData/Cell_Line_Data/cl_5cv_mat.Rdata', 
               'none', 'none',
               'HEISERresults/test/',
               'cores=1', 
               'quantile=F', 'center=F',
               'scale=none', 'squash=F', 
 	             'warm_per=0',
               'layers=1', 'hidden=1500',
               'activation=RectifierWithDropout',
               'input_dropout=0.2', 'hidden_dropout=0.5',
               'epochs=5000',
               'stop_metric=MSE', 'stop_toler=0.01',
               'stop_rounds=5', 'cv_folds=0',
               'fold=0')

# Report arguments
print(unlist(args))

# Read in options
m1.file <-       args[[1]]
m2.file <-       args[[2]]
m3.file <-       args[[3]]
resp.file <-     args[[4]]
m1.fold.file <-  args[[5]]
m2.fold.file <-  args[[6]]
m3.fold.file <-  args[[7]]
out.dir <-       paste0(args[[8]], '/')
seed <-          as.numeric(gsub('seed=', '', args[grepl('seed=', args)]))
warm.per <-      as.numeric(gsub('warm_per=', '', args[grepl('warm_per=', args)]))
quantile <-      as.logical(gsub('quantile=', '', args[grepl('quantile=', args)]))
scale <-         gsub('scale=', '', args[grepl('scale=', args)])
center <-        as.logical(gsub('center=', '', args[grepl('center=', args)]))
squash <-        as.logical(gsub('squash=', '', args[grepl('squash=', args)]))
cores <-         as.numeric(gsub('cores=', '', args[grepl('cores=', args)]))
layers <-         as.numeric(gsub('layers=', '', args[grepl('layers=', args)]))
hidden <-         as.numeric(gsub('hidden=', '', args[grepl('hidden=', args)]))
activation <-     gsub('activation=', '', args[grepl('activation=', args)])
input_dropout <-  as.numeric(gsub('input_dropout=', '', args[grepl('input_dropout=', args)]))
hidden_dropout <- as.numeric(gsub('hidden_dropout=', '', args[grepl('hidden_dropout=', args)]))
epochs <-         as.numeric(gsub('epochs=', '', args[grepl('epochs=', args)]))
stop_metric <-    gsub('stop_metric=', '', args[grepl('stop_metric=', args)])
stop_toler <-     as.numeric(gsub('stop_toler=', '', args[grepl('stop_toler=', args)]))
stop_rounds <-    as.numeric(gsub('stop_rounds=', '', args[grepl('stop_rounds=', args)]))
cv_folds <-       as.numeric(gsub('cv_folds=', '', args[grepl('cv_folds=', args)]))
fold <-          as.numeric(gsub('fold=', '', args[grepl('fold=', args)])) +1
                 # +1 because HT condor is zero based.

# Set defaults
if(length(seed)==0)     seed       <- NA       # Default to a random seed
if(length(warm.per)==0) warm.per   <- 0.01     # Default to removing 1% responses for warm prediction
if(length(quantile)==0) quantile <- F          # Default to quantile normalizing non-binary predictors
if(length(scale)==0) scale <- 'all'            # Default to scaling predictors
if(length(center)==0) center <- T              # Default to centering non-binary predictors
if(length(squash)==0) squash <- F              # Default to not squashing down outliers
if(length(cores)==0) cores         <- 1        # Default to using one core
if(length(layers)==0)    layers      <- 1                      # 1 hidden layer
if(length(hidden)==0)    hidden      <- 512                    # 512 hidden nodes per layer
if(length(activation)==0) activation <- 'RectifierWithDropout' # Relu activation w/ dropout
if(length(input_dropout)==0) input_dropout <- 0                # No dropout
if(length(hidden_dropout)==0) hidden_dropout <- 0              # No dropout
if(length(epochs)==0) epochs <- 5000                           # Max of 5000 epochs
if(length(stop_metric)==0) stop_metric <- 'MSE' # Default to mean squared error
if(length(stop_toler)==0) stop_toler   <- 0.01  # Default to stopping when change < 0.01
if(length(stop_rounds)==0) stop_rounds <- 5     # Default to looking at 5 rounds when stopping
if(length(cv_folds)==0) cv_folds         <- 1   # Default to no cross-validation
if(length(fold)==0) fold           <- 1        # Default to fold 1

# Functions
###########################################################
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Predict function
predict <- function(net, predictors) {
  if(ncol(predictors) == P+Q+S)
    res <- h2o.predict(net, as.h2o(predictors))
  if(ncol(predictors) == P+Q+S+1)
    res <- h2o.predict(net, as.h2o(predictors[,1:(P+Q+S)]))
  return(as.vector(res))
}

p_cor <- function(x,y)
  if((sum(!is.na(x)) > 0) && length(y) && sum(!is.na(y)>0)) {
    return(stats::cor(as.vector(x),as.vector(y),use='complete.obs'))
  } else return(NA)

s_cor <- function(x,y)
  if((sum(!is.na(x)) > 0) && length(y) && sum(!is.na(y)>0)) {
    return(stats::cor(as.vector(x),as.vector(y),use='complete.obs', method='spearman'))
  } else return(NA)

##################################################
# If the output directory doesn't exist, create it
print(paste('Output directory:', out.dir))
try(dir.create(file.path(out.dir)))

# Load data
###################################################

# Load the response data
resp.tens <- loadRData(resp.file)

# Load and rename input data if supplied otherwise make empty matrices
if(m1.file != 'none') m1.mat <- loadRData(m1.file) else {
  m1.mat <- matrix(nrow=dim(resp.tens)[1], ncol=0)
  dimnames(m1.mat)[[1]] <- dimnames(resp.tens)[[1]]
}
if(m2.file != 'none') m2.mat <- loadRData(m2.file) else {
  m2.mat <- matrix(nrow=dim(resp.tens)[2], ncol=0)
  dimnames(m2.mat)[[1]] <- dimnames(resp.tens)[[2]]
}
if(m3.file != 'none') m3.mat <- loadRData(m3.file) else {
  m3.mat <- matrix(nrow=dim(resp.tens)[3], ncol=0)
  dimnames(m3.mat)[[1]] <- dimnames(resp.tens)[[3]]
}

# Check that the input are the same size and in the same order
m1.mat <- m1.mat[match(dimnames(resp.tens)[[1]], rownames(m1.mat), nomatch=0),]
m2.mat <- m2.mat[match(dimnames(resp.tens)[[2]], rownames(m2.mat), nomatch=0),]
m3.mat <- m3.mat[match(dimnames(resp.tens)[[3]], rownames(m3.mat), nomatch=0),]

resp.tens <- resp.tens[match(rownames(m1.mat), dimnames(resp.tens)[[1]], nomatch=0),
                       match(rownames(m2.mat), dimnames(resp.tens)[[2]], nomatch=0),
                       match(rownames(m3.mat), dimnames(resp.tens)[[3]], nomatch=0)]

# Load the matrices with names for cross validation
if(m1.fold.file != 'none') {
  m1.cv.fold <- loadRData(m1.fold.file)
  all.m1 <- unique(as.vector(m1.cv.fold))
  m1.mat <- m1.mat[rownames(m1.mat) %in% all.m1,]
  resp.tens <- resp.tens[dimnames(resp.tens)[[1]] %in% all.m1,,]
} else 
  m1.cv.fold <- matrix('', nrow=fold, ncol=0)
if(m2.fold.file != 'none') {
  m2.cv.fold <- loadRData(m2.fold.file)
  all.m2 <- unique(as.vector(m2.cv.fold))
  m2.mat <- m2.mat[rownames(m2.mat) %in% all.m2,]
  resp.tens <- resp.tens[,dimnames(resp.tens)[[2]] %in% all.m2,]
} else 
  m2.cv.fold <- matrix('', nrow=fold, ncol=0)

# Remove predictors with no variance
m1.mat <- m1.mat[,apply(m1.mat, 2, sd)!=0]
m2.mat <- m2.mat[,apply(m2.mat, 2, sd)!=0]

# Determine which predictors are binary
m1.bin <- apply(m1.mat, 2, function(x) identical(range(x), c(0,1)))
m2.bin <- apply(m2.mat, 2, function(x) identical(range(x), c(0,1)))
m3.bin <- apply(m3.mat, 2, function(x) identical(range(x), c(0,1)))

# Center non-binary columns
if(center) {
  # Center non-binary columns of the predictor matrices
  m1.mat[,!m1.bin] <- scale(m1.mat[,!m1.bin], scale=F, center=T)
  m2.mat[,!m2.bin] <- scale(m2.mat[,!m2.bin], scale=F, center=T)
  m3.mat[,!m3.bin] <- scale(m3.mat[,!m3.bin], scale=F, center=T)
}

# Scale non-binary columns
if(scale=='all') {
  m1.mat <- scale(m1.mat, center=F, scale=apply(m1.mat, 2, sd, na.rm=T))
  m2.mat <- scale(m2.mat, center=F, scale=apply(m2.mat, 2, sd, na.rm=T))
  m3.mat <- scale(m3.mat, center=F, scale=apply(m3.mat, 2, sd, na.rm=T))
  } else if(scale=='cont') {
  # scale columns of the predictor matrices
  m1.mat[,!m1.bin] <- scale(m1.mat[,!m1.bin], center=F,
                            scale=apply(m1.mat[,!m1.bin], 2, sd, na.rm=T))
  m2.mat[,!m2.bin] <- scale(m2.mat[,!m2.bin], center=F,
                            scale=apply(m2.mat[,!m2.bin], 2, sd, na.rm=T))
  m3.mat[,!m3.bin] <- scale(m3.mat[,!m3.bin], center=F,
                            scale=apply(m3.mat[,!m3.bin], 2, sd, na.rm=T))
}

# squash non-binary outliers more than <sqash> s.d. from their mean
if(squash) {
  m1.mat[,!m1.bin][m1.mat[,!m1.bin] < -squash] <- -squash
  m1.mat[,!m1.bin][m1.mat[,!m1.bin] > squash] <- squash
  m2.mat[,!m2.bin][m2.mat[,!m2.bin] < -squash] <- -squash
  m2.mat[,!m2.bin][m2.mat[,!m2.bin] > squash] <- squash
  m3.mat[,!m3.bin][m2.mat[,!m3.bin] < -squash] <- -squash
  m3.mat[,!m3.bin][m2.mat[,!m3.bin] > squash] <- squash
}

# quantile normalize non-binary columns (treat gex and cnv separately)
if(quantile) {
  m1.cnv <- grepl('cnv', colnames(m1.mat))
  m1.gex <- grepl('gex', colnames(m1.mat))

  if(ncol(m1.mat)) m1.mat[,m1.cnv] <- normalize.quantiles(t(m1.mat[,m1.cnv]))
  if(ncol(m1.mat)) m1.mat[,m1.gex] <- normalize.quantiles(t(m1.mat[,m1.gex]))

  # if(ncol(m1.mat)) m1.mat[,!m1.bin] <- normalize.quantiles(t(m1.mat[,!m1.bin]))
  if(ncol(m2.mat)) m2.mat[,!m2.bin] <- normalize.quantiles(t(m2.mat[,!m2.bin]))
  if(ncol(m3.mat)) m3.mat[,!m3.bin] <- normalize.quantiles(t(m3.mat[,!m3.bin]))
}

# Remove the data for this cross validation fold 
# (removing columns if these are kernel matrices)
###############################################################
test.m1 <- m1.cv.fold[fold,]
if(length(test.m1)) {
  test.m1.mat <- m1.mat[row.names(m1.mat) %in% test.m1,
    !(colnames(m1.mat) %in% test.m1), drop=F]
} else test.m1.mat <- m1.mat[0,]
train.m1.mat <- m1.mat[!(row.names(m1.mat) %in% test.m1),
  !(colnames(m1.mat) %in% test.m1),drop=F]

test.m2 <- m2.cv.fold[fold,]
if(length(test.m2)) {
  test.m2.mat <- m2.mat[row.names(m2.mat) %in% test.m2,
    !(colnames(m2.mat) %in% test.m2),drop=F]
} else test.m2.mat <- m2.mat[0,]
train.m2.mat <- m2.mat[!(row.names(m2.mat) %in% test.m2),
  !(colnames(m2.mat) %in% test.m2),drop=F]

# Subset and reorder responses to match predictor set
resp.train <- resp.tens[dimnames(resp.tens)[[1]] %in% row.names(train.m1.mat),
                        dimnames(resp.tens)[[2]] %in% row.names(train.m2.mat),,drop=F]
resp.test.m1 <- resp.tens[dimnames(resp.tens)[[1]] %in% test.m1,
                          dimnames(resp.tens)[[2]] %in% row.names(train.m2.mat),,drop=F]
resp.test.m2 <- resp.tens[dimnames(resp.tens)[[1]] %in% row.names(train.m1.mat),
                          dimnames(resp.tens)[[2]] %in% test.m2,,drop=F]
resp.test.m1m2 <- resp.tens[dimnames(resp.tens)[[1]] %in% test.m1,
                            dimnames(resp.tens)[[2]] %in% test.m2,,drop=F]

# Remove warm.per percent of training responses for warm prediction
# Note: Some training responses are already missing due to incomplete data
all.resp.train <- resp.train   # Stores responses before warms removed
if(warm.per > 0) {
  mask <- sample(prod(dim(resp.train)[1:2]), round(warm.per*prod(dim(resp.train)[1:2])))
  warm.test.resp <- list()
  for(i in 1:dim(resp.train)[3]) {
    resp.train[,,i][mask] <- NA
    warm.test.resp[[i]] <- all.resp.train[,,i][is.na(resp.train[,,i])]
  }
}

###############################################################
I <- dim(resp.train)[1]
J <- dim(resp.train)[2]
K <- dim(resp.train)[3]
P <- ncol(train.m1.mat)
Q <- ncol(train.m2.mat)
S <- ncol(m3.mat)

print(sprintf("Random forest is being trained with %d mode1 samples, %d mode2 samples, and %d mode3 samples",
              I,J,K))
print(sprintf("Cold start predictions will be made for %d mode1 samples, %d mode2 samples, and %d mode3 samples",
              dim(resp.test.m1)[[1]], dim(resp.test.m2)[[2]], 0))
print(sprintf("Using %d mode 1 predictors, %d mode 2 predictors and %d mode 3 predictors",
              ncol(m1.mat), ncol(m2.mat), 0))
print(paste0(round(sum(is.na(resp.train))/prod(dim(resp.train))*100, 2),
            "% of responses are missing, ",
            round(sum(is.na(resp.train) & !is.na(all.resp.train))/prod(dim(resp.train))*100,2),
            "% are used for warm predictions"))

# Get the responses for training and test data if we just predict the mean
# response for that sample given others in the training data.
m1.means <- apply(resp.train, c(1,3), mean, na.rm=T)
m1.sds <- apply(resp.train, c(1,3), sd, na.rm=T)
m2.means <- apply(resp.train, c(2,3), mean, na.rm=T)
m2.sds <- apply(resp.train, c(2,3), sd, na.rm=T)
m1m2.means <- apply(resp.train, 3, mean, na.rm=T)
m1m2.sds <- apply(resp.train, 3, sd, na.rm=T)

# Get RMSE, exp. var., Pearson correlation & Spearman correlation for predicting 
# the mean response at each mode3
train.m1.mean.pred <- aperm(array(m2.means, dim=c(dim(m2.means), I)), c(3,1,2))
train.m2.mean.pred <- aperm(array(m1.means, dim=c(dim(m1.means), J)), c(1,3,2))
train.m1m2.mean.pred <- aperm(array(m1m2.means, dim=c(K, I, J)), c(2,3,1))

warm.m1.mean.pred <- train.m1.mean.pred[which(is.na(resp.train))]
warm.m2.mean.pred <- train.m2.mean.pred[which(is.na(resp.train))]
warm.m1m2.mean.pred <- train.m1m2.mean.pred[which(is.na(resp.train))]

m1.mean.pred <- aperm(array(m2.means, dim=c(dim(m2.means), nrow(test.m1.mat))), c(3,1,2))
m2.mean.pred <- aperm(array(m1.means, dim=c(dim(m1.means), nrow(test.m2.mat))), c(1,3,2))
m1m2.mean.pred <- aperm(array(m1m2.means, dim=c(K, nrow(test.m1.mat), nrow(test.m2.mat))), c(2,3,1))

train.m1.mean.rmse <- nrmse(resp.train, train.m1.mean.pred)
train.m2.mean.rmse <- nrmse(resp.train, train.m2.mean.pred)
train.m1m2.mean.rmse <- nrmse(resp.train, train.m1m2.mean.pred)
if(warm.per > 0) {
  warm.m1.mean.rmse <- nrmse(unlist(warm.test.resp), warm.m1.mean.pred)
  warm.m2.mean.rmse <- nrmse(unlist(warm.test.resp), warm.m2.mean.pred)
  warm.m1m2.mean.rmse <- nrmse(unlist(warm.test.resp), warm.m1m2.mean.pred)
}
m1.mean.rmse <- nrmse(resp.test.m1, m1.mean.pred)
m2.mean.rmse <- nrmse(resp.test.m2, m2.mean.pred)
m1m2.mean.rmse <- nrmse(resp.test.m1m2, m1m2.mean.pred)

train.m1.mean.exp.var <- exp_var(resp.train, train.m1.mean.pred)
train.m2.mean.exp.var <- exp_var(resp.train, train.m2.mean.pred)
train.m1m2.mean.exp.var <- exp_var(resp.train, train.m1m2.mean.pred)
if(warm.per > 0) {
  warm.m1.mean.exp.var <- exp_var(unlist(warm.test.resp), warm.m1.mean.pred)
  warm.m2.mean.exp.var <- exp_var(unlist(warm.test.resp), warm.m2.mean.pred)
  warm.m1m2.mean.exp.var <- exp_var(unlist(warm.test.resp), warm.m1m2.mean.pred)
}
m1.mean.exp.var <- exp_var(resp.test.m1, m1.mean.pred)
m2.mean.exp.var <- exp_var(resp.test.m2, m2.mean.pred)
m1m2.mean.exp.var <- exp_var(resp.test.m1m2, m1m2.mean.pred)

train.m1.mean.p.cor <- p_cor(resp.train, train.m1.mean.pred)
train.m2.mean.p.cor <- p_cor(resp.train, train.m2.mean.pred)
train.m1m2.mean.p.cor <- p_cor(resp.train, train.m1m2.mean.pred)
if(warm.per > 0) {
  warm.m1.mean.p.cor <- p_cor(unlist(warm.test.resp), warm.m1.mean.pred)
  warm.m2.mean.p.cor <- p_cor(unlist(warm.test.resp), warm.m2.mean.pred)
  warm.m1m2.mean.p.cor <- p_cor(unlist(warm.test.resp), warm.m1m2.mean.pred)
}
if(dim(resp.test.m1)[[1]])
  m1.mean.p.cor <- p_cor(resp.test.m1, m1.mean.pred)
if(dim(resp.test.m2)[[2]])
  m2.mean.p.cor <- p_cor(resp.test.m2, m2.mean.pred)
if(dim(resp.test.m1m2)[[1]] & dim(resp.test.m1m2)[[2]])
  m1m2.mean.p.cor <- p_cor(resp.test.m1m2, m1m2.mean.pred)

train.m1.mean.s.cor <- s_cor(resp.train, train.m1.mean.pred)
train.m2.mean.s.cor <- s_cor(resp.train, train.m2.mean.pred)
train.m1m2.mean.s.cor <- s_cor(resp.train, train.m1m2.mean.pred)
if(warm.per > 0) {
  warm.m1.mean.s.cor <- s_cor(unlist(warm.test.resp), warm.m1.mean.pred)
  warm.m2.mean.s.cor <- s_cor(unlist(warm.test.resp), warm.m2.mean.pred)
  warm.m1m2.mean.s.cor <- s_cor(unlist(warm.test.resp), warm.m1m2.mean.pred)
}
if(dim(resp.test.m1)[[1]])
  m1.mean.s.cor <- s_cor(resp.test.m1, m1.mean.pred)
if(dim(resp.test.m2)[[2]])
  m2.mean.s.cor <- s_cor(resp.test.m2, m2.mean.pred)
if(dim(resp.test.m1m2)[[1]] & dim(resp.test.m1m2)[[2]])
  m1m2.mean.s.cor <- s_cor(resp.test.m1m2, m1m2.mean.pred)

print('Predicting mean responses per mode3 from training data')
print('RMSEs relative to training response std. dev.:')
print(sprintf("RMSE for training data predicting means across mode2: %.3f", m1.mean.rmse))
print(sprintf("RMSE for training data predicting means across mode1: %.3f", m1.mean.rmse))
print(sprintf("RMSE for training data predicting means across mode1&2: %.3f", m1.mean.rmse))

print(sprintf("RMSE for test mode1: %.3f", m1.mean.rmse))
print(sprintf('RMSE for test mode2: %.3f', m2.mean.rmse))
print(sprintf('RMSE for testing both mode1 and mode2: %.3f', m1m2.mean.rmse))

# Reshape responses into dataframes with the mode1 name,
# mode2 name and response.
resp.melt <- melt(resp.train, varnames=c('m1', 'm2', 'm3'))
train.resp.df <- dcast(resp.melt, m1 + m2 ~ m3)

resp.melt <- melt(all.resp.train, varnames=c('m1', 'm2', 'm3'))
all.train.resp.df <- dcast(resp.melt, m1 + m2 ~ m3)

if(dim(resp.test.m1)[[1]]) {
  resp.melt <- melt(resp.test.m1, varnames=c('m1', 'm2', 'm3'))
  m1.test.resp.df <- dcast(resp.melt, m1 + m2 ~ m3)
}

if(dim(resp.test.m2)[[2]]) {
  resp.melt <- melt(resp.test.m2, varnames=c('m1', 'm2', 'm3'))
  m2.test.resp.df <- dcast(resp.melt, m1 + m2 ~ m3)
}

if(dim(resp.test.m1m2)[[1]] & dim(resp.test.m1m2)[[2]]) {
  resp.melt <- melt(resp.test.m1m2, varnames=c('m1', 'm2', 'm3'))
  m1m2.test.resp.df <- dcast(resp.melt, m1 + m2 ~ m3)
}

# predictors for both mode1 and mode2 are used to predict each combination
input.train <- matrix(0, nrow=nrow(train.resp.df), ncol=P+Q,
  dimnames = list(rownames(train.resp.df), 
  	          c(colnames(train.m1.mat), colnames(train.m2.mat))))
input.train[,1:P] <- train.m1.mat[match(train.resp.df$m1, row.names(train.m1.mat)),]
input.train[,(P+1):(P+Q)] <- train.m2.mat[match(train.resp.df$m2, row.names(train.m2.mat)),]

if(P & dim(resp.test.m1)[[1]]) {
  input.test.m1 <- matrix(0, nrow=nrow(m1.test.resp.df), ncol=P+Q,
    dimnames = list(rownames(m1.test.resp.df),
                    c(colnames(train.m1.mat), colnames(train.m2.mat))))
  input.test.m1[,1:P] <- test.m1.mat[match(m1.test.resp.df$m1, row.names(test.m1.mat)),]
  if(Q & dim(resp.test.m2)[[2]])
    input.test.m1[,(P+1):(P+Q)] <- train.m2.mat[match(m1.test.resp.df$m2, row.names(train.m2.mat)),]
}

if(Q & dim(resp.test.m2)[[2]]) {
  input.test.m2 <- matrix(0, nrow=nrow(m2.test.resp.df), ncol=P+Q,
    dimnames = list(rownames(m2.test.resp.df),
                    c(colnames(train.m1.mat), colnames(train.m2.mat))))
  if(P & dim(resp.test.m1)[[1]])
    input.test.m2[,1:P] <- train.m1.mat[match(m2.test.resp.df$m1, row.names(train.m1.mat)),]
  input.test.m2[,(P+1):(P+Q)] <- test.m2.mat[match(m2.test.resp.df$m2, row.names(test.m2.mat)),]
}

if(P & dim(resp.test.m1)[[1]] & Q & dim(resp.test.m2)[[2]]) {
  input.test.m1m2 <- matrix(0, nrow=nrow(m1m2.test.resp.df), ncol=P+Q,
    dimnames = list(rownames(m1m2.test.resp.df),
                    c(colnames(train.m1.mat), colnames(train.m2.mat))))
  input.test.m1m2[,1:P] <- test.m1.mat[match(m1m2.test.resp.df$m1, row.names(test.m1.mat)),]
  input.test.m1m2[,(P+1):(P+Q)] <- test.m2.mat[match(m1m2.test.resp.df$m2, row.names(test.m2.mat)),]
}

# Preallocate storage
train.preds <- list()
warm.test.preds <- list()
if(P & dim(resp.test.m1)[[1]])
  m1.test.preds <- matrix(0, nrow(input.test.m1), K)
if(Q & dim(resp.test.m2)[[2]])
  m2.test.preds <- matrix(0, nrow(input.test.m2), K)
if(P & dim(resp.test.m1)[[1]] & Q & dim(resp.test.m2)[[2]])
  m1m2.test.preds <- matrix(0, nrow(input.test.m1m2), K)
  models <- list()  # Don't save models for space?

# Store the coefficients weights for prediction at each value of mode3
coefs <- matrix(0, P+Q, K,
  dimnames=list(c(colnames(train.m1.mat), colnames(train.m2.mat)),
                  dimnames(resp.tens)[[3]]))

# Construct a vector of the number of nodes per hidden layer
# and the dropout ratios for each hidden layer
hid_vec=rep(hidden, layers)
hid_drop_vec=rep(hidden_dropout, layers)

# Start up h2o instance
localH2O = h2o.init(nthreads=cores)

# Run a separate neural net for each mode3 sample
for(k in 1:K) {
  # Get response vector for combinations w/ responses
  train.resp <- train.resp.df[!is.na(train.resp.df[,k+2]),k+2]
  train.input <- input.train[!is.na(train.resp.df[,k+2]),]
  warm.test.input <- input.train[is.na(train.resp.df[,k+2]),]

  # Paste together responses and predictors
  tmp.train <- cbind(train.resp, train.input)

  # Run the random forest model
  model <- h2o.deeplearning(y = 1,
                            training_frame=as.h2o(tmp.train),
                            activation = activation,
                            hidden = hid_vec,
                            input_dropout_ratio = input_dropout,
                            hidden_dropout_ratios = hid_drop_vec,
                            stopping_metric=stop_metric,
                            stopping_tolerance=stop_toler,
                            stopping_rounds=stop_rounds,
                            epochs = epochs,
                            nfolds = cv_folds,
                            variable_importances=T)

  # Store the weights of predictors used for each mode3
  infl <- model@model$variable_importances$scaled_importance
  names(infl) <- model@model$variable_importances$variable
  infl <- infl[match(colnames(train.input), names(infl), nomatch=0)]
  coefs[match(names(infl), rownames(coefs), nomatch=0),k]  <- infl

  # Predict 
  train.preds[[k]] <- predict(model, tmp.train[,-1])
  warm.test.preds[[k]] <- predict(model, warm.test.input)
  if(P & dim(resp.test.m1)[[1]])
    m1.test.preds[,k]   <- predict(model, input.test.m1)
  if(Q & dim(resp.test.m2)[[2]])
    m2.test.preds[,k]   <- predict(model, input.test.m2)
  if(P & dim(resp.test.m1)[[1]] & Q & dim(resp.test.m2)[[2]])
    m1m2.test.preds[,k] <- predict(model, input.test.m1m2)

  models[[k]] <- model
  # rm(fit) # Don't save models for space?
}

h2o.shutdown(prompt=F)

save.image(paste0(out.dir, 'image.Rdata'))

warm.test.resp <- list()
train.exp.var.mode3 <- rep(0, K)
warm.exp.var.mode3 <- rep(0, K); m1.exp.var.mode3 <- rep(0, K)
m2.exp.var.mode3 <- rep(0, K); m1m2.exp.var.mode3 <- rep(0, K)
train.p.cor.mode3 <- rep(0, K)
warm.p.cor.mode3 <- rep(0, K); m1.p.cor.mode3 <- rep(0, K)
m2.p.cor.mode3 <- rep(0, K); m1m2.p.cor.mode3 <- rep(0, K)
train.s.cor.mode3 <- rep(0, K)
warm.s.cor.mode3 <- rep(0, K); m1.s.cor.mode3 <- rep(0, K)
m2.s.cor.mode3 <- rep(0, K); m1m2.s.cor.mode3 <- rep(0, K)
for(k in 1:K) {
  warm.test.resp[[k]] <- all.train.resp.df[is.na(train.resp.df[,k+2]) ,k+2]

  train.exp.var.mode3[k] <- exp_var(train.resp.df[,k+2][!is.na(train.resp.df[,k+2])], train.preds[[k]])
  warm.exp.var.mode3[k] <- exp_var(warm.test.resp[[k]], warm.test.preds[[k]])
  if(P & dim(resp.test.m1)[[1]])
    m1.exp.var.mode3[k]   <- exp_var(m1.test.resp.df[,k+2], m1.test.preds[,k])
  if(Q & dim(resp.test.m2)[[2]])
    m2.exp.var.mode3[k]   <- exp_var(m2.test.resp.df[,k+2], m2.test.preds[,k])
  if(P & dim(resp.test.m1)[[1]] & Q & dim(resp.test.m2)[[2]])
    m1m2.exp.var.mode3[k] <- exp_var(m1m2.test.resp.df[,k+2], m1m2.test.preds[,k])

  train.p.cor.mode3[k] <- p_cor(train.resp.df[,k+2][!is.na(train.resp.df[,k+2])], train.preds[[k]])
  warm.p.cor.mode3[k] <- p_cor(warm.test.resp[[k]], warm.test.preds[[k]])
  if(P & dim(resp.test.m1)[[1]])
    m1.p.cor.mode3[k]   <- p_cor(m1.test.resp.df[,k+2], m1.test.preds[,k])
  if(Q & dim(resp.test.m2)[[2]])
    m2.p.cor.mode3[k]   <- p_cor(m2.test.resp.df[,k+2], m2.test.preds[,k])
  if(P & dim(resp.test.m1)[[1]] & Q & dim(resp.test.m2)[[2]])
    m1m2.p.cor.mode3[k] <- p_cor(m1m2.test.resp.df[,k+2], m1m2.test.preds[,k])

  train.s.cor.mode3[k] <- s_cor(train.resp.df[,k+2][!is.na(train.resp.df[,k+2])], train.preds[[k]])
  warm.s.cor.mode3[k] <- s_cor(warm.test.resp[[k]], warm.test.preds[[k]])
  if(P & dim(resp.test.m1)[[1]])
    m1.s.cor.mode3[k]   <- s_cor(m1.test.resp.df[,k+2], m1.test.preds[,k])
  if(Q & dim(resp.test.m2)[[2]])
    m2.s.cor.mode3[k]   <- s_cor(m2.test.resp.df[,k+2], m2.test.preds[,k])
  if(P & dim(resp.test.m1)[[1]] & Q & dim(resp.test.m2)[[2]])
    m1m2.s.cor.mode3[k] <- s_cor(m1m2.test.resp.df[,k+2], m1m2.test.preds[,k])
}

# Get the explained variances predictions
train.exp.var <- exp_var(train.resp.df[,3:(K+2)][!is.na(train.resp.df[,3:(K+2)])],
                             unlist(train.preds))
warm.exp.var <- exp_var(unlist(warm.test.resp), unlist(warm.test.preds))
if(P & dim(resp.test.m1)[[1]])
  m1.exp.var <- exp_var(m1.test.resp.df[,3:(K+2)], m1.test.preds)
if(Q & dim(resp.test.m2)[[2]])
  m2.exp.var <- exp_var(m2.test.resp.df[,3:(K+2)], m2.test.preds)
if(P & dim(resp.test.m1)[[1]] & Q & dim(resp.test.m2)[[2]])
  m1m2.exp.var <- exp_var(m1m2.test.resp.df[,3:(K+2)], m1m2.test.preds)
 
# Get correlations
train.p.cor <- p_cor(train.resp.df[,3:(K+2)][!is.na(train.resp.df[,3:(K+2)])],
                             unlist(train.preds))
warm.p.cor <- p_cor(unlist(warm.test.resp), unlist(warm.test.preds))
if(P & dim(resp.test.m1)[[1]])
  m1.p.cor <- p_cor(as.vector(unlist(m1.test.resp.df[,3:(K+2)])), as.vector(m1.test.preds))
if(Q & dim(resp.test.m2)[[2]])
  m2.p.cor <- p_cor(as.vector(unlist(m2.test.resp.df[,3:(K+2)])), as.vector(m2.test.preds))
if(P & dim(resp.test.m1)[[1]] & Q & dim(resp.test.m2)[[2]])
  m1m2.p.cor <- p_cor(as.vector(unlist(m1m2.test.resp.df[,3:(K+2)])), as.vector(m1m2.test.preds))

train.s.cor <- s_cor(train.resp.df[,3:(K+2)][!is.na(train.resp.df[,3:(K+2)])],
                             unlist(train.preds))
warm.s.cor <- s_cor(unlist(warm.test.resp), unlist(warm.test.preds))
if(P & dim(resp.test.m1)[[1]])
  m1.s.cor <- s_cor(as.vector(unlist(m1.test.resp.df[,3:(K+2)])), as.vector(m1.test.preds))
if(Q & dim(resp.test.m2)[[2]])
  m2.s.cor <- s_cor(as.vector(unlist(m2.test.resp.df[,3:(K+2)])), as.vector(m2.test.preds))
if(P & dim(resp.test.m1)[[1]] & Q & dim(resp.test.m2)[[2]])
  m1m2.s.cor <- s_cor(as.vector(unlist(m1m2.test.resp.df[,3:(K+2)])), as.vector(m1m2.test.preds))

# Print results
print('################## Explained variance ######################')
print(sprintf('Explained variance for training data: %.4f', train.exp.var))
print('Explained variance for training data by mode3'); print(train.exp.var.mode3)
print(sprintf('Explained variance for warm start data: %.4f', warm.exp.var))
print('Explained variance for warm start data by mode3'); print(warm.exp.var.mode3)
if(P & dim(resp.test.m1)[[1]]) {
  print(sprintf('Explained variance for for cold start mode1: %.4f', m1.exp.var))
  print('Explained variance for cold start mode1 by mode3'); print(m1.exp.var.mode3)
}
if(Q & dim(resp.test.m2)[[2]]) {
  print(sprintf('Explained variance for cold start mode2: %.4f', m2.exp.var))
  print('Explained variance for cold start mode2 by mode3'); print(m2.exp.var.mode3)
}
if(P & dim(resp.test.m1)[[1]] & Q & dim(resp.test.m2)[[2]]) {
  print(sprintf('Explained variance for cold start mode1 and mode2: %.4f', m1m2.exp.var))
  print('Explained variance for cold start mode1 and mode2 by mode3'); print(m1m2.exp.var.mode3)
}

print('################## Pearson correlation ######################')
print(sprintf('Pearson correlation for training data: %.4f', train.p.cor))
print('Pearson correlation for training data by mode3'); print(train.p.cor.mode3)

print(sprintf('Pearson correlation for warm start data: %.4f', warm.p.cor))
print('Pearson correlation for warm start data by mode3'); print(warm.p.cor.mode3)
if(P & dim(resp.test.m1)[[1]]) {
  print(sprintf('Pearson correlation for for cold start mode1: %.4f', m1.p.cor))
  print('Pearson correlation for cold start mode1 by mode3'); print(m1.p.cor.mode3)
}
if(Q & dim(resp.test.m2)[[2]]) {
  print(sprintf('Pearson correlation for cold start mode2: %.4f', m2.p.cor))
  print('Pearson correlation for cold start mode2 by mode3'); print(m2.p.cor.mode3)
}
if(P & dim(resp.test.m1)[[1]] & Q & dim(resp.test.m2)[[2]]) {
  print(sprintf('Pearson correlation for cold start mode1 and mode2: %.4f', m1m2.p.cor))
  print('Pearson correlation for cold start mode1 and mode2 by mode3'); print(m1m2.p.cor.mode3)
}

print('################## Spearman correlation ######################')
print(sprintf('Spearman correlation for training data: %.4f', train.p.cor))
print('Spearman correlation for training data by mode3'); print(train.p.cor.mode3)

print(sprintf('Spearman correlation for warm start data: %.4f', warm.s.cor))
print('Spearman correlation for warm start data by mode3'); print(warm.s.cor.mode3)
if(P & dim(resp.test.m1)[[1]]) {
  print(sprintf('Spearman correlation for for cold start mode1: %.4f', m1.s.cor))
  print('Spearman correlation for cold start mode1 by mode3'); print(m1.s.cor.mode3)
}
if(Q & dim(resp.test.m2)[[2]]) {
  print(sprintf('Spearman correlation for cold start mode2: %.4f', m2.s.cor))
  print('Spearman correlation for cold start mode2 by mode3'); print(m2.s.cor.mode3)
}
if(P & dim(resp.test.m1)[[1]] & Q & dim(resp.test.m2)[[2]]) {
  print(sprintf('Spearman correlation for cold start mode1 and mode2: %.4f', m1m2.s.cor))
  print('Spearman correlation for cold start mode1 and mode2 by mode3'); print(m1m2.s.cor.mode3)
}

# Get RMSEs for all mode3 relative to the standard deviation of responses
train.rmse <- nrmse(train.resp.df[,3:(K+2)][!is.na(train.resp.df[,3:(K+2)])], unlist(train.preds))
warm.rmse <- nrmse(unlist(warm.test.resp), unlist(warm.test.preds))
if(P & dim(resp.test.m1)[[1]])
  m1.rmse <- nrmse(m1.test.preds, m1.test.resp.df[,3:(K+2)])
if(Q & dim(resp.test.m2)[[2]])
  m2.rmse <- nrmse(m2.test.preds, m2.test.resp.df[,3:(K+2)])
if(P & dim(resp.test.m1)[[1]] & Q & dim(resp.test.m2)[[2]])
  m1m2.rmse <- nrmse(m1m2.test.preds, m1m2.test.resp.df[,3:(K+2)])

# Get RMSEs for individual mode3 relative to the standard deviation of responses
train.rmse.mode3 <- rep(0, K)
warm.rmse.mode3 <- rep(0, K)
if(P & dim(resp.test.m1)[[1]])
  m1.rmse.mode3 <- rep(0, K)
if(Q & dim(resp.test.m2)[[2]])
  m2.rmse.mode3 <- rep(0, K)
if(P & dim(resp.test.m1)[[1]] & Q & dim(resp.test.m2)[[2]])
  m1m2.rmse.mode3 <- rep(0, K)
for(k in 1:K) {
  train.rmse.mode3[k] <- nrmse(train.resp.df[!is.na(train.resp.df[,k+2]),k+2], train.preds[[k]])
  warm.rmse.mode3[k] <- nrmse(warm.test.resp[[k]], warm.test.preds[[k]])
  if(P & dim(resp.test.m1)[[1]])
    m1.rmse.mode3[k] <- nrmse(m1.test.preds[,k], m1.test.resp.df[,k+2])
  if(Q & dim(resp.test.m2)[[2]])
    m2.rmse.mode3[k] <- nrmse(m2.test.preds[,k], m2.test.resp.df[,k+2])
  if(P & dim(resp.test.m1)[[1]] & Q & dim(resp.test.m2)[[2]])
    m1m2.rmse.mode3[k] <- nrmse(m1m2.test.preds[,k], m1m2.test.resp.df[,k+2])
}

print('############ RMSEs #############')
print(sprintf('Training RMSE relative to sd of training responses: %.4f', train.rmse))
print('Relative training RMSE by mode3'); print(train.rmse.mode3)
print(sprintf('Warm start RMSE relative to sd of training responses: %.4f', warm.rmse))
print('Relative warm start RMSE by mode3'); print(warm.rmse.mode3)
if(P & dim(resp.test.m1)[[1]]) {
  print(sprintf('Relative RMSE for cold start mode1: %.4f', m1.rmse))
  print('Relative RMSE for cold start mode1 by mode3'); print(m1.rmse.mode3)
}
if(Q & dim(resp.test.m2)[[2]]) {
  print(sprintf('Relative RMSE for cold start mode2: %.4f', m2.rmse))
  print('Relative RMSE for cold start mode2 by mode3'); print(m2.rmse.mode3)
}
if(P & dim(resp.test.m1)[[1]] & Q & dim(resp.test.m2)[[2]]) {
  print(sprintf('Relative RMSE for cold start mode1 and mode2: %.4f', m1m2.rmse))
  print('Relative RMSE for cold start mode1 and mode2 by mode3'); print(m1m2.rmse.mode3)
}

# Plots 
###############################################

# Plot heatmap of the coefficients for each mode3
pdf(paste0(out.dir, 'coef_mat_', fold, '.pdf'))
im_mat(coefs)
dev.off()

# Barplot of the number of predictors used for each mode3
pdf(paste0(out.dir, 'preds_by_mode3_', fold, '.pdf'))
n.preds <- apply(coefs, 2, function(x) sum(x!=0))
names(n.preds) <- 1:K
barplot(n.preds, xlab='Mode3', main='Predictors used')
dev.off()

# Barplots of the explained variance per mode3
pdf(paste0(out.dir, 'exp_var_barplots', fold, '.pdf'), height=14, width=14)
par(mfrow=c(3,2))
barplot(train.exp.var.mode3, xlab='Mode3', main='Explained variance for training data')
barplot(warm.exp.var.mode3, xlab='Mode3', main='Explained variance for warm test data')
barplot(m1.exp.var.mode3, xlab='Mode3', main='Explained variance for cold test mode1')
barplot(m2.exp.var.mode3, xlab='Mode3', main='Explained variance for cold test mode2')
barplot(m1m2.exp.var.mode3, xlab='Mode3', main='Explained variance for cold test \nmode1 and mode2')
dev.off()

# Scatter plots of predicted vs. observed responses
pdf(paste0(out.dir, "pred_plots_", fold, ".pdf"), height=28, width=28)
par(mfrow=c(4,4))
# Plot predictions for training data
for(k in 1:K) {
  plot(train.resp.df[!is.na(train.resp.df[,k+2]),k+2], train.preds[[k]],
       pch=20, cex=0.1, col=rgb(0,0,0,0.1), #xlim=c(0, 1.75),
       main=sprintf("Training predictions for mode3 %.0f \n RMSE = %.4f",
                    k, train.rmse[k]),
       xlab="Observed responses",
       ylab="Predicted responses")
  points(warm.test.resp[[k]], warm.test.preds[[k]], pch=20, cex=0.1, col='red')
  if(P & dim(resp.test.m1)[[1]])
    points(m1.test.resp.df[,k+2], m1.test.preds[,k], pch=20, cex=0.1, col='green')
  if(Q & dim(resp.test.m2)[[2]])
    points(m2.test.resp.df[,k+2], m2.test.preds[,k], pch=20, cex=0.1, col='blue')
  if(P & dim(resp.test.m1)[[1]] & Q & dim(resp.test.m2)[[2]])
    points(m1m2.test.resp.df[,k+2], m1m2.test.preds[,k], col='purple')
  abline(a=0,b=1, lwd=2)
  legend('bottomright', legend=c('train', 'warm', 'mode1', 'mode2', 'm1 and m2'),
         col=c('black', 'red', 'green', 'blue', 'purple'),
         pch=c(rep(20,4), 1))
}
dev.off()

# Summarize results
###################

print('########################################################')

# Print the number of mode1 and mode2 predictors used at each mode3
preds.df <- data.frame(mode3=1:K, m1=0, m2=0)
preds.df$m1 <- apply(coefs, 2, function(x) sum(x[2:(P+1)]!=0))
preds.df$m2 <- apply(coefs, 2, function(x) sum(x[(P+2):(P+Q+1)]!=0))
preds.df$m1.per <- round(preds.df$m1/P * 100, 2)
preds.df$m2.per <- round(preds.df$m2/Q * 100, 2)
print('Number and percentage of predictors used for each mode3:')
print(preds.df)

# Count total number of predictors chosen
shared.preds <- list()
shared.pred.counts <- apply(coefs, 1, function(x) sum(x!=0))
for(k in 1:K) shared.preds[[k]] <- names(shared.pred.counts)[shared.pred.counts>k]

print("Number of variables used for prediction in at least n mode3:")
n.shared.preds <- sapply(shared.preds, length)
names(n.shared.preds) <- 1:K
print(n.shared.preds)
print("Percentage of predictors shared out of the total number of predictors")
print(round(n.shared.preds/(P+Q+1) * 100, 2))

print(sprintf("%d (%.3f%%) of mode1 predictors were used at any mode3",
              sum(shared.pred.counts[1:P]!=0), 
              sum(shared.pred.counts[1:P]!=0)/P * 100)) 

print(sprintf("%d (%.3f%%) of mode2 predictors were used at any mode3",
              sum(shared.pred.counts[(P+1):(P+Q)]!=0), 
              sum(shared.pred.counts[(P+1):(P+Q)]!=0)/Q * 100)) 

# Barplot of the number of shared predictors
pdf(paste0(out.dir, 'shared_preds_barplot', fold, '.pdf'))
barplot(n.preds, xlab='Mode3', main='Shared predictors')
dev.off()

# Get RMSEs for each mode1 for predicting the mean or using glmnet
# m1.mean.diff.df <- m1.test.resp.df
# m1.pred.diff.df <- m1.test.resp.df
# m1.mean.diff.df[,3:(K+2)] <- m1.test.resp.df[,3:(K+2)] -
#   m2.means[match(m1.test.resp.df$m2, row.names(m2.means)),]
# m1.pred.diff.df[,3:(K+2)] <- m1.test.resp.df[,3:(K+2)] - m1.test.preds
# by.m1.rmse <- data.frame(m1=unique(m1.mean.diff.df$m1), mean=0, pred=0)
# for(i in 1:nrow(by.m1.rmse)) {
#   by.m1.rmse$mean[i] <- nrmse(m1.mean.diff.df[m1.mean.diff.df$m1==by.m1.rmse$m1[i],3:(K+2)])
#   by.m1.rmse$pred[i] <- nrmse(m1.pred.diff.df[m1.pred.diff.df$m1==by.m1.rmse$m1[i],
#                                                    3:(K+2)])^2, na.rm=T))
# }
# print('RMSEs for left out mode1')
# print(by.m1.rmse)
# print(sprintf("Random forest predicts %d of %d (%.2f%%) left out mode1 better than just predicting the mean response",
#               sum(by.m1.rmse$pred < by.m1.rmse$mean), nrow(by.m1.rmse),
#               sum(by.m1.rmse$pred < by.m1.rmse$mean)/nrow(by.m1.rmse) * 100))
# 
# # Get RMSEs for each mode2 for predicting the mean or using glmnet
# if(Q & dim(resp.test.m2)[[2]]) {
#   m2.mean.diff.df <- m2.test.resp.df
#   m2.pred.diff.df <- m2.test.resp.df
#   m2.mean.diff.df[,3:(K+2)] <- m2.test.resp.df[,3:(K+2)] -
#     m1.means[match(m2.test.resp.df$m1, row.names(m1.means)),]
#   m2.pred.diff.df[,3:(K+2)] <- m2.test.resp.df[,3:(K+2)] - m2.test.preds
#   by.m2.rmse <- data.frame(m2=unique(m2.mean.diff.df$m2), mean=0, pred=0)
#   for(i in 1:nrow(by.m2.rmse)) {
#     by.m2.rmse$mean[i] <- nrmse(m2.mean.diff.df[m2.mean.diff.df$m2==by.m2.rmse$m2[i],
#                                                      3:(K+2)])^2, na.rm=T))
#     by.m2.rmse$pred[i] <- nrmse(m2.pred.diff.df[m2.pred.diff.df$m2==by.m2.rmse$m2[i],
#                                                      3:(K+2)])^2, na.rm=T))
#   }
#   print('RMSEs for left out mode2')
#   print(by.m2.rmse)
#   print(sprintf("Random forest predicts %d of %d (%.2f%%) left out mode2 better than just predicting the mean response",
#                 sum(by.m2.rmse$pred < by.m2.rmse$mean, na.rm=T), nrow(by.m2.rmse),
#                 sum(by.m2.rmse$pred < by.m2.rmse$mean, na.rm=T)/nrow(by.m2.rmse) * 100))
# }
# 
# # Get RMSEs for left out mode1/mode2 combos for predicting the mean and using glmnet
# if(P & dim(resp.test.m1)[[1]] & Q & dim(resp.test.m2)[[2]]) {
#   m1m2.mean.diff.df <- m1m2.test.resp.df
#   m1m2.pred.diff.df <- m1m2.test.resp.df
#   m1m2.mean.diff.df[,3:(K+2)] <- m1m2.test.resp.df[,3:(K+2)] -
#     matrix(m1m2.means, nrow(m1m2.mean.diff.df), K, byrow=T)
#   m1m2.pred.diff.df[,3:(K+2)] <- m1m2.test.resp.df[,3:(K+2)] - m1m2.test.preds
#   by.m1m2.rmse <- m1m2.pred.diff.df[,1:2]
#   by.m1m2.rmse$mean <- apply(m1m2.mean.diff.df[,3:(K+2)], 1, function(x) nrmsex, na.rm=T)))
#   by.m1m2.rmse$pred <- apply(m1m2.pred.diff.df[,3:(K+2)], 1, function(x) nrmsex, na.rm=T)))
#   by.m1.m1m2.rmse <- data.frame(m1=unique(m1m2.mean.diff.df$m1), mean=0, pred=0)
#   by.m2.m1m2.rmse <- data.frame(m2=unique(m1m2.mean.diff.df$m2), mean=0, pred=0)
#   for(i in 1:nrow(by.m1.m1m2.rmse)) {
#     by.m1.m1m2.rmse$mean[i] <- nrmse(m1m2.mean.diff.df[m1m2.mean.diff.df$m1==by.m1.m1m2.rmse$m1[i],
#                                                      3:(K+2)])^2, na.rm=T))
#     by.m1.m1m2.rmse$pred[i] <- nrmse(m1m2.pred.diff.df[m1m2.pred.diff.df$m1==by.m1.m1m2.rmse$m1[i],
#                                                      3:(K+2)])^2, na.rm=T))
#   }
#   for(i in 1:nrow(by.m2.m1m2.rmse)) {
#     by.m2.m1m2.rmse$mean[i] <- nrmse(m1m2.mean.diff.df[m1m2.mean.diff.df$m2==by.m2.m1m2.rmse$m2[i],
#                                                             3:(K+2)])^2, na.rm=T))
#     by.m2.m1m2.rmse$pred[i] <- nrmse(m1m2.pred.diff.df[m1m2.pred.diff.df$m2==by.m2.m1m2.rmse$m2[i],
#                                                             3:(K+2)])^2, na.rm=T))
#   }
#   print('RMSEs for left out mode1/mode2 combinations')
#   # print(by.m1m2.rmse)
#   print(by.m1.m1m2.rmse)
#   print(by.m2.m1m2.rmse)
#   print(sprintf("Random forest predicts %d of %d (%.2f%%) left out mode1/mode2 combinations better than just predicting the mean response",
#                 sum(by.m1m2.rmse$pred < by.m1m2.rmse$mean, na.rm=T), nrow(by.m1m2.rmse),
#                 sum(by.m1m2.rmse$pred < by.m1m2.rmse$mean, na.rm=T)/nrow(by.m1m2.rmse) * 100))
# 
#   print(sprintf("Random forest predicts %d of %d (%.2f%%) left out mode1 better than just predicting the mean response on mode2 that were not included in the training set",
#                 sum(by.m1.m1m2.rmse$pred < by.m1.m1m2.rmse$mean, na.rm=T), nrow(by.m1.m1m2.rmse),
#                 sum(by.m1.m1m2.rmse$pred < by.m1.m1m2.rmse$mean, na.rm=T)/nrow(by.m1.m1m2.rmse) * 100))
# 
# 
#   print(sprintf("Random forest predicts %d of %d (%.2f%%) left out mode2 better than just predicting the mean response on mode1 that were not included in the training set",
#                 sum(by.m2.m1m2.rmse$pred < by.m2.m1m2.rmse$mean, na.rm=T), nrow(by.m2.m1m2.rmse),
#                 sum(by.m2.m1m2.rmse$pred < by.m2.m1m2.rmse$mean, na.rm=T)/nrow(by.m2.m1m2.rmse) * 100))
# }

save.image(paste0(out.dir, 'image.Rdata'))

