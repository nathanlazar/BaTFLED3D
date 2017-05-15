# Testing out the idea of fitting a linear model to the responses 
# given the row and column values and subtracting this effect

library(BaTFLED3D)
library(reshape2) # Look into bypassing this
library(ggplot2)  # For plotting histograms

load("D:/Box Sync/Thesis/NewHeiserDataset/Responses/all_norm_tens.Rdata")
# Remove cell lines with no data
norm.tens <- norm.tens[apply(norm.tens, 1, function(x) sum(!is.na(x))!=0),,]

slice <- function(tens, dim, n) {
  # Slice a 3D tensor in dimension 'dim' at element 'n'
  if(dim==1) return(tens[n,,])
  if(dim==2) return(tens[,n,])
  if(dim==3) return(tens[,,n])
}

# Make tensors of fitted values and residuals for each mode slice
mode <- 1
tens <- norm.tens
resid.tens <- norm.tens
resid.tens[,,] <- NA
fit.tens <- resid.tens
for(n in 1:dim(tens)[[mode]]) {
  fit.mat <- slice(tens, mode, n)
  fit.mat[,] <- NA
  resid.mat <- fit.mat

  melted <- melt(slice(tens, mode, n))
  anov <- aov(value ~ Var1 + Var2, melted)
  melted$resid <- NA
  melted$resid[!is.na(melted$value)] <- anov$residuals

  # Reconstruct tensors using the fit
  coef <- anov$coefficients
  names(coef) <- sub('Var1', '', sub('Var2', '', names(coef)))
  for(m1 in rownames(fit.mat)) for(m2 in colnames(fit.mat)) 
    fit.mat[m1, m2] <- coef['(Intercept)'] + coef[m1] + coef[m2]
  fit.mat[1,1] <- coef['(Intercept)']
  fit.mat[1,colnames(fit.mat) %in% names(coef)] <- coef[names(coef)  %in% colnames(fit.mat)] + coef['(Intercept)']
  fit.mat[rownames(fit.mat) %in% names(coef),1] <- coef[names(coef) %in% rownames(fit.mat)] + coef['(Intercept)']
  
  for(i in 1:nrow(melted)) 
    resid.mat[melted$Var1[i], melted$Var2[i]] <- melted$resid[i]
  
  if(mode ==1) {
    fit.tens[n,,] <- fit.mat
    resid.tens[n,,] <- resid.mat
  }
  if(mode ==2) {
    fit.tens[,n,] <- fit.mat
    resid.tens[,n,] <- resid.mat
  }
  if(mode ==3) {
    fit.tens[,,n] <- fit.mat
    resid.tens[,,n] <- resid.mat
  }
}

# Look at slices of these tensors 
for(i in 1:dim(norm.tens)[[1]]) {
  im_mat(norm.tens[i,,], zlim=range(norm.tens, na.rm=T), main=i)
  im_mat(fit.tens[i,,], zlim=range(norm.tens, na.rm=T), main=i)
  im_mat(resid.tens[i,,], zlim=range(norm.tens, na.rm=T), main=i)
  im_mat(fit.tens[i,,] + resid.tens[i,,], zlim=range(norm.tens, na.rm=T), main=i)
}

for(j in 1:dim(norm.tens)[[2]]) {
  im_mat(norm.tens[,j,], zlim=range(norm.tens, na.rm=T), main=j)
  im_mat(fit.tens[,j,], zlim=range(norm.tens, na.rm=T), main=j)
  im_mat(resid.tens[,j,], zlim=range(norm.tens, na.rm=T), main=j)
  im_mat(fit.tens[,j,] + resid.tens[,j,], zlim=range(norm.tens, na.rm=T), main=j)
}

for(k in 1:dim(norm.tens)[[3]]) {
  im_mat(norm.tens[,,k], zlim=range(norm.tens, na.rm=T), main=k)
  im_mat(fit.tens[,,k], zlim=range(norm.tens, na.rm=T), main=k)
  im_mat(resid.tens[,,k], zlim=range(norm.tens, na.rm=T), main=k)
  im_mat(fit.tens[,,k] + resid.tens[,,k], zlim=range(norm.tens, na.rm=T), main=k)
}

# Look at histograms of data, fits and residuals
par(mfrow=c(3,1))
hist(norm.tens)
hist(fit.tens)
hist(resid.tens)

df <- melt(norm.tens)
df$fit <- melt(fit.tens)$value
df$resid <- melt(resid.tens)$value
df.melt <- melt(df)

ggplot(df.melt, aes(x=value, fill=variable)) + geom_density(alpha=0.25)
# ggplot(df.melt, aes(x=value, fill=variable)) + geom_histogram(alpha=0.25)
# ggplot(df.melt, aes(x=variable, y=value, fill=variable)) + geom_boxplot()

# Save tensors
save(fit.tens, file="D:/Box Sync/Thesis/NewHeiserDataset/Responses/norm_fit_tens.Rdata")
save(resid.tens, file="D:/Box Sync/Thesis/NewHeiserDataset/Responses/norm_resid_tens.Rdata")

####################################################################################
# Look at results from a test run using the residuals as the responses

load('Results/Testing_ANOVA/image.Rdata')
load('Responses/norm_fit_tens.Rdata')

train.fit.tens <- fit.tens[dimnames(fit.tens)[[1]] %in% dimnames(trained$resp)[[1]],
                           dimnames(fit.tens)[[2]] %in% dimnames(trained$resp)[[2]],
                           dimnames(fit.tens)[[3]] %in% dimnames(trained$resp)[[3]]]

test.m1.fit.tens <- fit.tens[dimnames(fit.tens)[[1]] %in% dimnames(test.data.m1$resp)[[1]],
                             dimnames(fit.tens)[[2]] %in% dimnames(test.data.m1$resp)[[2]],
                             dimnames(fit.tens)[[3]] %in% dimnames(test.data.m1$resp)[[3]]]

m1.resp.untrans <- m1.resp + test.m1.fit.tens
m1.mean.untrans <- mean.tens.list$m1 + test.m1.fit.tens
m1.pred.resp.untrans <- m1.pred.resp + test.m1.fit.tens

plot(m1.resp.untrans, m1.pred.resp.untrans)
abline(a=0,b=1, col='blue', lwd=2)
nrmse(m1.resp.untrans, m1.pred.resp.untrans)

plot(m1.resp.untrans, m1.mean.untrans)
abline(a=0,b=1, col='blue', lwd=2)
nrmse(m1.resp.untrans, m1.mean.untrans)


