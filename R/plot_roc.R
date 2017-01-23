#' Plot reciever operating characteristic (ROC) curves for two projection (A) matrices
#' 
#' This is a little different than a typical ROC curve since any rows of the true matrix
#' that are non-zero are treated as equal true positives.
#' 
#' @importFrom graphics plot abline
#' 
#' @export
#' @param true projection matrix where rows of true predictors have non-zero values
#' @param pred projection matrix where rows of learned predictors have larger values
#' @param main title of the ROC plot

plot_roc <- function(true, pred, main=character(0)) {
  # Plots reciever operating characteristic curves given a 
  # known projection matrix and a learned projection matrix.
  
  if(is.null(rownames(pred)) & nrow(true) == nrow(pred))
    rownames(pred) <- rownames(true)
  
  true.ss <- apply(true^2, 1, sum)
  pred.ss <- apply(pred^2, 1, sum)
  true.pos <- names(which(true.ss > min(true.ss)))
  true.neg <- names(which(true.ss == min(true.ss)))

    # Scale columns of both matrices to avoid scaling problems between latent dimensions
  true <- scale(true)
  pred <- scale(pred)
  
  cuts <- c(exp(seq(log(max(pred.ss)), -50, -.005)),0) + min(pred.ss)
  n <- length(cuts)
  tpr <- rep(0, n)
  fpr <- rep(0, n)
  
  for(i in 1:n) {
    pos <- names(which(pred.ss >= cuts[i]))
    tpr[i] <- mean(true.pos %in% pos)
    fpr[i] <- mean(true.neg %in% pos)
  }
  
  # Add ones on the end for predictors that were removed
  tpr <- c(tpr, 1)
  fpr <- c(fpr, 1)
  
  auroc <- sum((fpr[-1]-fpr[-(n+1)]) * (tpr[-1] + tpr[-(n+1)]) /2)
  
  plot(fpr, tpr, type = 'l', xlim=c(0,1), ylim=c(0,1), 
       xlab="False positive rate", ylab="True positive rate",
       main=paste(main, 'AUROC =', round(auroc, 3)))
  
  abline(a=0,b=1, lty=3, col='red')
  return(auroc)
}
