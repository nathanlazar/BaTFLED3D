#' Plot heatmaps of two matrices in red and blue
#' 
#' Displays two heatmaps of matrices using red and blue colors. Options to scale
#' and sort as well as any other graphical parameters with ... Sorting attempts to match
#' columns between the two matrices using their correlation over rows. If \code{sort==TRUE}
#' then the new ordering for the second matrix is returned.
#' 
#' @export
#' @param x1 matrix
#' @param x2 matrix
#' @param high string of either 'red' or 'blue' used to show higher values
#' @param xaxt string indicating how to display the x axis. Suppress x axis with 'n'
#' @param yaxt string indicating how to display the y axis. Suppress y axis with 'n'
#' @param scale logical indicating whether the matrices should be z scaled to have 
#' columns with norm zero and standard deviation one.
#' @param absol logical indicating whether to take absolute value of the entries 
#' before plotting
#' @param sort logical indicating whether the columns of the matrix should 
#' be sorted in decreasing order of their means
#' @param center logical indicating wether to center ranges for x and y around zero
#' @param main1 string to be used as the main title for the first matrix image
#' @param main2 string to be used as the main title for the second matrix image
#' @param ... other graphical parameters passed to image
#' 
#' @return If \code{sort==TRUE} the ordering of the second matrix used to match columns.
#' 
#' @examples
#' par(mfrow=c(1,2))
#' im_2_mat(matrix(1:12, nrow=3, ncol=4),  matrix(13:24, nrow=3, ncol=4), sort=FALSE, scale=FALSE)
#' im_2_mat(matrix(1:12, nrow=3, ncol=4),  matrix(13:24, nrow=3, ncol=4), sort=TRUE, scale=FALSE)
#' im_2_mat(matrix(1:12, nrow=3, ncol=4),  matrix(13:24, nrow=3, ncol=4), sort=TRUE, scale=TRUE)
#' im_2_mat(matrix(1:12, nrow=3, ncol=4),  matrix(13:24, nrow=3, ncol=4), sort=FALSE, 
#'          scale=FALSE, center=TRUE)

im_2_mat <- function(x1, x2, high='red', xaxt='n', yaxt='n', scale='col', 
                     absol=FALSE, sort=TRUE, center=FALSE, main1='', main2='', ...) {
  # display images of rotated matrices with the option to sort columns
  # by absolute magnitude, zlimits are set to the max and min of both matrices
  
  rwb <- colorRampPalette(c('red', 'white', 'blue'), space='rgb')
  
  # Remove constant columns if present
  if(sum(apply(x1, 2, function(x) all(x == rep(1,nrow(x1)))))>0) {
    x1 <- x1[,colnames(x1)!='const',drop=FALSE]
    const1 <-TRUE
  } else const1 <- FALSE
  if(sum(apply(x2, 2, function(x) all(x == rep(1,nrow(x2)))))>0) {
    x2 <- x2[,colnames(x2)!='const',drop=FALSE]
    const2 <- TRUE
  } else const2 <- FALSE
  
  if(scale=='col') {
    # Safe version of scale (if sd of columns = 0, set to minimum value)
    x1.col.means <- apply(x1, 2, mean)
    x1.col.sds <- apply(x1, 2, sd)
    x1.col.sds[x1.col.sds==0] <- 1e-300
    x1 <- sweep(x1, MARGIN=2, STATS=x1.col.means, FUN='-')
    x1 <- sweep(x1, MARGIN=2, STATS=x1.col.sds, FUN='/')

    x2.col.means <- apply(x2, 2, mean)
    x2.col.sds <- apply(x2, 2, sd)
    x2.col.sds[x2.col.sds==0] <- 1e-300
    x2 <- sweep(x2, MARGIN=2, STATS=x2.col.means, FUN='-')
    x2 <- sweep(x2, MARGIN=2, STATS=x2.col.sds, FUN='/')
  } else if(scale=='all') {
    x1 <- (x1 - mean(x1))/sd(x1)
    x2 <- (x2 - mean(x2))/sd(x2)
  }

  if(absol) {x1 <- abs(x1); x2 <- abs(x2)}
  
  # Rearrange columns 
  if(sort) {
    if(!all(dim(x1)==dim(x2))) {
      if(sum(grepl('const', rownames(x1))) + sum(grepl('const', rownames(x2)))==1) {
        x1 <- x1[!grepl('const', rownames(x1)),]
        x2 <- x2[!grepl('const', rownames(x1)),]
      } else {
        print('Can not sort, matrices are different sizes')
        sort <- FALSE
      }
    } else {
      # Make a correlation matrix between the columns of the two matrices
      cor.mat <- cor(x2,x1)
      reorder <- apply(abs(cor.mat), 2, which.max)
      reorder[which(duplicated(reorder))] <- (1:ncol(x1))[!(1:ncol(x1) %in% reorder)]
      
      x2 <- x2[,reorder]
      cor.mat <- cor(x2, x1)
      cor <- diag(cor.mat)
      for(i in 1:ncol(x2)) 
        if(cor[i] < 0) {
          x2[,i] <- -x2[,i]
          reorder[i] <- -reorder[i]
        }
    }
  }
  
  if(center) {
    zlim=c(-max(c(abs(x1), abs(x2))), max(c(abs(x1), abs(x2))))
  } else {
    zlim <- range(x1, x2, na.rm=TRUE)
  }
  
  # Add back in constant columns
  if(const1)
    x1 <- cbind(const=1, x1)
  if(const2) {
    x2 <- cbind(const=1, x2)
    if(sort) {
      reorder[reorder > 0] <- reorder[reorder > 0] + 1
      reorder[reorder < 0] <- reorder[reorder < 0] - 1
      reorder <- c(const=1, reorder)
    }
  }
  
  if(high=='blue') {
    image(rot(x1), col=rwb(256), zlim=zlim, xaxt=xaxt, yaxt=xaxt, main=main1, ...)
    image(rot(x2), col=rwb(256), zlim=zlim, xaxt=xaxt, yaxt=xaxt, main=main2, ...)
  } else if(high=='red') {
    image(rot(x1), col=rev(rwb(256)), zlim=zlim, xaxt=xaxt, yaxt=xaxt, main=main1, ...)
    image(rot(x2), col=rev(rwb(256)), zlim=zlim, xaxt=xaxt, yaxt=xaxt, main=main2, ...)
  }
  
  if(sort) {
    return(reorder)
  } else return(1:ncol(x2))
} 
