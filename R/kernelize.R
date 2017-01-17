#' Transform a matrix of input data into a matrix of concatenated kernel matrices
#' 
#' The input matrix should have samples as the rows and features as columns. A kernel will 
#' computed across all training data for each unique column category. The column names
#' should begin with the category. E.g. for a cell line input matrix the column names
#' may be snp.<gene>, meth.<loci>, etc. and this function will produce similarity kernels
#' for all snp data separately from the methylation data.
#' 
#' @export
#' @param m matrix on which to compute kernels
#' @param sigma2 numeric bandwidth (variance) for Gaussian kernels (default: 100)
#' 
#' @return m.kern matrix with the kernels concatenated column-wise
#' 
#' @examples
#' m <- matrix(rnorm(800), 8, 100, dimnames=list(paste0('sample.', 1:8), 
#'                           c(paste0('rna.', 1:20), paste0('meth.', 1:30), paste0('mut.', 1:50))))
#' m[,51:100] <- rbinom(8*50, 1, 0.5)
#' m.kern <- kernelize(m)

kernelize <- function(m, sigma2=100) {
  if(ncol(m)==0) return(m)
  
  col.names <- colnames(m)
  cat <- unique(sapply(strsplit(col.names, '.', fixed=T), '[', 1))
  kerns <- list()
  
  for(i in 1:length(cat)) {
    sub <- m[,grepl(paste0('^', cat[[i]]), colnames(m))]
    # If the submatrix is binary, compute Jaccard kernel
    # otherwise compute the Gaussian kernel Pearson correlations
    if(min(sub, na.rm=T)==0 & max(sub, na.rm=T)==1 & length(unique(as.vector(sub)))==2) {
      # compute Jaccard kernel
      kerns[[i]] <- matrix(NA, nrow(m), nrow(m), dimnames=list(row.names(m),
                                               paste0(cat[[i]], '.', rownames(m))))
      for(j in 1:nrow(m)) for (k in 1:nrow(m)) {
        if(j == k) {
          kerns[[i]][j,k] <- 1 
        } else {
          union <- sum(sub[j,] | sub[k,])
          if(union != 0) {
            kerns[[i]][j,k] <- sum(sub[j,] & sub[k,])/union
          } else kerns[[i]][j,k] <- 0
        }
      }
    } else {
      kerns[[i]] <- exp(as.matrix(-(dist(sub, diag=T, upper=T)^2)/sigma2))
        #cor(t(sub))
      colnames(kerns[[i]]) <- paste0(cat[[i]], '.', rownames(m))
    }
  }
  # Combine the kernels
  m.kern <- matrix(NA, nrow(m), 0, dimnames=list(rownames(m), c()))
  for(i in 1:length(kerns))
    m.kern <- cbind(m.kern, kerns[[i]])
  return(m.kern)
}
