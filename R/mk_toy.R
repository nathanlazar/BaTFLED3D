#' Make a toy dataset to test the 3d BaTFLED model.
#' 
#' Returns a toy model with the specified size, sparsity and noise generated
#' either with a CP or Tucker factorization model.
#' Values in predictor matrices (X1, X2, X3) are pulled from a standard normal 
#' distribuion. Dummy names are given to the predictors.
#' 
#' @export
#' @param params list of parameters created with \code{get_data_params}
#' @return a list containing elements of the model
#' \describe{
#' \item{mode1.X}{Input data for mode 1}
#' \item{mode2.X}{Input data for mode 2}
#' \item{mode3.X}{Input data for mode 3}
#' \item{mode1.A}{Projection matrix for mode 1}
#' \item{mode2.A}{Projection matrix for mode 2}
#' \item{mode3.A}{Projection matrix for mode 3}
#' \item{mode1.H}{Latent matrix for mode 1}
#' \item{mode2.H}{Latent matrix for mode 2}
#' \item{mode3.H}{Latent matrix for mode 3}
#' \item{core}{Core tensor if \code{params$decomp=='Tucker'}}
#' \item{resp}{Response tensor}
#' }
#' @examples
#' data.params <- get_data_params(c('decomp=Tucker'))
#' toy <- mk_toy(data.params)
#' 
#' data.params <- get_data_params(c('decomp=CP'))
#' toy <- mk_toy(data.params)

mk_toy <- function(params) {
  
  # Make all param variables available locally
  for(i in 1:length(params)) assign(names(params)[i], params[i][[1]])

  # TODO: Put in checks on arguments (Errors)

  if(!(decomp %in% c('CP', 'Tucker'))) stop("decomp must be 'CP' or 'Tucker'")
  
  if(decomp=='CP') {R1 <- R; R2 <- R; R3 <- R}

  if(!is.na(seed)) set.seed(seed)
  toy <- list()

  # The input X matrices are pulled from standard normal dist or given no columns
  toy$mode1.X <- matrix(rnorm(m1.rows*m1.cols), nrow=m1.rows, ncol=m1.cols)    
  toy$mode2.X <- matrix(rnorm(m2.rows*m2.cols), nrow=m2.rows, ncol=m2.cols)    
  toy$mode3.X <- matrix(rnorm(m3.rows*m3.cols), nrow=m3.rows, ncol=m3.cols)
  
  # If A intercepts, add a column of 1s
  if(A1.intercept) toy$mode1.X <- cbind(1, toy$mode1.X)
  if(A2.intercept) toy$mode2.X <- cbind(1, toy$mode2.X)
  if(A3.intercept) toy$mode3.X <- cbind(1, toy$mode3.X)
  
  # Add dummy names
  rownames(toy$mode1.X) <- paste0('m1.samp', 1:m1.rows)
  rownames(toy$mode2.X) <- paste0('m2.samp', 1:m2.rows)
  rownames(toy$mode3.X) <- paste0('m3.samp', 1:m3.rows)
  if(m1.cols != 0) if(A1.intercept) {
    colnames(toy$mode1.X) <- c('const', paste0('m1.pred', 1:m1.cols))
  } else colnames(toy$mode1.X) <- paste0('m1.pred', 1:m1.cols)
  if(m2.cols != 0) if(A2.intercept) {
    colnames(toy$mode2.X) <- c('const', paste0('m2.pred', 1:m2.cols))
  } else colnames(toy$mode2.X) <- paste0('m2.pred', 1:m2.cols)
  if(m3.cols != 0) if(A3.intercept) {
    colnames(toy$mode3.X) <- c('const', paste0('m3.pred', 1:m3.cols))
  } else colnames(toy$mode3.X) <- paste0('m3.pred', 1:m3.cols)
  
  # Scale the input data by columns so that each feature has 
  # the same chance of being selected during training
  if(scale & nrow(toy$mode1.X)>1) if(A1.intercept) {
    toy$mode1.X[,-1] <- scale(toy$mode1.X[,-1])
  } else toy$mode1.X <- scale(toy$mode1.X)
  if(scale & nrow(toy$mode2.X)>1) if(A2.intercept) {
    toy$mode2.X[,-1] <- scale(toy$mode2.X[,-1])
  } else toy$mode2.X <- scale(toy$mode2.X)
  if(scale & nrow(toy$mode3.X)>1) if(A3.intercept) {
    toy$mode3.X[,-1] <- scale(toy$mode3.X[,-1])
  } else toy$mode3.X <- scale(toy$mode3.X)

  # Projection (A) matrices are sparse w/ non-zeros pulled from the standard normal distribution.
  toy$mode1.A <- matrix(0, nrow=ncol(toy$mode1.X), ncol=R1)
  toy$mode2.A <- matrix(0, nrow=ncol(toy$mode2.X), ncol=R2)
  toy$mode3.A <- matrix(0, nrow=ncol(toy$mode3.X), ncol=R3)
  
  if(!row.share) {
    # Make the A projection matrices sparse separately for each latent factor
    if(A1.intercept) {
      toy$mode1.A[1,sample(1:R1, A1.const.prob*R1)] <- 1
      for(r1 in 1:R1) toy$mode1.A[-1,][sample(1:m1.cols, m1.true),r1] <- 1
    } else for(r1 in 1:R1) toy$mode1.A[sample(1:m1.cols, m1.true),r1] <- 1

    if(A2.intercept) {
      toy$mode2.A[1,sample(1:R2, A2.const.prob*R2)] <- 1
      for(r2 in 1:R2) toy$mode2.A[-1,][sample(1:m2.cols, m2.true),r2] <- 1
    } else for(r2 in 1:R2) toy$mode2.A[sample(1:m2.cols, m2.true),r2] <- 1

    if(A3.intercept) {
      toy$mode3.A[1,sample(1:R3, A3.const.prob*R3)] <- 1
      for(r3 in 1:R3) toy$mode3.A[-1,][sample(1:m3.cols, m3.true),r3] <- 1
    } else for(r3 in 1:R3) toy$mode3.A[sample(1:m3.cols, m3.true),r3] <- 1
    
  } else {
    # Make the A projection matrices sparse row-wise to simulate a few genes or 
    # drug properties influencing the response.
    if(A1.intercept) {
      toy$mode1.A[1, sample(1:R1, A1.const.prob*R1)] <- 1
      toy$mode1.A[-1,][sample(1:m1.cols, m1.true),] <- 1
    } else toy$mode1.A[sample(1:m1.cols, m1.true),] <- 1
    if(A2.intercept) {
      toy$mode2.A[1, sample(1:R2, A2.const.prob*R2)] <- 1
      toy$mode2.A[-1,][sample(1:m2.cols, m2.true),] <- 1
    } else toy$mode2.A[sample(1:m2.cols, m2.true),] <- 1
    if(A3.intercept) {
      toy$mode3.A[1, sample(1:R3, A3.const.prob*R3)] <- 1
      toy$mode3.A[-1,][sample(1:m3.cols, m3.true),] <- 1
    } else toy$mode3.A[sample(1:m3.cols, m3.true),] <- 1
  } 
  
  # Multiply by randomly sampled matrices 
  toy$mode1.A <- toy$mode1.A * matrix(rnorm(prod(dim(toy$mode1.A))), nrow=nrow(toy$mode1.A), ncol=R1)
  rownames(toy$mode1.A) <- colnames(toy$mode1.X)
  toy$mode2.A <- toy$mode2.A * matrix(rnorm(prod(dim(toy$mode2.A))), nrow=nrow(toy$mode2.A), ncol=R2)
  rownames(toy$mode2.A) <- colnames(toy$mode2.X)
  toy$mode3.A <- toy$mode3.A * matrix(rnorm(prod(dim(toy$mode3.A))), nrow=nrow(toy$mode3.A), ncol=R3)
  rownames(toy$mode3.A) <- colnames(toy$mode3.X)
  
  colnames(toy$mode1.A) <- paste0('m1.lat', 1:R1)
  colnames(toy$mode2.A) <- paste0('m2.lat', 1:R2)
  colnames(toy$mode3.A) <- paste0('m3.lat', 1:R3)

  # Calculate the latent factor matrices (H) which multiply to generate
  # the response matrix if predictors are given. 
  # Otherwise, generate random H matrices
  if(m1.cols > 0){
    toy$mode1.H <- toy$mode1.X %*% toy$mode1.A
  } else {
    toy$mode1.H <- matrix(rnorm(m1.rows*R1), m1.rows, R1,
                          dimnames=list(rownames(toy$mode1.X), colnames(toy$mode1.A)))
  }
  if(m2.cols > 0){
    toy$mode2.H <- toy$mode2.X %*% toy$mode2.A
  } else {
    toy$mode2.H <- matrix(rnorm(m2.rows*R2), m2.rows, R2,
                          dimnames=list(rownames(toy$mode2.X), colnames(toy$mode2.A)))
  }
  if(m3.cols > 0){
    toy$mode3.H <- toy$mode3.X %*% toy$mode3.A
  } else {
    toy$mode3.H <- matrix(rnorm(m3.rows*R3), m3.rows, R3,
                          dimnames=list(rownames(toy$mode3.X), colnames(toy$mode3.A)))
  }
  
  # Make the core tensor
  if(decomp=='CP') {
    toy$core <- array(0, dim=c(ncol(toy$mode1.H), ncol(toy$mode2.H), ncol(toy$mode3.H)))
    for(i in 1:R) toy$core[i,i,i] <- 1
  } else if(decomp=='Tucker') {
    if(H1.intercept) toy$mode1.H <- cbind(const=1, toy$mode1.H)
    if(H2.intercept) toy$mode2.H <- cbind(const=1, toy$mode2.H)
    if(H3.intercept) toy$mode3.H <- cbind(const=1, toy$mode3.H)
    
    toy$core <- array(0, dim=c(ncol(toy$mode1.H), ncol(toy$mode2.H), ncol(toy$mode3.H)))
    dimnames(toy$core) <- list(colnames(toy$mode1.H), 
                               colnames(toy$mode2.H), 
                               colnames(toy$mode3.H))
    if(H1.intercept & H2.intercept & H3.intercept) toy$core[1,1,1] <- true.0D
    if(H2.intercept & H3.intercept) 
      toy$core[!grepl('const', colnames(toy$mode1.H)),1,1][sample(R1, true.1D.m1)] <- 1
    if(H1.intercept & H3.intercept) 
      toy$core[1,!grepl('const', colnames(toy$mode2.H)),1][sample(R2, true.1D.m2)] <- 1
    if(H1.intercept & H2.intercept) 
      toy$core[1,1,!grepl('const', colnames(toy$mode3.H))][sample(R3, true.1D.m3)] <- 1

    if(H3.intercept)
      toy$core[!grepl('const', colnames(toy$mode1.H)),
               !grepl('const', colnames(toy$mode2.H)),1][sample(R1*R2, true.2D.m1m2)] <- 1
    if(H2.intercept) 
      toy$core[!grepl('const', colnames(toy$mode1.H)),1,
               !grepl('const', colnames(toy$mode3.H))][sample(R1*R3, true.2D.m1m3)] <- 1
    if(H1.intercept) 
      toy$core[1,!grepl('const', colnames(toy$mode2.H)),
               !grepl('const', colnames(toy$mode3.H))][sample(R2*R3, true.2D.m2m3)] <- 1
    toy$core[!grepl('const', colnames(toy$mode1.H)),
             !grepl('const', colnames(toy$mode2.H)),
             !grepl('const', colnames(toy$mode3.H))][sample(R1*R2*R3, true.3D)] <- 1
    
    toy$core <- toy$core * array(rnorm(prod(dim(toy$core))), dim=dim(toy$core))
  }

  toy$resp <- mult_3d(toy$core, toy$mode1.H, toy$mode2.H, toy$mode3.H)
  if(noise.sd > 0) {
    toy$resp <- toy$resp + array(rnorm(prod(dim(toy$resp)), mean=0, 
                                       sd=noise.sd * sd(toy$resp)), 
                                 dim=dim(toy$resp))
  }
  
  return(toy)
}