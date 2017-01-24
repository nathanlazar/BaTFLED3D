#' Factorization object for 3D Tucker models.
#' 
#' \code{Tucker_model} objects are 'R6' objects so that their values can be updated in place. 
#' The object is treated like an environment and components are accessed using the \code{$} 
#' operator. 
#' When creating a new Tucker_model object it will be populated with default values and empty matrices.
#' To initialize a \code{Tucker_model} call the \code{initialize()} method. 
#' 
#' @docType class
#' @importFrom R6 R6Class
#' @exportClass Tucker_model
#' @export Tucker_model
#' @format An \code{\link{R6Class}} generator object
#' @section Members:
#' \describe{
#'   \item{iter}{Integer showing the number of iterations that have been run on this object.}
#'   \item{early.stop}{Stop if the lower bound increases by less than this value.}
#'   \item{lower.bnd}{Vector storing the lower bound values during training.}
#'   \item{RMSE}{Vector of the root mean squared error of the predictions during training.}
#'   \item{H.RMSE}{vector of the root mean squared error of predictions made by multiplying the H matrices.}
#'   \item{exp.var}{Vector of the explained variance of predictions during training.}
#'   \item{p.cor}{Vector of the Pearson correlation of predictions during training.}
#'   \item{s.cor}{Vector of the Spearman correlation of predictions during training.}
#'   \item{times}{Vector of the time taken for each update iteration.}
#'   \item{core.mean}{Mean parameters of the q Gaussian distributions in the core tensor.}
#'   \item{core.var}{Variance parameters of the q Gaussian distributions in the core tensor.}
#'   \item{core.lambda.shape}{Prior for the shape parameter of the gamma distribution on the core precision.}
#'   \item{core.lambda.scale}{Prior for the scale parameter of the gamma distribution on the core precision.}
#'   \item{resp}{array storing the predicted response tensor.}
#'   \item{delta}{binary array indicating whether the response is observed.}
#'   \item{core.var}{variance parameters of the q Gaussian distributions in the core tensor.}
#'   \item{m1Xm1X}{Product of mode1.X with itself stored to avoid recalculating.}
#'   \item{m2Xm2X}{Product of mode2.X with itself stored to avoid recalculating.}
#'   \item{m3Xm3X}{Product of mode3.X with itself stored to avoid recalculating.}
#'   \item{mode1.lambda.shape}{Matrix storing the shape parameters for the gamma distributions on the mode 1 projection (A) matrix.}
#'   \item{mode1.lambda.scale}{Matrix storing the scale parameters for the gamma distributions on the mode 1 projection (A) matrix.}
#'   \item{mode2.lambda.shape}{Matrix storing the shape parameters for the gamma distributions on the mode 2 projection (A) matrix.}
#'   \item{mode2.lambda.scale}{Matrix storing the scale parameters for the gamma distributions on the mode 2 projection (A) matrix.}
#'   \item{mode3.lambda.shape}{Matrix storing the shape parameters for the gamma distributions on the mode 3 projection (A) matrix.}
#'   \item{mode3.lambda.scale}{Matrix storing the scale parameters for the gamma distributions on the mode 3 projection (A) matrix.}
#'   \item{mode1.A.mean}{Matrix storing the mean parameters for the normal distributions on the mode 1 projection (A) matrix.}
#'   \item{mode1.A.cov}{Array storing the covariance parameters for the normal distributions on the mode 1 projection (A) matrix.}
#'   \item{mode2.A.mean}{Matrix storing the mean parameters for the normal distributions on the mode 2 projection (A) matrix.}
#'   \item{mode2.A.cov}{Array storing the covariance parameters for the normal distributions on the mode 2 projection (A) matrix.}
#'   \item{mode3.A.mean}{Matrix storing the mean parameters for the normal distributions on the mode 3 projection (A) matrix.}
#'   \item{mode3.A.cov}{Array storing the covariance parameters for the normal distributions on the mode 3 projection (A) matrix.}
#'   \item{mode1.H.mean}{Matrix storing the mean parameters for the normal distributions on the mode 1 latent (H) matrix.}
#'   \item{mode1.H.var}{Matrix storing the variance parameters for the normal distributions on the mode 1 latent (H) matrix.}
#'   \item{mode2.H.mean}{Matrix storing the mean parameters for the normal distributions on the mode 2 latent (H) matrix.}
#'   \item{mode2.H.var}{Matrix storing the variance parameters for the normal distributions on the mode 2 latent (H) matrix.}
#'   \item{mode3.H.mean}{Matrix storing the mean parameters for the normal distributions on the mode 3 latent (H) matrix.}
#'   \item{mode3.H.var}{Matrix storing the variance parameters for the normal distributions on the mode 3 latent (H) matrix.}
#'   \item{sigma2}{Variance for the response tensor.}
#'   \item{m1.sigma}{Variance for the mode 1 latent (H) matrix.}
#'   \item{m2.sigma}{Variance for the mode 2 latent (H) matrix.}
#'   \item{m3.sigma}{Variance for the mode 3 latent (H) matrix.}
#'   \item{m1.alpha}{Prior shape parameter for the gamma distribution on the precision of the mode 1 projection (A) matrix.}
#'   \item{m1.beta}{Prior scale paramet for the gamma distribution on the precision of the mode 1 projection (A) matrix.}
#'   \item{m2.alpha}{Prior shape parameter for the gamma distribution on the precision of the mode 2 projection (A) matrix.}
#'   \item{m2.beta}{Prior scale paramet for the gamma distribution on the precision of the mode 2 projection (A) matrix.}
#'   \item{m3.alpha}{Prior shape parameter for the gamma distribution on the precision of the mode 3 projection (A) matrix.}
#'   \item{m3.beta}{Prior scale paramet for the gamma distribution on the precision of the mode 3 projection (A) matrix.}
#'   \item{core.alpha}{Prior shape parameter for the gamma distribution on the precision of the core tensor.}
#'   \item{core.beta}{Prior scale parameter for the gamma distribution on the precision of the core tensor.}
#'   \item{core.0D.alpha}{Prior shape parameter for the gamma distribution on the precision of the 0D subset of the core tensor.}
#'   \item{core.0D.beta}{Prior scale parameter for the gamma distribution on the precision of the 0D subset of the core tensor.}
#'   \item{core.1D.alpha}{Prior shape parameter for the gamma distribution on the precision of the 1D subset of the core tensor.}
#'   \item{core.1D.beta}{Prior scale parameter for the gamma distribution on the precision of the 1D subset of the core tensor.}
#'   \item{core.2D.alpha}{Prior shape parameter for the gamma distribution on the precision of the 2D subset of the core tensor.}
#'   \item{core.2D.beta}{Prior scale parameter for the gamma distribution on the precision of the 2D subset of the core tensor.}
#'   \item{core.3D.alpha}{Prior shape parameter for the gamma distribution on the precision of the 3D subset of the core tensor.}
#'   \item{core.3D.beta}{Prior scale parameter for the gamma distribution on the precision of the 3D subset of the core tensor.}
#' }
#' @section Methods:
#' \describe{
#'   \item{\code{new(data, params)}}{Creates a new \code{Tucker_model} object with
#'   matrices sized accoring to the matrices in \code{data}.}
#'   \item{\code{rand_init(params)}}{Initializes the \code{Tucker_model} with
#'   random values accoring to \code{params}.}
#' }
#'
#' @examples
#' data.params <- get_data_params(c('decomp=Tucker'))
#' toy <- mk_toy(data.params)
#' train.data <- input_data$new(mode1.X=toy$mode1.X[,-1],
#'                              mode2.X=toy$mode2.X[,-1],
#'                              mode3.X=toy$mode3.X[,-1],
#'                              resp=toy$resp)
#' model.params <- get_model_params(c('decomp=Tucker'))
#' toy.model <- mk_model(train.data, model.params)
#' toy.model$rand_init(model.params)

Tucker_model <- R6Class("Tucker_model",
   portable=F, 
   class=F, 
   # cloneable=F,
   public=list(
     iter = 0,                              # Number of iterations of training
     early.stop = 0,                        # If the lower bound increases by less than this percentage, stop training.
     lower.bnd = numeric(),                 # Vector of lower bound values
     RMSE = numeric(),                      # Vector of training root mean squared errors for A matrices
     H.RMSE = numeric(),                    # Vector of training root mean squared errors for H matrices
     exp.var = numeric(),                   # Vector of explained variance for training data
     p.cor = numeric(),                     # Vector of Pearson correlations for training data
     s.cor = numeric(),                     # Vector of Spearman correlations for training data
     times = numeric(),                     # Vector of time to execute each loop
     core.mean = array(0,c(0,0,0)),         # Parameters for the core tensor
     core.var = array(0,c(0,0,0)),
     core.lambda.shape = array(0,c(0,0,0)),
     core.lambda.scale = array(0,c(0,0,0)),
     resp = array(0,c(0,0,0)),              # Predicted response tensor
     m1Xm1X = matrix(0,0,0),
     m2Xm2X = matrix(0,0,0),
     m3Xm3X = matrix(0,0,0),
     mode1.lambda.shape = matrix(0,0,0),    # Parameters for q Gamma distributions
     mode1.lambda.scale = matrix(0,0,0),
     mode2.lambda.shape = matrix(0,0,0),
     mode2.lambda.scale = matrix(0,0,0),
     mode3.lambda.shape = matrix(0,0,0),
     mode3.lambda.scale = matrix(0,0,0),
     mode1.A.mean = matrix(0,0,0),          # Parameters for the q A projection matrix distributions
     mode1.A.cov = array(0,c(0,0,0)),
     mode2.A.mean = matrix(0,0,0),
     mode2.A.cov = array(0,c(0,0,0)),
     mode3.A.mean = matrix(0,0,0),
     mode3.A.cov = array(0,c(0,0,0)),
     mode1.H.mean = matrix(0,0,0),          # Paramters for the q H latent matrices distributions
     mode1.H.var = matrix(0,0,0),
     mode2.H.mean = matrix(0,0,0),
     mode2.H.var = matrix(0,0,0),
     mode3.H.mean = matrix(0,0,0),
     mode3.H.var = matrix(0,0,0),
     sigma2 = 1,                            # Prior variance for response tensor
     m1.sigma2 = 1,                      # Prior variance for H matrices
     m2.sigma2 = 1,
     m3.sigma2 = 1,
     m1.alpha = 1,                       # Prior shape parameter for lambdas
     m1.beta = 1,                        # Prior scale parameter for lambdas
     m2.alpha = 1,
     m2.beta = 1,
     m3.alpha = 1,
     m3.beta = 1,
     core.alpha = array(0, c(0,0,0)),
     core.beta = array(0, c(0,0,0)),
     core.0D.alpha = 1,
     core.0D.beta = 1,
     core.1D.alpha = 1,
     core.1D.beta = 1,
     core.2D.alpha = 1,
     core.2D.beta = 1,
     core.3D.alpha = 1,
     core.3D.beta = 1,
     initialize = function(d, params) {
       # Make all variables in params accessable 
       for(i in 1:length(params)) {
         assign(names(params)[i], params[i][[1]])
       }
       
       # If sigma2=='auto' set it to the square root of the variance of the response data
       if(sigma2=='auto') sigma2 <- sqrt(var(d$resp, na.rm=T))
       
       I <- nrow(d$mode1.X); P <- ncol(d$mode1.X)
       J <- nrow(d$mode2.X); Q <- ncol(d$mode2.X)
       K <- nrow(d$mode3.X); S <- ncol(d$mode3.X)
       
       # These cross products are used in training
       if(A1.intercept) {
         m1Xm1X <<- crossprod(cbind(const=1,d$mode1.X), cbind(const=1,d$mode1.X))
       } else m1Xm1X <<- crossprod(d$mode1.X, d$mode1.X)
       if(A2.intercept) {
         m2Xm2X <<- crossprod(cbind(const=1,d$mode2.X), cbind(const=1,d$mode2.X))
       } else m2Xm2X <<- crossprod(d$mode2.X, d$mode2.X)
       if(A3.intercept) {
         m3Xm3X <<- crossprod(cbind(const=1,d$mode3.X), cbind(const=1,d$mode3.X))
       } else m3Xm3X <<- crossprod(d$mode3.X, d$mode3.X)

       # Initialize parameters for gamma priors on A matrices
       if(row.share) {
         if(A1.intercept) {
           mode1.lambda.shape <<- rep(m1.alpha, P+1)
           mode1.lambda.scale <<- rep(m1.beta, P+1)
           names(mode1.lambda.shape) <<- c('const', colnames(d$mode1.X))
           names(mode1.lambda.scale) <<- c('const', colnames(d$mode1.X))
         } else {
           mode1.lambda.shape <<- rep(m1.alpha, P)
           mode1.lambda.scale <<- rep(m1.beta, P)
           names(mode1.lambda.shape) <<- colnames(d$mode1.X)
           names(mode1.lambda.scale) <<- colnames(d$mode1.X)
         }
         if(A2.intercept) {
           mode2.lambda.shape <<- rep(m2.alpha, Q+1)
           mode2.lambda.scale <<- rep(m2.beta, Q+1)
           names(mode2.lambda.shape) <<- c('const', colnames(d$mode2.X))
           names(mode2.lambda.scale) <<- c('const', colnames(d$mode2.X))
         } else {
           mode2.lambda.shape <<- rep(m2.alpha, Q)
           mode2.lambda.scale <<- rep(m2.beta, Q)
           names(mode2.lambda.shape) <<- colnames(d$mode2.X)
           names(mode2.lambda.scale) <<- colnames(d$mode2.X)
         }
         if(A3.intercept) {
           mode3.lambda.shape <<- rep(m3.alpha, S+1)
           mode3.lambda.scale <<- rep(m3.beta, S+1)
           names(mode3.lambda.shape) <<- c('const', colnames(d$mode3.X))
           names(mode3.lambda.scale) <<- c('const', colnames(d$mode3.X))
         } else {
           mode3.lambda.shape <<- rep(m3.alpha, S)
           mode3.lambda.scale <<- rep(m3.beta, S)
           names(mode3.lambda.shape) <<- colnames(d$mode3.X)
           names(mode3.lambda.scale) <<- colnames(d$mode3.X)
         }
       } else {
         if(A1.intercept) {
           mode1.lambda.shape <<- matrix(m1.alpha, P+1, R1)
           mode1.lambda.scale <<- matrix(m1.beta, P+1, R1)
           dimnames(mode1.lambda.shape) <<- list(c('const', colnames(d$mode1.X)),
                                                 paste0('m1.r', 1:R1))
           dimnames(mode1.lambda.scale) <<- dimnames(mode1.lambda.shape)
         } else {
           mode1.lambda.shape <<- matrix(m1.alpha, P, R1)
           mode1.lambda.scale <<- matrix(m1.beta, P, R1)
           dimnames(mode1.lambda.shape) <<- list(colnames(d$mode1.X), paste0('m1.r', 1:R1))
           dimnames(mode1.lambda.scale) <<- dimnames(mode1.lambda.shape)
         }
         if(A2.intercept) {
           mode2.lambda.shape <<- matrix(m2.alpha, Q+1, R2)
           mode2.lambda.scale <<- matrix(m2.beta, Q+1, R2)
           dimnames(mode2.lambda.shape) <<- list(c('const', colnames(d$mode2.X)),
                                                 paste0('m2.r', 1:R2))
           dimnames(mode2.lambda.scale) <<- dimnames(mode2.lambda.shape)
         } else {
           mode2.lambda.shape <<- matrix(m2.alpha, Q, R2)
           mode2.lambda.scale <<- matrix(m2.beta, Q, R2)
           dimnames(mode2.lambda.shape) <<- list(colnames(d$mode2.X), paste0('m2.r', 1:R2))
           dimnames(mode2.lambda.scale) <<- dimnames(mode2.lambda.shape)
         }
         if(A3.intercept) {
           mode3.lambda.shape <<- matrix(m3.alpha, S+1, R3)
           mode3.lambda.scale <<- matrix(m3.beta, S+1, R3)
           dimnames(mode3.lambda.shape) <<- list(c('const', colnames(d$mode3.X)),
                                                 paste0('m3.r', 1:R3))
           dimnames(mode3.lambda.scale) <<- dimnames(mode3.lambda.shape)
         } else {
           mode3.lambda.shape <<- matrix(m3.alpha, S, R3)
           mode3.lambda.scale <<- matrix(m3.beta, S, R3)
           dimnames(mode3.lambda.shape) <<- list(colnames(d$mode3.X), paste0('m3.r', 1:R3))
           dimnames(mode3.lambda.scale) <<- dimnames(mode3.lambda.shape)
         }
       }
       
       # Initialize the sizes of the projection (A) matrices
       if(A1.intercept) {
         mode1.A.mean <<- matrix(0,P+1,R1)
         mode1.A.cov <<- array(0,dim=c(P+1,P+1,R1))
         dimnames(mode1.A.mean) <<- list(c('const', colnames(d$mode1.X)),
                                         paste0('m1.r', 1:R1))
         dimnames(mode1.A.cov) <<- list(c('const', colnames(d$mode1.X)),
                                        c('const', colnames(d$mode1.X)),
                                        paste0('m1.r', 1:R1))
       } else {
         mode1.A.mean <<- matrix(0,P,R1)
         mode1.A.cov <<- array(0,dim=c(P,P,R1))
         dimnames(mode1.A.mean) <<- list(colnames(d$mode1.X),
                                         paste0('m1.r', 1:R1))
         dimnames(mode1.A.cov) <<- list(colnames(d$mode1.X),
                                        colnames(d$mode1.X),
                                        paste0('m1.r', 1:R1))
       }
       if(A2.intercept) {
         mode2.A.mean <<- matrix(0,Q+1,R2)
         mode2.A.cov <<- array(0,dim=c(Q+1,Q+1,R2))
         dimnames(mode2.A.mean) <<- list(c('const', colnames(d$mode2.X)),
                                         paste0('m2.r', 1:R2))
         dimnames(mode2.A.cov) <<- list(c('const', colnames(d$mode2.X)),
                                        c('const', colnames(d$mode2.X)),
                                        paste0('m2.r', 1:R2))
       } else {
         mode2.A.mean <<- matrix(0,Q,R2)
         mode2.A.cov <<- array(0,dim=c(Q,Q,R2))
         dimnames(mode2.A.mean) <<- list(colnames(d$mode2.X),
                                         paste0('m2.r', 1:R2))
         dimnames(mode2.A.cov) <<- list(colnames(d$mode2.X),
                                        colnames(d$mode2.X),
                                        paste0('m2.r', 1:R2))
       }
       if(A3.intercept) {
         mode3.A.mean <<- matrix(0,S+1,R3)
         mode3.A.cov <<- array(0,dim=c(S+1,S+1,R3))
         dimnames(mode3.A.mean) <<- list(c('const', colnames(d$mode3.X)),
                                         paste0('m3.r', 1:R3))
         dimnames(mode3.A.cov) <<- list(c('const', colnames(d$mode3.X)),
                                        c('const', colnames(d$mode3.X)),
                                        paste0('m3.r', 1:R3))
       } else {
         mode3.A.mean <<- matrix(0,S,R3)
         mode3.A.cov <<- array(0,dim=c(S,S,R3))
         dimnames(mode3.A.mean) <<- list(colnames(d$mode3.X),
                                         paste0('m3.r', 1:R3))
         dimnames(mode3.A.cov) <<- list(colnames(d$mode3.X),
                                        colnames(d$mode3.X),
                                        paste0('m3.r', 1:R3))
       }

       # Initialize the sizes of the latent (H) matrices
       if(H1.intercept) {
         mode1.H.mean <<- matrix(0,I,R1+1)
         mode1.H.var <<- matrix(0,I,R1+1)
         dimnames(mode1.H.mean) <<- list(dimnames(d$mode1.X)[[1]],
                                         c('const', paste0('m1.r', 1:R1)))
         dimnames(mode1.H.var) <<- dimnames(mode1.H.mean)
       } else {
         mode1.H.mean <<- matrix(0,I,R1)
         mode1.H.var <<- matrix(0,I,R1)
         dimnames(mode1.H.mean) <<- list(dimnames(d$mode1.X)[[1]],
                                         paste0('m1.r', 1:R1))
         dimnames(mode1.H.var) <<- dimnames(mode1.H.mean)
       }
       if(H2.intercept) {
         mode2.H.mean <<- matrix(0,J,R2+1)
         mode2.H.var <<- matrix(0,J,R2+1)
         dimnames(mode2.H.mean) <<- list(dimnames(d$mode2.X)[[1]],
                                         c('const', paste0('m2.r', 1:R2)))
         dimnames(mode2.H.var) <<- dimnames(mode2.H.mean)
       } else {
         mode2.H.mean <<- matrix(0,J,R2)
         mode2.H.var <<- matrix(0,J,R2)
         dimnames(mode2.H.mean) <<- list(dimnames(d$mode2.X)[[1]],
                                         paste0('m2.r', 1:R2))
         dimnames(mode2.H.var) <<- dimnames(mode2.H.mean)
       }
       if(H3.intercept) {
         mode3.H.mean <<- matrix(0,K,R3+1)
         mode3.H.var <<- matrix(0,K,R3+1)
         dimnames(mode3.H.mean) <<- list(dimnames(d$mode3.X)[[1]],
                                         c('const', paste0('m3.r', 1:R3)))
         dimnames(mode3.H.var) <<- dimnames(mode3.H.mean)
       } else {
         mode3.H.mean <<- matrix(0,K,R3)
         mode3.H.var <<- matrix(0,K,R3)
         dimnames(mode3.H.mean) <<- list(dimnames(d$mode3.X)[[1]],
                                         paste0('m3.r', 1:R3))
         dimnames(mode3.H.var) <<- dimnames(mode3.H.mean)
       }
       
       # Initialize the core to the correct size
       if(H1.intercept) core1 <- R1+1 else core1 <- R1
       if(H2.intercept) core2 <- R2+1 else core2 <- R2
       if(H3.intercept) core3 <- R3+1 else core3 <- R3
       core.mean <<- array(0, dim=c(core1,core2,core3))
       if(H1.intercept) {
         dimnames(core.mean)[[1]] <<- c('const', paste0('m1.r', 1:R1))
       } else dimnames(core.mean)[[1]] <<- paste0('m1.r', 1:R1)
       if(H2.intercept) {
         dimnames(core.mean)[[2]] <<- c('const', paste0('m2.r', 1:R2))
       } else dimnames(core.mean)[[2]] <<- paste0('m2.r', 1:R2)
       if(H3.intercept) {
         dimnames(core.mean)[[3]] <<- c('const', paste0('m3.r', 1:R3))
       } else dimnames(core.mean)[[3]] <<- paste0('m3.r', 1:R3)
       core.var <<- array(0, dim=dim(core.mean), 
                          dimnames=dimnames(core.mean))
       
       # initialize the parameters for the gamma prior on the core
       core.alpha <<- array(core.3D.alpha, dim=dim(core.mean),
                                   dimnames=dimnames(core.mean))
       core.beta <<- array(core.3D.beta, dim=dim(core.mean),
                                   dimnames=dimnames(core.mean))
       if(H1.intercept) {
         core.alpha[1,,] <<- core.2D.alpha
         core.beta[1,,] <<- core.2D.beta
       }
       if(H2.intercept) {
         core.alpha[,1,] <<- core.2D.alpha
         core.beta[,1,] <<- core.2D.beta
       }
       if(H3.intercept) {
         core.alpha[,,1] <<- core.2D.alpha
         core.beta[,,1] <<- core.2D.beta
       }
       if(H1.intercept & H2.intercept) {
         core.alpha[1,1,] <<- core.1D.alpha
         core.beta[1,1,] <<- core.1D.beta
       }
       if(H1.intercept & H3.intercept) {
         core.alpha[1,,1] <<- core.1D.alpha
         core.beta[1,,1] <<- core.1D.beta
       }
       if(H2.intercept & H3.intercept) {
         core.alpha[,1,1] <<- core.1D.alpha
         core.beta[,1,1] <<- core.1D.beta
       }
       if(H1.intercept & H2.intercept & H3.intercept) {
         core.alpha[1,1,1] <<- core.0D.alpha
         core.beta[1,1,1] <<- core.0D.beta
       }
       core.lambda.shape <<- core.alpha
       core.lambda.scale <<- core.beta

       # Initialize the array of predicted responses
       resp <<- array(0,dim=dim(d$resp))
       
       # Store hyperparameters
       sigma2 <<- sigma2
       m1.sigma2 <<- m1.sigma2
       m2.sigma2 <<- m2.sigma2
       m3.sigma2 <<- m3.sigma2
       m1.alpha  <<- m1.alpha
       m1.beta   <<- m1.beta
       m2.alpha  <<- m2.alpha
       m2.beta   <<- m2.beta
       m3.alpha  <<- m3.alpha
       m3.beta   <<- m3.beta
     },   
     rand_init = function(params) {
       # Make all variables in params accessable 
       for(i in 1:length(params)) {
         assign(names(params)[i], params[i][[1]])
       }
       
       # Use random seed if given
       if(!is.na(seed)) set.seed(seed)
       
       # Get sizes for easy reference (A sizes include A.intercept)
       I <- nrow(mode1.H.mean); P <- nrow(mode1.A.mean)
       J <- nrow(mode2.H.mean); Q <- nrow(mode2.A.mean)
       K <- nrow(mode3.H.mean); S <- nrow(mode3.A.mean)
       R1 <- ncol(mode1.A.mean); R2 <- ncol(mode2.A.mean); R3 <- ncol(mode3.A.mean)
       core1 <- ncol(mode1.H.mean); core2 <- ncol(mode2.H.mean); core3 <- ncol(mode3.H.mean)
       
       # Initialize the A mean matrices with mean zero and standard dev. given by A.samp.sd
       # Normally want these to be over dispersed initially
       if(P) mode1.A.mean <<- matrix(rnorm(P*R1, sd=A.samp.sd), P, R1, 
                                     dimnames=dimnames(mode1.A.mean))
       if(Q) mode2.A.mean <<- matrix(rnorm(Q*R2, sd=A.samp.sd), Q, R2, 
                                     dimnames=dimnames(mode2.A.mean))
       if(S) mode3.A.mean <<- matrix(rnorm(S*R3, sd=A.samp.sd), S, R3, 
                                     dimnames=dimnames(mode3.A.mean))
       
       # Initialize the A variance arrays with the value A.var
       if(P) {
         for(r1 in 1:R1) mode1.A.cov[,,r1] <<- diag(A.var , P)
         dimnames(mode1.A.cov) <<- list(rownames(mode1.A.mean), 
                                        rownames(mode1.A.mean),
                                        paste0('m1.r', 1:R1))
       }
       if(Q) {
         for(r2 in 1:R2) mode2.A.cov[,,r2] <<- diag(A.var , Q)
         dimnames(mode2.A.cov) <<- list(rownames(mode2.A.mean), 
                                        rownames(mode2.A.mean),
                                        paste0('m2.r', 1:R2))
       }
       if(S) {
         for(r3 in 1:R3) mode3.A.cov[,,r3] <<- diag(A.var , S)
         dimnames(mode3.A.cov) <<- list(rownames(mode3.A.mean), 
                                        rownames(mode3.A.mean),
                                        paste0('m3.r', 1:R3))
       }
       
       # Initialize the H variance matrices with H.var
       mode1.H.var <<- matrix(H.var, I, ncol=core1, dimnames=dimnames(mode1.H.var))
       mode2.H.var <<- matrix(H.var, J, ncol=core2, dimnames=dimnames(mode2.H.var))
       mode3.H.var <<- matrix(H.var, K, ncol=core3, dimnames=dimnames(mode3.H.var))
       
       # If random.H is true, initialize latent (H) mean matrices with mean 0 and std. 
       # dev. given by H.samp.sd. The first constant column are all 1.
       # If random.H is false and A matrices exist, initialize as the product 
       # of the X and A matrices with added constant column
       if(random.H | P==0) {
         if(H1.intercept) {
           mode1.H.mean <<- matrix(c(rep(1,I), rnorm(I*R1, sd=H.samp.sd)),
                                   I, R1+1, dimnames=dimnames(mode1.H.mean))
         } else {
           mode1.H.mean <<- matrix(rnorm(I*R1, sd=H.samp.sd),
                                   I, R1, dimnames=dimnames(mode1.H.mean))
         }
       } else {
         if(H1.intercept) {
           dn <- dimnames(mode1.H.mean)
           mode1.H.mean <<- cbind(1, d$mode1.X %*% mode1.A.mean)
           dimnames(mode1.H.mean) <- dn
         } else {
           dn <- dimnames(mode1.H.mean)
           mode1.H.mean <<- d$mode1.X %*% mode1.A.mean
           dimnames(mode1.H.mean) <- dn
         }
       }
       if(random.H | Q==0) {
         if(H2.intercept) {
           mode2.H.mean <<- matrix(c(rep(1,J), rnorm(J*R2, sd=H.samp.sd)),
                                   J, R2+1, dimnames=dimnames(mode2.H.mean))
         } else {
           mode2.H.mean <<- matrix(rnorm(J*R2, sd=H.samp.sd),
                                   J, R2, dimnames=dimnames(mode2.H.mean))
         }
       } else {
         if(H2.intercept) {
           dn <- dimnames(mode2.H.mean)
           mode2.H.mean <<- cbind(1, d$mode2.X %*% mode2.A.mean)
           dimnames(mode2.H.mean) <- dn
         } else {
           dn <- dimnames(mode2.H.mean)
           mode2.H.mean <<- d$mode2.X %*% mode2.A.mean
           dimnames(mode2.H.mean) <- dn
         }
       }
       if(random.H | S==0) {
         if(H3.intercept) {
           mode3.H.mean <<- matrix(c(rep(1,K), rnorm(K*R3, sd=H.samp.sd)),
                                   K, R3+1, dimnames=dimnames(mode3.H.mean))
         } else {
           mode3.H.mean <<- matrix(rnorm(K*R3, sd=H.samp.sd),
                                   K, R3, dimnames=dimnames(mode3.H.mean))
         }
       } else {
         if(H3.intercept) {
           dn <- dimnames(mode3.H.mean)
           mode3.H.mean <<- cbind(1, d$mode3.X %*% mode3.A.mean)
           dimnames(mode3.H.mean) <- dn
         } else {
           dn <- dimnames(mode3.H.mean)
           mode3.H.mean <<- d$mode3.X %*% mode3.A.mean
           dimnames(mode3.H.mean) <- dn
         }
       }
       
       # Initialize the core tensor with mean 0 and variance R.var
       core.mean <<- array(rnorm((core1)*(core2)*(core3), sd=R.samp.sd), 
                           dim = c(core1, core2, core3),
                           dimnames=dimnames(core.mean))
       core.var <<- array(R.var, dim = c(core1, core2, core3), 
                          dimnames=dimnames(core.mean))
       
       print('** Initialization complete **')
     }
   )
)
