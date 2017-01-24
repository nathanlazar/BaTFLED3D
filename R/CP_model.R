#' BaTFLED model object for 3-D response tensor with CP decomposition.
#' 
#' \code{CP_model} objects are 'R6' objects so that their values can be updated in place. 
#' The object is treated like an environment and components are accessed using the \code{$} 
#' operator. 
#' When creating a new CP_model object it will be populated with default values and empty matrices.
#' To initialize a \code{CP_model} call the \code{initialize()} method. 
#' 
#' @docType class
#' @importFrom R6 R6Class
#' @importFrom stats var
#' @exportClass CP_model
#' @export CP_model
#' @format An \code{\link{R6Class}} generator object
#' @keywords data
#' @section Methods:
#' \describe{
#'   \item{\code{new(data, params)}}{Creates a new \code{CP_model} object with
#'   matrices sized accoring to the matrices in \code{data}.}
#'   \item{\code{rand_init(params)}}{Initializes the \code{CP_model} with
#'   random values accoring to \code{params}.}
#' }
#' 
#' 
#' @method clone() Creates a copy of the object with a new memory space.
#' @method initialize() Initializes a model with random values.
#' 
#' @seealso \code{get_model_params}, \code{input_data}, \code{Tucker_model}
#' 
#' @examples
#' data.params <- get_data_params(c('decomp=CP'))
#' toy <- mk_toy(data.params)
#' train.data <- input_data$new(mode1.X=toy$mode1.X[,-1],
#'                              mode2.X=toy$mode2.X[,-1],
#'                              mode3.X=toy$mode3.X[,-1],
#'                              resp=toy$resp)
#' model.params <- get_model_params(c('decomp=CP'))
#' toy.model <- mk_model(train.data, model.params)
#' toy.model$rand_init(model.params)

CP_model <- R6Class("CP_model",
    portable=F,
    class=F,
    public=list(
      iter = 0,                         # Number of iterations of training
      early.stop = 0,                   # Limit to stop training for norm of A matrices
      lower.bnd = numeric(),            # Vector of lower bound values
      RMSE = numeric(),                 # Vector of training root mean squared errors for A matrices
      H.RMSE = numeric(),               # Vector of training root mean squared errors for H matrices
      exp.var = numeric(),              # Vector of explained variance for training data
      p.cor = numeric(),                # Vector of Pearson correlations for training data
      s.cor = numeric(),                # Vector of Spearman correlations for training data
      times = numeric(),                # Vector of time to execute each loop
      resp = array(),                   # Predicted response tensor
      m1Xm1X = matrix(),
      m2Xm2X = matrix(),
      m3Xm3X = matrix(),
      mode1.lambda.shape = matrix(),    # Parameters for q Gamma distributions
      mode1.lambda.scale = matrix(),
      mode2.lambda.shape = matrix(),
      mode2.lambda.scale = matrix(),
      mode3.lambda.shape = matrix(),
      mode3.lambda.scale = matrix(),
      mode1.A.mean = matrix(),          # Parameters for the q A projection matrix distributions
      mode1.A.cov = array(),
      mode2.A.mean = matrix(),
      mode2.A.cov = array(),
      mode3.A.mean = matrix(),
      mode3.A.cov = array(),
      mode1.H.mean = matrix(),          # Paramters for the q H latent matrices distributions
      mode1.H.var = matrix(),
      mode2.H.mean = matrix(),
      mode2.H.var = matrix(),
      mode3.H.mean = matrix(),
      mode3.H.var = matrix(),
      sigma2 = 1,                       # Prior variance for response tensor
      m1.sigma2 = 1,                 # Prior variance for H matrices
      m2.sigma2 = 1,
      m3.sigma2 = 1,
      m1.alpha = 1,                  # Prior shape parameter for lambdas
      m1.beta = 1,                   # Prior scale parameter for lambdas
      m2.alpha = 1,
      m2.beta = 1,
      m3.alpha = 1,
      m3.beta = 1,
      initialize = function(d, params) {
        # Make all variables in params accessable 
        for(i in 1:length(params)) {
          assign(names(params)[i], params[i][[1]])
        }
        
        I <- nrow(d$mode1.X); P <- ncol(d$mode1.X)
        J <- nrow(d$mode2.X); Q <- ncol(d$mode2.X)
        K <- nrow(d$mode3.X); S <- ncol(d$mode3.X)
        
        # If sigma2=='auto' set it to the square root of the variance of the response data
        if(sigma2=='auto') sigma2 <- sqrt(var(d$resp, na.rm=T))
  
        # If no names are given for the input matrices, use defaults
        if(is.null(rownames(d$mode1.X))) rownames(d$mode1.X) <- paste0('m1.samp.', 1:I)
        if(is.null(colnames(d$mode1.X)) & P) colnames(d$mode1.X) <- paste0('m1.feat.', 1:P)
        if(is.null(rownames(d$mode2.X))) rownames(d$mode2.X) <- paste0('m2.samp.', 1:J)
        if(is.null(colnames(d$mode2.X)) & Q) colnames(d$mode2.X) <- paste0('m2.feat.', 1:Q)
        if(is.null(rownames(d$mode3.X))) rownames(d$mode3.X) <- paste0('m3.samp.', 1:K)
        if(is.null(colnames(d$mode3.X)) & S) colnames(d$mode3.X) <- paste0('m3.feat.', 1:S)
        if(is.null(dimnames(d$resp)[[1]])) dimnames(d$resp)[[1]] <- rownames(d$mode1.X)
        if(is.null(dimnames(d$resp)[[2]])) dimnames(d$resp)[[2]] <- rownames(d$mode2.X)
        if(is.null(dimnames(d$resp)[[3]])) dimnames(d$resp)[[3]] <- rownames(d$mode3.X)
        
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
            mode1.lambda.shape <<- matrix(m1.alpha, P+1, R)
            mode1.lambda.scale <<- matrix(m1.beta, P+1, R)
            dimnames(mode1.lambda.shape) <<- list(c('const', colnames(d$mode1.X)),
                                                  paste0('m1.r', 1:R))
            dimnames(mode1.lambda.scale) <<- dimnames(mode1.lambda.shape)
          } else {
            mode1.lambda.shape <<- matrix(m1.alpha, P, R)
            mode1.lambda.scale <<- matrix(m1.beta, P, R)
            dimnames(mode1.lambda.shape) <<- list(colnames(d$mode1.X), paste0('m1.r', 1:R))
            dimnames(mode1.lambda.scale) <<- dimnames(mode1.lambda.shape)
          }
          if(A2.intercept) {
            mode2.lambda.shape <<- matrix(m2.alpha, Q+1, R)
            mode2.lambda.scale <<- matrix(m2.beta, Q+1, R)
            dimnames(mode2.lambda.shape) <<- list(c('const', colnames(d$mode2.X)),
                                                  paste0('m2.r', 1:R))
            dimnames(mode2.lambda.scale) <<- dimnames(mode2.lambda.shape)
          } else {
            mode2.lambda.shape <<- matrix(m2.alpha, Q, R)
            mode2.lambda.scale <<- matrix(m2.beta, Q, R)
            dimnames(mode2.lambda.shape) <<- list(colnames(d$mode2.X), paste0('m2.r', 1:R))
            dimnames(mode2.lambda.scale) <<- dimnames(mode2.lambda.shape)
          }
          if(A3.intercept) {
            mode3.lambda.shape <<- matrix(m3.alpha, S+1, R)
            mode3.lambda.scale <<- matrix(m3.beta, S+1, R)
            dimnames(mode3.lambda.shape) <<- list(c('const', colnames(d$mode3.X)),
                                                  paste0('m3.r', 1:R))
            dimnames(mode3.lambda.scale) <<- dimnames(mode3.lambda.shape)
          } else {
            mode3.lambda.shape <<- matrix(m3.alpha, S, R)
            mode3.lambda.scale <<- matrix(m3.beta, S, R)
            dimnames(mode3.lambda.shape) <<- list(colnames(d$mode3.X), paste0('m3.r', 1:R))
            dimnames(mode3.lambda.scale) <<- dimnames(mode3.lambda.shape)
          }
        }
        
        # Initialize the sizes of the projection (A) matrices
        if(A1.intercept) {
          mode1.A.mean <<- matrix(0,P+1,R)
          mode1.A.cov <<- array(0,dim=c(P+1,P+1,R))
          dimnames(mode1.A.mean) <<- list(c('const', colnames(d$mode1.X)),
                                          paste0('m1.r', 1:R))
          dimnames(mode1.A.cov) <<- list(c('const', colnames(d$mode1.X)),
                                         c('const', colnames(d$mode1.X)),
                                         paste0('m1.r', 1:R))
        } else {
          mode1.A.mean <<- matrix(0,P,R)
          mode1.A.cov <<- array(0,dim=c(P,P,R))
          dimnames(mode1.A.mean) <<- list(colnames(d$mode1.X),
                                          paste0('m1.r', 1:R))
          dimnames(mode1.A.cov) <<- list(colnames(d$mode1.X),
                                         colnames(d$mode1.X),
                                         paste0('m1.r', 1:R))
        }
        if(A2.intercept) {
          mode2.A.mean <<- matrix(0,Q+1,R)
          mode2.A.cov <<- array(0,dim=c(Q+1,Q+1,R))
          dimnames(mode2.A.mean) <<- list(c('const', colnames(d$mode2.X)),
                                          paste0('m2.r', 1:R))
          dimnames(mode2.A.cov) <<- list(c('const', colnames(d$mode2.X)),
                                         c('const', colnames(d$mode2.X)),
                                         paste0('m2.r', 1:R))
        } else {
          mode2.A.mean <<- matrix(0,Q,R)
          mode2.A.cov <<- array(0,dim=c(Q,Q,R))
          dimnames(mode2.A.mean) <<- list(colnames(d$mode2.X),
                                          paste0('m2.r', 1:R))
          dimnames(mode2.A.cov) <<- list(colnames(d$mode2.X),
                                         colnames(d$mode2.X),
                                         paste0('m2.r', 1:R))
        }
        if(A3.intercept) {
          mode3.A.mean <<- matrix(0,S+1,R)
          mode3.A.cov <<- array(0,dim=c(S+1,S+1,R))
          dimnames(mode3.A.mean) <<- list(c('const', colnames(d$mode3.X)),
                                          paste0('m3.r', 1:R))
          dimnames(mode3.A.cov) <<- list(c('const', colnames(d$mode3.X)),
                                         c('const', colnames(d$mode3.X)),
                                         paste0('m3.r', 1:R))
        } else {
          mode3.A.mean <<- matrix(0,S,R)
          mode3.A.cov <<- array(0,dim=c(S,S,R))
          dimnames(mode3.A.mean) <<- list(colnames(d$mode3.X),
                                          paste0('m3.r', 1:R))
          dimnames(mode3.A.cov) <<- list(colnames(d$mode3.X),
                                         colnames(d$mode3.X),
                                         paste0('m3.r', 1:R))
        }

        # Initialize the sizes of the latent (H) matrices
        mode1.H.mean <<- matrix(0,I,R)
        mode1.H.var <<- matrix(0,I,R)
        dimnames(mode1.H.mean) <<- list(dimnames(d$mode1.X)[[1]], paste0('m1.r', 1:R))
        dimnames(mode1.H.var) <<- dimnames(mode1.H.mean)
        mode2.H.mean <<- matrix(0,J,R)
        mode2.H.var <<- matrix(0,J,R)
        dimnames(mode2.H.mean) <<- list(dimnames(d$mode2.X)[[1]], paste0('m2.r', 1:R))
        dimnames(mode2.H.var) <<- dimnames(mode2.H.mean)
        mode3.H.mean <<- matrix(0,K,R)
        mode3.H.var <<- matrix(0,K,R)
        dimnames(mode3.H.mean) <<- list(dimnames(d$mode3.X)[[1]], paste0('m3.r', 1:R))
        dimnames(mode3.H.var) <<- dimnames(mode3.H.mean)
        
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
          
        # Get sizes for easy reference
        R <- ncol(mode1.H.mean)
        I <- nrow(mode1.H.mean); P <- nrow(mode1.A.mean)
        J <- nrow(mode2.H.mean); Q <- nrow(mode2.A.mean)
        K <- nrow(mode3.H.mean); S <- nrow(mode3.A.mean)

        # Initialize the A mean matrices with mean zero and standard dev. given by A.samp.sd
        # Normally want these to be over dispersed initially
        if(P) mode1.A.mean <<- matrix(rnorm(P*R, sd=A.samp.sd), P, R, dimnames=dimnames(mode1.A.mean))
        if(Q) mode2.A.mean <<- matrix(rnorm(Q*R, sd=A.samp.sd), Q, R, dimnames=dimnames(mode2.A.mean))
        if(S) mode3.A.mean <<- matrix(rnorm(S*R, sd=A.samp.sd), S, R, dimnames=dimnames(mode3.A.mean))
        
        # Initialize the A variance arrays with the value A.var
        if(P) {
          for(r in 1:R) mode1.A.cov[,,r] <<- diag(A.var , P)
          dimnames(mode1.A.cov) <<- list(rownames(mode1.A.mean), rownames(mode1.A.mean), paste0('m1.r', 1:R))
        }
        if(Q) {
          for(r in 1:R) mode2.A.cov[,,r] <<- diag(A.var , Q)
          dimnames(mode2.A.cov) <<- list(rownames(mode2.A.mean), rownames(mode2.A.mean), paste0('m2.r', 1:R))
        }
        if(S) {
          for(r in 1:R) mode3.A.cov[,,r] <<- diag(A.var , S)
          dimnames(mode3.A.cov) <<- list(rownames(mode3.A.mean), rownames(mode3.A.mean), paste0('m3.r', 1:R))
        }

        # Initialize the H variance matrices with H.var
        mode1.H.var <- matrix(H.var, I, ncol=R, dimnames=dimnames(mode1.H.var))
        mode2.H.var <- matrix(H.var, J, ncol=R, dimnames=dimnames(mode2.H.var))
        mode3.H.var <- matrix(H.var, K, ncol=R, dimnames=dimnames(mode3.H.var))
          
        # Initialize latent (H) mean matrices with mean 0 and std. dev. given by H.samp.sd
        # if the random.H parameter is true otherwise initialize as the product of the X
        # and A matrices.
        if(random.H | P==0) {
          mode1.H.mean <<- matrix(rnorm(I*R, sd=H.samp.sd), I, R,
                                  dimnames=dimnames(mode1.H.mean))
        } else mode1.H.mean <<- d$mode1.X %*% mode1.A.mean
        if(random.H | Q==0) {
          mode2.H.mean <<- matrix(rnorm(J*R, sd=H.samp.sd), J, R,
                                  dimnames=dimnames(mode2.H.mean))
        } else mode2.H.mean <<- d$mode2.X %*% mode2.A.mean
        if(random.H | S==0) {
          mode3.H.mean <<- matrix(rnorm(K*R, sd=H.samp.sd), K, R,
                                  dimnames=dimnames(mode3.H.mean))
        } else mode3.H.mean <<- d$mode3.X %*% mode3.A.mean
          
        print('** Initialization complete **')
      }
    )
)
