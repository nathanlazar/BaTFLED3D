#' Update the first mode in a CP model.
#'
#' Update is performed in place to avoid memory issues. There is no return value.
#' 
#' @export
#' @param m A \code{CP_model} object created with \code{mk_model} 
#' @param d Input data object created with \code{input_data}
#' @param params List of parameters created with \code{get_model_params()}
#' @examples
#' data.params <- get_data_params(c('decomp=CP'))
#' toy <- mk_toy(data.params)
#' train.data <- input_data$new(mode1.X=toy$mode1.X[,-1],
#'                              mode2.X=toy$mode2.X[,-1],
#'                              mode3.X=toy$mode3.X,
#'                              resp=toy$resp)
#' model.params <- get_model_params(c('decomp=CP'))
#' toy.model <- mk_model(train.data, model.params)
#' toy.model$rand_init(model.params)
#'
#' update_mode1_CP(m=toy.model, d=train.data, params=model.params)

update_mode1_CP <- function(m, d, params) {
  I <- dim(d$resp)[1]
  P <- nrow(m$mode1.A.mean)
  R <- ncol(m$mode1.H.mean)
  
  # If the intercept term is removed change A1.intercept  
  A1.intercept <- ifelse('const' %in% rownames(m$mode1.A.mean), T, F)

  if(P != 0) { # If there is no input data, skip updates for lambda and A
    if(params$verbose) print("Updating prior lambda vector for mode 1")
    
    m1.A.var <- matrix(0, P, R)
    for(r in 1:R) m1.A.var[,r] <- diag(m$mode1.A.cov[,,r])
    if(params$row.share) {
      m$mode1.lambda.scale <- 1/(.5*(rowSums(m$mode1.A.mean^2 + m1.A.var)) + 1/m$m1.beta)
    } else m$mode1.lambda.scale <- 1/(.5*(m$mode1.A.mean^2 + m1.A.var) + 1/m$m1.beta)
    
    if(params$verbose) print("Updating projection (A) matrix for mode 1")
    # Update mode1.A covariance parameters. They only rely on X and lambdas
    lambda.exp <- m$mode1.lambda.shape * m$mode1.lambda.scale
    for(r in 1:R) {
      if(params$row.share) {
        m$mode1.A.cov[,,r] <- chol2inv(chol(diag(lambda.exp) + (1/m$m1.sigma2) * m$m1Xm1X))
      } else 
        m$mode1.A.cov[,,r] <- chol2inv(chol(diag(lambda.exp[,r]) + (1/m$m1.sigma2) * m$m1Xm1X))
    }

    # Update A means
    if(A1.intercept) {
      for(r in 1:R) m$mode1.A.mean[,r] <- (1/m$m1.sigma2) * 
          (m$mode1.A.cov[,,r] %*% t(cbind(1,d$mode1.X)) %*% m$mode1.H.mean[,r])
    } else {
      for(r in 1:R) m$mode1.A.mean[,r] <- (1/m$m1.sigma2) * 
          (m$mode1.A.cov[,,r] %*% t(d$mode1.X) %*% m$mode1.H.mean[,r])
    }
  }
  
  # Update the variance and mean for the H factor matrices
  if(params$verbose) print("Updating latent (H) matrix for mode 1")
  
  # Update the variance first. sapply vectorizes the updates for each row
  for(r in 1:R) {
    m$mode1.H.var[,r] <- sapply(1:I, function(i) 1/((1/m$sigma2) *
      sum(d$delta[i,,] * outer((m$mode2.H.mean[,r]^2 + m$mode2.H.var[,r]),
          (m$mode3.H.mean[,r]^2 + m$mode3.H.var[,r]))) + 1/m$m1.sigma2))
  }
  
  # Next update the mean
  for(i in 1:I) {
    for(r in 1:R) {
      if(P == 0) {x_times_a <- 0} else {
        x_times_a <- safe_prod(d$mode1.X[i,,drop=F], m$mode1.A.mean[,r,drop=F])
      }
      m$mode1.H.mean[i,r] <- m$mode1.H.var[i,r] * ((1/m$sigma2) *
        sum(outer(m$mode2.H.mean[,r], m$mode3.H.mean[,r]) *
            (d$resp[i,,] - (sweep(m$mode2.H.mean[,-r,drop=F], MARGIN=2, m$mode1.H.mean[i,-r], '*') %*%
            t(m$mode3.H.mean[,-r,drop=F]))), na.rm=T) + 1/m$m1.sigma2 * x_times_a)
    }
  }
}
