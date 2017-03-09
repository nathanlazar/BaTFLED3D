#' Update the second mode in a CP model.
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
#' update_mode2_CP(m=toy.model, d=train.data, params=model.params)

update_mode2_CP <- function(m, d, params) {
  J <- dim(d$resp)[2]
  Q <- nrow(m$mode2.A.mean)
  R <- ncol(m$mode2.H.mean)
  
  # If the intercept term is removed change A1.intercept  
  A2.intercept <- ifelse('const' %in% rownames(m$mode2.A.mean), T, F)

  if(Q != 0) { # If there is no input data, skip updates for lambda and A
    if(params$verbose) print("Updating prior lambda vector for mode 2")
    m2.A.var <- matrix(0, Q, R)
    for(r in 1:R) m2.A.var[,r] <- diagonal(m$mode2.A.cov[,,r])
    if(params$row.share) {
      m$mode2.lambda.scale <- 1/(.5*(rowSums(m$mode2.A.mean^2 + m2.A.var)) + 1/m$m2.beta)
    } else m$mode2.lambda.scale <- 1/(.5*(m$mode2.A.mean^2 + m2.A.var) + 1/m$m2.beta)
    
    if(params$verbose) print("Updating projection (A) matrix for mode 2")
    # Update mode2.A covariance parameters. They only rely on X and lambdas
    lambda.exp <- m$mode2.lambda.shape * m$mode2.lambda.scale
    for(r in 1:R) {
      if(params$row.share) {
        m$mode2.A.cov[,,r] <- chol2inv(chol(diagonal(lambda.exp) + (1/m$m2.sigma2) * m$m2Xm2X))
      } else
        m$mode2.A.cov[,,r] <- chol2inv(chol(diagonal(lambda.exp[,r]) + (1/m$m2.sigma2) * m$m2Xm2X))
    }
    
    # Update A means
    if(A2.intercept) {
      for(r in 1:R) m$mode2.A.mean[,r] <- (1/m$m2.sigma2) * 
          (m$mode2.A.cov[,,r] %*% t(cbind(1,d$mode2.X)) %*% m$mode2.H.mean[,r])
    } else {
      for(r in 1:R) m$mode2.A.mean[,r] <- (1/m$m2.sigma2) * 
          (m$mode2.A.cov[,,r] %*% t(d$mode2.X) %*% m$mode2.H.mean[,r])
    }
  }
  
  if(params$verbose) print("Updating latent (H) matrix for mode 2")
  # Update the variance first. sapply vectorizes the updates for each row
  for(r in 1:R) {
    m$mode2.H.var[,r] <- sapply(1:J, function(j) 1/((1/m$sigma2) *
      sum(d$delta[,j,] * outer((m$mode1.H.mean[,r]^2 + m$mode1.H.var[,r]),
          (m$mode3.H.mean[,r]^2 + m$mode3.H.var[,r]))) + 1/m$m2.sigma2))
  }
  
  # Next update the mean
  for(j in 1:J) {
    for(r in 1:R) {
      if(Q == 0) {x_times_a <- 0} else {
        x_times_a <- safe_prod(d$mode2.X[j,,drop=F], m$mode2.A.mean[,r,drop=F])
      }
      m$mode2.H.mean[j,r] <- m$mode2.H.var[j,r] * ((1/m$sigma2) *
        sum(outer(m$mode1.H.mean[,r], m$mode3.H.mean[,r]) *
            (d$resp[,j,] - (sweep(m$mode1.H.mean[,-r,drop=F], MARGIN=2, m$mode2.H.mean[j,-r], '*') %*%
            t(m$mode3.H.mean[,-r,drop=F]))), na.rm=T) + 1/m$m2.sigma2 * x_times_a)
    }
  }
}
