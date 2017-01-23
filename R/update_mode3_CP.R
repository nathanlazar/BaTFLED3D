#' Update the third mode in a CP model.
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
#' update_mode3_CP(m=toy.model, d=train.data, params=model.params)

update_mode3_CP <- function(m, d, params) {
  # Make all param variables available locally
  for(i in 1:length(params)) {
    assign(names(params)[i], params[i][[1]])
  }
  
  K <- dim(d$resp)[3]
  S <- nrow(m$mode3.A.mean)
  R <- ncol(m$mode3.H.mean)
  
  # If the intercept term is removed change A1.intercept  
  A3.intercept <- ifelse('const' %in% rownames(m$mode3.A.mean), T, F)
  
  if(S != 0) { # If there is no input data, skip updates for lambda and A
    if(params$verbose) print("Updating prior lambda vector for mode 3")
    
    m3.A.var <- matrix(0, S, R)
    for(r in 1:R) m3.A.var[,r] <- diag(m$mode3.A.cov[,,r])
    if(params$row.share) {
      m$mode3.lambda.scale <- 1/(.5*(rowSums(m$mode3.A.mean^2 + m3.A.var)) + 1/m$m3.beta)
    } else m$mode3.lambda.scale <- 1/(.5*(m$mode3.A.mean^2 + m3.A.var) + 1/m$m3.beta)
    
    if(params$verbose) print("Updating projection (A) matrix for mode 3")
    # Update mode3.A covariance parameters. They only rely on X and lambdas
    lambda.exp <- m$mode3.lambda.shape * m$mode3.lambda.scale
    for(r in 1:R) {
      if(params$row.share) {
        m$mode3.A.cov[,,r] <- chol2inv(chol(diag(lambda.exp) + (1/m$m3.sigma2) * m$m3Xm3X))
      } else
        m$mode3.A.cov[,,r] <- chol2inv(chol(diag(lambda.exp[,r]) + (1/m$m3.sigma2) * m$m3Xm3X))
    }

    # Update A means
    if(A3.intercept) {
      for(r in 1:R) m$mode3.A.mean[,r] <- (1/m$m3.sigma2) * 
          (m$mode3.A.cov[,,r] %*% t(cbind(1,d$mode3.X)) %*% m$mode3.H.mean[,r])
    } else {
      for(r in 1:R) m$mode3.A.mean[,r] <- (1/m$m3.sigma2) * 
          (m$mode3.A.cov[,,r] %*% t(d$mode3.X) %*% m$mode3.H.mean[,r])
    }
  }
  
  if(params$verbose) print("Updating latent (H) matrix for mode 3")
  # Update the variance first. sapply vectorizes the updates for each row
  for(r in 1:R) {
    m$mode3.H.var[,r] <- sapply(1:K, function(k) 1/((1/m$sigma2) *
      sum(d$delta[,,k] * outer((m$mode1.H.mean[,r]^2 + m$mode1.H.var[,r]),
          (m$mode2.H.mean[,r]^2 + m$mode2.H.var[,r]))) + 1/m$m3.sigma2))
  }
  
  # Next update the mean
  for(k in 1:K) {
    for(r in 1:R) {
      if(S == 0) {x_times_a <- 0} else {
        x_times_a <- safe_prod(d$mode3.X[k,,drop=F], m$mode3.A.mean[,r,drop=F])
      }
      m$mode3.H.mean[k,r] <- m$mode3.H.var[k,r] * ((1/m$sigma2) *
        sum(outer(m$mode2.H.mean[,r], m$mode1.H.mean[,r]) *
            (t(d$resp[,,k]) - (sweep(m$mode2.H.mean[,-r], MARGIN=2, m$mode3.H.mean[k,-r], '*') %*%
            t(m$mode1.H.mean[,-r,drop=F]))), na.rm=T) + 1/m$m3.sigma2 * x_times_a)
    }
  }
}
