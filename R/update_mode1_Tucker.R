#' Update the first mode in a Tucker model.
#'
#' Update is performed in place to avoid memory issues. There is no return value.
#' 
#' @importFrom foreach %dopar%
#' @importFrom foreach %:%
#' @export
#' @param m A \code{Tucker_model} object created with \code{mk_model} 
#' @param d Input data object created with \code{input_data}
#' @param params List of parameters created with \code{get_model_params()}
#' @examples
#' data.params <- get_data_params(c('decomp=Tucker'))
#' toy <- mk_toy(data.params)
#' train.data <- input_data$new(mode1.X=toy$mode1.X[,-1],
#'                              mode2.X=toy$mode2.X[,-1],
#'                              mode3.X=toy$mode3.X,
#'                              resp=toy$resp)
#' model.params <- get_model_params(c('decomp=Tucker'))
#' toy.model <- mk_model(train.data, model.params)
#' toy.model$rand_init(model.params)
#'
#' update_mode1_Tucker(m=toy.model, d=train.data, params=model.params)

update_mode1_Tucker <- function(m, d, params) {
  # Make all param variables available locally
  for(i in 1:length(params)) {
    assign(names(params)[i], params[i][[1]])
  }
  
  I <- dim(d$resp)[1]; J <- dim(d$resp)[2]; K <- dim(d$resp)[3]
  R1 <- ncol(m$mode1.A.mean); R2 <- ncol(m$mode2.A.mean); R3 <- ncol(m$mode3.A.mean)
  core1 <- ncol(m$mode1.H.mean); core2 <- ncol(m$mode2.H.mean); core3 <- ncol(m$mode3.H.mean)
  P <- nrow(m$mode1.A.mean)

  # If the intercept term is removed change A1.intercept  
  A1.intercept <- ifelse('const' %in% rownames(m$mode1.A.mean), T, F)

  if(P != 0) { # If there is no input data, skip updates for lambda and A
    if(verbose) print("Updating prior lambda vector for mode 1")

    m1.A.var <- matrix(0, P, R1)
    for(r1 in 1:R1) m1.A.var[,r1] <- diag(m$mode1.A.cov[,,r1])
    if(row.share) {
      m$mode1.lambda.scale <- 1/(.5*(rowSums(m$mode1.A.mean^2 + m1.A.var)) + 1/m$m1.beta)
    } else m$mode1.lambda.scale <- 1/(.5*(m$mode1.A.mean^2 + m1.A.var) + 1/m$m1.beta)

    if(verbose) print("Updating projection (A) matrix for mode 1")
    # Update mode1.A covariance parameters. They only rely on X and lambdas
    lambda.exp <- m$mode1.lambda.shape * m$mode1.lambda.scale
    for(r1 in 1:R1) {
      if(row.share) {
        m$mode1.A.cov[,,r1] <- chol2inv(chol(diag(lambda.exp) + (1/m$m1.sigma2) * m$m1Xm1X))
      } else 
        m$mode1.A.cov[,,r1] <- chol2inv(chol(diag(lambda.exp[,r1]) + (1/m$m1.sigma2) * m$m1Xm1X))
    }
    
    # Update each column of A
    if(A1.intercept) {
      if(H1.intercept) {
        for(r1 in 1:R1) m$mode1.A.mean[,r1] <- (1/m$m1.sigma2) * 
            (m$mode1.A.cov[,,r1] %*% t(cbind(1,d$mode1.X)) %*% m$mode1.H.mean[,r1+1])
      } else {
        for(r1 in 1:R1) m$mode1.A.mean[,r1] <- (1/m$m1.sigma2) * 
            (m$mode1.A.cov[,,r1] %*% t(cbind(1,d$mode1.X)) %*% m$mode1.H.mean[,r1])
      }
    } else {
      if(H1.intercept) {
        for(r1 in 1:R1) m$mode1.A.mean[,r1] <- (1/m$m1.sigma2) * 
            (m$mode1.A.cov[,,r1] %*% t(d$mode1.X) %*% m$mode1.H.mean[,r1+1])
      } else {
        for(r1 in 1:R1) m$mode1.A.mean[,r1] <- (1/m$m1.sigma2) * 
            (m$mode1.A.cov[,,r1] %*% t(d$mode1.X) %*% m$mode1.H.mean[,r1])
      }
    }
  }
  
  # Update the variance and mean for the H factor matrices
  if(verbose) print("Updating latent (H) matrix for mode 1")
  
  # Updating variance matrix (m$mode1.H.var)
  # Copy data so all of m and d aren't sent out to worker nodes
  mode2.H.mean <- m$mode2.H.mean
  mode2.H.var <- m$mode2.H.var
  mode3.H.mean <- m$mode3.H.mean
  mode3.H.var <- m$mode3.H.var
  sigma2 <- m$sigma2
  m1.sigma2 <- m$m1.sigma2
  
  if(H1.intercept) {
    m$mode1.H.var[,-1] <- foreach(delta=iterators::iapply(d$delta, 1), .combine='rbind') %:%
      foreach(core.mean=iterators::iapply(m$core.mean[-1,,,drop=F], 1), 
        core.var=iterators::iapply(m$core.var[-1,,,drop=F], 1), .combine='c') %dopar% {
          sum1 <- matrix(0, J, K); sum2 <- matrix(0, J, K)
          sum3 <- matrix(0, J, K); sum4 <- matrix(0, J, K)
          for(r2 in 1:core2) for(r3 in 1:core3) {
            sum1 <- sum1 + core.mean[r2,r3] * outer(mode2.H.mean[,r2], mode3.H.mean[,r3]) *
              (mode2.H.mean[,-r2,drop=F] %*% core.mean[-r2,-r3]) %*% t(mode3.H.mean[,-r3,drop=F])
            sum2 <- sum2 + core.mean[r2,r3] *
              outer((mode2.H.mean[,r2]^2 + mode2.H.var[,r2]), mode3.H.mean[,r3]) *
              matrix(mode3.H.mean[,-r3,drop=F] %*% core.mean[r2,-r3], J, K, byrow=T)
            sum3 <- sum3 + core.mean[r2,r3] *
              outer(mode2.H.mean[,r2], (mode3.H.mean[,r3]^2 + mode3.H.var[,r3])) *
              matrix(mode2.H.mean[,-r2,drop=F] %*% core.mean[-r2,r3], J, K)
            sum4 <- sum4 + (core.mean[r2,r3]^2 + core.var[r2,r3]) *
              outer((mode2.H.mean[,r2]^2 + mode2.H.var[,r2]),
                    (mode3.H.mean[,r3]^2 + mode3.H.var[,r3]))
          }
          1/((1/sigma2) * sum(delta * (sum1 + sum2 + sum3 + sum4)) + (1/m1.sigma2))
      }
  } else {
    dm <- dimnames(m$mode1.H.var)
    m$mode1.H.var <- foreach(delta=iterators::iapply(d$delta, 1), .combine='rbind') %:%
      foreach(core.mean=iterators::iapply(m$core.mean[,,], 1), 
              core.var=iterators::iapply(m$core.var[,,], 1), .combine='c') %dopar% {
        sum1 <- matrix(0, J, K); sum2 <- matrix(0, J, K)
        sum3 <- matrix(0, J, K); sum4 <- matrix(0, J, K)
        if(is.null(dim(delta))) delta <- matrix(delta, ncol=1)
        if(is.null(dim(core.mean))) core.mean <- matrix(core.mean, ncol=1)
        if(is.null(dim(core.var))) core.var <- matrix(core.var, ncol=1)
        for(r2 in 1:core2) for(r3 in 1:core3) {
          sum1 <- sum1 + core.mean[r2,r3] * outer(mode2.H.mean[,r2], mode3.H.mean[,r3]) *
            (mode2.H.mean[,-r2,drop=F] %*% core.mean[-r2,-r3]) %*% t(mode3.H.mean[,-r3,drop=F])
          sum2 <- sum2 + core.mean[r2,r3] *
            outer((mode2.H.mean[,r2]^2 + mode2.H.var[,r2]), mode3.H.mean[,r3]) *
            matrix(mode3.H.mean[,-r3,drop=F] %*% core.mean[r2,-r3], J, K, byrow=T)
          sum3 <- sum3 + core.mean[r2,r3] *
            outer(mode2.H.mean[,r2], (mode3.H.mean[,r3]^2 + mode3.H.var[,r3])) *
            matrix(mode2.H.mean[,-r2,drop=F] %*% core.mean[-r2,r3], J, K)
          sum4 <- sum4 + (core.mean[r2,r3]^2 + core.var[r2,r3]) *
            outer((mode2.H.mean[,r2]^2 + mode2.H.var[,r2]),
                  (mode3.H.mean[,r3]^2 + mode3.H.var[,r3]))
        }
        1/((1/sigma2) * sum(delta * (sum1 + sum2 + sum3 + sum4)) + (1/m1.sigma2))
    }
    if(is.null(dim(m$mode1.H.var)))
      m$mode1.H.var <- matrix(m$mode1.H.var, 1,1)
    dimnames(m$mode1.H.var) <- dm
  }
  
  # Compute as much as possible outside of loops
  if(P == 0) {
    x_times_a <- matrix(0, I, R1)
  } else x_times_a <- safe_prod(d$mode1.X, m$mode1.A.mean)
  if(H1.intercept) x_times_a <- cbind(1, x_times_a)
  sum0 <- rTensor::ttl(rTensor::as.tensor(m$core.mean), list(m$mode2.H.mean, m$mode3.H.mean), c(2,3))@data
  
  # Update the mean parameters (m$mode1.H.mean)
  core.mean <- m$core.mean
  dm <- dimnames(m$mode1.H.mean)
  
  # Loop is over samples (I)
  if(H1.intercept) R1.rng <- 2:core1 else R1.rng <- 1:core1 # Don't update the constant column
  m$mode1.H.mean <- foreach(mode1.H.mean = iterators::iter(m$mode1.H.mean, by='row'), 
                            mode1.H.var = iterators::iter(m$mode1.H.var, by='row'),
                            resp = iterators::iapply(d$resp, 1),
                            x_t_a = iterators::iter(x_times_a, by='row'), .combine='rbind') %dopar% {
    for(r1 in R1.rng) { 
      big_sum <- matrix(0,J,K)
      for(r1. in (1:core1)[-r1]) {
        sum1 <- matrix(0, J, K); sum2 <- matrix(0, J, K)
        sum3 <- matrix(0, J, K); sum4 <- matrix(0, J, K)
        for(r2 in 1:core2) for(r3 in 1:core3) {
          sum1 <- sum1 + core.mean[r1,r2,r3] * outer(mode2.H.mean[,r2], mode3.H.mean[,r3]) *
            (mode2.H.mean[,-r2,drop=F] %*% core.mean[r1.,-r2,-r3]) %*% t(mode3.H.mean[,-r3,drop=F])
          sum2 <- sum2 + core.mean[r1,r2,r3] *
            outer((mode2.H.mean[,r2]^2 + mode2.H.var[,r2]), mode3.H.mean[,r3]) *
            matrix(mode3.H.mean[,-r3,drop=F] %*% matrix(core.mean[r1.,r2,-r3],core3-1,1) ,J,K,byrow=T)
          sum3 <- sum3 + core.mean[r1,r2,r3] *
            outer(mode2.H.mean[,r2], (mode3.H.mean[,r3]^2 + mode3.H.var[,r3])) *
            matrix(mode2.H.mean[,-r2,drop=F] %*% core.mean[r1.,-r2,r3], J, K)
          sum4 <- sum4 + core.mean[r1,r2,r3] * core.mean[r1.,r2,r3] *
            outer((mode2.H.mean[,r2]^2 + mode2.H.var[,r2]),
                  (mode3.H.mean[,r3]^2 + mode3.H.var[,r3]))
        }
        big_sum <- big_sum + mode1.H.mean[r1.] * (sum1 + sum2 + sum3 + sum4)
      }
      # move out of loop?
      mode1.H.mean[r1] <- mode1.H.var[r1] * ((1/sigma2) * sum((resp * sum0[r1,,]) - big_sum, na.rm=T) +
                                               (1/m1.sigma2) * x_t_a[r1])
    }
    mode1.H.mean
  }
  dimnames(m$mode1.H.mean) <- dm
}