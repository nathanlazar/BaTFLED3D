#' Update the second mode in a Tucker model.
#'
#' Update is performed in place to avoid memory issues. There is no return value.
#' 
#' @importFrom foreach foreach
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
#' update_mode2_Tucker(m=toy.model, d=train.data, params=model.params)

update_mode2_Tucker <- function(m, d, params) {
  # Make all param variables available locally
  for(i in 1:length(params)) {
    assign(names(params)[i], params[i][[1]])
  }
  
  I <- dim(d$resp)[1]; J <- dim(d$resp)[2]; K <- dim(d$resp)[3]
  R1 <- ncol(m$mode1.A.mean); R2 <- ncol(m$mode2.A.mean); R3 <- ncol(m$mode3.A.mean)
  core1 <- ncol(m$mode1.H.mean); core2 <- ncol(m$mode2.H.mean); core3 <- ncol(m$mode3.H.mean)
  Q <- nrow(m$mode2.A.mean)
  
  # If the intercept term is removed change A2.intercept  
  A2.intercept <- ifelse('const' %in% rownames(m$mode2.A.mean), T, F)

  if(Q != 0) { # If there is no input data, skip updates for lambda and A
    if(verbose) print("Updating prior lambda vector for mode 2")
    
    m2.A.var <- matrix(0, Q, R2)
    for(r2 in 1:R2) m2.A.var[,r2] <- diag(m$mode2.A.cov[,,r2])
    if(row.share) {
      m$mode2.lambda.scale <- 1/(.5*(rowSums(m$mode2.A.mean^2 + m2.A.var)) + 1/m$m2.beta)
    } else m$mode2.lambda.scale <- 1/(.5*(m$mode2.A.mean^2 + m2.A.var) + 1/m$m2.beta)
    
    if(verbose) print("Updating projection (A) matrix for mode 2")
    # Update mode2.A covariance parameters. They only rely on X and lambdas
    lambda.exp <- m$mode2.lambda.shape * m$mode2.lambda.scale
    for(r2 in 1:R2) {
      if(row.share) {
        m$mode2.A.cov[,,r2] <- chol2inv(chol(diag(lambda.exp) + (1/m$m2.sigma2) * m$m2Xm2X))
      } else
        m$mode2.A.cov[,,r2] <- chol2inv(chol(diag(lambda.exp[,r2]) + (1/m$m2.sigma2) * m$m2Xm2X))
    }

    # Update each column of A
    if(A2.intercept) {
      if(H2.intercept) {
        for(r2 in 1:R2) m$mode2.A.mean[,r2] <- (1/m$m2.sigma2) * 
            (m$mode2.A.cov[,,r2] %*% t(cbind(1,d$mode2.X)) %*% m$mode2.H.mean[,r2+1])
      } else {
        for(r2 in 1:R2) m$mode2.A.mean[,r2] <- (1/m$m2.sigma2) * 
            (m$mode2.A.cov[,,r2] %*% t(cbind(1,d$mode2.X)) %*% m$mode2.H.mean[,r2])
      }
    } else {
      if(H2.intercept) {
        for(r2 in 1:R2) m$mode2.A.mean[,r2] <- (1/m$m2.sigma2) * 
            (m$mode2.A.cov[,,r2] %*% t(d$mode2.X) %*% m$mode2.H.mean[,r2+1])
      } else {
        for(r2 in 1:R2) m$mode2.A.mean[,r2] <- (1/m$m2.sigma2) * 
            (m$mode2.A.cov[,,r2] %*% t(d$mode2.X) %*% m$mode2.H.mean[,r2])
      }
    }
  }
  
  if(verbose) print("Updating latent (H) matrix for mode 2")
  # Update the variance first
  # Copy data so all of m and d aren't sent out to worker nodes
  mode1.H.mean <- m$mode1.H.mean
  mode1.H.var  <- m$mode1.H.var
  mode3.H.mean <- m$mode3.H.mean
  mode3.H.var  <- m$mode3.H.var
  sigma2 <- m$sigma2
  m2.sigma2 <- m$m2.sigma2
  
  if(H2.intercept) {
    m$mode2.H.var[,-1] <- foreach(delta=iterators::iapply(d$delta, 2), .combine='rbind') %:%
      foreach(core.mean=iterators::iapply(m$core.mean[,-1,,drop=F], 2), 
              core.var=iterators::iapply(m$core.var[,-1,,drop=F], 2), .combine='c') %dopar% {
        sum1 <- matrix(0, I, K); sum2 <- matrix(0, I, K)
        sum3 <- matrix(0, I, K); sum4 <- matrix(0, I, K)
        for(r1 in 1:core1) for(r3 in 1:core3) {
          sum1 <- sum1 + core.mean[r1,r3] * outer(mode1.H.mean[,r1], mode3.H.mean[,r3]) *
            (mode1.H.mean[,-r1,drop=F] %*% core.mean[-r1,-r3]) %*% t(mode3.H.mean[,-r3,drop=F])
          sum2 <- sum2 + core.mean[r1,r3] *
            outer((mode1.H.mean[,r1]^2 + mode1.H.var[,r1]), mode3.H.mean[,r3]) *
            matrix(mode3.H.mean[,-r3,drop=F] %*% core.mean[r1,-r3], I, K, byrow=T)
          sum3 <- sum3 + core.mean[r1,r3] *
            outer(mode1.H.mean[,r1], (mode3.H.mean[,r3]^2 + mode3.H.var[,r3])) *
            matrix(mode1.H.mean[,-r1,drop=F] %*% core.mean[-r1,r3], I, K)
          sum4 <- sum4 + (core.mean[r1,r3]^2 + core.var[r1,r3]) *
            outer((mode1.H.mean[,r1]^2 + mode1.H.var[,r1]),
                  (mode3.H.mean[,r3]^2 + mode3.H.var[,r3]))
        }
        1/((1/sigma2) * sum(delta * (sum1 + sum2 + sum3 + sum4)) + (1/m2.sigma2))
      }
  } else {
    dm <- dimnames(m$mode2.H.var)
    m$mode2.H.var <- foreach(delta=iterators::iapply(d$delta, 2), .combine='rbind') %:%
      foreach(core.mean=iterators::iapply(m$core.mean, 2), 
              core.var=iterators::iapply(m$core.var, 2), .combine='c') %dopar% {
        sum1 <- matrix(0, I, K); sum2 <- matrix(0, I, K)
        sum3 <- matrix(0, I, K); sum4 <- matrix(0, I, K)
        if(is.null(dim(delta))) delta <- matrix(delta, ncol=1)
        if(is.null(dim(core.mean))) core.mean <- matrix(core.mean, ncol=1)
        if(is.null(dim(core.var))) core.var <- matrix(core.var, ncol=1)
        for(r1 in 1:core1) for(r3 in 1:core3) {
          sum1 <- sum1 + core.mean[r1,r3] * outer(mode1.H.mean[,r1], mode3.H.mean[,r3]) *
            (mode1.H.mean[,-r1,drop=F] %*% core.mean[-r1,-r3,drop=F]) %*% t(mode3.H.mean[,-r3,drop=F])
          sum2 <- sum2 + core.mean[r1,r3] *
            outer((mode1.H.mean[,r1]^2 + mode1.H.var[,r1]), mode3.H.mean[,r3]) *
            matrix(mode3.H.mean[,-r3] %*% core.mean[r1,-r3], I, K, byrow=T)
          sum3 <- sum3 + core.mean[r1,r3] *
            outer(mode1.H.mean[,r1], (mode3.H.mean[,r3]^2 + mode3.H.var[,r3])) *
            matrix(mode1.H.mean[,-r1] %*% core.mean[-r1,r3], I, K)
          sum4 <- sum4 + (core.mean[r1,r3]^2 + core.var[r1,r3]) *
            outer((mode1.H.mean[,r1]^2 + mode1.H.var[,r1]),
                  (mode3.H.mean[,r3]^2 + mode3.H.var[,r3]))
        }
        1/((1/sigma2) * sum(delta * (sum1 + sum2 + sum3 + sum4)) + (1/m2.sigma2))
    }
    if(is.null(dim(m$mode2.H.var)))
      m$mode2.H.var <- matrix(m$mode2.H.var, 1,1)
    dimnames(m$mode2.H.var) <- dm
  }
  
  # Compute as much as possible outside of loops
  if(Q == 0) {
    x_times_a <- matrix(0, J, R2)
  } else x_times_a <- safe_prod(d$mode2.X, m$mode2.A.mean)
  if(H2.intercept) x_times_a <- cbind(1, x_times_a)
  sum0 <- rTensor::ttl(rTensor::as.tensor(m$core.mean), list(m$mode1.H.mean, m$mode3.H.mean), c(1,3))@data
  
  # Update the mean parameters (m$mode2.H.mean)
  core.mean <- m$core.mean

  dm <- dimnames(m$mode2.H.mean)
  # Loop is over samples (J)
  if(H2.intercept) R2.rng <- 2:core2 else R2.rng <- 1:core2 # Don't update the constant column
  m$mode2.H.mean <- foreach(mode2.H.mean = iterators::iter(m$mode2.H.mean, by='row'), 
                            mode2.H.var = iterators::iter(m$mode2.H.var, by='row'),
                            resp = iterators::iapply(d$resp, 2),
                            x_t_a = iterators::iter(x_times_a, by='row'), .combine='rbind') %dopar% {
    for(r2 in R2.rng) { 
      big_sum <- matrix(0,I,K)
      for(r2. in (1:core2)[-r2]) {
        sum1 <- matrix(0,I,K); sum2 <- matrix(0,I,K)
        sum3 <- matrix(0,I,K); sum4 <- matrix(0,I,K)
        for(r1 in 1:core1) for(r3 in 1:core3) {
          sum1 <- sum1 + core.mean[r1,r2,r3] * outer(mode1.H.mean[,r1], mode3.H.mean[,r3]) *
            (mode1.H.mean[,-r1,drop=F] %*% core.mean[-r1,r2.,-r3]) %*% t(mode3.H.mean[,-r3,drop=F])
          sum2 <- sum2 + core.mean[r1,r2,r3] *
            outer((mode1.H.mean[,r1]^2 + mode1.H.var[,r1]), mode3.H.mean[,r3]) *
            matrix(mode3.H.mean[,-r3] %*% matrix(core.mean[r1,r2.,-r3],core3-1,1), I, K, byrow=T)
          sum3 <- sum3 + core.mean[r1,r2,r3] *
            outer(mode1.H.mean[,r1], (mode3.H.mean[,r3]^2 + mode3.H.var[,r3])) *
            matrix(mode1.H.mean[,-r1,drop=F] %*% core.mean[-r1,r2.,r3], I, K)
          sum4 <- sum4 + core.mean[r1,r2,r3] * core.mean[r1,r2.,r3] *
            outer((mode1.H.mean[,r1]^2 + mode1.H.var[,r1]),
                  (mode3.H.mean[,r3]^2 + mode3.H.var[,r3]))
        }
        big_sum <- big_sum + mode2.H.mean[r2.] * (sum1 + sum2 + sum3 + sum4)
      }
      mode2.H.mean[r2] <- mode2.H.var[r2] * ((1/sigma2) * sum((resp * sum0[,r2,]) - big_sum, na.rm=T) +
                                               (1/m2.sigma2) * x_t_a[r2])
    }
    mode2.H.mean
  }
  dimnames(m$mode2.H.mean) <- dm
}