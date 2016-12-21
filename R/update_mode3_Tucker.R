#' Update the third mode in a Tucker model.
#'
#' Update is performed in place to avoid memory issues. There is no return value.
#' 
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
#' update_mode3_Tucker(m=toy.model, d=train.data, params=model.params)

update_mode3_Tucker <- function(m, d, params, batch.samps) {
  library(foreach) # namespace must be loaded for %do% and %dopar%
  
  # Make all param variables available locally
  for(i in 1:length(params)) {
    assign(names(params)[i], params[i][[1]])
  }

    # Set the stepsize rho
  tau <- 1        # Make this a parameter
  kappa <- 0.9    # Make this a parameter (or adaptive)
  rho <- (m$iter + 1)^(-kappa)
  
  I <- dim(d$resp)[1]; J <- dim(d$resp)[2]; K <- dim(d$resp)[3]
  R1 <- ncol(m$mode1.A.mean); R2 <- ncol(m$mode2.A.mean); R3 <- ncol(m$mode3.A.mean)
  core1 <- ncol(m$mode1.H.mean); core2 <- ncol(m$mode2.H.mean); core3 <- ncol(m$mode3.H.mean)
  S <- nrow(m$mode3.A.mean)
  
  if(verbose) print("Updating latent (H) matrix for mode 3")
  # Update the variance first.
  # Copy data so all of m and d aren't sent out to worker nodes
  mode1.H.mean <- m$mode1.H.mean
  mode1.H.var  <- m$mode1.H.var
  mode2.H.mean <- m$mode2.H.mean
  mode2.H.var  <- m$mode2.H.var
  sigma2 <- m$sigma2
  m3.sigma2 <- m$m3.sigma2
  
  if(H3.intercept) {
    m$mode3.H.var[,-1] <- foreach(delta=iterators::iapply(d$delta, 3), .combine='rbind') %:%
      foreach(core.mean=iterators::iapply(m$core.mean[,,-1,drop=F], 3), 
              core.var=iterators::iapply(m$core.var[,,-1,drop=F], 3), .combine='c') %do% {
                sum1 <- matrix(0, I, J); sum2 <- matrix(0, I, J)
                sum3 <- matrix(0, I, J); sum4 <- matrix(0, I, J)
                for(r1 in 1:core1) for(r2 in 1:core2) {
                  sum1 <- sum1 + core.mean[r1,r2] * outer(mode1.H.mean[,r1], mode2.H.mean[,r2]) *
                    (mode1.H.mean[,-r1,drop=F] %*% core.mean[-r1,-r2,drop=F]) %*% t(mode2.H.mean[,-r2,drop=F])
                  sum2 <- sum2 + core.mean[r1,r2] *
                    outer((mode1.H.mean[,r1]^2 + mode1.H.var[,r1]), mode2.H.mean[,r2]) *
                    matrix(mode2.H.mean[,-r2,drop=F] %*% core.mean[r1,-r2], I, J, byrow=T)
                  sum3 <- sum3 + core.mean[r1,r2] *
                    outer(mode1.H.mean[,r1], (mode2.H.mean[,r2]^2 + mode2.H.var[,r2])) *
                    matrix(mode1.H.mean[,-r1,drop=F] %*% core.mean[-r1,r2], I, J)
                  sum4 <- sum4 + (core.mean[r1,r2]^2 + core.var[r1,r2]) *
                    outer((mode1.H.mean[,r1]^2 + mode1.H.var[,r1]),
                          (mode2.H.mean[,r2]^2 + mode2.H.var[,r2]))
                }
                1/((1/sigma2) * sum(delta * (sum1 + sum2 + sum3 + sum4)) + (1/m3.sigma2))
              }
  } else {
    m$mode3.H.var[,] <- foreach(delta=iterators::iapply(d$delta, 3), .combine='rbind') %:%
      foreach(core.mean=iterators::iapply(m$core.mean[,,], 3), 
              core.var=iterators::iapply(m$core.var[,,], 3), .combine='c') %do% {
                sum1 <- matrix(0, I, J); sum2 <- matrix(0, I, J)
                sum3 <- matrix(0, I, J); sum4 <- matrix(0, I, J)
                for(r1 in 1:core1) for(r2 in 1:core2) {
                  sum1 <- sum1 + core.mean[r1,r2] * outer(mode1.H.mean[,r1], mode2.H.mean[,r2]) *
                    (mode1.H.mean[,-r1,drop=F] %*% core.mean[-r1,-r2]) %*% t(mode2.H.mean[,-r2,drop=F])
                  sum2 <- sum2 + core.mean[r1,r2] *
                    outer((mode1.H.mean[,r1]^2 + mode1.H.var[,r1]), mode2.H.mean[,r2]) *
                    matrix(mode2.H.mean[,-r2,drop=F] %*% core.mean[r1,-r2], I, J, byrow=T)
                  sum3 <- sum3 + core.mean[r1,r2] *
                    outer(mode1.H.mean[,r1], (mode2.H.mean[,r2]^2 + mode2.H.var[,r2])) *
                    matrix(mode1.H.mean[,-r1,drop=F] %*% core.mean[-r1,r2], I, J)
                  sum4 <- sum4 + (core.mean[r1,r2]^2 + core.var[r1,r2]) *
                    outer((mode1.H.mean[,r1]^2 + mode1.H.var[,r1]),
                          (mode2.H.mean[,r2]^2 + mode2.H.var[,r2]))
                }
                1/((1/sigma2) * sum(delta * (sum1 + sum2 + sum3 + sum4)) + (1/m3.sigma2))
              }
  }
  
  # Update means
  
  # Compute as much as possible outside of loops
  if(S == 0) {
    x_times_a <- matrix(0, K, R3)
  } else x_times_a <- safe_prod(d$mode3.X, m$mode3.A.mean)
  if(H3.intercept) x_times_a <- cbind(1, x_times_a)
  
  sum0 <- rTensor::ttl(rTensor::as.tensor(m$core.mean), list(m$mode1.H.mean, m$mode2.H.mean), c(1,2))@data
  
  # Update the mean parameters (m$mode3.H.mean)
  core.mean <- m$core.mean
  
  dm <- dimnames(m$mode3.H.mean)
  # Loop is over samples (K)
  if(H3.intercept) R3.rng <- 2:core3 else R3.rng <- 1:core3 # Don't update the constant column
  m$mode3.H.mean <- foreach(mode3.H.mean = iterators::iter(m$mode3.H.mean, by='row'), 
                            mode3.H.var = iterators::iter(m$mode3.H.var, by='row'),
                            resp = iterators::iapply(d$resp, 3),
                            x_t_a = iterators::iter(x_times_a, by='row'), .combine='rbind') %do% {
                              for(r3 in R3.rng) { 
                                big_sum <- matrix(0,I,J)
                                for(r3. in (1:core3)[-r3]) {
                                  sum1 <- matrix(0,I,J); sum2 <- matrix(0,I,J)
                                  sum3 <- matrix(0,I,J); sum4 <- matrix(0,I,J)
                                  for(r1 in 1:core1) for(r2 in 1:core2) {
                                    sum1 <- sum1 + core.mean[r1,r2,r3] * outer(mode1.H.mean[,r1], mode2.H.mean[,r2]) *
                                      (mode1.H.mean[,-r1,drop=F] %*% core.mean[-r1,-r2,r3.]) %*% t(mode2.H.mean[,-r2,drop=F])
                                    sum2 <- sum2 + core.mean[r1,r2,r3] *
                                      outer((mode1.H.mean[,r1]^2 + mode1.H.var[,r1]), mode2.H.mean[,r2]) *
                                      matrix(mode2.H.mean[,-r2,drop=F] %*% core.mean[r1,-r2,r3.], I, J, byrow=T)
                                    sum3 <- sum3 + core.mean[r1,r2,r3] *
                                      outer(mode1.H.mean[,r1], (mode2.H.mean[,r2]^2 + mode2.H.var[,r2])) *
                                      matrix(mode1.H.mean[,-r1,drop=F] %*% core.mean[-r1,r2,r3.], I, J)
                                    sum4 <- sum4 + core.mean[r1,r2,r3] * core.mean[r1,r2,r3.] *
                                      outer((mode1.H.mean[,r1]^2 + mode1.H.var[,r1]),
                                            (mode2.H.mean[,r2]^2 + mode2.H.var[,r2]))
                                  }
                                  big_sum <- big_sum + mode3.H.mean[r3.] * (sum1 + sum2 + sum3 + sum4)
                                }
                                mode3.H.mean[r3] <- mode3.H.var[r3] * ((1/sigma2) * sum((resp * sum0[,,r3]) - big_sum, na.rm=T) +
                                                                         (1/m3.sigma2) * x_t_a[r3])
                              }
                              mode3.H.mean
                            }
  dimnames(m$mode3.H.mean) <- dm
  
  # If the intercept term is removed change A3.intercept  
  A3.intercept <- ifelse('const' %in% rownames(m$mode3.A.mean), T, F)
  
  # TODO: check this
  if(S != 0) { # If there is no input data, skip updates for lambda and A
    if(verbose) print("Updating prior lambda vector for mode 3")
    
    m3.A.var <- matrix(0, S, R3)
    for(r3 in 1:R3) m3.A.var[,r3] <- diag(m$mode3.A.cov[,,r3])
    if(row.share) {
      m$mode3.lambda.scale <- 1/(.5*(rowSums(m$mode3.A.mean^2 + m3.A.var)) + 1/m$m3.beta)
    } else m$mode3.lambda.scale <- 1/(.5*(m$mode3.A.mean^2 + m3.A.var) + 1/m$m3.beta)
    
    if(verbose) print("Updating projection (A) matrix for mode 3")
    # Update mode3.A covariance parameters. They only rely on X and lambdas
    lambda.exp <- m$mode3.lambda.shape * m$mode3.lambda.scale
    for(r3 in 1:R3) {
      if(row.share) {
        m$mode3.A.cov[,,r3] <- chol2inv(chol(diag(lambda.exp) + (1/m$m3.sigma2) * m$m3Xm3X))
      } else
        m$mode3.A.cov[,,r3] <- chol2inv(chol(diag(lambda.exp[,r3]) + (1/m$m3.sigma2) * m$m3Xm3X))
    }

    # Update each column of A means
    if(A3.intercept) {
      if(H3.intercept) {
        for(r3 in 1:R3) m$mode3.A.mean[,r3] <- (1/m$m3.sigma2) * 
            (m$mode3.A.cov[,,r3] %*% t(cbind(1,d$mode3.X)) %*% m$mode3.H.mean[,r3+1])
      } else {
        for(r3 in 1:R3) m$mode3.A.mean[,r3] <- (1/m$m3.sigma2) * 
            (m$mode3.A.cov[,,r3] %*% t(cbind(1,d$mode3.X)) %*% m$mode3.H.mean[,r3])
      }
    } else {
      if(H3.intercept) {
        for(r3 in 1:R3) m$mode3.A.mean[,r3] <- (1/m$m3.sigma2) * 
            (m$mode3.A.cov[,,r3] %*% t(d$mode3.X) %*% m$mode3.H.mean[,r3+1])
      } else {
        for(r3 in 1:R3) m$mode3.A.mean[,r3] <- (1/m$m3.sigma2) * 
            (m$mode3.A.cov[,,r3] %*% t(d$mode3.X) %*% m$mode3.H.mean[,r3])
      }
    }
  }
}