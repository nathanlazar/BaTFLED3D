#' Calculate the lower bound of the log likelihood for a trained CP model
#'
#' @export
#' @param m object of the class \code{CP_model}
#' @param d object of the class \code{input_data}
#' @return Returns a numerical value (should be negative)
#'
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
#' train(d=train.data, m=toy.model, new.iter=1, params=model.params)
#' 
#' lower_bnd_CP(toy.model, train.data)

lower_bnd_CP <- function(m, d) {
  # Calculate the lower bound
  I <- dim(d$resp)[1]; J <- dim(d$resp)[2]; K <- dim(d$resp)[3]
  P <- nrow(m$mode1.A.mean); Q <- nrow(m$mode2.A.mean); S <- nrow(m$mode3.A.mean)
  R <- ncol(m$mode1.H.mean)

  log.2.pi <- log(2*pi)
  
  if(P != 0) {
    # Add constant column to X if necessary and rename so d isn't changed
    if(ncol(d$mode1.X) == (nrow(m$mode1.A.mean)-1)) {mode1.X <- cbind(1, d$mode1.X)} else mode1.X <- d$mode1.X
    
    exp.mode1.lambda <- m$mode1.lambda.shape * m$mode1.lambda.scale
    exp.log.mode1.lambda <- digamma(m$mode1.lambda.shape) + safe_log(m$mode1.lambda.scale)
    
    exp.p.mode1.lambda <- sum((m$m1.alpha-1)*exp.log.mode1.lambda - (1/m$m1.beta)*exp.mode1.lambda -
                                m$m1.alpha * safe_log(m$m1.beta) - lgamma(m$m1.alpha))

    m1.A.var <- matrix(0, P, R) # This could be stored earlier
    for(r in 1:R) m1.A.var[,r] <- diag(m$mode1.A.cov[,,r])
    
    exp.p.mode1.A <- -.5*P*R*log.2.pi + .5*sum(exp.log.mode1.lambda - exp.mode1.lambda * 
                                               (m$mode1.A.mean^2 + m1.A.var))

    Xm1Am1AX <- matrix(0, I, R)
    for(i in 1:I) for(r in 1:R) 
      Xm1Am1AX[i,r] <- mode1.X[i,,drop=F] %*% 
        (outer(m$mode1.A.mean[,r], m$mode1.A.mean[,r]) + m$mode1.A.cov[,,r]) %*% t(mode1.X[i,,drop=F])
    
    exp.p.mode1.H <- -.5*I*R * safe_log(2*pi*m$m1.sigma2) - 
      1/(2*m$m1.sigma2) * sum(m$mode1.H.mean^2 + m$mode1.H.var) +
      1/(m$m1.sigma2) * sum(m$mode1.H.mean * (mode1.X %*% m$mode1.A.mean)) -
      1/(2*m$m1.sigma2) * sum(Xm1Am1AX)
    
    exp.q.mode1.lambda <- sum(-m$mode1.lambda.shape - safe_log(m$mode1.lambda.scale) -
                              safe_log(gamma(m$mode1.lambda.shape)) -
                              (1-m$mode1.lambda.shape) * digamma(m$mode1.lambda.shape))

    m1.A.dets <- rep(0, R)
    for(r in 1:R) m1.A.dets[r] <- determinant(m$mode1.A.cov[,,r], logarithm=T)$modulus
    exp.q.mode1.A <- -sum(P/2*(log.2.pi+1) -.5*m1.A.dets)
  } else {
    exp.p.mode1.lambda <- 0
    exp.p.mode1.A <- 0
    exp.p.mode1.H <- -.5*I*R * safe_log(2*pi*m$m1.sigma2) - 1/(2*m$m1.sigma2) * sum(m$mode1.H.mean^2 + m$mode1.H.var)
    exp.q.mode1.lambda <- 0
    exp.q.mode1.A <- 0
  }
  exp.q.mode1.H <- -.5*I*R*(log.2.pi+1) -.5*sum(safe_log(m$mode1.H.var))

  if(Q != 0) {
    # Add constant column to X if necessary and rename so d isn't changed
    if(ncol(d$mode2.X) == (nrow(m$mode2.A.mean)-1)) {mode2.X <- cbind(1, d$mode2.X)} else mode2.X <- d$mode2.X

    exp.mode2.lambda <- m$mode2.lambda.shape * m$mode2.lambda.scale
    exp.log.mode2.lambda <- digamma(m$mode2.lambda.shape) + safe_log(m$mode2.lambda.scale)
    
    exp.p.mode2.lambda <- sum((m$m2.alpha-1)*exp.log.mode2.lambda - (1/m$m2.beta)*exp.mode2.lambda -
                                m$m2.alpha * safe_log(m$m2.beta) - lgamma(m$m2.alpha))
    
    m2.A.var <- matrix(0, Q, R) # This could be stored earlier
    for(r in 1:R) m2.A.var[,r] <- diag(m$mode2.A.cov[,,r])
    
    exp.p.mode2.A <- -.5*Q*R*log.2.pi + .5*sum(exp.log.mode2.lambda - exp.mode2.lambda * 
                                                 (m$mode2.A.mean^2 + m2.A.var))
    
    Xm2Am2AX <- matrix(0, J, R)
    for(j in 1:J) for(r in 1:R) 
      Xm2Am2AX[j,r] <- mode2.X[j,,drop=F] %*% 
        (outer(m$mode2.A.mean[,r], m$mode2.A.mean[,r]) + m$mode2.A.cov[,,r]) %*% t(mode2.X[j,,drop=F])
    
      
    exp.p.mode2.H <- -.5*J*R * safe_log(2*pi*m$m2.sigma2) - 
      1/(2*m$m2.sigma2) * sum(m$mode2.H.mean^2 + m$mode2.H.var) +
      1/(m$m2.sigma2) * sum(m$mode2.H.mean * (mode2.X %*% m$mode2.A.mean)) -
      1/(2*m$m2.sigma2) * sum(Xm2Am2AX)
    
    exp.q.mode2.lambda <- sum(-m$mode2.lambda.shape - safe_log(m$mode2.lambda.scale) -
                                safe_log(gamma(m$mode2.lambda.shape)) -
                                (1-m$mode2.lambda.shape) * digamma(m$mode2.lambda.shape))
    
    m2.A.dets <- rep(0, R)
    for(r in 1:R) m2.A.dets[r] <- determinant(m$mode2.A.cov[,,r], logarithm=T)$modulus
    exp.q.mode2.A <- -sum(Q/2*(log.2.pi+1) -.5*m2.A.dets)
  } else {
    exp.p.mode2.lambda <- 0
    exp.p.mode2.A <- 0
    exp.p.mode2.H <- -.5*J*R * safe_log(2*pi*m$m2.sigma2) - 1/(2*m$m2.sigma2) * sum(m$mode2.H.mean^2 + m$mode2.H.var)
    exp.q.mode2.lambda <- 0
    exp.q.mode2.A <- 0
  }
  exp.q.mode2.H <- -.5*J*R*(log.2.pi+1) -.5*sum(safe_log(m$mode2.H.var))
  
  if(S != 0) {
    # Add constant column to X if necessary and rename so d isn't changed
    if(ncol(d$mode3.X) == (nrow(m$mode3.A.mean)-1)) {mode3.X <- cbind(1, d$mode3.X)} else mode3.X <- d$mode3.X
    
    exp.mode3.lambda <- m$mode3.lambda.shape * m$mode3.lambda.scale
    exp.log.mode3.lambda <- digamma(m$mode3.lambda.shape) + safe_log(m$mode3.lambda.scale)
    
    exp.p.mode3.lambda <- sum((m$m3.alpha-1)*exp.log.mode3.lambda - (1/m$m3.beta)*exp.mode3.lambda -
                               m$m3.alpha * safe_log(m$m3.beta) - lgamma(m$m3.alpha))
    
    m3.A.var <- matrix(0, S, R) # This could be stored earlier
    for(r in 1:R) m3.A.var[,r] <- diag(m$mode3.A.cov[,,r])
    
    exp.p.mode3.A <- -.5*S*R*log.2.pi + .5*sum(exp.log.mode3.lambda - exp.mode3.lambda * 
                                                 (m$mode3.A.mean^2 + m3.A.var))
    
    Xm3Am3AX <- matrix(0, K, R)
    for(k in 1:K) for(r in 1:R) 
      Xm3Am3AX[j,r] <- mode3.X[k,,drop=F] %*% 
        (outer(m$mode3.A.mean[,r], m$mode3.A.mean[,r]) + m$mode3.A.cov[,,r]) %*% t(mode3.X[k,,drop=F])
    
    exp.p.mode3.H <- -.5*K*R * safe_log(2*pi*m$m3.sigma2) - 
      1/(2*m$m3.sigma2) * sum(m$mode3.H.mean^2 + m$mode3.H.var) +
      1/(m$m3.sigma2) * sum(m$mode3.H.mean * (mode3.X %*% m$mode3.A.mean)) -
      1/(2*m$m3.sigma2) * sum(Xm3Am3AX)
    
    exp.q.mode3.lambda <- sum(-m$mode3.lambda.shape - safe_log(m$mode3.lambda.scale) -
                                safe_log(gamma(m$mode3.lambda.shape)) -
                                (1-m$mode3.lambda.shape) * digamma(m$mode3.lambda.shape))
    
    m3.A.dets <- rep(0, R)
    for(r in 1:R) m3.A.dets[r] <- determinant(m$mode3.A.cov[,,r], logarithm=T)$modulus
    exp.q.mode3.A <- -sum(S/2*(log.2.pi+1) -.5*m3.A.dets)
  } else {
    exp.p.mode3.lambda <- 0
    exp.p.mode3.A <- 0
    exp.p.mode3.H <- -.5*K*R * safe_log(2*pi*m$m3.sigma2) - 1/(2*m$m3.sigma2) * sum(m$mode3.H.mean^2 + m$mode3.H.var)
    exp.q.mode3.lambda <- 0
    exp.q.mode3.A <- 0
  }
  exp.q.mode3.H <- -.5*K*R*(log.2.pi+1) -.5*sum(safe_log(m$mode3.H.var))
  
  core <- array(0, dim=c(R,R,R))
  for(r in 1:R) core[r,r,r] <- 1
  
  minus.r <- array(0, dim=c(I,J,K,R))
  for(r in 1:R) 
    minus.r[,,,r] <- mult_3d(core[r,r,r,drop=F], m$mode1.H.mean[,r,drop=F], m$mode2.H.mean[,r,drop=F], m$mode3.H.mean[,r,drop=F]) *
                     mult_3d(core[-r,-r,-r], m$mode1.H.mean[,-r], m$mode2.H.mean[,-r], m$mode3.H.mean[,-r])
  
  exp.p.y <- -.5 * safe_log(2*pi*m$sigma2) * sum(d$delta) -
    (.5/m$sigma2) * sum(d$resp^2, na.rm=T) -
    (.5/m$sigma2) * sum(-2*d$resp * mult_3d(core, m$mode1.H.mean,  m$mode2.H.mean, m$mode3.H.mean) +
                        mult_3d(core, m$mode1.H.mean^2 + m$mode1.H.var, m$mode2.H.mean^2 + m$mode2.H.var, m$mode3.H.mean^2 + m$mode3.H.var) +
                        apply(minus.r, c(1,2,3), sum), na.rm=T)

  lower.bnd <- exp.p.mode1.lambda + exp.p.mode2.lambda + exp.p.mode3.lambda +
    exp.p.mode1.A + exp.p.mode2.A + exp.p.mode3.A +
    exp.p.mode1.H + exp.p.mode2.H + exp.p.mode3.H + exp.p.y -
    exp.q.mode1.lambda - exp.q.mode2.lambda - exp.q.mode3.lambda -
    exp.q.mode1.A - exp.q.mode2.A - exp.q.mode3.A -
    exp.q.mode1.H - exp.q.mode2.H - exp.q.mode3.H

  return(lower.bnd)
}
