#' Calculate the lower bound of the log likelihood for a trained Tucker model
#'
#' @export
#' @param m object of the class \code{Tucker_model}
#' @param d object of the class \code{input_data}
#' @return Returns a numerical value (should be negative)
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
#'
#' lower_bnd_Tucker(toy.model, train.data)

lower_bnd_Tucker <- function(m, d) {
  # Calculate the lower bound
  I <- dim(d$resp)[1]; J <- dim(d$resp)[2]; K <- dim(d$resp)[3]
  P <- nrow(m$mode1.A.mean); Q <- nrow(m$mode2.A.mean); S <- nrow(m$mode3.A.mean)
  R1 <- ncol(m$mode1.A.mean); R2 <- ncol(m$mode2.A.mean); R3 <- ncol(m$mode3.A.mean)
  core1 <- dim(m$core.mean)[1]; core2 <- dim(m$core.mean)[2]; core3 <- dim(m$core.mean)[3]

  # Remove the constant column from the H matrices if present and rename them so that 
  # m isn't modified in place.
  mode1.H.mean <- m$mode1.H.mean[,!grepl('const', colnames(m$mode1.H.mean)), drop=F]
  mode2.H.mean <- m$mode2.H.mean[,!grepl('const', colnames(m$mode2.H.mean)), drop=F]
  mode3.H.mean <- m$mode3.H.mean[,!grepl('const', colnames(m$mode3.H.mean)), drop=F]
  mode1.H.var <- m$mode1.H.var[,!grepl('const', colnames(m$mode1.H.var)), drop=F]
  mode2.H.var <- m$mode2.H.var[,!grepl('const', colnames(m$mode2.H.var)), drop=F]
  mode3.H.var <- m$mode3.H.var[,!grepl('const', colnames(m$mode3.H.var)), drop=F]
  
  log.2.pi <- log(2*pi)
  
  if(P != 0) {
    exp.mode1.lambda <- m$mode1.lambda.shape * m$mode1.lambda.scale
    exp.log.mode1.lambda <- digamma(m$mode1.lambda.shape) + safe_log(m$mode1.lambda.scale)
    #TODO: check whether this is necessary...
    # if(class(m$mode1.lambda.shape)=='numeric') { # variance is shared across rows
    exp.p.mode1.lambda <- sum((m$m1.alpha-1)*exp.log.mode1.lambda -
                                (1/m$m1.beta)*exp.mode1.lambda -
                                m$m1.alpha * safe_log(m$m1.beta) -
                                safe_log(gamma(m$m1.alpha)))
    # } else {

    m1.A.var <- matrix(0, P, R1) # This could be stored earlier
    for(r1 in 1:R1) m1.A.var[,r1] <- diagonal(m$mode1.A.cov[,,r1])
    
    exp.p.mode1.A <- -.5*P*R1*log.2.pi + .5*sum(exp.log.mode1.lambda -
      exp.mode1.lambda * (m$mode1.A.mean^2 + m1.A.var))

    # Sum used in for the expected value of the log probability for H
    cov.sum <- 0
    for(i in 1:I) for(r1 in 1:R1) {
      if(ncol(d$mode1.X) == (nrow(m$mode1.A.mean)-1)) {
        X_out_X <- outer(c(const=1, d$mode1.X[i,]), c(const=1, d$mode1.X[i,]))
      } else X_out_X <- outer(d$mode1.X[i,], d$mode1.X[i,])
      cov.sum <- cov.sum + sum(X_out_X *
        (outer(m$mode1.A.mean[,r1], m$mode1.A.mean[,r1]) + m$mode1.A.cov[,,r1]))
    }

    exp.p.mode1.H <- -.5*I*R1 * safe_log(2*pi*m$m1.sigma2) -
      1/(2*m$m1.sigma2) * sum(mode1.H.mean^2 + mode1.H.var) +
      1/(m$m1.sigma2) * sum(mode1.H.mean * safe_prod(d$mode1.X, m$mode1.A.mean)) -
      1/(2*m$m1.sigma2) * cov.sum

    exp.q.mode1.lambda <- sum(-m$mode1.lambda.shape - safe_log(m$mode1.lambda.scale) -
                              safe_log(gamma(m$mode1.lambda.shape)) -
                              (1-m$mode1.lambda.shape) * digamma(m$mode1.lambda.shape))

    m1.A.dets <- rep(0, R1)
    for(r1 in 1:R1) m1.A.dets[r1] <- determinant(m$mode1.A.cov[,,r1], logarithm=T)$modulus
    exp.q.mode1.A <- -sum(P/2*(log.2.pi+1) -.5*m1.A.dets)
  } else {
    exp.p.mode1.lambda <- 0
    exp.p.mode1.A <- 0
    exp.p.mode1.H <- -.5*I*R1 * safe_log(2*pi*m$m1.sigma2) - 1/(2*m$m1.sigma2) * sum(mode1.H.mean^2 + mode1.H.var)
    exp.q.mode1.lambda <- 0
    exp.q.mode1.A <- 0
  }
  exp.q.mode1.H <- -.5*I*R1*(log.2.pi+1) -.5*sum(safe_log(mode1.H.var))

  if(Q != 0) {
    exp.mode2.lambda <- m$mode2.lambda.shape * m$mode2.lambda.scale
    exp.log.mode2.lambda <- digamma(m$mode2.lambda.shape) + safe_log(m$mode2.lambda.scale)

    exp.p.mode2.lambda <- sum((m$m2.alpha-1)*exp.log.mode2.lambda -
                                (1/m$m2.beta)*exp.mode2.lambda -
                                m$m2.alpha * safe_log(m$m2.beta) -
                                safe_log(gamma(m$m2.alpha)))

    m2.A.var <- matrix(0, Q, R2) # This could be stored earlier
    for(r2 in 1:R2) m2.A.var[,r2] <- diagonal(m$mode2.A.cov[,,r2])

    exp.p.mode2.A <- -.5 * Q*R2 * log.2.pi + 1/2 * sum(exp.log.mode2.lambda -
      exp.mode2.lambda * (m$mode2.A.mean^2 + m2.A.var))

    # Sum used in for the expected value of the log probability for H
    cov.sum <- 0
    for(j in 1:J) for(r2 in 1:R2) {
      if(ncol(d$mode2.X) == (nrow(m$mode2.A.mean)-1)) {
        X_out_X <- outer(c(const=1, d$mode2.X[j,]), c(const=1, d$mode2.X[j,]))
      } else X_out_X <- outer(d$mode2.X[j,], d$mode2.X[j,])
      cov.sum <- cov.sum + sum(X_out_X *
        (outer(m$mode2.A.mean[,r2], m$mode2.A.mean[,r2]) + m$mode2.A.cov[,,r2]))
    }
    
    exp.p.mode2.H <- -.5*J*R2 * safe_log(2*pi*m$m2.sigma2) -
      1/(2*m$m2.sigma2) * sum(mode2.H.mean^2 + mode2.H.var) +
      1/(m$m2.sigma2) * sum(mode2.H.mean * safe_prod(d$mode2.X, m$mode2.A.mean)) -
      1/(2*m$m2.sigma2) * cov.sum

    exp.q.mode2.lambda <- sum(-m$mode2.lambda.shape - safe_log(m$mode2.lambda.scale) -
                                safe_log(gamma(m$mode2.lambda.shape)) -
                                (1-m$mode2.lambda.shape) * digamma(m$mode2.lambda.shape))

    m2.A.dets <- rep(0, R2)
    for(r2 in 1:R2) m2.A.dets[r2] <- determinant(m$mode2.A.cov[,,r2], logarithm=T)$modulus
    exp.q.mode2.A <- -sum(Q/2*(log.2.pi+1) -.5*m2.A.dets)
  } else {
    exp.p.mode2.lambda <- 0
    exp.p.mode2.A <- 0
    exp.p.mode2.H <- -.5*J*R2 * safe_log(2*pi*m$m2.sigma2) - 1/(2*m$m2.sigma2) * sum(mode2.H.mean^2 + mode2.H.var)
    exp.q.mode2.lambda <- 0
    exp.q.mode2.A <- 0
  }
  exp.q.mode2.H <- -.5*J*R2*(log.2.pi+1) -.5*sum(safe_log(mode2.H.var))

  if(S != 0) {
    exp.mode3.lambda <- m$mode3.lambda.shape * m$mode3.lambda.scale
    exp.log.mode3.lambda <- digamma(m$mode3.lambda.shape) + safe_log(m$mode3.lambda.scale)

    exp.p.mode3.lambda <- sum((m$m3.alpha-1)*exp.log.mode3.lambda -
                                (1/m$m3.beta)*exp.mode3.lambda -
                                m$m3.alpha * safe_log(m$m3.beta) -
                                safe_log(gamma(m$m3.alpha)))

    m3.A.var <- matrix(0, S, R3) # This could be stored earlier
    for(r3 in 1:R3) m3.A.var[,r3] <- diagonal(m$mode3.A.cov[,,r3])

    exp.p.mode3.A <- -.5 * S*R3 * log.2.pi + 1/2 * sum(exp.log.mode3.lambda -
      exp.mode3.lambda * (m$mode3.A.mean^2 + m3.A.var))
    
    # Sum used in for the expected value of the log probability for H
    cov.sum <- 0
    for(k in 1:K) for(r3 in 1:R3) {
      if(ncol(d$mode3.X) == (nrow(m$mode3.A.mean)-1)) {
        X_out_X <- outer(c(const=1, d$mode3.X[k,]), c(const=1, d$mode3.X[k,]))
      } else X_out_X <- outer(d$mode3.X[k,], d$mode3.X[k,])
      cov.sum <- cov.sum + sum(X_out_X *
        (outer(m$mode3.A.mean[,r3], m$mode3.A.mean[,r3]) + m$mode3.A.cov[,,r3]))
    }
 
    exp.p.mode3.H <- -.5*K*R3 * safe_log(2*pi*m$m3.sigma2) -
      1/(2*m$m3.sigma2) * sum(mode3.H.mean^2 + mode3.H.var) +
      1/(m$m3.sigma2) * sum(mode3.H.mean * safe_prod(d$mode3.X, m$mode3.A.mean)) -
      1/(2*m$m3.sigma2) * cov.sum

    exp.q.mode3.lambda <- sum(-m$mode3.lambda.shape - safe_log(m$mode3.lambda.scale) -
                                safe_log(gamma(m$mode3.lambda.shape)) -
                                (1-m$mode3.lambda.shape) * digamma(m$mode3.lambda.shape))

    m3.A.dets <- rep(0, R3)
    for(r3 in 1:R3) m3.A.dets[r3] <- determinant(m$mode3.A.cov[,,r3], logarithm=T)$modulus
    exp.q.mode3.A <- -sum(S/2*(log.2.pi+1) -.5*m3.A.dets)
  } else {
    exp.p.mode3.lambda <- 0
    exp.p.mode3.A <- 0
    exp.p.mode3.H <- -.5*K*R3 * safe_log(2*pi*m$m3.sigma2) - 1/(2*m$m3.sigma2) * sum(mode3.H.mean^2 + mode3.H.var)
    exp.q.mode3.lambda <- 0
    exp.q.mode3.A <- 0
  }
  exp.q.mode3.H <- -.5*K*R3*(log.2.pi+1) -.5*sum(safe_log(mode3.H.var))

  exp.core.lambda <- m$core.lambda.shape * m$core.lambda.scale
  exp.log.core.lambda <- digamma(m$core.lambda.shape) + safe_log(m$core.lambda.scale)
  
  exp.p.core.lambda <- sum((m$core.alpha-1) * (exp.log.core.lambda) - 1/m$core.beta * exp.core.lambda - 
                           m$core.alpha * log(m$core.beta) - log(gamma(m$core.alpha)))

  exp.q.core.lambda <- sum(-m$core.lambda.shape - log(m$core.lambda.scale) -
                           log(gamma(m$core.lambda.shape)) -
                           (1-m$core.lambda.shape) * digamma(m$core.lambda.shape))

  exp.p.core <- -.5 * core1 * core2 * core3 * log.2.pi +
    .5 * sum(exp.log.core.lambda - exp.core.lambda * (m$core.mean^2 + m$core.var))

  exp.q.core <- -.5*core1*core2*core3*(log.2.pi+1) -.5*sum(safe_log(m$core.var))

  sums <- array(0, dim=dim(d$resp))
    for(r1 in 1:core1) for(r2 in 1:core2) for(r3 in 1:core3) {
      subsums <- array(0, dim=c(I,J,K))
      subsums <- subsums + mult_3d(m$core.mean[r1,r2,r3,drop=F], m$mode1.H.mean[,r1,drop=F], m$mode2.H.mean[,r2,drop=F], m$mode3.H.mean[,r3,drop=F]) *
        mult_3d(m$core.mean[-r1,-r2,-r3,drop=F], m$mode1.H.mean[,-r1,drop=F], m$mode2.H.mean[,-r2,drop=F], m$mode3.H.mean[,-r3,drop=F])
      subsums <- subsums + sweep(mult_3d(m$core.mean[r1,r2,r3,drop=F],
                                   (m$mode1.H.mean[,r1,drop=F]^2 + m$mode1.H.var[,r1,drop=F]),
                                   m$mode2.H.mean[,r2,drop=F],
                                   m$mode3.H.mean[,r3,drop=F]),
        c(2,3), m$mode2.H.mean[,-r2,drop=F] %*% m$core.mean[r1,-r2,-r3] %*% t(m$mode3.H.mean[,-r3,drop=F]), FUN='*')
      subsums <- subsums + sweep(mult_3d(m$core.mean[r1,r2,r3,drop=F],
                                         m$mode1.H.mean[,r1,drop=F],
                                   (m$mode2.H.mean[,r2,drop=F]^2 + m$mode2.H.var[,r2,drop=F]),
                                   m$mode3.H.mean[,r3,drop=F]),
        c(1,3), m$mode1.H.mean[,-r1,drop=F] %*% m$core.mean[-r1,r2,-r3] %*% t(m$mode3.H.mean[,-r3,drop=F]), FUN='*')
      subsums <- subsums + sweep(mult_3d(m$core.mean[r1,r2,r3,drop=F],
                                         m$mode1.H.mean[,r1,drop=F],
                                         m$mode2.H.mean[,r2,drop=F],
                                   (m$mode3.H.mean[,r3,drop=F]^2 + m$mode3.H.var[,r3,drop=F])),
        c(1,2), m$mode1.H.mean[,-r1,drop=F] %*% m$core.mean[-r1,-r2,r3,drop=F][,,1] %*% t(m$mode2.H.mean[,-r2,drop=F]), FUN='*')
      subsums <- subsums + sweep(mult_3d(m$core.mean[r1,r2,r3,drop=F],
                                   (m$mode1.H.mean[,r1,drop=F]^2 + m$mode1.H.var[,r1,drop=F]),
                                   (m$mode2.H.mean[,r2,drop=F]^2 + m$mode2.H.var[,r2,drop=F]),
                                   m$mode3.H.mean[,r3,drop=F]),
        3, tcrossprod(m$core.mean[r1,r2,-r3,drop=F], m$mode3.H.mean[,-r3,drop=F]), FUN='*')
      subsums <- subsums + sweep(mult_3d(m$core.mean[r1,r2,r3,drop=F],
                                   (m$mode1.H.mean[,r1,drop=F]^2 + m$mode1.H.var[,r1,drop=F]),
                                   m$mode2.H.mean[,r2,drop=F],
                                   (m$mode3.H.mean[,r3,drop=F]^2 + m$mode3.H.var[,r3,drop=F])),
        2, tcrossprod(m$core.mean[r1,-r2,r3,drop=F], m$mode2.H.mean[,-r2,drop=F]), FUN='*')
      subsums <- subsums + sweep(mult_3d(m$core.mean[r1,r2,r3,drop=F],
                                         m$mode1.H.mean[,r1,drop=F],
                                   (m$mode2.H.mean[,r2,drop=F]^2 + m$mode2.H.var[,r2,drop=F]),
                                   (m$mode3.H.mean[,r3,drop=F]^2 + m$mode3.H.var[,r3,drop=F])),
        1, tcrossprod(m$core.mean[-r1,r2,r3,drop=F], m$mode1.H.mean[,-r1,drop=F]), FUN='*')
      subsums <- subsums + (m$core.mean[r1,r2,r3]^2 + m$core.var[r1,r2,r3]) *
        outer(outer(m$mode1.H.mean[,r1]^2+m$mode1.H.var[,r1],
                    m$mode2.H.mean[,r2]^2+m$mode2.H.var[,r2]),
              m$mode3.H.mean[,r3]^2+m$mode3.H.var[,r3])
      sums <- sums + subsums
  }

  exp.p.y <- -.5 * safe_log(2 * pi * m$sigma2) * sum(d$delta) -
    (1/(2*m$sigma2)) * sum(d$resp^2, na.rm=T) +
    (1/m$sigma2) * sum(d$resp * mult_3d(m$core.mean, m$mode1.H.mean,
                                        m$mode2.H.mean, m$mode3.H.mean), na.rm=T) -
    (1/(2*m$sigma2)) * sum(d$delta * sums)

  lower.bnd <- exp.p.mode1.lambda + exp.p.mode2.lambda + exp.p.mode3.lambda +
    exp.p.mode1.A + exp.p.mode2.A + exp.p.mode3.A +
    exp.p.mode1.H + exp.p.mode2.H + exp.p.mode3.H +
    exp.p.core.lambda + exp.p.core + exp.p.y -
    exp.q.mode1.lambda - exp.q.mode2.lambda - exp.q.mode3.lambda -
    exp.q.mode1.A - exp.q.mode2.A - exp.q.mode3.A -
    exp.q.mode1.H - exp.q.mode2.H - exp.q.mode3.H -
    exp.q.core.lambda - exp.q.core

  return(lower.bnd)
}
