#' Update values in the core tensor for a Tucker model.
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
#' update_core_Tucker(m=toy.model, d=train.data, params=model.params)

update_core_Tucker <- function(m, d, params) {
  # Make all param variables available locally
  for(i in 1:length(params)) {
    assign(names(params)[i], params[i][[1]])
  }
  
  # Number of core updates to perform per iteration.
  if(core.updates == Inf) core.updates <- prod(dim(m$core.mean))
  
  if(verbose) print("Updating core tensor")
  
  I <- dim(d$resp)[1]; J <- dim(d$resp)[2]; K <- dim(d$resp)[3]
  R1 <- dim(m$core.mean)[1]; R2 <- dim(m$core.mean)[2]; R3 <- dim(m$core.mean)[3]
  # Note: these R's include constant slices
  
  # Update scale parameter for core lambda
  m$core.lambda.scale <- 1/(.5*((m$core.mean)^2 + m$core.var) + 1/m$core.beta)
  
  # Store some values to avoid repetition
  exp.m1.H.sq <- m$mode1.H.mean^2 + m$mode1.H.var
  exp.m2.H.sq <- m$mode2.H.mean^2 + m$mode2.H.var
  exp.m3.H.sq <- m$mode3.H.mean^2 + m$mode3.H.var
  
  # Update coefficients of constant terms first
  eg.mix <- expand.grid(1:R1, 1:R2, 1:R3)
  # eg.mix <- eg.mix[order(apply(eg.mix, 1, function(x) sum(x!=1))), ]
  eg.mix <- eg.mix[sample(nrow(eg.mix), core.updates),]
  
  # TODO: speed up this loop (it can be parallelized)
  # May need to restrict it to the subset of core elements that are to be updated each time randomly
  for(n in 1:core.updates) {
    r1 <- eg.mix[n,1]; r2 <- eg.mix[n,2]; r3 <- eg.mix[n,3]
    m$core.var[r1,r2,r3] <- 1/((1/m$sigma2) *
       sum(d$delta * (exp.m1.H.sq[,r1] %o% exp.m2.H.sq[,r2] %o% exp.m3.H.sq[,r3])) +
       m$core.lambda.shape[r1,r2,r3] * m$core.lambda.scale[r1,r2,r3])
  }
  
  # Update core mean
  # for(r3 in 1:R3) for(r2 in 1:R2) for(r1 in 1:R1) {
  # for(n in 1:nrow(eg.mix)) {
  for(n in 1:core.updates) {
    r1 <- eg.mix[n,1]; r2 <- eg.mix[n,2]; r3 <- eg.mix[n,3]
    sum0 <- m$mode1.H.mean[,r1] %o% m$mode2.H.mean[,r2] %o% m$mode3.H.mean[,r3]

    # Sum when all r's are not equal
    big_sum <- (m$mode1.H.mean[,r1] %o% m$mode2.H.mean[,r2] %o% m$mode3.H.mean[,r3]) *
<<<<<<< HEAD
      mult_3d(m$core.mean[-r1,-r2,-r3,drop=F], m$mode1.H.mean[,-r1,drop=F],
              m$mode2.H.mean[,-r2,drop=F], m$mode3.H.mean[,-r3,drop=F])

    # Sums when two r's are not equal
    big_sum <- big_sum + (exp.m1.H.sq[,r1] %o% ((m$mode2.H.mean[,r2] %o% m$mode3.H.mean[,r3]) *
                            (m$mode2.H.mean[,-r2,drop=F] %*% m$core.mean[r1,-r2,-r3] %*%
                             t(m$mode3.H.mean[,-r3,drop=F]))))
                               # rTensor::ttl(rTensor::as.tensor(m$core.mean[r1,-r2,-r3,drop=F]),
                               #     list(m$mode2.H.mean[,-r2,drop=F], m$mode3.H.mean[,-r3,drop=F]),
                               #     c(2,3))@data[1,,]))
    # Rotate this one with aperm
    big_sum <- big_sum + aperm((exp.m2.H.sq[,r2] %o% ((m$mode1.H.mean[,r1] %o% m$mode3.H.mean[,r3]) *
                               (m$mode1.H.mean[,-r1,drop=F] %*% m$core.mean[-r1,r2,-r3] %*%
                                 t(m$mode3.H.mean[,-r3,drop=F])))), c(2,1,3))
                         # (rTensor::ttl(rTensor::as.tensor(m$core.mean[-r1,r2,-r3,drop=F]),
                         #      list(m$mode1.H.mean[,-r1,drop=F], m$mode3.H.mean[,-r3,drop=F]),
                         #      c(1,3))@data[,1,]))), c(2,1,3))
    big_sum <- big_sum + ((m$mode1.H.mean[,r1] %o% m$mode2.H.mean[,r2]) *
                            (m$mode1.H.mean[,-r1,drop=F] %*% m$core.mean[-r1,-r2,r3] %*%
                            t(m$mode2.H.mean[,-r2,drop=F]))) %o% exp.m3.H.sq[,r3]
                            # rTensor::ttl(rTensor::as.tensor(m$core.mean[-r1,-r2,r3,drop=F]),
                            #     list(m$mode1.H.mean[,-r1,drop=F], m$mode2.H.mean[,-r2,drop=F]),
                            #     c(1,2))@data[,,1]) %o% exp.m3.H.sq[,r3]
=======
      rTensor::ttl(rTensor::as.tensor(m$core.mean[-r1,-r2,-r3,drop=F]),
          list(m$mode1.H.mean[,-r1,drop=F],
               m$mode2.H.mean[,-r2,drop=F],
               m$mode3.H.mean[,-r3,drop=F]), c(1,2,3))@data
    
    # Sums when two r's are not equal
    big_sum <- big_sum + (exp.m1.H.sq[,r1] %o%
                            ((m$mode2.H.mean[,r2] %o% m$mode3.H.mean[,r3]) *
                               rTensor::ttl(rTensor::as.tensor(m$core.mean[r1,-r2,-r3,drop=F]),
                                   list(m$mode2.H.mean[,-r2,drop=F], m$mode3.H.mean[,-r3,drop=F]),
                                   c(2,3))@data[1,,]))
    # Rotate this one with aperm
    big_sum <- big_sum + aperm((exp.m2.H.sq[,r2] %o%
                                  ((m$mode1.H.mean[,r1] %o% m$mode3.H.mean[,r3]) *
                                     (rTensor::ttl(rTensor::as.tensor(m$core.mean[-r1,r2,-r3,drop=F]),
                                          list(m$mode1.H.mean[,-r1,drop=F], m$mode3.H.mean[,-r3,drop=F]),
                                          c(1,3))@data[,1,]))), c(2,1,3))
    big_sum <- big_sum + ((m$mode1.H.mean[,r1] %o% m$mode2.H.mean[,r2]) *
                            rTensor::ttl(rTensor::as.tensor(m$core.mean[-r1,-r2,r3,drop=F]),
                                list(m$mode1.H.mean[,-r1,drop=F], m$mode2.H.mean[,-r2,drop=F]),
                                c(1,2))@data[,,1]) %o% exp.m3.H.sq[,r3]
>>>>>>> b452907b1fbd5d5edc1cccb97e72c2c1d3cebb36
    
    # Sums when 1 r is not equal
    big_sum <- big_sum + ((exp.m1.H.sq[,r1] %o% exp.m2.H.sq[,r2]) %o%
                            (m$mode3.H.mean[,r3] *
                               drop(tcrossprod(matrix(m$core.mean[r1,r2,-r3],1,R3-1), m$mode3.H.mean[,-r3,drop=F]))))
    big_sum <- big_sum + exp.m1.H.sq[,r1] %o%
      (m$mode2.H.mean[,r2] *
         drop(tcrossprod(matrix(m$core.mean[r1,-r2,r3],1,R2-1), m$mode2.H.mean[,-r2,drop=F]))) %o%
      exp.m3.H.sq[,r3]
    big_sum <- big_sum + (m$mode1.H.mean[,r1] *
                            drop(tcrossprod(matrix(m$core.mean[-r1,r2,r3],1, R1-1), m$mode1.H.mean[,-r1,drop=F]))) %o%
      exp.m2.H.sq[,r2] %o% exp.m3.H.sq[,r3]
    
    m$core.mean[r1,r2,r3] <- m$core.var[r1,r2,r3] * (1/m$sigma2) * sum(d$resp * sum0 - big_sum, na.rm=T)
  }
}
