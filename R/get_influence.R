#' Given a \code{model} object, rank the input predictors (and combinations thereof)
#' by thier influence on the output
#' 
#' If \code{method} is \code{'add'} then the baseline prediction is made using just the
#' constant coefficients (if used) and the mean squared error (MSE) is measured between 
#' the baseline and predictions made with each predictor added alone (univariate analysis).
#' 
#' If \code{method} is \code{'sub'} then the baseline is made using all predictors and
#' MSE measured for predictions made with each predictor removed.
#' 
#' If \code{interactions==TRUE} then MSE for predictions made with predictors for each mode
#' interacting are measured
#' 
#' @export
#' @param m \code{Tucker_model} or \code{CP_model} object
#' @param d \code{input_data} object
#' @param method string 'sub' or 'add' indicating whether to start with a full or empty 
#' feature vector and remove or add features to judge their influence.
#' @param interactions logical indicating whether to get influence for two-way interactions 
#' between predictors (def: sub)
 
get_influence <- function(m, d, method='sub', interactions=TRUE) {
  
  P <- nrow(m$mode1.A.mean)
  Q <- nrow(m$mode2.A.mean)
  S <- nrow(m$mode3.A.mean)

  m1.inf <- rep(NA, P);  m2.inf <- rep(NA, Q);  m3.inf <- rep(NA, S)  
  names(m1.inf) <- rownames(m$mode1.A.mean)
  names(m2.inf) <- rownames(m$mode2.A.mean)
  names(m3.inf) <- rownames(m$mode3.A.mean)

  if(method=='add') base.value <- 0
  if(method=='sub') base.value <- 1
  
  base.d <- input_data$new(mode1.X=matrix(base.value,1,P), 
                           mode2.X=matrix(base.value,1,Q), 
                           mode3.X=matrix(base.value,1,S))
  colnames(base.d$mode1.X) <- rownames(m$mode1.A.mean)
  colnames(base.d$mode2.X) <- rownames(m$mode2.A.mean)
  colnames(base.d$mode3.X) <- rownames(m$mode3.A.mean)
  
  baseline <- test(base.d, m)

  # Get univariate influence for mode 1
  if(P) for(p in 1:P) {
    new.d <- base.d$clone()
    new.d$mode1.X[p] <- -new.d$mode1.X[p] + 1 # 0 -> 1 and 1 -> 0
    m1.inf[p] <- mean((test(new.d, m) - baseline)^2, na.rm=T)
  }

  # Get univariate influence for mode 2
  if(Q) for(q in 1:Q) {
    new.d <- base.d$clone()
    new.d$mode2.X[q] <- -new.d$mode2.X[q] + 1 # 0 -> 1 and 1 -> 0
    m2.inf[q] <- mean((test(new.d, m) - baseline)^2, na.rm=T)
  }

  # Get univariate influence for mode 3
  if(S) for(s in 1:S) {
    new.d <- base.d$clone()
    new.d$mode3.X[s] <- -new.d$mode3.X[s] + 1 # 0 -> 1 and 1 -> 0
    m3.inf[s] <- mean((test(new.d, m) - baseline)^2, na.rm=T)
  }

  # Get multivariate influences  
  if(interactions) {
    if(length(m1.inf)) m1.names <- names(m1.inf) else m1.names <- 'none'
    if(length(m2.inf)) m2.names <- names(m2.inf) else m2.names <- 'none'
    if(length(m3.inf)) m3.names <- names(m3.inf) else m3.names <- 'none'
    inter.inf <- array(NA, dim=c(ifelse(P, P, 1),
                                 ifelse(Q, Q, 1),
                                 ifelse(S, S, 1)),
                       dimnames=list(m1.names, m2.names, m3.names))
    
    for(p in 1:ifelse(P,P,1)) for(q in 1:ifelse(Q,Q,1)) for(s in 1:ifelse(S,S,1)) {
      new.d <- base.d$clone()
      new.d$mode1.X[p] <- -new.d$mode1.X[p] + 1
      new.d$mode2.X[q] <- -new.d$mode2.X[q] + 1
      new.d$mode3.X[s] <- -new.d$mode3.X[s] + 1
      if(is.na(new.d$mode1.X[1])) new.d$mode1.X <- matrix(NA, 0, 0)
      if(is.na(new.d$mode2.X[1])) new.d$mode2.X <- matrix(NA, 0, 0)
      if(is.na(new.d$mode3.X[1])) new.d$mode3.X <- matrix(NA, 0, 0)
      inter.inf[p,q,s] <- mean((test(new.d, m) - baseline)^2, na.rm=T)
    }
  } else inter.inf <- NA
  
  return(list(m1.inf=m1.inf, m2.inf=m2.inf, m3.inf=m3.inf, 
              inter.inf=inter.inf))
}