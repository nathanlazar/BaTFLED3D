#' Remove predictors from A and X matrices if indicated by remove.lmt  
#' 
#' \code{test} 
#'
#' Laurem ispsum... This is a generic function: methods can be defined for it
#' directly or via the \code{\link{Summary}} group generic. For this to work
#' properly, the arguments \code{...} should be unnamed, and dispatch is on the
#' first argument.
#'
#' @export
#' @param An object 
#' @return Laurem ipsum

remove_preds <- function(m, d, params) {
  # Make all param variables available locally
  for(i in 1:length(params)) {
    assign(names(params)[i], params[i][[1]])
  }

  if(m$iter > remove.start) {
    m1.rem <- integer(0)
    # if(row.share) {
    #   m1.exp.A.var <- 1/(m$mode1.lambda.shape * m$mode1.lambda.scale)
    # } else {
    #   m1.exp.A.var <- rowMeans(1/(m$mode1.lambda.shape * m$mode1.lambda.scale))
    # }
    # if(m1.remove.lmt > 0) m1.rem <- which(m1.exp.A.var < m1.remove.lmt)
    m1.ab.means <- apply(abs(m$mode1.A.mean), 1, mean)
    if(m1.remove.lmt > 0) m1.rem <- which(m1.ab.means < m1.remove.lmt)
    if(remove.per > 0) {
      # m1.rem.cut <- sort(m1.exp.A.var)[length(m1.exp.A.var)*remove.per]
      m1.rem.cut <- sort(m1.ab.means)[length(m1.ab.means)*remove.per]
      m1.rem <- which(m1.exp.A.var < m1.rem.cut)
    }
    if(length(m1.rem) > 0) {
      print(sprintf('Removing %.0f predictors from mode 1', length(m1.rem)))
      if(row.share) {
        m$mode1.lambda.shape <- m$mode1.lambda.shape[-m1.rem]
        m$mode1.lambda.scale <- m$mode1.lambda.scale[-m1.rem]
        m$mode1.A.cov <- m$mode1.A.cov[-m1.rem,-m1.rem]
      } else {
        m$mode1.lambda.shape <- m$mode1.lambda.shape[-m1.rem,]
        m$mode1.lambda.scale <- m$mode1.lambda.scale[-m1.rem,]
        m$mode1.A.cov <- m$mode1.A.cov[-m1.rem,-m1.rem,]
      }
      m$mode1.A.mean <- m$mode1.A.mean[-m1.rem,]
    }
    # Mode 2
    m2.rem <- integer(0)
    # if(row.share) {
    #   m2.exp.A.var <- 1/(m$mode2.lambda.shape * m$mode2.lambda.scale)
    # } else {
    #   m2.exp.A.var <- rowMeans(1/(m$mode2.lambda.shape * m$mode2.lambda.scale))
    # }
    # if(m2.remove.lmt > 0) m2.rem <- which(m2.exp.A.var < m2.remove.lmt)
    m2.ab.means <- apply(abs(m$mode2.A.mean), 1, mean)
    if(m2.remove.lmt > 0) m2.rem <- which(m2.ab.means < m2.remove.lmt)
    if(remove.per > 0) {
      # m2.rem.cut <- sort(m2.exp.A.var)[length(m2.exp.A.var)*remove.per]
      m2.rem.cut <- sort(m2.ab.means)[length(m2.ab.means)*remove.per]
      m2.rem <- which(m2.exp.A.var < m2.rem.cut)
    }
    if(length(m2.rem) > 0) {
      print(sprintf('Removing %.0f predictors from mode 2', length(m2.rem)))
      if(row.share) {
        m$mode2.lambda.shape <- m$mode2.lambda.shape[-m2.rem]
        m$mode2.lambda.scale <- m$mode2.lambda.scale[-m2.rem]
        m$mode2.A.cov <- m$mode2.A.cov[-m2.rem,-m2.rem]
      } else {
        m$mode2.lambda.shape <- m$mode2.lambda.shape[-m2.rem,]
        m$mode2.lambda.scale <- m$mode2.lambda.scale[-m2.rem,]
        m$mode2.A.cov <- m$mode2.A.cov[-m2.rem,-m2.rem,]
      }
      m$mode2.A.mean <- m$mode2.A.mean[-m2.rem,]
    }
    # Mode 3
    m3.rem <- integer(0)
    # if(row.share) {
    #   m3.exp.A.var <- 1/(m$mode3.lambda.shape * m$mode3.lambda.scale)
    # } else {
    #   m3.exp.A.var <- rowMeans(1/(m$mode3.lambda.shape * m$mode3.lambda.scale))
    # }
    # if(m3.remove.lmt > 0) m3.rem <- which(m3.exp.A.var < m3.remove.lmt)
    m3.ab.means <- apply(abs(m$mode3.A.mean), 1, mean)
    if(m3.remove.lmt > 0) m3.rem <- which(m3.ab.means < m3.remove.lmt)
    if(remove.per > 0) {
      # m3.rem.cut <- sort(m3.exp.A.var)[length(m3.exp.A.var)*remove.per]
      m3.rem.cut <- sort(m3.ab.means)[length(m3.ab.means)*remove.per]
      m3.rem <- which(m3.exp.A.var < m3.rem.cut)
    }
    if(length(m3.rem) > 0) {
      print(sprintf('Removing %.0f predictors from mode 3', length(m3.rem)))
      if(row.share) {
        m$mode3.lambda.shape <- m$mode3.lambda.shape[-m3.rem]
        m$mode3.lambda.scale <- m$mode3.lambda.scale[-m3.rem]
        m$mode3.A.cov <- m$mode3.A.cov[-m3.rem,-m3.rem]
      } else {
        m$mode3.lambda.shape <- m$mode3.lambda.shape[-m3.rem,]
        m$mode3.lambda.scale <- m$mode3.lambda.scale[-m3.rem,]
        m$mode3.A.cov <- m$mode3.A.cov[-m3.rem,-m3.rem,]
      }
      m$mode3.A.mean <- m$mode3.A.mean[-m3.rem,]
    }
  }
  
  # Check dimensions of input data d and tensor factorization model m
  # and subset d if necessary
  if(length(rownames(m$mode1.A.mean))) {
    if(ncol(d$mode1.X) != sum(rownames(m$mode1.A.mean)!='const')) {
      d$mode1.X <- d$mode1.X[,dimnames(d$mode1.X)[[2]] %in% dimnames(m$mode1.A.mean)[[1]], drop=F]
      d$mode1.X <- d$mode1.X[,match(dimnames(m$mode1.A.mean)[[1]], dimnames(d$mode1.X)[[2]]), drop=F]
      if('const' %in% rownames(m$mode1.A.mean)) {
        m$m1Xm1X <- crossprod(cbind(1,d$mode1.X), cbind(1,d$mode1.X))
      } else m$m1Xm1X <- crossprod(d$mode1.X, d$mode1.X)
    }
  }
  if(length(rownames(m$mode2.A.mean))) {
    if(ncol(d$mode2.X) != sum(rownames(m$mode2.A.mean) != 'const')) {
      d$mode2.X <- d$mode2.X[,dimnames(d$mode2.X)[[2]] %in% dimnames(m$mode2.A.mean)[[1]], drop=F]
      d$mode2.X <- d$mode2.X[,match(dimnames(m$mode2.A.mean)[[1]], dimnames(d$mode2.X)[[2]]), drop=F]
      if('const' %in% rownames(m$mode2.A.mean)) {
        m$m2Xm2X <- crossprod(cbind(1,d$mode2.X), cbind(1,d$mode2.X))
      } else m$m2Xm2X <- crossprod(d$mode2.X, d$mode2.X)
    }
  }
  if(length(rownames(m$mode3.A.mean))) {
    if(ncol(d$mode3.X) != sum(rownames(m$mode3.A.mean) != 'const')) {
      d$mode3.X <- d$mode3.X[,dimnames(d$mode3.X)[[2]] %in% dimnames(m$mode3.A.mean)[[1]], drop=F]
      d$mode3.X <- d$mode3.X[,match(dimnames(m$mode3.A.mean)[[1]], dimnames(d$mode3.X)[[2]]), drop=F]
      if('const' %in% rownames(m$mode2.A.mean)) {
        m$m3Xm3X <- crossprod(cbind(1,d$mode3.X), cbind(1,d$mode3.X))
      } else m$m3Xm3X <- crossprod(d$mode3.X, d$mode3.X)
    }
  }
}