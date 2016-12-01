#' Get parameters for building a model with known relationships
#' 
#' Read in vector of arguments, check their types and add them to a list \code{params} for 
#' building a model of input and response data with known relationships. If a parameter 
#' isn't in the given list the default is used.
#' 
#' @export
#' @param args A character vector of arguments (character strings) of the form "<name>=<value>". 
#' Values will be converted to logical or numeric when necessary.
#' Accepted <names> are below. Defaults in parenthesis:
#' \describe{
#' \item{decomp}{Either 'CP' or 'Tucker'. (Tucker)}
#' \item{row.share}{Logical. Should the variance be shared across rows of the projection matrices? 
#' This will cause predictors to be or excluded for the whole model, instead of just for particular 
#' latent factors. (T)}
#' \item{seed}{Numeric. Seed used for random initialization. (NA)}
#' \item{scale}{Logical. Should the input data columns should be scaled to have mean 0 and 
#' standard deviation 1. (TRUE)}
#' \item{m1.rows}{Numeric. Number of rows (samples) for mode 1. (20)}
#' \item{m2.rows}{Numeric. Number of rows (samples) for mode 2. (25)}
#' \item{m3.rows}{Numeric. Number of rows (samples) for mode 3. (10)}
#' \item{m1.cols}{Numeric. Number of columns (predictors) for mode 1. (100)}
#' \item{m2.cols}{Numeric. Number of columns (predictors) for mode 2. (150)}
#' \item{m3.cols}{Numeric. Number of columns (predictors) for mode 3. (0)}
#' \item{R}{Numeric. If \code{decomp=='CP'} the dimension of the latent space for all modes. (4)}
#' \item{R1}{Numeric. If \code{decomp=='Tucker'} the dimension of the core (latent space) for
#' mode 1. (4)}
#' \item{R2}{Numeric. If \code{decomp=='Tucker'} the dimension of the core (latent space) for
#' mode 2. (4)}
#' \item{R3}{Numeric. If \code{decomp=='Tucker'} the dimension of the core (latent space) for
#' mode 3. (3)}
#' \item{A1.intercept}{Logical. Should a column of 1s be added to the input data for mode 1. (TRUE)}
#' \item{A2.intercept}{Logical. Should a column of 1s be added to the input data for mode 2. (TRUE)}
#' \item{A3.intercept}{Logical. Should a column of 1s be added to the input data for mode 3. (TRUE)}
#' \item{H1.intercept}{Logical. Should a column of 1s be added to the latent (H) matrix for mode 1. (TRUE)}
#' \item{H2.intercept}{Logical. Should a column of 1s be added to the latent (H) matrix for mode 2. (TRUE)}
#' \item{H3.intercept}{Logical. Should a column of 1s be added to the latent (H) matrix for mode 3. (TRUE)}
#' \item{m1.true}{Numeric. Number of predictors for mode 1 (not counting the constant) 
#' contributing to the response. (15)}
#' \item{m2.true}{Numeric. Number of predictors for mode 2 (not counting the constant) 
#' contributing to the response. (20)}
#' \item{m3.true}{Numeric. Number of predictors for mode 3 (not counting the constant) 
#' contributing to the response. (0)}
#' \item{A1.const.prob}{Numeric. Probability (0-1) of the constant term for mode 1 contributing 
#' to the response for mode 1. (1)}
#' \item{A2.const.prob}{Numeric. Probability (0-1) of the constant term for mode 2 contributing 
#' to the response. (1)}
#' \item{A3.const.prob}{Numeric. Probability (0-1) of the constant term for mode 3 contributing 
#' to the response. (1)}
#' \item{A.samp.sd}{Numeric. Standard deviation for sampling values for the projection (A) matrices. (1)}
#' \item{H.samp.sd}{Numeric. Standard deviation for sampling values for the latent (H) matrices. (1)}
#' \item{R.samp.sd}{Numeric. Standard deviation for sampling values for the core tensor. (1)}
#' \item{true.0D}{Numeric. 0 or 1, should a global intercept (0 dimensional intercept) be added
#' to all responses? Only possible if \code{H1.intercept==H2.intercept==H3.intercept==TRUE}.
#' \code{core.spar} is used if equal to \code{NA}. (NA)}
#' \item{true.1D.m[1-3]}{Numeric. Number of interactions of 1 dimension in the core tensor (non-zero elements
#' on the edges of the core tensor if \code{H#.intercept==TRUE}). \code{core.spar} is used if 
#' equal to \code{NA}. (NA)}
#' \item{true.2D.m[1-3]m[1-3]}{Numeric. Number of interactions of 2 dimensions in the core tensor (non-zero elements
#' of the faces of the core tensor if \code{H#.intercept==TRUE}). \code{core.spar} is used if 
#' equal to \code{NA}. (NA)}
#' \item{true.3D}{Numeric. Number of interactions of 3 dimensions in the core tensor (non-zero elements
#' internal to the core tensor). \code{core.spar} is used if equal to \code{NA}. (NA)}
#' \item{core.spar}{Numeric. Fraction of core elements that are non-zero. (1)}
#' \item{noise.sd}{Numeric. Relative standard deviation of noise added to response tensor. (0.1)}
#' }
#' @return list of parameters used by \code{mk_toy} function. Values in \code{args} that
#' are not accepted parameters will be excluded and a warning displayed.
#' @examples
#' args <- c('decomp=Tucker', 'row.share=F',
#'           'A1.intercept=T', 'A2.intercept=T', 'A3.intercept=F',
#'           'H1.intercept=T', 'H2.intercept=T', 'H3.intercept=T',
#'           'R1=4', 'R2=4', 'R3=2')
#' data.params <- get_data_params(args)
#' @seealso \code{\link{mk_toy}}

get_data_params <- function(args) {
  
  get_name <- function(arg) {
    strsplit(arg, split = '=')[[1]][1]
  }

  convert <- function(args, name) {
    value <- gsub(name, '', args[grepl(name, args)])
    value <- gsub('=', '', value)
    if(is.na(as.logical(value))) {
      if(is.na(suppressWarnings(as.numeric(value)))) {
        return(value)
      } else return(as.numeric(value))
    } else return(as.logical(value))
  }

  # Set defaults
  params <- list(decomp='Tucker', row.share=T, seed=NA, scale=T,
                 m1.rows=20, m2.rows=25, m3.rows=10,
                 m1.cols=100, m2.cols=150, m3.cols=0,
                 R=4, R1=4, R2=4, R3=3,
                 A1.const.prob=1, A2.const.prob=1, A3.const.prob=1,
                 m1.true=15, m2.true=20, m3.true=0,
                 A1.intercept=T, A2.intercept=T, A3.intercept=F,
                 H1.intercept=T, H2.intercept=T, H3.intercept=T,
                 A.samp.sd=1, H.samp.sd=1, R.samp.sd=1,
                 true.0D=NA, true.1D.m1=NA, true.1D.m2=NA, true.1D.m3=NA,
                 true.2D.m1m2=NA, true.2D.m1m3=NA, true.2D.m2m3=NA,
                 true.3D=NA, core.spar=1, noise.sd=0.1)
                 
  accepted <- names(params)

  # Over ride defaults if arguments are supplied
  for(arg in args) {
    name <- get_name(arg)
    params[name] <- convert(args, name)
  }

  # If any of true.<>D are 'NA' then use core.spar to fill them in
  if(is.na(params$true.0D)) 
    params$true.0D <- rbinom(1, 1, params$core.spar)
  if(is.na(params$true.1D.m1)) 
    params$true.1D.m1 <- sum(rbinom(params$R1, 1, params$core.spar))
  if(is.na(params$true.1D.m2)) 
    params$true.1D.m2 <- sum(rbinom(params$R2, 1, params$core.spar))
  if(is.na(params$true.1D.m3)) 
    params$true.1D.m3 <- sum(rbinom(params$R3, 1, params$core.spar))
  if(is.na(params$true.2D.m1m2)) 
    params$true.2D.m1m2 <- sum(rbinom(params$R1*params$R2, 1, params$core.spar))
  if(is.na(params$true.2D.m1m3)) 
    params$true.2D.m1m3 <- sum(rbinom(params$R1*params$R3, 1, params$core.spar))
  if(is.na(params$true.2D.m2m3)) 
    params$true.2D.m2m3 <- sum(rbinom(params$R2*params$R3, 1, params$core.spar))
  if(is.na(params$true.3D)) 
    params$true.3D <- sum(rbinom(params$R1*params$R2*params$R3, 1, params$core.spar))
  
  # Drop parameters if not needed
  if(params$decomp=='CP') {
    params[c('R1','R2','R3')] <- NULL
    params[c('H1.intercept','H2.intercept','H3.intercept')] <- NULL
    params[c('H1.const.prob','H2.const.prob','H3.const.prob')] <- NULL
    params[c('core.spar', 'true.0D', 'true.3D')] <- NULL
    params['R.samp.sd'] <- NULL
    params[grepl('true.1D', names(params))] <- NULL
    params[grepl('true.2D', names(params))] <- NULL
  }
  
  if(params$decomp=='Tucker') {
    params$R <- NULL
    if(!params$A1.intercept) params$A1.const.prob <- NULL
    if(!params$A2.intercept) params$A2.const.prob <- NULL
    if(!params$A3.intercept) params$A3.const.prob <- NULL
    if(!(params$H1.intercept & params$H2.intercept & params$H3.intercept)) 
      params$true.0D <- NULL
    if(!(params$H1.intercept & params$H2.intercept))
      params$true.1D.m3 <- NULL
    if(!(params$H1.intercept & params$H3.intercept))
      params$true.1D.m2 <- NULL
    if(!(params$H2.intercept & params$H3.intercept))
      params$true.1D.m1 <- NULL
    if(!params$H1.intercept) params$true.2D.m2m3 <- NULL
    if(!params$H2.intercept) params$true.2D.m1m3 <- NULL
    if(!params$H3.intercept) params$true.2D.m1m2 <- NULL
  }
  
  wrong <- params[!(names(params) %in% accepted)]
  if(length(wrong) > 0) {
    print('The following parameters are not accepted')
    print(wrong)
  }
  params <- params[names(params) %in% accepted]

  return(params)
}