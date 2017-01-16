#' Get parameters to build a BaTFLED model
#' 
#' Read in vector of arguments, check their types and add them to a list \code{params} for 
#' model training. If a parameter isn't in the given list the default is used.
#' 
#' @export
#' @param args A character vector of arguments (character strings) of the form "<name>=<value>". 
#' Values will be converted to logical or numeric when possible.
#' Accepted <names> are below. Defaults in parenthesis:
#' \describe{
#' \item{decomp}{Either 'CP' or 'Tucker'. (Tucker)}
#' \item{row.share}{Logical. Should the variance be shared across rows of the projection matrices? 
#' This will cause predictors to be or excluded for the whole model, instead of just for particular 
#' latent factors. (F)}
#' \item{seed}{Numeric. Seed used for random initialization. (NA)}
#' \item{verbose}{Logical. Display more messages during training. (F)}
#' \item{parallel}{Logical. Perform operations in parallel when possible. (T)}
#' \item{cores}{Numeric. The number of parallel threads to use. (2)}
#' \item{lower.bnd}{Logical. Should the lower bound be calculated during training. Setting to FALSE
#' saves time (F)}
#' \item{RMSE}{Logical. Should the root mean squared error for the training data be calculated during 
#' training. (T)}
#' \item{cor}{Logical. Should the Pearson correlation for the training data be calculated during 
#' training. (T)}
#' \item{A1.intercept}{Logical. Add a constant column to the mode 1 predictors. (T)}
#' \item{A2.intercept}{Logical. Add a constant column to the mode 2 predictors. (T)}
#' \item{A3.intercept}{Logical. Add a constant column to the mode 3 predictors. (F)}
#' \item{H1.intercept}{Logical. Add a constant column to the mode 1 latent (H) matrix. (F)}
#' \item{H2.intercept}{Logical. Add a constant column to the mode 2 latent (H) matrix. (F)}
#' \item{H3.intercept}{Logical. Add a constant column to the mode 3 latent (H) matrix. (F)}
#' \item{R}{Numeric. Number of latent factors used in a CP model. (4)}
#' \item{R1}{Numeric. Number of latent factors used for mode 1 in a Tucker decomposition. (4)}
#' \item{R2}{Numeric. Number of latent factors used for mode 2 in a Tucker decomposition. (4)}
#' \item{R3}{Numeric. Number of latent factors used for mode 3 in a Tucker decomposition. (3)}
#' \item{core.updates}{Numeric. Number of core elements to update each round for stochastic training. (all)}
#' \item{m1.alpha}{Numeric. Prior for the 'shape' parameter of the gamma distribution on the 
#' precision values in the mode 1 projection (A) matrix. Set this to a small value (ex. 1e-10)
#' to encourage sparsity in mode 1 predictors. (1e-10)}
#' \item{m2.alpha}{Numeric. Same as above for mode 2. (1e-10)}
#' \item{m3.alpha}{Numeric. Same as above for mode 3. (1)}
#' \item{m1.beta}{Numeric. Prior for the 'scale' parameter of the gamma distribution on the 
#' precision values in the mode 1 projection (A) matrix. Set this to a large value (ex. 1e10)
#' to encourage sparsity in mode 1 predictors. Note this should stay balanced with m1.alpha 
#' so thir product is 1. (1e10)}
#' \item{m2.beta}{Numeric. Same as above for mode 2. (1e10)}
#' \item{m3.beta}{Numeric. Same as above for mode 3. (1)}
#' \item{A.samp.sd}{Numeric. Standard deviation used when initializing values in the projection
#' (A) matrices. (1)}
#' \item{H.samp.sd}{Numeric. Standard deviation used when initializing values in the latent
#' (H) matrices. (1)}
#' \item{R.samp.sd}{Numeric. Standard deviation used when initializing values in the core
#' tensor for Tucker models. (1)}
#' \item{A.var}{Numeric. Initial variance for projection (A) matrices. (1)}
#' \item{H.var}{Numeric. Initial variance for latent (H) matrices. (1)}
#' \item{R.var}{Numeric. Initial variance for the core tensor in Tucker models. (1)}
#' \item{random.H}{Logical. Should the latent matrices be initialized randomly or be the result
#' of multiplying the input data by the projection matrices. (T)}
#' \item{core.0D.alpha}{Numeric. Prior for the 'scale' parameter of the gamma distribution on the 
#' precision value in the element of the core tensor corresponding to the intercept for all
#' three modes (core.mean[1,1,1]). Only used for Tucker models when all H intercepts are true. 
#' Set this to a small value (ex. 1e-10) to encourage sparsity in core predictor. (1e-10)}
#' \item{core.1D.alpha}{Numeric. As above for values corresponding to the intercepts for 
#' two modes (core.mean[1,1,], core.mean[1,,1] and core.mean[,1,1]). (1e-10)}
#' \item{core.2D.alpha}{Numeric. As above for values corresponding to the intercepts for 
#' one mode (core.mean[1,,], core.mean[,1,] and core.mean[,,1]). (1e-10)}
#' \item{core.3D.alpha}{Numeric. As above for values not corresponding to intercepts. (1e-10)}
#' \item{core.0D.beta}{Numeric. As above but a prior for the 'scale' parameter. (1e10)}
#' \item{core.1D.beta}{Numeric. As above but a prior for the 'scale' parameter. (1e10)}
#' \item{core.2D.beta}{Numeric. As above but a prior for the 'scale' parameter. (1e10)}
#' \item{core.3D.beta}{Numeric. As above but a prior for the 'scale' parameter. (1e10)}
#' \item{m1.sigma2}{Numeric. Variance for the mode 1 latent (H) matrix. Set small to link the
#' values in the latent matrices to the product of the input and projection matrices. If there 
#' is no input data, set to one or larger. (0.01)}
#' \item{m2.sigma2}{Numeric. As above for mode 2. (0.01)}
#' \item{m3.sigma2}{Numeric. As above for mode 3. (1)}
#' \item{Sigma2}{Numeric. Variance for the response tensor. (1)}
#' \item{remove.start}{Numeric. The iteration to begin removing predictors if any of 
#' \code{m1.remove.lmt}, \code{m2.remove.lmt}, \code{m3.remove.lmt} or \code{remove.per}
#' are set. (Inf)}
#' \item{remove.per}{Numeric. Percentage of predictors to remove with the lowest mean of 
#' squared values across rows of the projection matrix. (0)}
#' \item{m1.remove.lmt}{Numeric. Remove a mode 1 predictor if the mean squared value of 
#' its row in the projection matrix drop below this value. (0)}
#' \item{m2.remove.lmt}{As above for mode 2. (0)}
#' \item{m3.remove.lmt}{As above for mode 3. (0)}
#' \item{early.stop}{Numeric. Stop training if the lower bound value changes by less than 
#' this value. (0)}
#' \item{plot}{Logical. Show plots while training}
#' \item{show.mode}{Numeric vector. Display images of the projection and latent matrices
#' for these modes while training. (c(1,2,3))}
#' \item{update.order}{Numeric vector. Update the modes in this order (c(3,2,1))}
#' }
#' @return list of parameters used by \code{train} function. Values in \code{args} that
#' are not model parameters will be excluded and a warning displayed.
#' @examples
#' args <- c('decomp=Tucker', 'row.share=F',
#'           'A1.intercept=T', 'A2.intercept=T', 'A3.intercept=F',
#'           'H1.intercept=T', 'H2.intercept=T', 'H3.intercept=T',
#'           'plot=T', 'verbose=F','R1=4', 'R2=4', 'R3=2',
#'           'm1.alpha=1e-10', 'm2.alpha=1e-10', 'm3.alpha=1',
#'           'm1.beta=1e10', 'm2.beta=1e10', 'm3.beta=1',
#'           'core.3D.alpha=1e-10', 'core.3D.beta=1e10',
#'           'parallel=T', 'cores=5', 'lower.bnd=T',
#'           'update.order=c(3,2,1)', 'show.mode=c(1,2,3)',
#'           'wrong=1')
#' model.params <- get_model_params(args)
#' @seealso \code{\link{CP_model}} \code{\link{Tucker_model}}

get_model_params <- function(args) {

  get_name <- function(arg) {
    strsplit(arg, split = '=')[[1]][1]
  }
  
  convert <- function(args, name) {
    value <- gsub(paste0('^', name, '='), '', args[grepl(paste0('^', name, '='), args)])
    if(length(value)) {
      if(value=='NA') return(NA)
      if(is.na(as.logical(value))) {
        if(is.na(suppressWarnings(as.numeric(value)))) {
          if(grepl('c\\(', value)) {
            v <- strsplit(sub(')', '', sub('c\\(', '', value)), split=',')[[1]]
            return(as.numeric(v))
          } else return(value)
        } else return(as.numeric(value))
      } else return(as.logical(value))
    } else return(NA)
  }
  
  # Set defaults
  params <- list(decomp='Tucker', row.share=F, 
                 seed=NA, verbose=F, parallel=T, cores=2, 
                 lower.bnd=T, RMSE=T, exp.var=T, cor=T, time=T,
                 A1.intercept=T, A2.intercept=T, A3.intercept=F,
                 H1.intercept=F, H2.intercept=F, H3.intercept=F,
                 R=4, R1=4, R2=4, R3=3,
                 m1.alpha=1, m1.beta=1, 
                 m2.alpha=1, m2.beta=1,
                 m3.alpha=1, m3.beta=1,
                 A.samp.sd=1, H.samp.sd=1, R.samp.sd=1,
                 A.var=1, H.var=1, R.var=1,
                 random.H=T,
                 core.0D.alpha=1, core.0D.beta=1,
                 core.1D.alpha=1, core.1D.beta=1,
                 core.2D.alpha=1, core.2D.beta=1,
                 core.3D.alpha=1, core.3D.beta=1,
                 core.updates=Inf,
<<<<<<< HEAD
                 m1.sigma2=.01, m2.sigma2=.01, m3.sigma2=0.01,
=======
                 m1.sigma2=.01, m2.sigma2=.01, m3.sigma2=1,
>>>>>>> b452907b1fbd5d5edc1cccb97e72c2c1d3cebb36
                 sigma2=1, 
                 remove.start=Inf, 
                 remove.per=0,
                 m1.remove.lmt=0,
                 m2.remove.lmt=0,
                 m3.remove.lmt=0,
                 plot=T, early.stop=NA,
                 update.order=c(3,2,1),
                 show.mode=c(1,2,3))

  accepted <- names(params)
  
  # Over ride defaults if arguments are supplied
  for(arg in args) {
    name <- get_name(arg)
    value <- convert(args,name)
    for(n in 1:length(value)) params[name][[1]][n] <- value[n]
    params[name][[1]] <- params[name][[1]][1:length(value)]
  }

  # Drop parameters if not needed
  if(params$decomp=='CP') {
    params[c('R1','R2','R3')] <- NULL
    params[c('H1.intercept','H2.intercept','H3.intercept')] <- NULL
    params[grepl('core', names(params))] <- NULL
    params[c('R.samp.sd', 'R.var')] <- NULL
  }
  
  if(params$decomp=='Tucker') {
    params['R'] <- NULL
    if(!(params$H1.intercept & params$H2.intercept & params$H3.intercept)) 
      {params$core.0D.alpha <- NULL; params$core.0D.beta <- NULL}
    if(sum(params$H1.intercept, params$H2.intercept, params$H3.intercept)<2)
      {params$core.1D.alpha <- NULL; params$core.1D.beta <- NULL}
    if(sum(params$H1.intercept, params$H2.intercept, params$H3.intercept)==0)
      {params$core.2D.alpha <- NULL; params$core.2D.beta <- NULL}
  }
  if(!params$random.H) params['H.samp.sd'] <- NULL
  
  wrong <- params[!(names(params) %in% accepted)]
  params <- params[names(params) %in% accepted]
  
  if(length(wrong) > 0) {
    warning(paste('The following parameter is not accepted:', 
                  names(wrong), collapse = '\n'), call.=F)
  }
  return(params)
}