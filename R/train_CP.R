#' Train a CP model.
#' 
#' Model objects are updated in place to avoid memory issues. Nothing is returned.
#'
#' @importFrom stats cor
#' @importFrom graphics plot par
#'
#' @export
#' @param d an input data object created with \code{input_data}
#' @param m a \code{CP_model} object created with \code{mk_model} 
#' @param params List of parameters created with \code{get_model_params()}
#' @param new.iter numeric number of iterations to run (def: 1)
#' 
#' @examples
#' data.params <- get_data_params(c('decomp=CP'))
#' toy <- mk_toy(data.params)
#' train.data <- input_data$new(mode1.X=toy$mode1.X[,-1],
#'                              mode2.X=toy$mode2.X[,-1],
#'                              mode3.X=toy$mode3.X[,-1],
#'                              resp=toy$resp)
#' model.params <- get_model_params(c('decomp=CP'))
#' toy.model <- mk_model(train.data, model.params)
#' toy.model$rand_init(model.params)
#'
#' train(d=train.data, m=toy.model, new.iter=1, params=model.params)

# TODO: Remove parallel, cores?
#       What is "scale"? Just scales printed matrices.

train_CP <- function(d, m, new.iter=1, params) {
  # Make a copy of d so it is not changed
  d <- d$clone()
  
  # list of functions for updating (allows update order to be modified)
  up.funs <- list(update_mode1_CP, 
                  update_mode2_CP, 
                  update_mode3_CP)
  
# If remove.lmt != 0 rows of the A (projection) matrices are removed if the sum
  # of their absolute value drops below remove.lmt.
  # Multiplier is used to avoid underflow when normalizing.

  st <- F # Stopping condition

  if((params$m1.remove.lmt | params$m2.remove.lmt | params$m3.remove.lmt) & params$remove.per)
    print("Warning: only one of remove.lmt and remove.per should be set.")
  m$early.stop <- params$early.stop

  # Vector to store The lower bound of the log likelihood
  m$lower.bnd <- c(m$lower.bnd, rep(0, new.iter))

  # Vectors to store response measures (RMSE, explained variation & Pearson correlation)
  # for training data
  if(params$RMSE) {
    m$RMSE <- c(m$RMSE, rep(0, new.iter))       # RMSE from A matrices
    m$H.RMSE <- c(m$H.RMSE, rep(0, new.iter))   # RMse from H matrices
    m$exp.var <- c(m$exp.var, rep(0, new.iter)) # Explained variance for training data
    m$p.cor <- c(m$p.cor, rep(0, new.iter))     # Pearson correlation for training data
    m$s.cor <- c(m$s.cor, rep(0, new.iter))     # Spearman correlation for training data
  }
  
  # Vector to store the elapsed times for each iteration of the model
  if(params$time) m$times <- c(m$times, rep(0, new.iter))

  if(params$verbose) print('** Begining updates **')
  
  # The shape parameters for lambdas are all the same  and the update is only done once
  if(params$row.share) {
    m$mode1.lambda.shape[] <- rep(m$m1.alpha + params$R/2, nrow(m$mode1.A.mean))
    m$mode2.lambda.shape[] <- rep(m$m2.alpha + params$R/2, nrow(m$mode2.A.mean))
    m$mode3.lambda.shape[] <- rep(m$m3.alpha + params$R/2, nrow(m$mode3.A.mean))
  } else {
    m$mode1.lambda.shape[,] <- m$m1.alpha + .5
    m$mode2.lambda.shape[,] <- m$m2.alpha + .5
    m$mode3.lambda.shape[,] <- m$m3.alpha + .5
  }
  
  if(params$verbose) print('** Begining updates **')

  tot.iter <- m$iter + new.iter
  while(m$iter < tot.iter) {
    # Increase the total number of iterations run on the tens_fac object
    m$iter <- m$iter + 1

    print(paste('*************** Updating q distributions: round', m$iter, '******************'))

    # Start timer
    if(params$time) start.time <- proc.time()

    # Remove predictors if the mean of the absolute values of entries in the 
    # row of the projection matrix drops below remove.lmt. 
    # This speeds calculations and may give better predictions
    # remove_preds(m, d, params)

    # Choose which mode to update 
    # update.order <- sample(3,3)
    for(mode in params$update.order) up.funs[[mode]](m, d, params)

    # update_mode1_CP(m, d, params)
    # update_mode2_CP(m, d, params)
    # update_mode3_CP(m, d, params)

    # Update response predictions. If predictors are present responses
    # are the result of multiplying the X matrices through. Else
    # they are the result of just multiplying the X matrices
    if(ncol(d$mode1.X) == 0) {
      mode1.X.times.A <- m$mode1.H.mean
    } else mode1.X.times.A <- safe_prod(d$mode1.X, m$mode1.A.mean)

    if(ncol(d$mode2.X) == 0) {
      mode2.X.times.A <- m$mode2.H.mean
    } else mode2.X.times.A <- safe_prod(d$mode2.X, m$mode2.A.mean)

    if(ncol(d$mode3.X) == 0) {
      mode3.X.times.A <- m$mode3.H.mean
    } else mode3.X.times.A <- safe_prod(d$mode3.X, m$mode3.A.mean)

    core <- array(0, dim=c(params$R,params$R,params$R))
    for(r in 1:params$R) core[r,r,r] <- 1
    m$resp <- mult_3d(core, mode1.X.times.A, mode2.X.times.A, mode3.X.times.A)

    # Update the lower bound
    if(params$lower.bnd==T) {
      m$lower.bnd[m$iter] <- lower_bnd_CP(m, d)
      print(paste('Lower bound =',
        format(m$lower.bnd[m$iter], decimal.mark=".", big.mark=",", nsmall = 1)))
    }

    # Update the training RMSEs
    # Update the training performance measures
    if(params$RMSE) rmse(m, d, verbose=T)                               # Prints
    if(params$exp.var) m$exp.var[m$iter] <- exp_var(d$resp, m$resp, verbose=T) # Prints
    if(params$cor) {
      m$p.cor[m$iter] <- cor(d$resp, m$resp, use='complete.obs')
      print(sprintf("Pearson correlation: %.4f", m$p.cor[m$iter]))
      m$s.cor[m$iter] <- cor(d$resp, m$resp, use='complete.obs', method='spearman')
      print(sprintf("Spearman correlation: %.4f", m$s.cor[m$iter]))
    }

    # Plot the training predictions and projection matrices
    if(params$plot) {
      panels <- 1
      if(1 %in% params$show.mode) {
        panels <- panels + 1
        if(ncol(d$mode1.X)!=0) panels <- panels + 1
      }
      if(2 %in% params$show.mode) {
        panels <- panels + 1
        if(ncol(d$mode2.X)!=0) panels <- panels + 1
      }
      if(3 %in% params$show.mode) {
        panels <- panels + 1
        if(ncol(d$mode3.X)!=0) panels <- panels + 1
      }
      par(mfrow=c(panels%/%3, ceiling(panels/panels%/%3)))

      plot_preds(m$resp, d$resp, main=paste('Round: ', m$iter))

      show_mat(m, d, params$show.mode)
    }

    # Check early stopping criterion.
    # Stop if the lower bound increases by less than early.stop
    if(!is.na(m$early.stop) & params$lower.bnd & (m$iter > 20)) {
      if(m$lower.bnd[m$iter] < (m$lower.bnd[m$iter-1] + m$early.stop))  st <- T
    }

    if(params$time) m$times[m$iter] <- (proc.time() - start.time)[3]

    if(st) {
      m$lower.bnd <- m$lower.bnd[1:m$iter]
      m$RMSE <- m$RMSE[1:m$iter]
      m$H.RMSE <- m$H.RMSE[1:m$iter]
      m$times <- m$times[1:m$iter]
      m$exp.var <- m$exp.var[1:m$iter]
      m$cor <- m$cor[1:m$iter]
      return()
    }
  }
  return()
}