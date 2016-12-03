#' Train a Tucker model using BaTFLED algorthm 
#' 
#' Model objects are updated in place to avoid memory issues. Nothing is returned.
#'
#' @export
#' @param d an input data object created with \code{input_data}
#' @param m a \code{Tucker_model} object created with \code{mk_model} 
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
#' train(d=train.data, m=toy.model, new.iter=1, params=model.params)

train_Tucker <- function(d, m, new.iter=1, params) {
  # Make all param variables available locally
  for(i in 1:length(params)) {
    assign(names(params)[i], params[i][[1]])
  }
  
  # Make a copy of d so it is not changed
  d <- d$clone()
  
  # list of functions for updating (allows update order to be modified)
  up.funs <- list(update_mode1 = update_mode1_Tucker, 
                  update_mode2 = update_mode2_Tucker, 
                  update_mode3 = update_mode3_Tucker)
  
  # If remove.lmt != 0 rows of the A (projection) matrices and corresponding parameters
  # are removed if the expected value of the row variance (1/lambda) falls below remove.lmt.
  # If remove.start is set, then this doesn't begin until that iteration.
  # If remove.per is set than a percentage of predictors are removed instead.
  # Multiplier is used to avoid underflow when normalizing.

  st <- F # Stopping condition

  if((m1.remove.lmt | m2.remove.lmt | m3.remove.lmt) & remove.per)
      print("Warning: only one of remove.lmt and remove.per should be set.")

  m$early.stop <- early.stop

  # Vector to store The lower bound of the log likelihood
  m$lower.bnd <- c(m$lower.bnd, rep(0, new.iter))

  # Vectors to store response measures (RMSE, explained variation & Pearson correlation)
  # for training data
  if(RMSE) {
    m$RMSE <- c(m$RMSE, rep(0, new.iter))       # RMSE from A matrices
    m$H.RMSE <- c(m$H.RMSE, rep(0, new.iter))   # RMse from H matrices
    m$exp.var <- c(m$exp.var, rep(0, new.iter)) # Explained variance for training data
    m$p.cor <- c(m$p.cor, rep(0, new.iter))     # Pearson correlation for training data
    m$s.cor <- c(m$s.cor, rep(0, new.iter))     # Spearman correlation for training data
  }

  # Vector to store the elapsed times for each iteration of the model
  if(time) m$times <- c(m$times, rep(0, new.iter))

  if(verbose) print('** Begining updates **')

  # The shape parameters for lambdas are all the same  and the update is only done once
  if(row.share) {
    m$mode1.lambda.shape[] <- rep(m$m1.alpha + R1/2, nrow(m$mode1.A.mean))
    m$mode2.lambda.shape[] <- rep(m$m2.alpha + R2/2, nrow(m$mode2.A.mean))
    m$mode3.lambda.shape[] <- rep(m$m3.alpha + R3/2, nrow(m$mode3.A.mean))
  } else {
    m$mode1.lambda.shape[,] <- m$m1.alpha + .5
    m$mode2.lambda.shape[,] <- m$m2.alpha + .5
    m$mode3.lambda.shape[,] <- m$m3.alpha + .5
  }
  
  # The shape parameters for the core lambda are all the same and
  # the update is only done once
  m$core.lambda.shape <- m$core.alpha + .5

  tot.iter <- m$iter + new.iter
  while(m$iter < tot.iter) {
    # Increase the total number of iterations run on the tens_fac object
    m$iter <- m$iter + 1

    print(paste('*************** Updating q distributions: round', m$iter, '******************'))

    # Start timer
    if(time) start.time <- proc.time()

    # Remove predictors if the mean of the absolute values of entries in the 
    # row of the projection matrix drops below remove.lmt. 
    # This speeds calculations and may give better predictions
    remove_preds(m, d, params)

    # Choose which mode to update 
    #up.mode <- sample(3,3)
    update_core_Tucker(m, d, params)
    for(mode in update.order) up.funs[[mode]](m, d, params)

    # update_mode3_Tucker(m, d, params)
    # update_mode2_Tucker(m, d, params)
    # update_mode1_Tucker(m, d, params)
    # update_core_Tucker(m, d, params)

    # Update response predictions. If predictors are present responses
    # are the result of multiplying the X matrices through. Else
    # they are the result of just multiplying the X matrices
    if(ncol(d$mode1.X) == 0) {
      mode1.X.times.A <- m$mode1.H.mean
    } else if(H1.intercept) {
      mode1.X.times.A <- cbind(1, safe_prod(d$mode1.X, m$mode1.A.mean))
    } else mode1.X.times.A <- safe_prod(d$mode1.X, m$mode1.A.mean)

    if(ncol(d$mode2.X) == 0) {
      mode2.X.times.A <- m$mode2.H.mean
    } else if(H2.intercept) {
      mode2.X.times.A <- cbind(1, safe_prod(d$mode2.X, m$mode2.A.mean))
    } else mode2.X.times.A <- safe_prod(d$mode2.X, m$mode2.A.mean)

    if(ncol(d$mode3.X) == 0) {
      mode3.X.times.A <- m$mode3.H.mean
    } else if(H3.intercept) {
      mode3.X.times.A <- cbind(1, safe_prod(d$mode3.X, m$mode3.A.mean))
    } else mode3.X.times.A <- safe_prod(d$mode3.X, m$mode3.A.mean)

    m$resp <- mult_3d(m$core.mean, mode1.X.times.A, mode2.X.times.A, mode3.X.times.A)

    # Update the lower bound
    if(lower.bnd) {
      if(verbose) print('Calculating lower bound')
      m$lower.bnd[m$iter] <- lower_bnd_Tucker(m, d)
      print(paste('Lower bound =',
        format(m$lower.bnd[m$iter], decimal.mark=".", big.mark=",", nsmall = 1)))
    }

    # Update the training performance measures
    if(RMSE) rmse(m, d, verbose=T)                                      # Prints
    if(exp.var) m$exp.var[m$iter] <- exp_var(d$resp, m$resp, verbose=T) # Prints
    if(cor) {
      m$p.cor[m$iter] <- cor(d$resp, m$resp, use='complete.obs')
      print(sprintf("Pearson correlation: %.4f", m$p.cor[m$iter]))
      m$s.cor[m$iter] <- cor(d$resp, m$resp, use='complete.obs', method='spearman')
      print(sprintf("Spearman correlation: %.4f", m$s.cor[m$iter]))
    }

    # Plot the training predictions and projection matrices
    if(plot) {
      panels <- 1
      if(1 %in% show.mode) {
        panels <- panels + 1
        if(!all(dim(d$mode1.X)==0)) panels <- panels + 1
      }
      if(2 %in% show.mode) {
        panels <- panels + 1
        if(!all(dim(d$mode2.X)==0)) panels <- panels + 1
      }
      if(3 %in% show.mode) {
        panels <- panels + 1
        if(!all(dim(d$mode3.X)==0)) panels <- panels + 1
      }
      par(mfrow=c(panels%/%3, ceiling(panels/panels%/%3)))

      plot_preds(m$resp, d$resp, main=paste('Round: ', m$iter))

      show_mat(m, d, show.mode)
    }

    # Check early stopping criterion.
    # Stop if the lower bound increases by less than early.stop
    if(!is.na(m$early.stop) & lower.bnd & (m$iter > 20)) {
      if(m$lower.bnd[m$iter] < (m$lower.bnd[m$iter-1] + early.stop))  st <- T
    }

    if(time) m$times[m$iter] <- (proc.time() - start.time)[3]

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
