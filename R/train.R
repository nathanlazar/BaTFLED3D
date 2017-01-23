#' Train model using BaTFLED algorthm 
#' 
#' Model objects are updated in place to avoid memory issues. Nothing is returned.
#'
#' @export
#' @param d an input data object created with \code{input_data}
#' @param m a \code{CP_model} or \code{Tucker_model} object created with \code{mk_model} 
#' @param ... extra arguments (params) passed to train_CP or train_Tucker
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

train <- function(d, m, ...) {
  if(length(m$core.mean)==0) {
    train_CP(d, m, ...)
  } else {
    train_Tucker(d, m, ...)
  }
}
