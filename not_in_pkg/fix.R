
pkgTest <- function(x) { # load packages and install if needed
  if (!require(x,character.only = TRUE)) {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

pkgTest('foreach')      # For parallel processing in loops
pkgTest('R6')           # For making memory efficent R6 objects
pkgTest('rTensor')      # Needed for multiplying matrices into tensors (could be removed)
pkgTest('dplyr')        # General data frame manipulation
library(BaTFLED3D)

args <- commandArgs(TRUE)

load(args[1])

all.resp.train <- train.data$resp
all.resp.train[is.na(all.resp.train)] <- warm.resp

save.image(file=args[1])