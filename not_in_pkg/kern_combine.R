# nathan dot lazat at gmail dot com

kern_combine <- function(train.m1.mat, train.m2.mat, train.m3.mat, 
  test.m1.mat, test.m2.mat, test.m3.mat, 
  kern.combine='mean', args) {
  # Apply kernels and combine kernelized input matrices

  ret.list <- list()

  if(ncol(m1.mat)) {
    if(kern.combine=='linear') {
      ret.list[['orig.train.m1.mat']] <- train.m1.mat
      ret.list[['orig.test.m1.mat']] <- test.m1.mat

      # Replace any NAs with 0
      train.m1.mat[is.na(train.m1.mat)] <- 0
      test.m1.mat[is.na(test.m1.mat)] <- 0

      # Take the inner product (linear kernel)
      if(nrow(test.m1.mat))
        test.m1.mat <- test.m1.mat %*% t(train.m1.mat)
      train.m1.mat <- train.m1.mat %*% t(train.m1.mat)

    } else {
      m1.types <- unique(sapply(strsplit(colnames(m1.mat), split='.', fixed=T), '[', 1))
      # Get the kernel width for each type from args
      m1.s <- rep(1, length(m1.types))
      names(m1.s) <- m1.types
      for(type in m1.types)
        if(sum(grepl(type, args[9:length(args)])))
          m1.s[[type]] <- as.numeric(
            strsplit(grep(type, args[9:length(args)], value=T), split='=')[[1]][2])
  
      m1.train.list <- list()
      if(nrow(test.m1.mat))
        m1.test.list <- list()
      for(i in 1:length(m1.s)) {
        ind <- grep(m1.types[[i]], colnames(train.m1.mat))
        m1.train.list[[i]] <- kernelize(train.m1.mat[,ind,drop=F],
          train.m1.mat[,ind,drop=F], m1.s[[i]])
        names(m1.train.list)[[i]] <- m1.types[[i]]
        if(nrow(test.m1.mat)) {
          m1.test.list[[i]] <- kernelize(test.m1.mat[,ind,drop=F],
            train.m1.mat[,ind,drop=F], m1.s[[i]])
        }
      }
  
      if(kern.combine == 'cat') {
      ########### CHANGE THIS? #############
      # Concatenate the unkernelized annotation data
      # and the other two kernels
      ######################################
        train.m1.mat <- cbind(cbind(
          train.m1.mat[,grep('anno.', colnames(train.m1.mat)), drop=F],
                             m1.train.list[[2]]), m1.train.list[[3]])
        # Add prefixes to column names
        #  colnames(m1.train.list[[i]]) <- paste(m1.types[[i]],
        #  colnames(m1.train.list[[i]]), sep='.')
        #  colnames(m1.test.list[[i]]) <- paste(m1.types[[i]],
        #  colnames(m1.test.list[[i]]), sep='.')
  
        if(nrow(test.m1.mat))
          test.m1.mat <- cbind(cbind(
            test.m1.mat[,grep('anno.', colnames(test.m1.mat)), drop=F],
            m1.test.list[[2]]), m1.test.list[[3]])
      }
      if(kern.combine == 'mean') {
      # Take the average of the kernels
        train.m1.mat <- matrix(0, nrow(m1.train.list[[1]]), nrow(m1.train.list[[1]]),
          dimnames=list(rownames(m1.train.list[[1]]), rownames(m1.train.list[[1]])))
        train.m1.counts <- train.m1.mat
        if(nrow(test.m1.mat)) {
          test.m1.mat <- matrix(0, nrow(m1.test.list[[1]]), nrow(m1.train.list[[1]]),
            dimnames=list(rownames(m1.test.list[[1]]), rownames(m1.train.list[[1]])))
          test.m1.counts <- test.m1.mat
        }
        for(n in 1:length(m1.train.list)) {
          m1.train.list[[n]][is.na(m1.train.list[[n]])] <- 0
  
          train.m1.mat[, colnames(m1.train.list[[n]])] <-
            train.m1.mat[, colnames(m1.train.list[[n]])] +
            m1.train.list[[n]][rownames(train.m1.mat), colnames(m1.train.list[[n]])]
  
          train.m1.counts[, colnames(m1.train.list[[n]])] <-
            train.m1.counts[, colnames(m1.train.list[[n]])] +
            (m1.train.list[[n]][rownames(train.m1.mat), colnames(m1.train.list[[n]])] != 0)
  
          if(nrow(test.m1.mat)) {
            m1.test.list[[n]][is.na(m1.test.list[[n]])] <- 0
  
            test.m1.mat[, colnames(m1.test.list[[n]])] <-
              test.m1.mat[, colnames(m1.test.list[[n]])] +
              m1.test.list[[n]][rownames(test.m1.mat), colnames(m1.test.list[[n]])]
  
            test.m1.counts[, colnames(m1.test.list[[n]])] <-
              test.m1.counts[, colnames(m1.test.list[[n]])] +
              (m1.test.list[[n]][rownames(test.m1.mat), colnames(m1.test.list[[n]])] != 0)
          }
        }
        train.m1.mat <- train.m1.mat/train.m1.counts
        if(nrow(test.m1.mat))
          test.m1.mat <- test.m1.mat/test.m1.counts
      }
  
      if(kern.combine == 'prod') {  ## Untested ##
      # Take the product of the kernels
        train.m1.mat <- m1.train.list[[1]]
        if(nrow(test.m1.mat)) test.m1.mat <- m1.test.list[[1]]
        if(length(m1.train.list) > 1) {
          for(n in 2:length(m1.train.list)) {
            train.m1.mat <- train.m1.mat * m1.train.list[[n]]
            if(nrow(test.m1.mat)) test.m1.mat <- test.m1.mat * m1.test.list[[n]]
          }
        }
      }
  
      if(!nrow(test.m1.mat)) test.m1.mat <- train.m1.mat[0,]
    }
  }
  
  if(ncol(m2.mat)) {
    if(kern.combine=='linear') {
      ret.list[['orig.train.m2.mat']] <- train.m2.mat
      ret.list[['orig.test.m2.mat']] <- test.m2.mat

      # Replace any NAs with 0
      train.m2.mat[is.na(train.m2.mat)] <- 0
      test.m2.mat[is.na(test.m2.mat)] <- 0

      # Take the inner product (linear kernel)
      if(nrow(test.m2.mat))
        test.m2.mat <- test.m2.mat %*% t(train.m2.mat)
      train.m2.mat <- train.m2.mat %*% t(train.m2.mat)

    } else {
      m2.types <- unique(sapply(strsplit(colnames(m2.mat), split='.', fixed=T), '[', 1))
      m2.s <- rep(1, length(m2.types))
      names(m2.s) <- m2.types
      for(type in m2.types)
        if(sum(grepl(type, args[9:length(args)])))
          m2.s[[type]] <- as.numeric(
            strsplit(grep(type, args[9:length(args)], value=T), split='=')[[1]][2])
  
      m2.train.list <- list()
      if(nrow(test.m2.mat)) m2.test.list <- list()
      for(i in 1:length(m2.s)) {
        ind <- grep(m2.types[[i]], colnames(train.m2.mat))
        m2.train.list[[i]] <- kernelize(train.m2.mat[,ind,drop=F],
          train.m2.mat[,ind,drop=F], m2.s[[i]])
        colnames(m2.train.list[[i]]) <- paste(m2.types[[i]],
          colnames(m2.train.list[[i]]), sep='.')
        if(nrow(test.m2.mat)) {
          m2.test.list[[i]] <- kernelize(test.m2.mat[,ind,drop=F],
            train.m2.mat[,ind,drop=F], m2.s[[i]])
          colnames(m2.test.list[[i]]) <- paste(m2.types[[i]],
            colnames(m2.test.list[[i]]), sep='.')
        }
      }
      if(kern.combine == 'cat') {
      ########### CHANGE THIS? #############
      # Concatenate the unkernelized target and class data
      # and the other two kernels
      #######################################
        train.m2.mat <- cbind(cbind(
          train.m2.mat[,grep('targ.|class.', colnames(train.m2.mat)), drop=F],
                             m2.train.list[[3]]), m2.train.list[[4]])
        if(nrow(test.m2.mat))
          test.m2.mat <- cbind(cbind(
            test.m2.mat[,grep('targ.|class.', colnames(test.m2.mat)), drop=F],
            m2.test.list[[3]]), m2.test.list[[4]])
      }
      if(kern.combine == 'mean') {
      # Take the average of the kernels
        train.m2.mat <- m2.train.list[[1]]
        if(nrow(test.m2.mat)) test.m2.mat <- m2.test.list[[1]]
        if(length(m2.train.list) > 1) {
          for(n in 2:length(m2.train.list)) {
            train.m2.mat <- train.m2.mat + m2.train.list[[n]]
            if(nrow(test.m2.mat)) test.m2.mat <- test.m2.mat + m2.test.list[[n]]
          }
          train.m2.mat <- train.m2.mat/length(m2.train.list)
          if(nrow(test.m2.mat)) test.m2.mat <- test.m2.mat/length(m2.test.list)
        }
      }
      if(kern.combine == 'prod') {
      # Take the product of the kernels
        train.m2.mat <- m2.train.list[[1]]
        if(nrow(test.m2.mat)) test.m2.mat <- m2.test.list[[1]]
        if(length(m2.train.list) > 1) {
          for(n in 2:length(m2.train.list)) {
            train.m2.mat <- train.m2.mat * m2.train.list[[n]]
            if(nrow(test.m2.mat)) test.m2.mat <- test.m2.mat * m2.test.list[[n]]
          }
        }
      }
  
      if(!nrow(test.m2.mat)) test.m2.mat <- train.m2.mat[0,]
    }
  }
  
  if(ncol(m3.mat)) {
    if(kern.combine=='linear') {
      ret.list[['orig.train.m3.mat']] <- train.m3.mat
      ret.list[['orig.test.m3.mat']] <- test.m3.mat

      # Replace any NAs with 0
      train.m3.mat[is.na(train.m3.mat)] <- 0
      test.m3.mat[is.na(test.m3.mat)] <- 0

      # Take the inner product (linear kernel)
      if(nrow(test.m3.mat))
        test.m3.mat <- test.m3.mat %*% t(train.m3.mat)
      train.m3.mat <- train.m3.mat %*% t(train.m3.mat)

    } else {
      m3.types <- unique(sapply(strsplit(colnames(m3.mat), split='.', fixed=T), '[', 1))
      m3.s <- rep(1, length(m3.types))
      names(m3.s) <- m3.types
      for(type in m3.types)
        if(sum(grepl(type, args[9:length(args)])))
          m3.s[[type]] <- as.numeric(strsplit(
            grep(type, args[9:length(args)], value=T), split='=')[[1]][2])
  
      m3.train.list <- list()
      if(nrow(test.m3.mat)) m3.test.list <- list()
      for(i in 1:length(m3.s)) {
        ind <- grep(m3.types[[i]], colnames(train.m3.mat))
        m3.train.list[[i]] <- kernelize(train.m3.mat[,ind,drop=F],
          train.m3.mat[,ind,drop=F], m3.s[[i]])
        colnames(m3.train.list[[i]]) <- paste(m3.types[[i]],
          colnames(m3.train.list[[i]]), sep='.')
        if(nrow(test.m3.mat)) {
          m3.test.list[[i]] <- kernelize(test.m3.mat[,ind,drop=F],
            train.m3.mat[,ind,drop=F], m3.s[[i]])
          colnames(m3.test.list[[i]]) <- paste(m3.types[[i]],
            colnames(m3.test.list[[i]]), sep='.')
        }
      }
      if(kern.combine == 'cat') {
      # Concatenate the kernels
        train.m3.mat <- m3.train.list[[1]]
        if(nrow(test.m3.mat)) test.m3.mat <- m3.test.list[[1]]
        if(length(m3.train.list) > 1) {
          for(n in 2:length(m3.train.list)) {
            train.m3.mat <- cbind(train.m3.mat, m3.train.list[[n]])
            if(nrow(test.m3.mat)) test.m3.mat <- cbind(test.m3.mat, m3.test.list[[n]])
          }
        }
      }
  
      if(kern.combine == 'mean') {
      # Take the average of the kernels
        train.m3.mat <- m3.train.list[[1]]
        if(nrow(test.m3.mat)) test.m3.mat <- m3.test.list[[1]]
        if(length(m3.train.list) > 1) {
          for(n in 2:length(m3.train.list)) {
            train.m3.mat <- train.m3.mat + m3.train.list[[n]]
            if(nrow(test.m3.mat)) test.m3.mat <- test.m3.mat + m3.test.list[[n]]
          }
          if(nrow(test.m3.mat)) train.m3.mat <- train.m3.mat/length(m3.train.list)
          test.m3.mat <- test.m3.mat/length(m3.test.list)
        }
      }
      if(kern.combine == 'prod') {
      # Take the product of the kernels
        train.m3.mat <- m3.train.list[[1]]
        if(nrow(test.m3.mat)) test.m3.mat <- m3.test.list[[1]]
        if(length(m3.train.list) > 1) {
          for(n in 2:length(m3.train.list)) {
            train.m3.mat <- train.m3.mat * m3.train.list[[n]]
            if(nrow(test.m3.mat)) test.m3.mat <- test.m3.mat * m3.test.list[[n]]
          }
        }
      }
  
      if(!nrow(test.m3.mat)) test.m3.mat <- train.m3.mat[0,]
    }
  }
  
  ret.list[['train.m1.mat']] <- train.m1.mat
  ret.list[['test.m1.mat']] <- test.m1.mat
  ret.list[['train.m2.mat']] <- train.m2.mat
  ret.list[['test.m2.mat']] <- test.m2.mat
  ret.list[['train.m3.mat']] <- train.m3.mat
  ret.list[['test.m3.mat']] <- test.m3.mat
  return(ret.list)
}
