diagnostics <- function(m) {
  tbl <- matrix('', 4, 4, dimnames=list(c('A.mean', 'A.cov', 'H.mean', 'H.var'),
                                        c('Mode 1', 'Mode 2', 'Mode 3', 'Core')))

  tbl[1,1] <- sprintf('(%.3f, %.3f)', min(m$mode1.A.mean), max(m$mode1.A.mean))
  tbl[2,1] <- sprintf('(%.3f, %.3f)', min(m$mode1.A.cov), max(m$mode1.A.cov))
  tbl[3,1] <- sprintf('(%.3f, %.3f)', min(m$mode1.H.mean), max(m$mode1.H.mean))
  tbl[4,1] <- sprintf('(%.3f, %.3f)', min(m$mode1.H.var), max(m$mode1.H.var))

  tbl[1,2] <- sprintf('(%.3f, %.3f)', min(m$mode2.A.mean), max(m$mode2.A.mean))
  tbl[2,2] <- sprintf('(%.3f, %.3f)', min(m$mode2.A.cov), max(m$mode2.A.cov))
  tbl[3,2] <- sprintf('(%.3f, %.3f)', min(m$mode2.H.mean), max(m$mode2.H.mean))
  tbl[4,2] <- sprintf('(%.3f, %.3f)', min(m$mode2.H.var), max(m$mode2.H.var))
  
  tbl[1,3] <- sprintf('(%.3f, %.3f)', min(m$mode3.A.mean), max(m$mode3.A.mean))
  tbl[2,3] <- sprintf('(%.3f, %.3f)', min(m$mode3.A.cov), max(m$mode3.A.cov))
  tbl[3,3] <- sprintf('(%.3f, %.3f)', min(m$mode3.H.mean), max(m$mode3.H.mean))
  tbl[4,3] <- sprintf('(%.3f, %.3f)', min(m$mode3.H.var), max(m$mode3.H.var))

  tbl[1,4] <- sprintf('(%.3f, %.3f)', min(m$core.mean), max(m$core.mean))
  tbl[2,4] <- sprintf('(%.3f, %.3f)', min(m$core.var), max(m$core.var))

  print(tbl)
}