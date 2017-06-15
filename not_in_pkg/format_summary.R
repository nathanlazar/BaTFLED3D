# Parse the output of summarize.R to load easily into excel
# when predicting modes 1, 2 & 3

# nathan dot lazar at gmail dot com

# Usage: Rscript format_summary.R <run.Rdata> <out.txt>

# Example: Rscript format_summary.R SIMresults/run_1222681_summary.Rdata \
#  <SIMresults/run_1222681_summary.txt>

args <- commandArgs(TRUE)

in.file <- args[1]
# out.file <- args[2]

load(in.file)

out.table <- data.frame(
  train.rmse.mean = mean(results$summaries['train', 'RMSE', ], na.rm=T),
  train.rmse.sd = sd(results$summaries['train', 'RMSE', ], na.rm=T),
  warm.rmse.mean = mean(results$summaries['warm', 'RMSE', ], na.rm=T),
  warm.rmse.sd = sd(results$summaries['warm', 'RMSE', ], na.rm=T),
  m1.rmse.mean = mean(results$summaries['m1', 'RMSE', ], na.rm=T),
  m1.rmse.sd = sd(results$summaries['m1', 'RMSE', ], na.rm=T),
#  m1.rmse.better = better['m1', 'RMSE'],
  m2.rmse.mean = mean(results$summaries['m2', 'RMSE', ], na.rm=T),
  m2.rmse.sd = sd(results$summaries['m2', 'RMSE', ], na.rm=T),
  m3.rmse.mean = mean(results$summaries['m3', 'RMSE', ], na.rm=T),
  m3.rmse.sd = sd(results$summaries['m3', 'RMSE', ], na.rm=T),
  m1m2.rmse.mean = mean(results$summaries['m1m2', 'RMSE', ], na.rm=T),
  m1m2.rmse.sd = sd(results$summaries['m1m2', 'RMSE', ], na.rm=T),
  m1m3.rmse.mean = mean(results$summaries['m1m3', 'RMSE', ], na.rm=T),
  m1m3.rmse.sd = sd(results$summaries['m1m3', 'RMSE', ], na.rm=T),
  m2m3.rmse.mean = mean(results$summaries['m2m3', 'RMSE', ], na.rm=T),
  m2m3.rmse.sd = sd(results$summaries['m2m3', 'RMSE', ], na.rm=T),
  m1m2m3.rmse.mean = mean(results$summaries['m1m2m3', 'RMSE', ], na.rm=T),
  m1m2m3.rmse.sd = sd(results$summaries['m1m2m3', 'RMSE', ], na.rm=T),

  train.exp.var.mean = mean(results$summaries['train', 'exp.var', ], na.rm=T),
  train.exp.var.sd = sd(results$summaries['train', 'exp.var', ], na.rm=T),
  warm.exp.var.mean = mean(results$summaries['warm', 'exp.var', ], na.rm=T),
  warm.exp.var.sd = sd(results$summaries['warm', 'exp.var', ], na.rm=T),
  m1.exp.var.mean = mean(results$summaries['m1', 'exp.var', ], na.rm=T),
  m1.exp.var.sd = sd(results$summaries['m1', 'exp.var', ], na.rm=T),
#  m1.exp.var.better = better['m1', 'exp.var'],
  m2.exp.var.mean = mean(results$summaries['m2', 'exp.var', ], na.rm=T),
  m2.exp.var.sd = sd(results$summaries['m2', 'exp.var', ], na.rm=T),
  m3.exp.var.mean = mean(results$summaries['m3', 'exp.var', ], na.rm=T),
  m3.exp.var.sd = sd(results$summaries['m3', 'exp.var', ], na.rm=T),
  m1m2.exp.var.mean = mean(results$summaries['m1m2', 'exp.var', ], na.rm=T),
  m1m2.exp.var.sd = sd(results$summaries['m1m2', 'exp.var', ], na.rm=T),
  m1m3.exp.var.mean = mean(results$summaries['m1m3', 'exp.var', ], na.rm=T),
  m1m3.exp.var.sd = sd(results$summaries['m1m3', 'exp.var', ], na.rm=T),
  m2m3.exp.var.mean = mean(results$summaries['m2m3', 'exp.var', ], na.rm=T),
  m2m3.exp.var.sd = sd(results$summaries['m2m3', 'exp.var', ], na.rm=T),
  m1m2m3.exp.var.mean = mean(results$summaries['m1m2m3', 'exp.var', ], na.rm=T),
  m1m2m3.exp.var.sd = sd(results$summaries['m1m2m3', 'exp.var', ], na.rm=T),

  train.p.cor.mean = mean(results$summaries['train', 'p.cor', ], na.rm=T),
  train.p.cor.sd = sd(results$summaries['train', 'p.cor', ], na.rm=T),
  warm.p.cor.mean = mean(results$summaries['warm', 'p.cor', ], na.rm=T),
  warm.p.cor.sd = sd(results$summaries['warm', 'p.cor', ], na.rm=T),
  m1.p.cor.mean = mean(results$summaries['m1', 'p.cor', ], na.rm=T),
  m1.p.cor.sd = sd(results$summaries['m1', 'p.cor', ], na.rm=T),
#  m1.p.cor.better = better['m1', 'p.cor'],
  m2.p.cor.mean = mean(results$summaries['m2', 'p.cor', ], na.rm=T),
  m2.p.cor.sd = sd(results$summaries['m2', 'p.cor', ], na.rm=T),
  m3.p.cor.mean = mean(results$summaries['m3', 'p.cor', ], na.rm=T),
  m3.p.cor.sd = sd(results$summaries['m3', 'p.cor', ], na.rm=T),
  m1m2.p.cor.mean = mean(results$summaries['m1m2', 'p.cor', ], na.rm=T),
  m1m2.p.cor.sd = sd(results$summaries['m1m2', 'p.cor', ], na.rm=T),
  m1m3.p.cor.mean = mean(results$summaries['m1m3', 'p.cor', ], na.rm=T),
  m1m3.p.cor.sd = sd(results$summaries['m1m3', 'p.cor', ], na.rm=T),
  m2m3.p.cor.mean = mean(results$summaries['m2m3', 'p.cor', ], na.rm=T),
  m2m3.p.cor.sd = sd(results$summaries['m2m3', 'p.cor', ], na.rm=T),
  m1m2m3.p.cor.mean = mean(results$summaries['m1m2m3', 'p.cor', ], na.rm=T),
  m1m2m3.p.cor.sd = sd(results$summaries['m1m2m3', 'p.cor', ], na.rm=T),

  train.s.cor.mean = mean(results$summaries['train', 's.cor', ], na.rm=T),
  train.s.cor.sd = sd(results$summaries['train', 's.cor', ], na.rm=T),
  warm.s.cor.mean = mean(results$summaries['warm', 's.cor', ], na.rm=T),
  warm.s.cor.sd = sd(results$summaries['warm', 's.cor', ], na.rm=T),
  m1.s.cor.mean = mean(results$summaries['m1', 's.cor', ], na.rm=T),
  m1.s.cor.sd = sd(results$summaries['m1', 's.cor', ], na.rm=T),
#  m1.s.cor.better = better['m1', 's.cor'],
  m2.s.cor.mean = mean(results$summaries['m2', 's.cor', ], na.rm=T),
  m2.s.cor.sd = sd(results$summaries['m2', 's.cor', ], na.rm=T),
  m3.s.cor.mean = mean(results$summaries['m3', 's.cor', ], na.rm=T),
  m3.s.cor.sd = sd(results$summaries['m3', 's.cor', ], na.rm=T),
  m1m2.s.cor.mean = mean(results$summaries['m1m2', 's.cor', ], na.rm=T),
  m1m2.s.cor.sd = sd(results$summaries['m1m2', 's.cor', ], na.rm=T),
  m1m3.s.cor.mean = mean(results$summaries['m1m3', 's.cor', ], na.rm=T),
  m1m3.s.cor.sd = sd(results$summaries['m1m3', 's.cor', ], na.rm=T),
  m2m3.s.cor.mean = mean(results$summaries['m2m3', 's.cor', ], na.rm=T),
  m2m3.s.cor.sd = sd(results$summaries['m2m3', 's.cor', ], na.rm=T),
  m1m2m3.s.cor.mean = mean(results$summaries['m1m2m3', 's.cor', ], na.rm=T),
  m1m2m3.s.cor.sd = sd(results$summaries['m1m2m3', 's.cor', ], na.rm=T)
)

if(binary) {
  out.table$train.acc.mean <- mean(results$summaries['train', 'acc',], na.rm=T)
  out.table$train.acc.sd <- sd(results$summaries['train', 'acc', ], na.rm=T)
  out.table$warm.acc.mean <- mean(results$summaries['warm', 'acc', ], na.rm=T)
  out.table$warm.acc.sd <- sd(results$summaries['warm', 'acc', ], na.rm=T)
  out.table$m1.acc.mean <- mean(results$summaries['m1', 'acc', ], na.rm=T)
  out.table$m1.acc.sd <- sd(results$summaries['m1', 'acc', ], na.rm=T)
#  out.table$m1.acc.better <- better['m1', 'acc']
  out.table$m2.acc.mean <- mean(results$summaries['m2', 'acc', ], na.rm=T)
  out.table$m2.acc.sd <- sd(results$summaries['m2', 'acc', ], na.rm=T)
  out.table$m3.acc.mean <- mean(results$summaries['m3', 'acc', ], na.rm=T)
  out.table$m3.acc.sd <- sd(results$summaries['m3', 'acc', ], na.rm=T)
  out.table$m1m2.acc.mean <- mean(results$summaries['m1m2', 'acc', ], na.rm=T)
  out.table$m1m2.acc.sd <- sd(results$summaries['m1m2', 'acc', ], na.rm=T)
  out.table$m1m3.acc.mean <- mean(results$summaries['m1m3', 'acc', ], na.rm=T)
  out.table$m1m3.acc.sd <- sd(results$summaries['m1m3', 'acc', ], na.rm=T)
  out.table$m2m3.acc.mean <- mean(results$summaries['m2m3', 'acc', ], na.rm=T)
  out.table$m2m3.acc.sd <- sd(results$summaries['m2m3', 'acc', ], na.rm=T)
  out.table$m1m2m3.acc.mean <- mean(results$summaries['m1m2m3', 'acc', ], na.rm=T)
  out.table$m1m2m3.acc.sd <- sd(results$summaries['m1m2m3', 'acc', ], na.rm=T)
}

# Warn about NAs
for(type in dimnames(results$summaries)[[1]]) 
  for(resp in dimnames(results$summaries)[[2]]) {
    nas <- sum(is.na(results$summaries[type, resp,]))
    if(nas > 0 & nas < dim(results$summaries)[3])
      print(sprintf('%s prediction, %s response has %d NA value(s) out of %d', 
                    type, resp, nas, nrow(out.table)))
}

out.table[is.na(out.table)] <- NA

print(paste0(names(out.table), sep='', collapse=','))
print(paste0(out.table[1,], sep='', collapse=','))

