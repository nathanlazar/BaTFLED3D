# Parse the output of summarize.R to load easily into excel
# when predicting modes 1, 2 & 3

# nathan dot lazar at gmail dot com

# Usage: Rscript format_summary_m1m2m3.R <run.Rdata> <out.txt>

# Example: Rscript format_summary_m1m2m3.R SIMresults/run_1222681_summary.Rdata \
#  <SIMresults/run_1222681_summary.txt>

args <- commandArgs(TRUE)

in.file <- args[1]
out.file <- args[2]

load(in.file)

out.table <- data.frame(
  train.rmse.mean = mean(results$summaries['train', 'RMSE', ]),
  train.rmse.sd = sd(results$summaries['train', 'RMSE', ]),
  warm.rmse.mean = mean(results$summaries['warm', 'RMSE', ]),
  warm.rmse.sd = sd(results$summaries['warm', 'RMSE', ]),
  m1.rmse.mean = mean(results$summaries['m1', 'RMSE', ]),
  m1.rmse.sd = sd(results$summaries['m1', 'RMSE', ]),
  m2.rmse.mean = mean(results$summaries['m2', 'RMSE', ]),
  m2.rmse.sd = sd(results$summaries['m2', 'RMSE', ]),
  m3.rmse.mean = mean(results$summaries['m3', 'RMSE', ]),
  m3.rmse.sd = sd(results$summaries['m3', 'RMSE', ]),
  m1m2.rmse.mean = mean(results$summaries['m1m2', 'RMSE', ]),
  m1m2.rmse.sd = sd(results$summaries['m1m2', 'RMSE', ]),
  m1m3.rmse.mean = mean(results$summaries['m1m3', 'RMSE', ]),
  m1m3.rmse.sd = sd(results$summaries['m1m3', 'RMSE', ]),
  m2m3.rmse.mean = mean(results$summaries['m2m3', 'RMSE', ]),
  m2m3.rmse.sd = sd(results$summaries['m2m3', 'RMSE', ]),
  m1m2m3.rmse.mean = mean(results$summaries['m1m2m3', 'RMSE', ]),
  m1m2m3.rmse.sd = sd(results$summaries['m1m2m3', 'RMSE', ]),

  train.exp.var.mean = mean(results$summaries['train', 'exp.var', ]),
  train.exp.var.sd = sd(results$summaries['train', 'exp.var', ]),
  warm.exp.var.mean = mean(results$summaries['warm', 'exp.var', ]),
  warm.exp.var.sd = sd(results$summaries['warm', 'exp.var', ]),
  m1.exp.var.mean = mean(results$summaries['m1', 'exp.var', ]),
  m1.exp.var.sd = sd(results$summaries['m1', 'exp.var', ]),
  m2.exp.var.mean = mean(results$summaries['m2', 'exp.var', ]),
  m2.exp.var.sd = sd(results$summaries['m2', 'exp.var', ]),
  m3.exp.var.mean = mean(results$summaries['m3', 'exp.var', ]),
  m3.exp.var.sd = sd(results$summaries['m3', 'exp.var', ]),
  m1m2.exp.var.mean = mean(results$summaries['m1m2', 'exp.var', ]),
  m1m2.exp.var.sd = sd(results$summaries['m1m2', 'exp.var', ]),
  m1m3.exp.var.mean = mean(results$summaries['m1m3', 'exp.var', ]),
  m1m3.exp.var.sd = sd(results$summaries['m1m3', 'exp.var', ]),
  m2m3.exp.var.mean = mean(results$summaries['m2m3', 'exp.var', ]),
  m2m3.exp.var.sd = sd(results$summaries['m2m3', 'exp.var', ]),
  m1m2m3.exp.var.mean = mean(results$summaries['m1m2m3', 'exp.var', ]),
  m1m2m3.exp.var.sd = sd(results$summaries['m1m2m3', 'exp.var', ]),

  train.p.cor.mean = mean(results$summaries['train', 'p.cor', ]),
  train.p.cor.sd = sd(results$summaries['train', 'p.cor', ]),
  warm.p.cor.mean = mean(results$summaries['warm', 'p.cor', ]),
  warm.p.cor.sd = sd(results$summaries['warm', 'p.cor', ]),
  m1.p.cor.mean = mean(results$summaries['m1', 'p.cor', ]),
  m1.p.cor.sd = sd(results$summaries['m1', 'p.cor', ]),
  m2.p.cor.mean = mean(results$summaries['m2', 'p.cor', ]),
  m2.p.cor.sd = sd(results$summaries['m2', 'p.cor', ]),
  m3.p.cor.mean = mean(results$summaries['m3', 'p.cor', ]),
  m3.p.cor.sd = sd(results$summaries['m3', 'p.cor', ]),
  m1m2.p.cor.mean = mean(results$summaries['m1m2', 'p.cor', ]),
  m1m2.p.cor.sd = sd(results$summaries['m1m2', 'p.cor', ]),
  m1m3.p.cor.mean = mean(results$summaries['m1m3', 'p.cor', ]),
  m1m3.p.cor.sd = sd(results$summaries['m1m3', 'p.cor', ]),
  m2m3.p.cor.mean = mean(results$summaries['m2m3', 'p.cor', ]),
  m2m3.p.cor.sd = sd(results$summaries['m2m3', 'p.cor', ]),
  m1m2m3.p.cor.mean = mean(results$summaries['m1m2m3', 'p.cor', ]),
  m1m2m3.p.cor.sd = sd(results$summaries['m1m2m3', 'p.cor', ]),

  train.s.cor.mean = mean(results$summaries['train', 's.cor', ]),
  train.s.cor.sd = sd(results$summaries['train', 's.cor', ]),
  warm.s.cor.mean = mean(results$summaries['warm', 's.cor', ]),
  warm.s.cor.sd = sd(results$summaries['warm', 's.cor', ]),
  m1.s.cor.mean = mean(results$summaries['m1', 's.cor', ]),
  m1.s.cor.sd = sd(results$summaries['m1', 's.cor', ]),
  m2.s.cor.mean = mean(results$summaries['m2', 's.cor', ]),
  m2.s.cor.sd = sd(results$summaries['m2', 's.cor', ]),
  m3.s.cor.mean = mean(results$summaries['m3', 's.cor', ]),
  m3.s.cor.sd = sd(results$summaries['m3', 's.cor', ]),
  m1m2.s.cor.mean = mean(results$summaries['m1m2', 's.cor', ]),
  m1m2.s.cor.sd = sd(results$summaries['m1m2', 's.cor', ]),
  m1m3.s.cor.mean = mean(results$summaries['m1m3', 's.cor', ]),
  m1m3.s.cor.sd = sd(results$summaries['m1m3', 's.cor', ]),
  m2m3.s.cor.mean = mean(results$summaries['m2m3', 's.cor', ]),
  m2m3.s.cor.sd = sd(results$summaries['m2m3', 's.cor', ]),
  m1m2m3.s.cor.mean = mean(results$summaries['m1m2m3', 's.cor', ]),
  m1m2m3.s.cor.sd = sd(results$summaries['m1m2m3', 's.cor', ])
)

write.table(out.table, file=out.file, quote=F, sep=',', row.names=F)