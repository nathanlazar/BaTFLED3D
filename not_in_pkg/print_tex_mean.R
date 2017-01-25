# nathan dot lazar at gmail dot com

# Print out response measures and standard deviations for TeX tables

args <- commandArgs(TRUE)
in.data <- args[1]

load(in.data)

all.means <- apply(results$mean, c(1,2), mean)
all.sds <- apply(results$mean, c(1,2), sd)

train.mean <- apply(all.means[,c('train.m1','train.m2','train.m3')], 1, mean)
train.sd <- apply(all.means[,c('train.m1','train.m2','train.m3')], 1, sd)
warm.mean <- apply(all.means[,c('warm.m1','warm.m2','warm.m3')], 1, mean)
warm.sd <- apply(all.means[,c('warm.m1','warm.m2','warm.m3')], 1, sd)

for(measure in c('RMSE', 'exp.var', 'p.cor', 's.cor')) {
  print(measure)
#   print(paste0(paste(round(c(train.mean[measure], warm.mean[measure],
#                              all.means[measure, c('m1', 'm2', 'm3', 
#                                  'm1m2', 'm1m3', 'm2m3', 'm1m2m3')]), 2),
#                      collapse=' & '), ' \\'))
  print(sprintf('%.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f',
                train.mean[measure], warm.mean[measure],
                all.means[measure, 'm1'], all.means[measure, 'm2'],
                all.means[measure, 'm3'], all.means[measure, 'm1m2'],
                all.means[measure, 'm1m3'], all.means[measure, 'm2m3'],
                all.means[measure, 'm1m2m3']))
}

for(measure in c('RMSE', 'exp.var', 'p.cor', 's.cor')) {
  print(measure)
#   print(paste0(paste(round(all.means[c('train', 'warm', 'm1', 'm2', 'm3', 
#                                        'm1m2', 'm1m3', 'm2m3', 'm1m2m3'),measure], 2),
#                      round(all.sds[c('train', 'warm', 'm1', 'm2', 'm3', 
#                                      'm1m2', 'm1m3', 'm2m3', 'm1m2m3'),measure], 2), 
#                      sep='(', collapse=') & '), ') \\'))
  print(sprintf('%.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f)',
                train.mean[measure], train.sd[measure],
                warm.mean[measure], warm.sd[measure],
                all.means[measure, 'm1'], all.sds[measure, 'm1'], 
                all.means[measure, 'm2'], all.sds[measure, 'm2'], 
                all.means[measure, 'm3'], all.sds[measure, 'm3'], 
                all.means[measure, 'm1m2'], all.sds[measure, 'm1m2'], 
                all.means[measure, 'm1m3'], all.sds[measure, 'm1m3'], 
                all.means[measure, 'm2m3'], all.sds[measure, 'm2m3'], 
                all.means[measure, 'm1m2m3'], all.sds[measure, 'm1m2m3']))
}
