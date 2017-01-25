# nathan dot lazar at gmail dot com

# Print out response measures and standard deviations for TeX tables

args <- commandArgs(TRUE)
in.data <- args[1]

load(in.data)

all.means <- apply(results$summaries, c(1,2), mean)
all.sds <- apply(results$summaries, c(1,2), sd)

for(measure in c('RMSE', 'exp.var', 'p.cor', 's.cor')) {
  print(measure)
#   print(paste0(paste(round(all.means[c('train', 'warm', 'm1', 'm2', 'm3', 
#                                  'm1m2', 'm1m3', 'm2m3', 'm1m2m3'),measure], 2),
#                      collapse=' & '), ' \\'))

  print(sprintf('%.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f //',
                all.means['train',measure],
                all.means['warm',measure],
                all.means['m1',measure],
                all.means['m2',measure],
                all.means['m3',measure],
                all.means['m1m2',measure],
                all.means['m1m3',measure],
                all.means['m2m3',measure],
                all.means['m1m2m3',measure]))
}

for(measure in c('RMSE', 'exp.var', 'p.cor', 's.cor')) {
  print(measure)
#   print(paste0(paste(round(all.means[c('train', 'warm', 'm1', 'm2', 'm3', 
#                                        'm1m2', 'm1m3', 'm2m3', 'm1m2m3'),measure], 2),
#                      round(all.sds[c('train', 'warm', 'm1', 'm2', 'm3', 
#                                      'm1m2', 'm1m3', 'm2m3', 'm1m2m3'),measure], 2), 
#                      sep='(', collapse=') & '), ') \\'))
  print(sprintf('%.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f)',
                all.means['train',measure], all.sds['train',measure],
                all.means['warm',measure], all.sds['warm',measure],
                all.means['m1',measure], all.sds['m1',measure],
                all.means['m2',measure], all.sds['m2',measure],
                all.means['m3',measure], all.sds['m3',measure],
                all.means['m1m2',measure], all.sds['m1m2',measure],
                all.means['m1m3',measure], all.sds['m1m3',measure],
                all.means['m2m3',measure], all.sds['m2m3',measure],
                all.means['m1m2m3',measure], all.sds['m1m2m3',measure]))
}
