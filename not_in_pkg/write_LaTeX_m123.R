library(gtools)
library(dplyr)
args <- commandArgs(TRUE)

test <- T
if(test) {
  args <- c('Results/design_matrix.txt', 'Results/')
  if(exists('results.df')) rm(results.df)  
}

design.df <- data.frame(read.table(args[[1]], sep='\t', header=T, stringsAsFactors=F))
dir <- args[2]

design.df$file <- sapply(design.df$run, function(x)
  grep('summary', grep(x, list.files(path=dir, recursive=T), value=T), value=T))

loadRData <- function(file) {
  load(file)
  return(results$summaries)
}

loadMean <- function(file) {
  load(file)
  return(results$mean)
}

# First get results for mean prediction
for(i in which(design.df$class=='Mean')) {
  file <- design.df$file[i]
  res <- loadMean(paste0(dir, file))
  res <- res[,!(dimnames(res)[[2]] %in% c('train.m1m2', 'train.m1m3', 'train.m2m3', 'train.m1m2m3',
                                          'warm.m1m2', 'warm.m1m3', 'warm.m2m3', 'warm.m1m2m3')),]
  res.mean <- apply(res, c(1,2), mean)
  res.sd <- apply(res, c(1,2), sd)
  colnames(res.sd) <- paste0(colnames(res.sd), '.sd')
  res.df <- as.data.frame(cbind(res.mean, res.sd))
  res.df <- res.df[,apply(!is.na(res.df), 2, sum)>0]
  # Split train.m1, train.m2, train.m3 results
  res.df$model <- design.df$model[i]
  if(design.df$task[i] != 'Warm')
    res.df <- res.df[,!grepl('train', colnames(res.df))]
  else {
    if(design.df$model[i] == 'Cl. mean')
      res.df <- res.df[,!grepl('m2|m3', colnames(res.df))]
    if(design.df$model[i] == 'Lig. mean')
      res.df <- res.df[,!grepl('m1|m3', colnames(res.df))]
    if(design.df$model[i] == 'ECMP mean')
      res.df <- res.df[,!grepl('m1|m2', colnames(res.df))]
    colnames(res.df) <- sub('.m3', '', sub('.m2', '', sub('.m1', '', colnames(res.df))))
  }
  res.df$task <- design.df$task[i]
  res.df$measure <- rownames(res.df)
  # Remove sub task results (e.g. predicting cl when predicting cl & ligand)
  if(design.df$task[i] != 'Cell lines')
    res.df <- res.df[,!(colnames(res.df) %in% c('m1', 'm1.sd'))]
  if(design.df$task[i] != 'Ligands')
    res.df <- res.df[,!(colnames(res.df) %in% c('m2', 'm2.sd'))]
  if(design.df$task[i] != 'ECMPs')
    res.df <- res.df[,!(colnames(res.df) %in% c('m3', 'm3.sd'))]
  if(design.df$task[i] != 'Cl/ligands')
    res.df <- res.df[,!(colnames(res.df) %in% c('m1m2', 'm1m2.sd'))]
  if(design.df$task[i] != 'Cl/ECMPs')
    res.df <- res.df[,!(colnames(res.df) %in% c('m1m3', 'm1m3.sd'))]
  if(design.df$task[i] != 'Ligand/ECMPs')
    res.df <- res.df[,!(colnames(res.df) %in% c('m2m3', 'm2m3.sd'))]
  if(!exists('results.df')) {
    results.df <- res.df
  } else {
    results.df <- smartbind(results.df, res.df)
  }
}

for(i in which(design.df$class!='Mean')) {
  res <- loadRData(paste0(dir, design.df$file[i]))
  res.mean <- t(apply(res, c(1,2), mean))
  res.sd <- t(apply(res, c(1,2), sd))
  colnames(res.sd) <- paste0(colnames(res.sd), '.sd')
  res.df <- as.data.frame(cbind(res.mean, res.sd))
  res.df <- res.df[!grepl('clip|min|max|iter', rownames(res.df)),
                   apply(res.df, 2, function(x) sum(!is.na(x)))>0]
  res.df <- res.df[,!grepl('^H', names(res.df))]
  res.df$model <- design.df$model[i]
  res.df$task <- design.df$task[i]
  res.df$measure <- rownames(res.df)
  if(design.df$task[i] != 'Warm') # Only use warm and train responses from warm runs
    res.df <- select(res.df, -train, -train.sd)
  # Remove sub task results (e.g. predicting cl when predicting cl & ligand)
  if(design.df$task[i] != 'Cell lines')
    res.df <- res.df[,!(colnames(res.df) %in% c('m1', 'm1.sd'))]
  if(design.df$task[i] != 'Ligands')
    res.df <- res.df[,!(colnames(res.df) %in% c('m2', 'm2.sd'))]
  if(design.df$task[i] != 'ECMPs')
    res.df <- res.df[,!(colnames(res.df) %in% c('m3', 'm3.sd'))]
  if(design.df$task[i] != 'Cl/ligands')
    res.df <- res.df[,!(colnames(res.df) %in% c('m1m2', 'm1m2.sd'))]
  if(design.df$task[i] != 'Cl/ECMPs')
    res.df <- res.df[,!(colnames(res.df) %in% c('m1m3', 'm1m3.sd'))]
  if(design.df$task[i] != 'Ligand/ECMPs')
    res.df <- res.df[,!(colnames(res.df) %in% c('m2m3', 'm2m3.sd'))]

  # if(design.df$model[i] %in% results.df$model) { # Fill in values 
  #    results.df[results.df$model==design.df$model[i],
  #               colnames(res.df)] <- res.df
  results.df <- smartbind(results.df, res.df)
}

results.df <- results.df[,c('measure', 'model', 'train', 'train.sd', 'warm', 'warm.sd',
                            'm1', 'm1.sd', 'm2', 'm2.sd', 'm3', 'm3.sd', 
                            'm1m2', 'm1m2.sd', 'm1m3', 'm1m3.sd', 'm2m3', 'm2m3.sd',
                            'm1m2m3', 'm1m2m3.sd')]

# Collapse results into rows for each model
results.tmp <- results.df
for(i in 1:nrow(results.df)) {
  for(j in 3:ncol(results.df)) {
    if(is.na(results.df[i,j]))
      values <- c(NA)
      if(grepl('mean', results.df$model[i])) {
        values <- results.df[results.df$measure==results.df$measure[i] & results.df$model=='Mean', j]
      } else {
        values <- results.df[results.df$measure==results.df$measure[i] & results.df$model==results.df$model[i], j]
      }
      if(sum(!is.na(values)))
        results.df[i,j] <- unique(values[!is.na(values)])
  }
}

results.df <- unique(results.df)

for(meas in c('RMSE', 'exp.var', 'p.cor', 's.cor')) {
  results.df[results.df$measure==meas & results.df$model=='Cl. mean',]
}

for(meas in c('RMSE', 'exp.var', 'p.cor', 's.cor')) {
  sub.df <- results.df[results.df$measure==meas,]
  sub.df <- sub.df[,-1]
  # sub.df <- sub.df[order(sub.df$model)]
  cat(sprintf('**** %s ****\n', meas))
  cat('\\resizebox*{\\textwidth}{!}{\n')
  cat('\\begin{tabular}{lccccccccc}\n')
  cat('\\toprule\n')
  cat('& Training & Warm & Cell lines & Ligands & ECMPs')
  cat('& \\multicolumn{1}{p{1.5cm}}{\\centering Cell lines \\& Ligands}')
  cat('& \\multicolumn{1}{p{1.5cm}}{\\centering Cell lines \\& ECMPs}')
  cat('& \\multicolumn{1}{p{1.5cm}}{\\centering Ligands \\& ECMPs} & All \\\\\n')
  cat('\\cmidrule{2-10}\n')
  for(i in 1:nrow(sub.df))
    cat(do.call(sprintf, 
      c(list('%s & %.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f) \\\\\n'), 
                             sub.df[i,])))
  cat('\\bottomrule\n')
  cat('\\end{tabular}}\n')
}

# OLD
################################################

# for(i in which(design.df$class=='Mean')) {
#   file <- design.df$file[i]
#   res <- loadMean(paste0(dir, file))
#   res.mean <- t(apply(res, c(1,2), mean))
#   res.sd <- t(apply(res, c(1,2), sd))
#   colnames(res.sd) <- paste0(colnames(res.sd), '.sd')
#   res.df <- as.data.frame(cbind(res.mean, res.sd))
#   res.df$model <- design.df$model[i]
#   res.df$task <- design.df$task[i]
#   res.df$measure <- rownames(res.df)
#   res.df <- res.df[,!grepl('m1m|m2m', colnames(res.df))]
#   if(!exists('results.df')) {
#     results.df <- res.df
#   } else {
#     results.df <- smartbind(results.df, res.df)
#   }
# }
# 
# # Remove mean results for side effect tasks (i.e. m1 when predicting m1m2)
# results.df <- results.df %>% 
#   filter((task=='Warm' & measure %in% c('train.m1', 'warm.m1')) |
#          (task=='Cell lines' & measure=='m1') |
#          (task=='Ligands' & measure=='m2') |
#          (task=='ECMPs' & measure=='m3') |
#          (task=='Cl/ligands' & measure=='m1m2') |
#          (task=='Cl/ECMPs' & measure=='m1m3') |
#          (task=='Ligand/ECMPs' & measure=='m2m3') |
#          (task=='Cl/ligand/ECMPs' & measure=='m1m2m3'))
# Rotate to get a models X tasks matrix
# results.df$measure <- sub('.m2', '', sub('.m1', '', results.df$measure, fixed=T), fixed=T)
# results.mat <- as.matrix(select(results.df, -task, -model, -measure))
# rownames(results.mat) <- results.df$measure
# results.mat <- t(results.mat)
# res.mean <- results.mat[!grepl('sd', rownames(results.mat)),]
# res.sd <- results.mat[grepl('sd', rownames(results.mat)),]
# colnames(res.sd) <- paste0(colnames(res.sd), '.sd')
# results.mat <- cbind(res.mean, res.sd)
# results.df <- data.frame(apply(results.mat, 2, as.numeric))
# results.df$measure <- rownames(results.mat)
# results.df$model <- 'Mean'
# results.df$task <- c('Warm', 'Cell lines', 'Drugs', 'Cl. & Dr.')
