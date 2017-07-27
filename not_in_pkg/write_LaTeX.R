library(gtools)
library(dplyr)
args <- commandArgs(TRUE)

test <- T
if(test) {
  # args <- c('Results/CV_design_matrix.txt', 'Results/CV/')
  # args <- c('Results/CV_kern_design_matrix.txt', 'Results/CV_kernels/')
  # args <- c('Results/Final_design_matrix.txt', 'Results/Final/')
  args <- c('Results/Final_kern_design_matrix.txt', 'Results/Final_kernels/')
  if(exists('results.df')) rm(results.df)  
}

design.df <- data.frame(read.table(args[[1]], sep='\t', header=T, stringsAsFactors=F))
dir <- args[2]

design.df$file <- sapply(design.df$run, function(x)
  grep(x, list.files(path=dir, recursive=T), value=T))

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
  res.mean <- t(apply(res, c(1,2), mean))
  res.sd <- t(apply(res, c(1,2), sd))
  colnames(res.sd) <- paste0(colnames(res.sd), '.sd')
  res.df <- as.data.frame(cbind(res.mean, res.sd))
  res.df$model <- 'Mean'
  res.df$task <- design.df$task[i]
  res.df$measure <- rownames(res.df)
  if(!exists('results.df')) {
    results.df <- res.df
  } else {
    results.df <- smartbind(results.df, res.df)
  }
}

# Remove rows of NAs 
results.df <- results.df[!is.na(results.df$RMSE),]

# Remove mean results for side effect tasks (i.e. m1 when predicting m1m2)
results.df <- results.df %>% 
  filter((task=='Warm' & measure %in% c('train.m1', 'warm.m1')) |
         (task=='Cell lines' & measure=='m1') |
         (task=='Drugs' & measure=='m2') |
         (task=='Cl. & Dr.' & measure=='m1m2'))

# Rotate to get a models X tasks matrix
results.df$measure <- sub('.m2', '', sub('.m1', '', results.df$measure, fixed=T), fixed=T)
results.mat <- as.matrix(select(results.df, -task, -model, -measure))
rownames(results.mat) <- results.df$measure
results.mat <- t(results.mat)
res.mean <- results.mat[!grepl('sd', rownames(results.mat)),]
res.sd <- results.mat[grepl('sd', rownames(results.mat)),]
colnames(res.sd) <- paste0(colnames(res.sd), '.sd')
results.mat <- cbind(res.mean, res.sd)
results.df <- data.frame(apply(results.mat, 2, as.numeric))
results.df$measure <- rownames(results.mat)
results.df$model <- 'Mean'
results.df$task <- c('Warm', 'Cell lines', 'Drugs', 'Cl. & Dr.')

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
  if(design.df$task[i] == 'Cl. & Dr.') # remove m1 and m2 responses from m1m2 runs
    res.df <- select(res.df, -m1, -m1.sd, -m2, -m2.sd)

  if(design.df$model[i] %in% results.df$model) { # Fill in values 
     results.df[results.df$model==design.df$model[i],
                colnames(res.df)] <- res.df
  } else {
    results.df <- smartbind(results.df, res.df)
  }
}

results.df <- results.df[,c('measure', 'model', 'train', 'train.sd', 'warm', 'warm.sd',
                            'm1', 'm1.sd', 'm2', 'm2.sd', 'm1m2', 'm1m2.sd')]

for(meas in c('RMSE', 'exp.var', 'p.cor', 's.cor')) {
  sub.df <- results.df[results.df$measure==meas,]
  sub.df <- sub.df[,-1]
  # sub.df <- sub.df[order(sub.df$model)]
  cat(sprintf('**** %s ****\n', meas))
  cat('\\resizebox*{\\textwidth}{!}{\n')
  cat('\\begin{tabular}{lccccc}\n')
  cat('\\toprule\n')
  cat('& Training & Warm & Cell lines & Drugs\n')
  cat('& \\multicolumn{1}{p{1.5cm}}{\\centering Cell lines \\& drugs}  \\\\\n')
  cat('\\cmidrule{2-6}\n')
  for(i in 1:nrow(sub.df))
    cat(do.call(sprintf, c(list('%s & %.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f) \\\\\n'), 
                             sub.df[i,])))
  cat('\\bottomrule\n')
  cat('\\end{tabular}}\n')
}
