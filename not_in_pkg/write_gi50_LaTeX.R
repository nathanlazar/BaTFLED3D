library(gtools)
library(dplyr)
library(reshape2)
args <- commandArgs(TRUE)

load('Results/gi50_predictions/cv_results.Rdata')
# load('Results/gi50_predictions/final_results.Rdata')

# Rotate to get a models X tasks matrix
mat.list <- list()
for(meas in c('NRMSE', 'Exp. var.', 'Pearson cor.', 'Spearman cor.')) {
  means <- filter(results.df, stat=='mean', measure==meas) %>%
    select(model, mode, value) %>%
    dcast(model~mode)
  sds <- filter(results.df, stat=='sd', measure==meas) %>%
    select(model, mode, value) %>%
    dcast(model~mode)
  colnames(sds)[-1] <- paste0(colnames(sds)[-1], '.sd')

  mat <- merge(means, sds, by='model')  
  mat$model <- factor(mat$model, levels=c(
    'Mean', 'LASSO a=1', 'E. net a=0.9', 'E. net a=0.5', 'Ridge a=0',
    'R.F. 1000x5', 'R.F. 5000x5', 'R.F. 1000x10', 
    'N.N. 1L', 'N.N. 2L', 'N.N. 3L',
    'CP', 'CP retrain', 'CP selected', 
    'Tucker', 'Tucker retrain', 'Tucker selected'))
  mat <- arrange(mat, model)
  col.order <- c('model', as.vector(rbind(colnames(means)[-1], colnames(sds)[-1])))
  mat <- mat[,col.order]
  mat.list[[meas]] <- mat
} 

for(meas in c('NRMSE', 'Exp. var.', 'Pearson cor.', 'Spearman cor.')) {
  cat(sprintf('**** %s ****\n', meas))
  cat('\\resizebox*{\\textwidth}{!}{\n')
  cat('\\begin{tabular}{lccccc}\n')
  cat('\\toprule\n')
  cat('& Training* & Warm & Cell lines & Drugs\n')
  cat('& \\multicolumn{1}{p{1.5cm}}{\\centering Cell lines \\& drugs}  \\\\\n')
  cat('\\cmidrule{2-6}\n')
  for(i in 1:nrow(mat.list[[meas]]))
    cat(do.call(sprintf, c(list('%s & %.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f) & %.2f(%.2f) \\\\\n'), 
                             mat.list[[meas]][i,])))
  cat('\\bottomrule\n')
  cat('\\end{tabular}}\n')
}
