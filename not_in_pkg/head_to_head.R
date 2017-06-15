# Make plots of predictions comparing BaTFLED to mean for each sample

# Usage1: head_to_head_mean.R <run_prefix1> <run_prefix2>
# Usage2: head_to_head.R <run_prefix> mean

args <- commandArgs(TRUE)

run1 <- args[1]
run2 <- args[2]

# Determine the number of runs with this prefix
n.files <- length(list.files(path = dirname(run1),
  pattern = paste0(basename(run1), '.[0-9]+.out')))

if(run2=='mean') {
# Compare to predicting the mean

  for(fld in 1:nfiles) {

    load(paste0(run1, '.', (fld-1), '/image.Rdata'))

    if(nrow(test.m1.mat)) {
      
    }


  }
}

load(sum_file)

prefix.split <- strsplit(run_prefix, split='/')
run_prefix <- prefix.split[[1]][length(prefix.split[[1]])]

for(resp in dimnames(results$summaries)[[2]])
  for(type in dimnames(results$summaries)[[1]]) 
    for(type.mean in dimnames(results$mean)[[2]])
      if(resp %in% dimnames(results$mean)[[1]]) {
        pred <- results$summaries[type, resp,]
        mean <- results$mean[resp, type.mean,]
        # print(pred)
        # print(mean)
        if(sum(is.na(pred))==0 & sum(is.na(mean))==0) {
          png(paste0(run_prefix, '_', type.mean, '_', resp, '.png'))
          plot(mean, pred, xlim=range(pred,mean), ylim=range(pred,mean),
               main=paste(type, resp), pch=19, cex = 2)
          abline(0,1, col='blue')
          # Write percentages on plots
          upper <- sum(pred > mean)/length(pred)
          lower <- sum(mean > pred)/length(pred)
          mtext(sprintf('%.2f%%', upper*100),side=3,line=-1.5, 
                at=par("usr")[1]+0.5*diff(par("usr")[1:2]), cex=1.2)
          mtext(sprintf('%.2f%%', lower*100),side=1,line=-1.5, 
                at=par("usr")[1]+.5*diff(par("usr")[1:2]), cex=1.2)
          dev.off()
        } else {
          # Reset these if only one was all NAs
          pred <- c(NA)
          mean <- c(NA)
        }
      }

