# Make plots of predictions comparing BaTFLED to mean for each sample

# Usage: head_to_head_mean.R <run_summary.Rdata>

args <- commandArgs(TRUE)

sum_file <- args[1]
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

