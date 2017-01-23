# I'm wondering if modeling the variance of responses better would help prediction
# Hints that this would work: 
# - poster by --- at ISMB
# - poor performance on summary measures (maybe these could just be z-scaled for each measure)

load("D:/Box Sync/Thesis/DREAM7/DREAM7_resp.Rdata")
dim(resp.tens)
par(mfrow=c(5,2))
apply(resp.tens, 3, hist, breaks=100)
