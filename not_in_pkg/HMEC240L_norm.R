
library(synapseClient)
library(dplyr)
library(BaTFLED3D)
library(rTensor)
library(ggplot2)
library(reshape2)
ifelse(.Platform$OS.type == "windows", library(doParallel), library(doMC))

# Functions
############
plot_by_well <- function(tens, name) {
  png(paste0(name, '.png'), height=480*1.5)
  rng <- range(tens, na.rm=T)
  par(mfrow=c(4,2), mar=c(1,6,2,1))
  for(i in 1:8) {
    im_mat(tens[i,,], main=dimnames(tens)[[1]][i], zlim=rng)
    axis(2, at = seq(0, 1,length.out=8), labels=rev(dimnames(tens)[[2]]), tick=F, las=2)
  }
  dev.off()
}
############

## Change this to prompt for password?
# Look in ignore.txt file
# synapseLogin('nathanlazar', **********)

# HMEC240L_SS1_Level3 <- synGet(id='syn7121242')
# HMEC240L_SS4_Level3 <- synGet(id='syn7121244')
# HMEC240L_SSC_Level3 <- synGet(id='syn7122627')

# Make an output directory for the cell line.
cl <- 'HMEC240L'
if (!file.exists(cl))
  dir.create(file.path(cl))

# MCF10A data 
# Nuclei_PA_gated_Edu_positiveProportion in staining set 2

# ss4 <- tbl_df(read.table(HMEC240L_SS4_Level3@filePath, sep='\t', header=T))
ss4 <- tbl_df(read.table('/home/users/lazar/.synapseCache/218/13087218/HMEC240L_SS4_Level3.txt',
                          sep='\t', header=T))

# ss4.nuclei <- ss4 %>% select(Well, Spot, ArrayRow, ArrayColumn, ECMp, Ligand, MEP, 
#                              BWL, k, matches('Nuclei_PA'))

ss4.edu <- select(ss4, Well, Spot, ArrayRow, ArrayColumn, ECMp, Ligand, MEP, 
                  BWL, k, matches('Nuclei_PA_Gated_EdU'))

rm(ss4)

# Form Nuclei_PA_Gated_EdUPositiveProportionLogit into a tensor 
# Array X plate X spots

unique(ss4.edu$Well)          # 8 wells
unique(ss4.edu$Spot)          # 698 spots (subset of 700 combos below?) should be 694?
unique(ss4.edu$ArrayRow)      # 35 rows of spots?
unique(ss4.edu$ArrayColumn)   # 20 columns of spots?
levels(ss4.edu$ECMp)          # 48 extra-cellular matrix proteins
levels(ss4.edu$MEP)           # 2736 combinations of 57 ligands and 48 ECM protens
levels(ss4.edu$BWL)           # 64 unique plate, well, ligand ids
unique(ss4.edu$k)             # k = 64 always

# Extract plate ids from BWL column
ss4.edu$Plate <- as.factor(sapply(as.character(ss4.edu$BWL), 
                                  function(x) strsplit(x, split='_')[[1]][1]))

# 8 wells X 8 plates X 698 spots (range from 2 to 699)
# Tensor of edu values
edu.tens <- array(NA, dim=c(8,8,698), dimnames=list(levels(ss4.edu$Well), levels(ss4.edu$Plate), 
                                                    paste0('spot.', 2:699)))
# Tensor of normalized edu values (RUVLoess)
edu.ruvloess.tens <- edu.tens
# Tensor of the MEP name at each position
mep.name.tens <- array('', dim=c(8,8,698), dimnames=dimnames(edu.tens))

# fill tensors
for(well in levels(ss4.edu$Well)) {
  for(plate in levels(ss4.edu$Plate)) {
    spots <- ss4.edu %>% filter(Well==well, Plate==plate) %>%
               select(Spot, MEP, Nuclei_PA_Gated_EdUPositiveProportionLogit,
                      Nuclei_PA_Gated_EdUPositiveProportionLogitRUVLoess)
    edu.vec <- rep(NA, 698)
    edu.vec[spots$Spot-1] <- spots$Nuclei_PA_Gated_EdUPositiveProportionLogit
    edu.tens[well, plate, ] <- edu.vec

    edu.ruvloess.vec <- rep(NA, 698)
    edu.ruvloess.vec[spots$Spot-1] <- spots$Nuclei_PA_Gated_EdUPositiveProportionLogitRUVLoess
    edu.ruvloess.tens[well, plate, ] <- edu.ruvloess.vec

    mep.vec <- rep('', 698)
    mep.vec[spots$Spot-1] <- as.character(spots$MEP)
    mep.name.tens[well, plate, ] <- mep.vec
  }
}

plot_by_well(edu.tens, paste0(cl, '/raw_EdU_by_well'))
# for(i in 1:8) im_mat(edu.tens[,i,], main=dimnames(edu.tens)[[2]][i])

plot_by_well(edu.ruvloess.tens, paste0(cl, '/ruvloess_EdU_by_well'))
# for(i in 1:8) im_mat(edu.ruvloess.tens[,i,], main=dimnames(edu.ruvloess.tens)[[2]][i])

# Make a dataframe of the median edu response for each MEP
mep.med.df <- ss4.edu %>% group_by(MEP) %>%
  summarise(med.edu = median(Nuclei_PA_Gated_EdUPositiveProportionLogit))

# Make a tensor of the median response across replicates for each MEP
mep.med.tens <- array(NA, dim=c(8,8,698), dimnames=dimnames(edu.tens))

for(well in levels(ss4.edu$Well)) 
  for(plate in levels(ss4.edu$Plate)) 
    mep.med.tens[well, plate, ] <- 
      mep.med.df$med.edu[match(mep.name.tens[well, plate, ], mep.med.df$MEP,  nomatch=NA)]

plot_by_well(mep.med.tens, paste0(cl, '/median_EdU_by_well'))
# for(i in 1:8) im_mat(mep.med.tens[,i,], main=dimnames(mep.med.tens)[[2]][i])

# Subtract median tensor from edu tensor 
resid.tens <- edu.tens - mep.med.tens

# Look at the residual tensor
plot_by_well(resid.tens, paste0(cl, '/resid_EdU_by_well'))

# Save the tensors
save(mep.name.tens,     file=paste0(cl, '/Nuclei_PA_Gated_EdUPositiveProportionLogit_median_names.Rdata'))
save(edu.tens,          file=paste0(cl, '/Nuclei_PA_Gated_EdUPositiveProportionLogit_tensor.Rdata'))
save(mep.med.tens,      file=paste0(cl, '/Nuclei_PA_Gated_EdUPositiveProportionLogit_median_tensor.Rdata'))
save(resid.tens,        file=paste0(cl, '/Nuclei_PA_Gated_EdUPositiveProportionLogit_residual_tensor.Rdata'))
save(edu.ruvloess.tens, file=paste0(cl, '/Nuclei_PA_Gated_EdUPositiveProportionLogit_RUVLoess_tensor.Rdata'))

# Factor a tensor of just the controls, reconstruct and remove efects.
#################################################
cores <- 20
if(.Platform$OS.type == "windows") {
  clust <- makeCluster(cores)
  registerDoParallel(clust)
} else registerDoMC(cores)

mep.name.tens[1:4,1:4,1:2]
names <- unique(as.vector(mep.name.tens))
names[grepl('COL1', names)]

col1.ind <- mep.med.tens
col1.ind[,,] <- grepl('COL1', mep.name.tens)
sum(col1.ind)/prod(dim(col1.ind))
# 17.18% of responses are these COL1 spots. They are not randomly distributed

plot_by_well(col1.ind, paste0(cl, '/col1_spots.png'))

edu.col1.tens <- edu.tens
edu.col1.tens[!col1.ind] <- NA

plot_by_well(edu.col1.tens, paste0(cl, '/col1_values.png'))

# Factor the colagen 1 tensor w/ BaTFLED3D
cores <- 6
if(.Platform$OS.type == "windows") {
  clust <- makeCluster(cores)
  registerDoParallel(clust)
} else registerDoMC(cores)

train.data <- input_data$new(
  mode1.X=matrix(0, 8, 0, dimnames=list(dimnames(edu.col1.tens)[[1]], c())),
  mode2.X=matrix(0, 8, 0, dimnames=list(dimnames(edu.col1.tens)[[2]], c())),
  mode3.X=matrix(0, 698, 0, dimnames=list(dimnames(edu.col1.tens)[[3]], c())),
  resp=edu.col1.tens)

# Initialize Tucker decomposition model
args <- list('decomp=Tucker', 
             'A1.intercept=F', 'A2.intercept=F', 'A3.intercept=F',
             'H1.intercept=T', 'H2.intercept=T', 'H3.intercept=T',
             'plot=T', 'verbose=F',
             'R1=8', 'R2=8', 'R3=8',
             paste0('sigma2=', var(edu.col1.tens, na.rm=T)),
             'm1.sigma2=1', 'm2.sigma2=1', 'm3.sigma2=1',
             'parallel=T', 'lower.bnd=T', paste0('cores=', cores),
             'update.order=c(3,1,2)', 'show.mode=c(1,2,3)')
model.params <- get_model_params(args)
tuck.model <- mk_model(train.data, model.params)
tuck.model$rand_init(model.params)

# Train the Tucker factorization model
reps <- 20
tuck.trained <- tuck.model$clone()
train(d=train.data, m=tuck.trained, params=model.params, new.iter=reps)

save.image(paste0(cl, '/image.Rdata'))

par(mfrow=c(4,2), mar=c(1,5,2,1))
rng <- range(edu.col1.tens, tuck.trained$resp, na.rm=T)
# Plots of each plate
for(i in 1:8) {
  # im_mat(edu.col1.tens[,i,], main=paste('true', dimnames(edu.col1.tens)[[2]][i]), zlim=rng)
  # axis(2, at = seq(0, 1,length.out=8), 
  #      labels=rev(dimnames(edu.col1.tens)[[1]]), tick=F, las=2)
  im_mat(edu.tens[,i,], main=paste('true', dimnames(edu.col1.tens)[[2]][i]), zlim=rng)
  axis(2, at = seq(0, 1,length.out=8),
       labels=rev(dimnames(edu.col1.tens)[[1]]), tick=F, las=2)
  im_mat(tuck.trained$resp[,i,], main=paste('learned', dimnames(tuck.trained$resp)[[2]][i]), zlim=rng)
  axis(2, at = seq(0, 1,length.out=8), labels=rev(dimnames(tuck.trained$resp)[[1]]), tick=F, las=2)
}

# Use the reconstructed tensor to normalize
norm.tens <- edu.tens - tuck.trained$resp

rng <- range(norm.tens, edu.tens, na.rm=T)
for(i in 1:8) {
  im_mat(edu.tens[,i,], main=paste('original', dimnames(edu.tens)[[2]][i]), zlim=rng)
  axis(2, at = seq(0, 1,length.out=8), labels=rev(dimnames(edu.tens)[[1]]), tick=F, las=2)
  im_mat(norm.tens[,i,], main=paste('normalized', dimnames(norm.tens)[[2]][i]), zlim=rng)
  axis(2, at = seq(0, 1,length.out=8), labels=rev(dimnames(norm.tens)[[1]]), tick=F, las=2)
} 

out.dir <- paste0(cl, '/')

# Compare this normalization to RUV-Loess normalization
png(paste0(out.dir, 'Tensor_vs_ruvloess_norm_R3_', model.params$R3, '.png'), height=480*3, width=480*3)
par(mfrow=c(8,3))
rng <- range(edu.tens, norm.tens, edu.ruvloess.tens, na.rm=T)
for(i in 1:8) {
  im_mat(edu.tens[,i,], main=paste('Original', dimnames(edu.tens)[[2]][i]), zlim=rng)
  axis(2, at = seq(0, 1,length.out=8), labels=rev(dimnames(edu.tens)[[1]]), tick=F, las=2)
  im_mat(norm.tens[,i,], main=paste('Tensor normalized', dimnames(norm.tens)[[2]][i]), zlim=rng)
  axis(2, at = seq(0, 1,length.out=8), labels=rev(dimnames(norm.tens)[[1]]), tick=F, las=2)
  im_mat(edu.ruvloess.tens[,i,], main=paste('RUV Loess norm', dimnames(edu.ruvloess.tens)[[2]][i]), zlim=rng)
  axis(2, at = seq(0, 1,length.out=8), labels=rev(dimnames(edu.ruvloess.tens)[[1]]), tick=F, las=2)
}
dev.off()

# Compare the original tensor to the smoothed version in a scatter plot
png(paste0(out.dir, 'org_vs_tensor_smooth_R3_', model.params$R3, '.png'))
smooth.cor <- cor(norm.tens, edu.tens, use='complete.obs')
plot(norm.tens, edu.tens, col=rgb(0,0,0,0.1),
     main=paste(sprintf('Pearson correlation: %.4f', smooth.cor)),
     xlab="Smoothed EdU positive proportion",
     ylab="EdU positive proportion")
abline(a=0,b=1, lwd=2, col='blue')
dev.off()

# Look at density plots for the two types of normalization
# Get the standard deviation and number of replicates for each MEP in the original data 
# and the smoothed data (the smoothed data fills in missing values, so may have more 'replicates')
mep.df <- data.frame(mep=unique(as.vector(mep.name.tens)), 
                     org.med=0, org.sd=0, org.count=0,
                     ruvloess.smooth.med=0, ruvloess.smooth.sd=0, ruvloess.smooth.count=0,
                     tensor.smooth.med=0, tensor.smooth.sd=0, tensor.smooth.count=0)

for(i in 1:nrow(mep.df)) {
  mep.idx <- which(mep.name.tens == mep.df$mep[i])
  org.vec <- edu.tens[mep.idx]
  ruvloess.smooth.vec <- edu.ruvloess.tens[mep.idx]
  tensor.smooth.vec <- norm.tens[mep.idx]
  mep.df$org.med[i] <- median(org.vec, na.rm=T)
  mep.df$org.sd[i] <- sd(org.vec, na.rm=T)
  mep.df$org.count[i] <- sum(!is.na(org.vec))
  mep.df$ruvloess.smooth.med[i] <- median(ruvloess.smooth.vec, na.rm=T)
  mep.df$ruvloess.smooth.sd[i] <- sd(ruvloess.smooth.vec, na.rm=T)
  mep.df$ruvloess.smooth.count[i] <- sum(!is.na(ruvloess.smooth.vec))
  mep.df$tensor.smooth.med[i] <- median(tensor.smooth.vec, na.rm=T)
  mep.df$tensor.smooth.sd[i] <- sd(tensor.smooth.vec, na.rm=T)
  mep.df$tensor.smooth.count[i] <- sum(!is.na(tensor.smooth.vec))
}

my_theme <- theme(plot.title = element_text(size = 24, face='bold'),
                  axis.title = element_text(size = 20),
                  axis.text.y = element_text(size=18),
                  axis.text.x = element_text(angle = 90, size=11,
                                             hjust = 1, vjust=0.3),
                  plot.margin = unit(c(.6, .3, .3, .3), "in"),
                  legend.title = element_text(size=18),
                  legend.text=element_text(size=12),
                  legend.justification = c(1, 1), 
                  legend.position = c(1, 1))

sd.melt <- melt(mep.df[,c('org.sd', 'ruvloess.smooth.sd', 'tensor.smooth.sd')])
levels(sd.melt$variable) <- c('Original', 'RUV loess', 'Tensor')

sd.dens.plot <- ggplot(sd.melt, aes(x = value, fill = variable)) + 
  geom_density(alpha=.6) +
  labs(x="log of sd across replicates", y="Frequency") +
  ggtitle('Dist of sd. across replicates\n') +
  guides(fill = guide_legend(title = NULL)) +
  coord_cartesian(xlim = c(0, 1)) +
  my_theme

png(paste(out.dir, 'hist_log_of_sd_all_3_R3_', model.params$R3, '.png'))
sd.dens.plot
dev.off()

med.melt <- melt(mep.df[,c('org.med', 'ruvloess.smooth.med', 'tensor.smooth.med')])
levels(sd.melt$variable) <- c('Original', 'RUV loess', 'Tensor')

med.dens.plot <- ggplot(med.melt, aes(x = value, fill = variable)) + 
  geom_density(alpha=.6) +
  labs(x="Medians across replicates", y="Frequency") +
  ggtitle('\n') +
  guides(fill = guide_legend(title = NULL)) +
  my_theme

png(paste(out.dir, 'hist_med_all_3_R3_', model.params$R3, '.png'))
med.dens.plot
dev.off()

# Save image
save.image(paste0(out.dir, 'image_8x8x8_', model.params$R3, '.Rdata'))

######
# Summary: tensor normalization doesn't change the distribution of the sd across replicates,
# but it does shift the medians to zero
# Seems to smooth out a lot of unwanted spatial effects
# Needs to be tested with different sized cores.
######

## HERE ##

# TODO: look at tensor of LowessRUV values to see if the same patterns are present

ss3 <- tbl_df(read.table(MCF10A_SS3_Level3@filePath, sep='\t', header=T))

ss3.nuclei <- ss3 %>% select(Well, Spot, ArrayRow, ArrayColumn, ECMp, Ligand, MEP, 
                             BWL, k, matches('Nuclei_PA'))


# HMEC1 data
ss1 <- read.table(MCF10A_SS1_Level3@filePath, sep='\t', header=T)
ss1 <- tbl_df(ss1)

ss1.edu <- ss1 %>% select(Well, Spot, ArrayRow, ArrayColumn, ECMp, Ligand, MEP, 
                          BWL, k, matches('Nuclei_PA'))
ss1.edu
