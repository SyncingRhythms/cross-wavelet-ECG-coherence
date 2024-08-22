library(WaveletComp)
library(dplyr)
library(reshape2)
library(abind)


##########################
##########################
##########################
# Input whether running on Linux or Mac
linux <- 1

# set working directory
if (linux == 1) {
  work_dir <- "/home/cglab/projects/dorry/cross-wavlet/"
} else {
  work_dir <- "/Users/cglab/projects/dorry/cross-wavlet/"
}

setwd(work_dir)

# load data
dat=read.csv('dorry_w1_p_FoodIS_matched_CWT.csv')
# list unique participant ids
dat.id.list<-unique(dat$ID)

dat.all=c()
### get subs CWT data, who have DERS Parent data
# append eac sub CWT matrix in to full list
matlist <- list()
for (i in 1:length(dat.id.list)){
  ### load participant coherence data object
  sub_num <- dat.id.list[i]
  sdat <- readRDS(paste0(paste(work_dir, sep = "/"), "/data_dt-1div3_dj5/rsa_stress_coh_s",
                         paste(sub_num, sep = "_"), "_method-AR_sim-250_dj-1div25.rds"))
  
  matlist[[i]] <- as.matrix(sdat$Coherence)
}
# convert matrix list to 3d matrix period x time x sub
coh_mats <- abind( matlist, along=3 )
# write.csv(coh_mats,'rsa_coh_subs_mats_stress.csv', row.names=FALSE)
# dimensions
nperiods <- dim(coh_mats)[1]
ntimepoints <- dim(coh_mats)[2]
nsubs <- dim(coh_mats)[3]

subscales <- c('Strategies', 'Nonaccpt', 'Impulse', 'Goal', 'Aware', 'Clarity')
# get scale to correlate
scale <- dat$FISPW1_T
# empty matrix (2d array) to input correlations below
corrs <-  array(0L, dim(sdat$Coherence)) 
pvals <-  array(0L, dim(sdat$Coherence)) 
for (p in 1:nperiods){
  for (t in 1:ntimepoints){
    corrs[p, t] <- cor.test(coh_mats[p, t, ], scale)$estimate
    pvals[p, t] <- cor.test(coh_mats[p, t, ], scale)$p.value
  }
}


# significant correlations
corsig <- corrs
# zero out nonsignificant (>.05) correlations
corsig[pvals>0.05] <- 0
# False Discovery Rate (FDR) adjust p-values for multiple comparisons 
padj.v <- p.adjust(pvals, method='BH')
# transfrom to matrix of appropriate size
padj <- matrix(padj.v, nrow = 176, ncol = 300)
corfdr <- corrs
corfdr[padj>0.05] <- 0
coh <- data.frame(sdat$Coherence)
# pval <- data.frame(sdat$Coherence.pval)
# pval.thresh <- 0.05
# mask <- pval < pval.thresh

#### replace regions outside of COI with NA,..
#### so that they aren't included in averages
# Extract the COI boundaries
coi.1 <- sdat$coi.1  # 
coi.2 <- sdat$coi.2  # 
timepoints <- ncol(coh)
scales <- ncol(coh)

# the coi.1 and coi.2 vectors are 3 items longer than the x-axis
# the first 3 and last 3 numbers below 0 and above the cone of influence
# so those are removed, leaving us with a vector 3 items too short
# thus, adding two 0's to the beginning and one 0 to the end
coi.1s <- coi.1[4:301] 
coi.2s <- coi.2[4:301]
# append -1 to the beginning, to match scale
pad <- c(-1)
coi.2p <- c(pad, coi.2s)
# append one zero to the end
coi.2p <- c(coi.2p, pad)
# append .5 and 299.5, to map to near-edges of timepoints
coi.1p <- c(0.5, coi.1s)
coi.1p <- c(coi.1p, 299.5)

xaxis <- sdat$axis.1
yaxis <- sdat$axis.2
# Iterate over each time point and NA outside of cone of influence
for (t in 1:ncol(coh)){
  if (coi.2p[t] <= 6) {
    # find index in the scale that is closest to the current COI value...
    # ie., find the value with the smallest absolute difference than coi val
    yidx = which.min(abs(yaxis - coi.2p[t])) 
    coh[yidx:dim(coh)[1], t] <- NA
    corrs[yidx:dim(corrs)[1], t] <- 0
    pvals[yidx:dim(pvals)[1], t] <- 1
    corfdr[yidx:dim(corfdr)[1], t] <- 0
    padj[yidx:dim(padj)[1], t] <- 1
  }
}

coi.mask <- !is.na(coh)
pvals.coi <- pvals[coi.mask]
padj.coi <- p.adjust(pvals.coi, method = 'BH')

# Standard P-values
# input corr and pvalues into WaveletComp object to take advantage of it's plotting capabilities
sdat$Coherence <- corrs
sdat$Coherence.pval <- pvals
wc.image(sdat, which.image = "wc", main = "Correlation of Coherence with Food Insecurity (n=65)",
         plot.arrow = FALSE,
         legend.params = list(label.digits = 2),
         periodlab = "period",
         timelab = "time")

# FDR-corrected p-values
# input corr and pvalues into WaveletComp object to take advantage of it's plotting capabilities
sdat$Coherence <- corrs
sdat$Coherence.pval <- padj
wc.image(sdat, which.image = "wc", main = "Correlation of Coherence with Food Insecurity (n=65)",
         plot.arrow = FALSE,
         legend.params = list(label.digits = 2),
         periodlab = "period",
         timelab = "time")
image(xaxis, sdat$Period, t(corsig))
