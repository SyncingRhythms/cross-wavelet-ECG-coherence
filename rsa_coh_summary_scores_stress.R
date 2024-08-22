library(WaveletComp)
library(dplyr)
library(reshape2)

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
dat=read.csv('ydi_rsa_stress_clean.csv')
# list unique participant ids
dat.id.list<-unique(dat$id)
remove(dat)

dat.all=c()

for (i in 1:length(dat.id.list)){
  ### load participant coherence data object
  sub_num <- dat.id.list[i]
  sdat <- readRDS(paste0(paste(work_dir, sep = "/"), "/data_dt-1div3/rsa_stress_coh_s",
                         paste(sub_num, sep = "_"), "_method-AR_sim-250_dj-1div25.rds"))
  
  coh <- data.frame(sdat$Coherence)
  pval <- data.frame(sdat$Coherence.pval)
  pval.thresh <- 0.05
  mask <- pval < pval.thresh
  
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
  
  yaxis <- sdat$axis.2
  # Iterate over each time point and NA outside of cone of influence
  for (t in 1:ncol(coh)){
    if (coi.2p[t] <= 6) {
      # find index in the scale that is closest to the current COI value...
      # ie., find the value with the smallest absolute difference than coi val
      yidx = which.min(abs(yaxis - coi.2p[t])) 
      #if (yidx <= 25){
      print(yidx)
      #}
      coh[yidx:dim(coh)[1], t] <- NA
    }
  }
  
  # average +- .2 around each period of interest, eg., 3.8 - 4.2. 
  # increasing at higher periods, due to freq smearing in wavelets at high freqs
  #get coherence at timescales
  coh.2=colMeans(coh[1:15,], na.rm = T) # 
  coh.4=colMeans(coh[94:108,], na.rm = T) # 3.8-4.2
  coh.8=colMeans(coh[194:208,], na.rm = T) # 7.6-8.4
  coh.16=colMeans(coh[295:306,], na.rm = T) # 15.4-16.6
  coh.32=colMeans(coh[397:401,], na.rm = T) # 31.2-32
  
  #get mean coherence at timescales (entire phase)
  m.coh.2=mean(colMeans(coh[1:15,], na.rm = T))
  m.coh.4=mean(colMeans(coh[94:108,], na.rm = T))
  m.coh.8=mean(colMeans(coh[194:208,], na.rm = T))
  m.coh.16=mean(colMeans(coh[295:306,], na.rm = T))
  m.coh.32=mean(colMeans(coh[397:401,], na.rm = T))
  
  #get mean coherence at timescales (each minute)
  
  m.coh.2.1m=mean(colMeans(coh[1:15,1:60], na.rm = T))
  m.coh.2.2m=mean(colMeans(coh[1:15,61:120], na.rm = T))
  m.coh.2.3m=mean(colMeans(coh[1:15,121:180], na.rm = T))
  m.coh.2.4m=mean(colMeans(coh[1:15,181:240], na.rm = T))
  m.coh.2.5m=mean(colMeans(coh[1:15,240:length(coh[2,])], na.rm = T))
  
  m.coh.4.1m=mean(colMeans(coh[94:108,1:60], na.rm = T))
  m.coh.4.2m=mean(colMeans(coh[94:108,61:120], na.rm = T))
  m.coh.4.3m=mean(colMeans(coh[94:108,121:180], na.rm = T))
  m.coh.4.4m=mean(colMeans(coh[94:108,181:240], na.rm = T))
  m.coh.4.5m=mean(colMeans(coh[94:108,240:length(coh[2,])], na.rm = T))
  
  m.coh.8.1m=mean(colMeans(coh[194:208,1:60], na.rm = T))
  m.coh.8.2m=mean(colMeans(coh[194:208,61:120], na.rm = T))
  m.coh.8.3m=mean(colMeans(coh[194:208,121:180], na.rm = T))
  m.coh.8.4m=mean(colMeans(coh[194:208,181:240], na.rm = T))
  m.coh.8.5m=mean(colMeans(coh[194:208,240:length(coh[2,])], na.rm = T))
  
  m.coh.16.1m=mean(colMeans(coh[295:306,1:60], na.rm = T))
  m.coh.16.2m=mean(colMeans(coh[295:306,61:120], na.rm = T))
  m.coh.16.3m=mean(colMeans(coh[295:306,121:180], na.rm = T))
  m.coh.16.4m=mean(colMeans(coh[295:306,181:240], na.rm = T))
  m.coh.16.5m=mean(colMeans(coh[295:306,240:length(coh[2,])], na.rm = T))
  
  m.coh.32.1m=mean(colMeans(coh[397:401,1:60], na.rm = T))
  m.coh.32.2m=mean(colMeans(coh[397:401,61:120], na.rm = T))
  m.coh.32.3m=mean(colMeans(coh[397:401,121:180], na.rm = T))
  m.coh.32.4m=mean(colMeans(coh[397:401,181:240], na.rm = T))
  m.coh.32.5m=mean(colMeans(coh[397:401,240:length(coh[2,])], na.rm = T))
  
  
  #bring it all together
  dat.d=data.frame(dat.id.list[i],
                   m.coh.2,
                   m.coh.4,
                   m.coh.8,
                   m.coh.16,
                   m.coh.32, 
                   m.coh.2.1m, m.coh.2.2m, m.coh.2.3m, m.coh.2.4m,m.coh.2.5m,
                   m.coh.4.1m, m.coh.4.2m, m.coh.4.3m, m.coh.4.4m,m.coh.4.5m,
                   m.coh.8.1m, m.coh.8.2m, m.coh.8.3m, m.coh.8.4m,m.coh.8.5m,
                   m.coh.16.1m, m.coh.16.2m, m.coh.16.3m, m.coh.16.4m,m.coh.16.5m,
                   m.coh.32.1m, m.coh.32.2m, m.coh.32.3m, m.coh.32.4m,m.coh.32.5m)
  
  colnames(dat.d)=c('dyad',
                    'coh2s','coh4s','coh8s','coh16s','coh32s',
                    'm.coh.2.1m', 'm.coh.2.2m', 'm.coh.2.3m', 'm.coh.2.4m', 'm.coh.2.5m',
                    'm.coh.4.1m', 'm.coh.4.2m', 'm.coh.4.3m', 'm.coh.4.4m', 'm.coh.4.5m',
                    'm.coh.8.1m', 'm.coh.8.2m', 'm.coh.8.3m', 'm.coh.8.4m', 'm.coh.8.5m',
                    'm.coh.16.1m', 'm.coh.16.2m', 'm.coh.16.3m', 'm.coh.16.4m', 'm.coh.16.5m',
                    'm.coh.32.1m', 'm.coh.32.2m', 'm.coh.32.3m', 'm.coh.32.4m', 'm.coh.32.5m'
  )
  
  dat.all=rbind(dat.all,dat.d)
  
}

#write csv
write.csv(dat.all,'rsa_coh_p2-64_scores_stress.csv', row.names=FALSE)

###############debug
# image(t(as.matrix(coh)))
