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
dat=read.csv('ydi_rsa_baseline_clean.csv')
# list unique participant ids
dat.id.list<-unique(dat$id)
remove(dat)

dat.all=c()

for (i in 1:length(dat.id.list)){
  ### load participant coherence data object
  sub_num <- dat.id.list[i]
  sdat <- readRDS(paste0(paste(work_dir, sep = "/"), "/data_dt-1div3/rsa_baseline_coh_s",
                         paste(sub_num, sep = "_"), "_method-AR_sim-250_dj-1div25.rds"))
  
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
  
  yaxis <- sdat$axis.2
  # Iterate over each time point and NA outside of cone of influence
  for (t in 1:ncol(coh)){
    if (coi.2p[t] <= 6) {
      # find index in the scale that is closest to the current COI value...
      # ie., find the value with the smallest absolute difference than coi val
      yidx = which.min(abs(yaxis - coi.2p[t]))
      coh[yidx:dim(coh)[1], t] <- NA
    }
  }
  
  # average +- .2 around each period of interest, eg., 3.8 - 4.2. 
  # increasing at higher periods, due to freq smearing in wavelets at high freqs
  #get coherence at timescales
  coh.f=colMeans(coh[1:16,], na.rm = T) # 0.5 - .75
  coh.1=colMeans(coh[18:32,], na.rm = T) # 0.8 -1.2
  coh.2=colMeans(coh[47:55,], na.rm = T) # 1.8 = 2.2
  coh.4=colMeans(coh[74:78,], na.rm = T) # 3.8-4.2
  coh.8=colMeans(coh[99:103,], na.rm = T) # 7.6-8.4
  coh.16=colMeans(coh[124:128,], na.rm = T) # 15.2-16.8
  coh.32=colMeans(coh[146:155,], na.rm = T) # 28-36
  coh.64=colMeans(coh[167:175,], na.rm = T) # 50-62
  
  #get mean coherence at timescales (entire phase)
  m.coh.f=mean(colMeans(coh[1:16,], na.rm = T), na.rm=T)
  m.coh.1=mean(colMeans(coh[18:32,], na.rm = T), na.rm=T)
  m.coh.2=mean(colMeans(coh[47:55,], na.rm = T), na.rm=T)
  m.coh.4=mean(colMeans(coh[74:78,], na.rm = T), na.rm=T)
  m.coh.8=mean(colMeans(coh[99:103,], na.rm = T), na.rm=T)
  m.coh.16=mean(colMeans(coh[124:128,], na.rm = T), na.rm=T)
  m.coh.32=mean(colMeans(coh[146:155,], na.rm = T), na.rm=T)
  m.coh.64=mean(colMeans(coh[167:175,], na.rm = T), na.rm=T)
  
  #get mean coherence at timescales (each minute)
  
  m.coh.f.1hm=mean(colMeans(coh[1:16,1:30], na.rm = T), na.rm=T)
  m.coh.f.2hm=mean(colMeans(coh[1:16,31:60], na.rm = T), na.rm=T)
  m.coh.f.3hm=mean(colMeans(coh[1:16,61:90], na.rm = T), na.rm=T)
  m.coh.f.4hm=mean(colMeans(coh[1:16,91:120], na.rm = T), na.rm=T)
  m.coh.f.5hm=mean(colMeans(coh[1:16,121:length(coh[2,])], na.rm = T), na.rm=T)
  
  m.coh.1.1hm=mean(colMeans(coh[18:32,1:30], na.rm = T), na.rm=T)
  m.coh.1.2hm=mean(colMeans(coh[18:32,31:60], na.rm = T), na.rm=T)
  m.coh.1.3hm=mean(colMeans(coh[18:32,61:90], na.rm = T), na.rm=T)
  m.coh.1.4hm=mean(colMeans(coh[18:32,91:120], na.rm = T), na.rm=T)
  m.coh.1.5hm=mean(colMeans(coh[18:32,121:length(coh[2,])], na.rm = T), na.rm=T)
  
  m.coh.2.1hm=mean(colMeans(coh[47:55,1:30], na.rm = T), na.rm=T)
  m.coh.2.2hm=mean(colMeans(coh[47:55,31:60], na.rm = T), na.rm=T)
  m.coh.2.3hm=mean(colMeans(coh[47:55,61:90], na.rm = T), na.rm=T)
  m.coh.2.4hm=mean(colMeans(coh[47:55,91:120], na.rm = T), na.rm=T)
  m.coh.2.5hm=mean(colMeans(coh[47:55,121:length(coh[2,])], na.rm = T), na.rm=T)
  
  m.coh.4.1hm=mean(colMeans(coh[74:78,1:30], na.rm = T), na.rm=T)
  m.coh.4.2hm=mean(colMeans(coh[74:78,31:60], na.rm = T), na.rm=T)
  m.coh.4.3hm=mean(colMeans(coh[74:78,61:90], na.rm = T), na.rm=T)
  m.coh.4.4hm=mean(colMeans(coh[74:78,91:120], na.rm = T), na.rm=T)
  m.coh.4.5hm=mean(colMeans(coh[74:78,121:length(coh[2,])], na.rm = T), na.rm=T)
  
  m.coh.8.1hm=mean(colMeans(coh[99:103,1:30], na.rm = T), na.rm=T)
  m.coh.8.2hm=mean(colMeans(coh[99:103,31:60], na.rm = T), na.rm=T)
  m.coh.8.3hm=mean(colMeans(coh[99:103,61:90], na.rm = T), na.rm=T)
  m.coh.8.4hm=mean(colMeans(coh[99:103,91:120], na.rm = T), na.rm=T)
  m.coh.8.5hm=mean(colMeans(coh[99:103,121:length(coh[2,])], na.rm = T), na.rm=T)
  
  m.coh.16.1hm=mean(colMeans(coh[124:128,1:30], na.rm = T), na.rm=T)
  m.coh.16.2hm=mean(colMeans(coh[124:128,31:60], na.rm = T), na.rm=T)
  m.coh.16.3hm=mean(colMeans(coh[124:128,61:90], na.rm = T), na.rm=T)
  m.coh.16.4hm=mean(colMeans(coh[124:128,91:120], na.rm = T), na.rm=T)
  m.coh.16.5hm=mean(colMeans(coh[124:128,121:length(coh[2,])], na.rm = T), na.rm=T)
  
  m.coh.32.1hm=mean(colMeans(coh[146:155,1:30], na.rm = T), na.rm=T)
  m.coh.32.2hm=mean(colMeans(coh[146:155,31:60], na.rm = T), na.rm=T)
  m.coh.32.3hm=mean(colMeans(coh[146:155,61:90], na.rm = T), na.rm=T)
  m.coh.32.4hm=mean(colMeans(coh[146:155,91:120], na.rm = T), na.rm=T)
  m.coh.32.5hm=mean(colMeans(coh[146:155,121:length(coh[2,])], na.rm = T), na.rm=T)
  
  m.coh.64.1hm=mean(colMeans(coh[167:175,1:30], na.rm = T), na.rm=T)
  m.coh.64.2hm=mean(colMeans(coh[167:175,31:60], na.rm = T), na.rm=T)
  m.coh.64.3hm=mean(colMeans(coh[167:175,61:90], na.rm = T), na.rm=T)
  m.coh.64.4hm=mean(colMeans(coh[167:175,91:120], na.rm = T), na.rm=T)
  m.coh.64.5hm=mean(colMeans(coh[167:175,121:length(coh[2,])], na.rm = T), na.rm=T)
  
  
  #bring it all together
  dat.d=data.frame(dat.id.list[i],
                   m.coh.f,
                   m.coh.1,
                   m.coh.2,
                   m.coh.4,
                   m.coh.8,
                   m.coh.16,
                   m.coh.32, 
                   m.coh.64, 
                   m.coh.f.1hm, m.coh.f.2hm, m.coh.f.3hm, m.coh.f.4hm,m.coh.f.5hm,
                   m.coh.1.1hm, m.coh.1.2hm, m.coh.1.3hm, m.coh.1.4hm,m.coh.1.5hm,
                   m.coh.2.1hm, m.coh.2.2hm, m.coh.2.3hm, m.coh.2.4hm,m.coh.2.5hm,
                   m.coh.4.1hm, m.coh.4.2hm, m.coh.4.3hm, m.coh.4.4hm,m.coh.4.5hm,
                   m.coh.8.1hm, m.coh.8.2hm, m.coh.8.3hm, m.coh.8.4hm,m.coh.8.5hm,
                   m.coh.16.1hm, m.coh.16.2hm, m.coh.16.3hm, m.coh.16.4hm,m.coh.16.5hm,
                   m.coh.32.1hm, m.coh.32.2hm, m.coh.32.3hm, m.coh.32.4hm,m.coh.32.5hm,
                   m.coh.64.1hm, m.coh.64.2hm, m.coh.64.3hm, m.coh.64.4hm,m.coh.64.5hm)
  
  colnames(dat.d)=c('dyad',
                    'cohfs','coh1s','coh2s','coh4s','coh8s','coh16s','coh32s','coh64s',
                    'm.coh.f.1hm', 'm.coh.f.2hm', 'm.coh.f.3hm', 'm.coh.f.4hm', 'm.coh.f.5hm',
                    'm.coh.1.1hm', 'm.coh.1.2hm', 'm.coh.1.3hm', 'm.coh.1.4hm', 'm.coh.1.5hm',
                    'm.coh.2.1hm', 'm.coh.2.2hm', 'm.coh.2.3hm', 'm.coh.2.4hm', 'm.coh.2.5hm',
                    'm.coh.4.1hm', 'm.coh.4.2hm', 'm.coh.4.3hm', 'm.coh.4.4hm', 'm.coh.4.5hm',
                    'm.coh.8.1hm', 'm.coh.8.2hm', 'm.coh.8.3hm', 'm.coh.8.4hm', 'm.coh.8.5hm',
                    'm.coh.16.1hm', 'm.coh.16.2hm', 'm.coh.16.3hm', 'm.coh.16.4hm', 'm.coh.16.5hm',
                    'm.coh.32.1hm', 'm.coh.32.2hm', 'm.coh.32.3hm', 'm.coh.32.4hm', 'm.coh.32.5hm',
                    'm.coh.64.1hm', 'm.coh.64.2hm', 'm.coh.64.3hm', 'm.coh.64.4hm', 'm.coh.64.5hm'
  )
  
  dat.all=rbind(dat.all,dat.d)
  
}

#write csv
write.csv(dat.all,'rsa_coh_p3_scores_baseline.csv', row.names=FALSE)