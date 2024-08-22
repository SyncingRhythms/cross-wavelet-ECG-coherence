library(WaveletComp)

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

num.plots <- length(dat.id.list)
my.plots <- vector(num.plots, mode='list')

for (i in 1:length(dat.id.list)){
  ### load participant coherence data object
  sub_num <- dat.id.list[i]
  sdat <- readRDS(paste0(paste(work_dir, sep = "/"), "/data_dt-1div3/rsa_baseline_coh_s",
                         paste(sub_num, sep = "_"), "_method-AR_sim-250_dj-1div25.rds"))
  
  wc.image(sdat, which.image = "wc", main = paste0("wavelet coherence, Parent over Youth - sub: ", paste(sub_num, sep = " ")),
           legend.params = list(label.digits = 3),
           periodlab = "period",
           timelab = "time")
  my.plots[[i]] <- recordPlot()
}
graphics.off()

pdf('coh_sub_plots_baseline_dt-1div3_method-AR_sim-250_dj-1div25.pdf', onefile=TRUE)
for (my.plot in my.plots) {
  replayPlot(my.plot)
}
graphics.off()
