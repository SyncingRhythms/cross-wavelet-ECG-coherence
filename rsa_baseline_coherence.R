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

for (i in 1:length(dat.id.list)){
  sub_num <- dat.id.list[i]
  dat.id <- dat[ which(dat$id==dat.id.list[i]), ]
  
  # parent = data.frame(dat.id$epoch,dat.id$rsa_p_ms)
  # youth = data.frame(dat.id$epoch,dat.id$rsa_y_ms)
  
  # get parent and youth rsa timeseries
  x <- dat.id$rsa_p_b # parent
  y <- dat.id$rsa_y_b # youth
  # asign to dataframe for input to coherence computation
  wave.data <- data.frame(x = x, y = y)
  
  AR <- list(AR = list(p = 1))
  # Compute Coherence -
  wave.wc <- analyze.coherency(wave.data, my.pair = c("x","y"),
                               loess.span = 0,
                               dt = 1/3, dj = 1/5, 
                               window.type.t = 1, window.type.s = 1,
                               window.size.t = 5, window.size.s = 1,
                               method = "AR",
                               params = AR,
                               lowerPeriod = 1/2,
                               upperPeriod = 64,
                               make.pval = TRUE,
                               n.sim = 250)
  
  # save participant's Coherence data object
  saveRDS(wave.wc, paste0(paste(work_dir, sep = "/"), "/data_dt-1div3_dj5/rsa_baseline_coh_s",
                          paste(sub_num, sep = "_"), "_method-AR_sim-250_dj-1div25.rds"))
}