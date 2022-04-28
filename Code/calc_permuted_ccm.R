### Cross Map System State with FD ###

## Preamble ##
require(tidyverse) # dplyr, ggplot etc.
require(pbmcapply) # paralled lapply 
require(data.table) # rbindlist function
source("Code/ccm_perm_fn.R") # custom function

###########################################################################
## Read in Data ##
###########################################################################
phyto.kin.fuzFDs.mth <- read.csv("Data/raw_FD/FD_kin_phyto_mth_raw.csv")
phyto.LZ.fuzFDs.mth <- read.csv("Data/raw_FD/FD_LZ_phyto_mth_raw.csv")
phyto.mad.fuzFDs.mth <- read.csv("Data/raw_FD/FD_mad_phyto_mth_raw.csv")
phyto.wind.fuzFDs.mth <- read.csv("Data/raw_FD/FD_wind_phyto_mth_raw.csv")
phyto.kas.fuzFDs.mth <- read.csv("Data/raw_FD/FD_kas_phyto_mth_raw.csv")

zoo.kin.fuzFDs.mth <- read.csv("Data/raw_FD/FD_kin_zoo_mth_raw.csv")
zoo.LZ.fuzFDs.mth <- read.csv("Data/raw_FD/FD_LZ_zoo_mth_raw.csv")
zoo.mad.fuzFDs.mth <- read.csv("Data/raw_FD/FD_mad_zoo_mth_raw.csv")
zoo.wind.fuzFDs.mth <- read.csv("Data/raw_FD/FD_wind_zoo_mth_raw.csv")
zoo.kas.fuzFDs.mth <- read.csv("Data/raw_FD/FD_kas_zoo_mth_raw.csv")

load("Data/all.system.states.RData")

###########################################################################
## Raw FD Prep ##
###########################################################################

kin.tot <- cbind(phyto.kin.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$kin.mth[,-c(9)])%>%
  mutate(zooFDis =  zoo.kin.fuzFDs.mth[,"FDis"],zooFEve = zoo.kin.fuzFDs.mth[,"FEve"],zooFRic = zoo.kin.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>% # log density, mvi and zo.ration to linearise
  mutate(across(-c(date,data.source,res),~scale(.x))) # center and scale to unit variance for plotting

mad.tot <- cbind(phyto.mad.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$mad.mth[,-c(9)])%>%
  mutate(zooFDis =  zoo.mad.fuzFDs.mth[,"FDis"],zooFEve = zoo.mad.fuzFDs.mth[,"FEve"],zooFRic = zoo.mad.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

LZ.tot <- cbind(phyto.LZ.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$LZ.mth[,-c(9)])%>%
  mutate(zooFDis =  zoo.LZ.fuzFDs.mth[,"FDis"],zooFEve = zoo.LZ.fuzFDs.mth[,"FEve"],zooFRic = zoo.LZ.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

wind.tot <- cbind(phyto.wind.fuzFDs.mth[c("FDis","FEve","FRic")],all.system.states$wind.mth[,-c(9)])%>%
  mutate(zooFDis =  zoo.wind.fuzFDs.mth[,"FDis"],zooFEve = zoo.wind.fuzFDs.mth[,"FEve"],zooFRic = zoo.wind.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

kas.tot <- cbind(phyto.kas.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$kas.mth[,-c(9)])%>%
  mutate(zooFDis =  zoo.kas.fuzFDs.mth[,"FDis"],zooFEve = zoo.kas.fuzFDs.mth[,"FEve"],zooFRic = zoo.kas.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

###########################################################################
## Estimate convergent cross map and permute ##
###########################################################################

## Kinneret CCM ##

kin.phytomth.ccm<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  bio <- suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"density")],
                                       iter = 10000,span =12*5,return.raw = T,
                                       detrend.method = "lm"))
  pc <-  suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"community")],
                                    iter = 10000,span =12*5,return.raw = T,
                                    detrend.method = "lm"))
  fi <-  suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"FI")],
                                   iter = 10000,span =12*5,return.raw = T,
                                   detrend.method = "lm"))
  mvi <-  suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"mvi")],
                                    iter = 10000,span =12*5,return.raw = T,
                                    detrend.method = "lm"))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"zp.ratio")],
                                         iter = 10000,span =12*5,return.raw = T,
                                         detrend.method = "lm"))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  # extract the observed correlation coefs for FD vs each system state
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  # extract all permuted correlation coefs for FD vs each system state
  out.raw.obs <- data.frame(rbind(pc$raw.obs,bio$raw.obs,fi$raw.obs,mvi$raw.obs,zp.ratio$raw.obs),
                            "state.metric" = c(rep("Community",nrow(pc$raw.obs)),rep("Density",nrow(bio$raw.obs)),rep("FI",nrow(fi$raw.obs)),rep("MVI",nrow(mvi$raw.obs)),rep("Z_P.ratio",nrow(zp.ratio$raw.obs))))
  # extract raw cross skill with lags for FD vs each system state
  out <- list("summary" = out.val,"perm.dens" = out.dens,"raw.obs" = out.raw.obs) 
  return(out)
},mc.cores = 1) 
# only single core specified as 'ccm.perm' already paralled within the function. 
# Provides opportunity for further parallelisation if desired
names(kin.phytomth.ccm) <- c("FDis","FEve","FRic") # name list elements
save(kin.phytomth.ccm,file="kin.phytomth.ccm.RData")
kin.phytomth.ccm.summary <- lapply(kin.phytomth.ccm, `[[`, 'summary')%>% # extract second level list elements (i.e. 'summary')
  data.table::rbindlist(idcol = "FD.metric") %>% # rbind the list with FD.metric id column
  mutate(system = "Kinneret", res = "Month",troph = "Phytoplankton") %>% # specify metadata for future plotting
  mutate(x_y.sig = ifelse(is.na(x_y.obs_value) | x_y.obs_value == "NaN","",
                          ifelse(x_y.quantile >= 0.95 ,"*","")))%>%
  mutate(y_x.sig = ifelse(is.na(y_x.obs_value) | y_x.obs_value == "NaN","",
                          ifelse(y_x.quantile >= 0.95 ,"*","")))# assess significance of observed cross correlation by comparing to 2.5 and 97.5 quartiles (two tailed)
kin.phytomth.ccm.raw <- lapply(kin.phytomth.ccm, `[[`, 'perm.dens')%>% # extract second level list elements (i.e. 'perm.dens')
  data.table::rbindlist(idcol = "FD.metric")%>% # rbind the list with FD.metric id column
  mutate(system = "Kinneret", res = "Month",troph = "Phytoplankton")# specify metadata for future plotting
kin.phytomth.ccm.lag <- lapply(kin.phytomth.ccm, `[[`, 'raw.obs')%>% # extract second level list elements (i.e. 'perm.dens')
  data.table::rbindlist(idcol = "FD.metric")%>% # rbind the list with FD.metric id column
  mutate(system = "Kinneret", res = "Month",troph = "Phytoplankton")# specify metadata for future plotting

kin.zoomth.ccm<- pbmcapply::pbmclapply(c("zooFDis","zooFEve","zooFRic"),function(x){
  bio <- suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"density")],
                                   detrend.method = "lm",
                                  iter = 10000,span =12*5,return.raw = T))
  pc <-  suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"community")],
                                   detrend.method = "lm",
                                    iter = 10000,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"FI")],
                                   detrend.method = "lm",
                                   iter = 10000,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"mvi")],
                                    detrend.method = "lm",
                                    iter = 10000,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"zp.ratio")],
                                         detrend.method = "lm",
                                         iter = 10000,span =12*5,return.raw = T))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out.raw.obs <- data.frame(rbind(pc$raw.obs,bio$raw.obs,fi$raw.obs,mvi$raw.obs,zp.ratio$raw.obs),
                            "state.metric" = c(rep("Community",nrow(pc$raw.obs)),rep("Density",nrow(bio$raw.obs)),rep("FI",nrow(fi$raw.obs)),rep("MVI",nrow(mvi$raw.obs)),rep("Z_P.ratio",nrow(zp.ratio$raw.obs))))
  # extract raw cross skill with lags for FD vs each system state
  out <- list("summary" = out.val,"perm.dens" = out.dens,"raw.obs" = out.raw.obs) 
  return(out)
},mc.cores = 1) 
names(kin.zoomth.ccm) <- c("FDis","FEve","FRic") # name list elements
kin.zoomth.ccm.summary <- lapply(kin.zoomth.ccm, `[[`, 'summary')%>% # extract second level list elements (i.e. 'summary')
  data.table::rbindlist(idcol = "FD.metric") %>% # rbind the list with FD.metric id column
  mutate(system = "Kinneret", res = "Month",troph = "Zooplankton") %>% # specify metadata for future plotting
  mutate(x_y.sig = ifelse(is.na(x_y.obs_value) | x_y.obs_value == "NaN","",
                          ifelse(x_y.quantile >= 0.95 ,"*","")))%>%
  mutate(y_x.sig = ifelse(is.na(y_x.obs_value) | y_x.obs_value == "NaN","",
                          ifelse(y_x.quantile >= 0.95 ,"*","")))# assess significance of observed cross correlation by comparing to 2.5 and 97.5 quartiles (two tailed)
kin.zoomth.ccm.raw <- lapply(kin.zoomth.ccm, `[[`, 'perm.dens')%>% # extract second level list elements (i.e. 'perm.dens')
  data.table::rbindlist(idcol = "FD.metric")%>% # rbind the list with FD.metric id column
  mutate(system = "Kinneret", res = "Month",troph = "Zooplankton") # specify metadata for future plotting
kin.zoomth.ccm.lag <- lapply(kin.zoomth.ccm, `[[`, 'raw.obs')%>% # extract second level list elements (i.e. 'perm.dens')
  data.table::rbindlist(idcol = "FD.metric")%>% # rbind the list with FD.metric id column
  mutate(system = "Kinneret", res = "Month",troph = "Zooplankton")# specify metadata for future plotting

kin.summary.ccm <- rbind(kin.phytomth.ccm.summary,kin.zoomth.ccm.summary)
kin.raw.ccm <- rbind(kin.phytomth.ccm.raw,kin.zoomth.ccm.raw)
kin.lag.ccm <- rbind(kin.phytomth.ccm.lag,kin.zoomth.ccm.lag)

write.csv(kin.summary.ccm,file ="Results/ccm/raw_data/kin_ccm_summary.csv",row.names = F)
kin.summary.ccm <- read.csv(file ="Results/ccm/raw_data/kin_ccm_summary.csv")
write.csv(kin.raw.ccm,file ="Results/ccm/raw_data/kin_ccm_raw.csv",row.names = F)
kin.raw.ccm <- read.csv(file ="Results/ccm/raw_data/kin_ccm_raw.csv")
write.csv(kin.lag.ccm,file ="Results/ccm/raw_data/kin_ccm_lag.csv",row.names = F)
kin.lag.ccm <- read.csv(file ="Results/ccm/raw_data/kin_ccm_lag.csv")

## Mendota CCM ##

mad.phytomth.ccm<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  bio <- suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"density")],
                                   detrend.method = "lm",
                                  iter = 10000,span =12*5,return.raw = T))
  pc <-  suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"community")],
                                   detrend.method = "lm",
                                    iter = 10000,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"FI")],
                                   detrend.method = "lm",
                                   iter = 10000,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"mvi")],
                                    detrend.method = "lm",
                                    iter = 10000,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"zp.ratio")],
                                         detrend.method = "lm",
                                         iter = 10000,span =12*5,return.raw = T))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out.raw.obs <- data.frame(rbind(pc$raw.obs,bio$raw.obs,fi$raw.obs,mvi$raw.obs,zp.ratio$raw.obs),
                            "state.metric" = c(rep("Community",nrow(pc$raw.obs)),rep("Density",nrow(bio$raw.obs)),rep("FI",nrow(fi$raw.obs)),rep("MVI",nrow(mvi$raw.obs)),rep("Z_P.ratio",nrow(zp.ratio$raw.obs))))
  out <- list("summary" = out.val,"perm.dens" = out.dens,"raw.obs" = out.raw.obs) 
  return(out)
},mc.cores = 1) 
names(mad.phytomth.ccm) <- c("FDis","FEve","FRic")
mad.phytomth.ccm.summary <- lapply(mad.phytomth.ccm, `[[`, 'summary')%>% 
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Mendota", res = "Month",troph = "Phytoplankton") %>% 
  mutate(x_y.sig = ifelse(is.na(x_y.obs_value) | x_y.obs_value == "NaN","",
                          ifelse(x_y.quantile >= 0.95,"*","")))%>%
  mutate(y_x.sig = ifelse(is.na(y_x.obs_value) | y_x.obs_value == "NaN","",
                          ifelse(y_x.quantile >= 0.95 | y_x.quantile <= 0.025,"*","")))
mad.phytomth.ccm.raw <- lapply(mad.phytomth.ccm, `[[`, 'perm.dens')%>% 
  data.table::rbindlist(idcol = "FD.metric")%>% 
  mutate(system = "Mendota", res = "Month",troph = "Phytoplankton")
mad.phytomth.ccm.lag <- lapply(mad.phytomth.ccm, `[[`, 'raw.obs')%>% 
  data.table::rbindlist(idcol = "FD.metric")%>% 
  mutate(system = "Mendota", res = "Month",troph = "Phytoplankton")

mad.zoomth.ccm<- pbmcapply::pbmclapply(c("zooFDis","zooFEve","zooFRic"),function(x){
  bio <- suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"density")],
                                   detrend.method = "lm",
                                  iter = 10000,span =12*5,return.raw = T))
  pc <-  suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"community")],
                                   detrend.method = "lm",
                                    iter = 10000,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"FI")],
                                   detrend.method = "lm",
                                   iter = 10000,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"mvi")],
                                    detrend.method = "lm",
                                    iter = 10000,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"zp.ratio")],
                                         detrend.method = "lm",
                                         iter = 10000,span =12*5,return.raw = T))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out.raw.obs <- data.frame(rbind(pc$raw.obs,bio$raw.obs,fi$raw.obs,mvi$raw.obs,zp.ratio$raw.obs),
                            "state.metric" = c(rep("Community",nrow(pc$raw.obs)),rep("Density",nrow(bio$raw.obs)),rep("FI",nrow(fi$raw.obs)),rep("MVI",nrow(mvi$raw.obs)),rep("Z_P.ratio",nrow(zp.ratio$raw.obs))))
  out <- list("summary" = out.val,"perm.dens" = out.dens,"raw.obs" = out.raw.obs) 
  return(out)
},mc.cores = 1) 
names(mad.zoomth.ccm) <- c("FDis","FEve","FRic") 
mad.zoomth.ccm.summary <- lapply(mad.zoomth.ccm, `[[`, 'summary')%>% 
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Mendota", res = "Month",troph = "Zooplankton") %>% 
  mutate(x_y.sig = ifelse(is.na(x_y.obs_value) | x_y.obs_value == "NaN","",
                          ifelse(x_y.quantile >= 0.95 ,"*","")))%>%
  mutate(y_x.sig = ifelse(is.na(y_x.obs_value) | y_x.obs_value == "NaN","",
                          ifelse(y_x.quantile >= 0.95,"*","")))
mad.zoomth.ccm.raw <- lapply(mad.zoomth.ccm, `[[`, 'perm.dens')%>% 
  data.table::rbindlist(idcol = "FD.metric")%>% 
  mutate(system = "Mendota", res = "Month",troph = "Zooplankton") 
mad.zoomth.ccm.lag <- lapply(mad.zoomth.ccm, `[[`, 'raw.obs')%>% 
  data.table::rbindlist(idcol = "FD.metric")%>% 
  mutate(system = "Mendota", res = "Month",troph = "Zooplankton")

mad.summary.ccm <- rbind(mad.phytomth.ccm.summary,mad.zoomth.ccm.summary)
mad.raw.ccm <- rbind(mad.phytomth.ccm.raw,mad.zoomth.ccm.raw)
mad.lag.ccm <- rbind(mad.phytomth.ccm.lag,mad.zoomth.ccm.lag)

write.csv(mad.summary.ccm,file ="Results/ccm/raw_data/mad_ccm_summary.csv",row.names = F)
mad.summary.ccm <- read.csv(file ="Results/ccm/raw_data/mad_ccm_summary.csv")
write.csv(mad.raw.ccm,file ="Results/ccm/raw_data/mad_ccm_raw.csv",row.names = F)
mad.raw.ccm <- read.csv(file ="Results/ccm/raw_data/mad_ccm_raw.csv")
write.csv(mad.lag.ccm,file ="Results/ccm/raw_data/mad_ccm_lag.csv",row.names = F)
mad.lag.ccm <- read.csv(file ="Results/ccm/raw_data/mad_ccm_lag.csv")

## Lower Zurich CCM ##

LZ.phytomth.ccm<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  bio <- suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"density")],
                                   detrend.method = "lm",
                                  iter = 10000,span =12*5,return.raw = T))
  pc <-  suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"community")],
                                   detrend.method = "lm",
                                    iter = 10000,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"FI")],
                                   detrend.method = "lm",
                                   iter = 10000,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"mvi")],
                                    detrend.method = "lm",
                                    iter = 10000,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"zp.ratio")],
                                         detrend.method = "lm",
                                         iter = 10000,span =12*5,return.raw = T))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out.raw.obs <- data.frame(rbind(pc$raw.obs,bio$raw.obs,fi$raw.obs,mvi$raw.obs,zp.ratio$raw.obs),
                            "state.metric" = c(rep("Community",nrow(pc$raw.obs)),rep("Density",nrow(bio$raw.obs)),rep("FI",nrow(fi$raw.obs)),rep("MVI",nrow(mvi$raw.obs)),rep("Z_P.ratio",nrow(zp.ratio$raw.obs))))
  out <- list("summary" = out.val,"perm.dens" = out.dens,"raw.obs" = out.raw.obs) 
  return(out)
},mc.cores = 1) 
names(LZ.phytomth.ccm) <- c("FDis","FEve","FRic")
LZ.phytomth.ccm.summary <- lapply(LZ.phytomth.ccm, `[[`, 'summary')%>% 
  data.table::rbindlist(idcol = "FD.metric")%>% 
  mutate(system = "Lower Zurich", res = "Month",troph = "Phytoplankton") %>% 
  mutate(x_y.sig = ifelse(is.na(x_y.obs_value) | x_y.obs_value == "NaN","",
                          ifelse(x_y.quantile >= 0.95 ,"*","")))%>%
  mutate(y_x.sig = ifelse(is.na(y_x.obs_value) | y_x.obs_value == "NaN","",
                          ifelse(y_x.quantile >= 0.95 ,"*","")))
LZ.phytomth.ccm.raw <- lapply(LZ.phytomth.ccm, `[[`, 'perm.dens')%>% 
  data.table::rbindlist(idcol = "FD.metric")%>% 
  mutate(system = "Lower Zurich", res = "Month",troph = "Phytoplankton")
LZ.phytomth.ccm.lag <- lapply(LZ.phytomth.ccm, `[[`, 'raw.obs')%>% 
  data.table::rbindlist(idcol = "FD.metric")%>% 
  mutate(system = "Lower Zurich", res = "Month",troph = "Phytoplankton")

LZ.zoomth.ccm<- pbmcapply::pbmclapply(c("zooFDis","zooFEve","zooFRic"),function(x){
  bio <- suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"density")],
                                   detrend.method = "lm",
                                  iter = 10000,span =12*5,return.raw = T))
  pc <-  suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"community")],
                                   detrend.method = "lm",
                                    iter = 10000,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"FI")],
                                   detrend.method = "lm",
                                   iter = 10000,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"mvi")],
                                    detrend.method = "lm",
                                    iter = 10000,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"zp.ratio")],
                                         detrend.method = "lm",
                                         iter = 10000,span =12*5,return.raw = T))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out.raw.obs <- data.frame(rbind(pc$raw.obs,bio$raw.obs,fi$raw.obs,mvi$raw.obs,zp.ratio$raw.obs),
                            "state.metric" = c(rep("Community",nrow(pc$raw.obs)),rep("Density",nrow(bio$raw.obs)),rep("FI",nrow(fi$raw.obs)),rep("MVI",nrow(mvi$raw.obs)),rep("Z_P.ratio",nrow(zp.ratio$raw.obs))))
  out <- list("summary" = out.val,"perm.dens" = out.dens,"raw.obs" = out.raw.obs) 
  return(out)
},mc.cores = 1) 
names(LZ.zoomth.ccm) <- c("FDis","FEve","FRic") 
LZ.zoomth.ccm.summary <- lapply(LZ.zoomth.ccm, `[[`, 'summary')%>% 
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Lower Zurich", res = "Month",troph = "Zooplankton") %>% 
  mutate(x_y.sig = ifelse(is.na(x_y.obs_value) | x_y.obs_value == "NaN","",
                          ifelse(x_y.quantile >= 0.95 ,"*","")))%>%
  mutate(y_x.sig = ifelse(is.na(y_x.obs_value) | y_x.obs_value == "NaN","",
                          ifelse(y_x.quantile >= 0.95 ,"*","")))
LZ.zoomth.ccm.raw <- lapply(LZ.zoomth.ccm, `[[`, 'perm.dens')%>% 
  data.table::rbindlist(idcol = "FD.metric")%>% 
  mutate(system = "Lower Zurich", res = "Month",troph = "Zooplankton") 
LZ.zoomth.ccm.lag <- lapply(LZ.zoomth.ccm, `[[`, 'raw.obs')%>% 
  data.table::rbindlist(idcol = "FD.metric")%>% 
  mutate(system = "Lower Zurich", res = "Month",troph = "Zooplankton")

LZ.summary.ccm <- rbind(LZ.phytomth.ccm.summary,LZ.zoomth.ccm.summary)
LZ.raw.ccm <- rbind(LZ.phytomth.ccm.raw,LZ.zoomth.ccm.raw)
LZ.lag.ccm <- rbind(LZ.phytomth.ccm.lag,LZ.zoomth.ccm.lag)

write.csv(LZ.summary.ccm,file ="Results/ccm/raw_data/LZ_ccm_summary.csv",row.names = F)
LZ.summary.ccm <- read.csv(file ="Results/ccm/raw_data/LZ_ccm_summary.csv")
write.csv(LZ.raw.ccm,file ="Results/ccm/raw_data/LZ_ccm_raw.csv",row.names = F)
LZ.raw.ccm <- read.csv(file ="Results/ccm/raw_data/LZ_ccm_raw.csv")
write.csv(LZ.lag.ccm,file ="Results/ccm/raw_data/LZ_ccm_lag.csv",row.names = F)
LZ.lag.ccm <- read.csv(file ="Results/ccm/raw_data/LZ_ccm_lag.csv")

## Windermere CCM ##

wind.phytomth.ccm<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  bio <- suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"density")],
                                   detrend.method = "lm",
                                  iter = 10000,span =12*5,return.raw = T))
  pc <-  suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"community")],
                                   detrend.method = "lm",
                                    iter = 10000,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"FI")],
                                   detrend.method = "lm",
                                   iter = 10000,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"mvi")],
                                    detrend.method = "lm",
                                    iter = 10000,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"zp.ratio")],
                                         detrend.method = "lm",
                                         iter = 10000,span =12*5,return.raw = T))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out.raw.obs <- data.frame(rbind(pc$raw.obs,bio$raw.obs,fi$raw.obs,mvi$raw.obs,zp.ratio$raw.obs),
                            "state.metric" = c(rep("Community",nrow(pc$raw.obs)),rep("Density",nrow(bio$raw.obs)),rep("FI",nrow(fi$raw.obs)),rep("MVI",nrow(mvi$raw.obs)),rep("Z_P.ratio",nrow(zp.ratio$raw.obs))))
  out <- list("summary" = out.val,"perm.dens" = out.dens,"raw.obs" = out.raw.obs) 
  return(out)
},mc.cores = 1) 
names(wind.phytomth.ccm) <- c("FDis","FEve","FRic")
wind.phytomth.ccm.summary <- lapply(wind.phytomth.ccm, `[[`, 'summary')%>% 
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Windermere", res = "Month",troph = "Phytoplankton") %>% 
  mutate(x_y.sig = ifelse(is.na(x_y.obs_value) | x_y.obs_value == "NaN","",
                          ifelse(x_y.quantile >= 0.95 ,"*","")))%>%
  mutate(y_x.sig = ifelse(is.na(y_x.obs_value) | y_x.obs_value == "NaN","",
                          ifelse(y_x.quantile >= 0.95 ,"*","")))
wind.phytomth.ccm.raw <- lapply(wind.phytomth.ccm, `[[`, 'perm.dens')%>% 
  data.table::rbindlist(idcol = "FD.metric")%>% 
  mutate(system = "Windermere", res = "Month",troph = "Phytoplankton")
wind.phytomth.ccm.lag <- lapply(wind.phytomth.ccm, `[[`, 'raw.obs')%>% 
  data.table::rbindlist(idcol = "FD.metric")%>% 
  mutate(system = "Windermere", res = "Month",troph = "Phytoplankton")

wind.zoomth.ccm<- pbmcapply::pbmclapply(c("zooFDis","zooFEve","zooFRic"),function(x){
  bio <- suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"density")],
                                   detrend.method = "lm",
                                  iter = 10000,span =12*5,return.raw = T))
  pc <-  suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"community")],
                                   detrend.method = "lm",
                                    iter = 10000,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"FI")],
                                   detrend.method = "lm",
                                   iter = 10000,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"mvi")],
                                    detrend.method = "lm",
                                    iter = 10000,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"zp.ratio")],
                                         detrend.method = "lm",
                                         iter = 10000,span =12*5,return.raw = T))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out.raw.obs <- data.frame(rbind(pc$raw.obs,bio$raw.obs,fi$raw.obs,mvi$raw.obs,zp.ratio$raw.obs),
                            "state.metric" = c(rep("Community",nrow(pc$raw.obs)),rep("Density",nrow(bio$raw.obs)),rep("FI",nrow(fi$raw.obs)),rep("MVI",nrow(mvi$raw.obs)),rep("Z_P.ratio",nrow(zp.ratio$raw.obs))))
  out <- list("summary" = out.val,"perm.dens" = out.dens,"raw.obs" = out.raw.obs) 
  return(out)
},mc.cores = 1) 
names(wind.zoomth.ccm) <- c("FDis","FEve","FRic") 
wind.zoomth.ccm.summary <- lapply(wind.zoomth.ccm, `[[`, 'summary')%>% 
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Windermere", res = "Month",troph = "Zooplankton") %>% 
  mutate(x_y.sig = ifelse(is.na(x_y.obs_value) | x_y.obs_value == "NaN","",
                          ifelse(x_y.quantile >= 0.95 ,"*","")))%>%
  mutate(y_x.sig = ifelse(is.na(y_x.obs_value) | y_x.obs_value == "NaN","",
                          ifelse(y_x.quantile >= 0.95 ,"*","")))
wind.zoomth.ccm.raw <- lapply(wind.zoomth.ccm, `[[`, 'perm.dens')%>% 
  data.table::rbindlist(idcol = "FD.metric")%>% 
  mutate(system = "Windermere", res = "Month",troph = "Zooplankton") 
wind.zoomth.ccm.lag <- lapply(wind.zoomth.ccm, `[[`, 'raw.obs')%>% 
  data.table::rbindlist(idcol = "FD.metric")%>% 
  mutate(system = "Windermere", res = "Month",troph = "Zooplankton")

wind.summary.ccm <- rbind(wind.phytomth.ccm.summary,wind.zoomth.ccm.summary)
wind.raw.ccm <- rbind(wind.phytomth.ccm.raw,wind.zoomth.ccm.raw)
wind.lag.ccm <- rbind(wind.phytomth.ccm.lag,wind.zoomth.ccm.lag)

write.csv(wind.summary.ccm,file ="Results/ccm/raw_data/wind_ccm_summary.csv",row.names = F)
wind.summary.ccm <- read.csv(file ="Results/ccm/raw_data/wind_ccm_summary.csv")
write.csv(wind.raw.ccm,file ="Results/ccm/raw_data/wind_ccm_raw.csv",row.names = F)
wind.raw.ccm <- read.csv(file ="Results/ccm/raw_data/wind_ccm_raw.csv")
write.csv(wind.lag.ccm,file ="Results/ccm/raw_data/wind_ccm_lag.csv",row.names = F)
wind.lag.ccm <- read.csv(file ="Results/ccm/raw_data/wind_ccm_lag.csv")

## Kasumigaura CCM ##

kas.phytomth.ccm<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  bio <- suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"density")],
                                   detrend.method = "lm",
                                  iter = 10000,span =12*5,return.raw = T))
  pc <-  suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"community")],
                                   detrend.method = "lm",
                                    iter = 10000,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"FI")],
                                   detrend.method = "lm",
                                   iter = 10000,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"mvi")],
                                    detrend.method = "lm",
                                    iter = 10000,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"zp.ratio")],
                                         detrend.method = "lm",
                                         iter = 10000,span =12*5,return.raw = T))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out.raw.obs <- data.frame(rbind(pc$raw.obs,bio$raw.obs,fi$raw.obs,mvi$raw.obs,zp.ratio$raw.obs),
                            "state.metric" = c(rep("Community",nrow(pc$raw.obs)),rep("Density",nrow(bio$raw.obs)),rep("FI",nrow(fi$raw.obs)),rep("MVI",nrow(mvi$raw.obs)),rep("Z_P.ratio",nrow(zp.ratio$raw.obs))))
  out <- list("summary" = out.val,"perm.dens" = out.dens,"raw.obs" = out.raw.obs) 
  return(out)
},mc.cores = 1) 
names(kas.phytomth.ccm) <- c("FDis","FEve","FRic")
kas.phytomth.ccm.summary <- lapply(kas.phytomth.ccm, `[[`, 'summary')%>% 
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Kasumigaura", res = "Month",troph = "Phytoplankton") %>% 
  mutate(x_y.sig = ifelse(is.na(x_y.obs_value) | x_y.obs_value == "NaN","",
                          ifelse(x_y.quantile >= 0.95,"*","")))%>%
  mutate(y_x.sig = ifelse(is.na(y_x.obs_value) | y_x.obs_value == "NaN","",
                          ifelse(y_x.quantile >= 0.95 ,"*","")))
kas.phytomth.ccm.raw <- lapply(kas.phytomth.ccm, `[[`, 'perm.dens')%>% 
  data.table::rbindlist(idcol = "FD.metric")%>% 
  mutate(system = "Kasumigaura", res = "Month",troph = "Phytoplankton")
kas.phytomth.ccm.lag <- lapply(kas.phytomth.ccm, `[[`, 'raw.obs')%>% 
  data.table::rbindlist(idcol = "FD.metric")%>% 
  mutate(system = "Kasumigaura", res = "Month",troph = "Phytoplankton")

kas.zoomth.ccm<- pbmcapply::pbmclapply(c("zooFDis","zooFEve","zooFRic"),function(x){
  bio <- suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"density")],
                                   detrend.method = "lm",
                                  iter = 10000,span =12*5,return.raw = T))
  pc <-  suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"community")],
                                   detrend.method = "lm",
                                    iter = 10000,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"FI")],
                                   detrend.method = "lm",
                                   iter = 10000,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"mvi")],
                                    detrend.method = "lm",
                                    iter = 10000,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"zp.ratio")],
                                         detrend.method = "lm",
                                         iter = 10000,span =12*5,return.raw = T))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out.raw.obs <- data.frame(rbind(pc$raw.obs,bio$raw.obs,fi$raw.obs,mvi$raw.obs,zp.ratio$raw.obs),
                            "state.metric" = c(rep("Community",nrow(pc$raw.obs)),rep("Density",nrow(bio$raw.obs)),rep("FI",nrow(fi$raw.obs)),rep("MVI",nrow(mvi$raw.obs)),rep("Z_P.ratio",nrow(zp.ratio$raw.obs))))
  out <- list("summary" = out.val,"perm.dens" = out.dens,"raw.obs" = out.raw.obs) 
  return(out)
},mc.cores = 1) 
names(kas.zoomth.ccm) <- c("FDis","FEve","FRic") 
kas.zoomth.ccm.summary <- lapply(kas.zoomth.ccm, `[[`, 'summary')%>% 
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Kasumigaura", res = "Month",troph = "Zooplankton") %>% 
  mutate(x_y.sig = ifelse(is.na(x_y.obs_value) | x_y.obs_value == "NaN","",
                          ifelse(x_y.quantile >= 0.95 ,"*","")))%>%
  mutate(y_x.sig = ifelse(is.na(y_x.obs_value) | y_x.obs_value == "NaN","",
                          ifelse(y_x.quantile >= 0.95 ,"*","")))
kas.zoomth.ccm.raw <- lapply(kas.zoomth.ccm, `[[`, 'perm.dens')%>% 
  data.table::rbindlist(idcol = "FD.metric")%>% 
  mutate(system = "Kasumigaura", res = "Month",troph = "Zooplankton") 
kas.zoomth.ccm.lag <- lapply(kas.zoomth.ccm, `[[`, 'raw.obs')%>% 
  data.table::rbindlist(idcol = "FD.metric")%>% 
  mutate(system = "Kasumigaura", res = "Month",troph = "Zooplankton")

kas.summary.ccm <- rbind(kas.phytomth.ccm.summary,kas.zoomth.ccm.summary)
kas.raw.ccm <- rbind(kas.phytomth.ccm.raw,kas.zoomth.ccm.raw)
kas.lag.ccm <- rbind(kas.phytomth.ccm.lag,kas.zoomth.ccm.lag)

write.csv(kas.summary.ccm,file ="Results/ccm/raw_data/kas_ccm_summary.csv",row.names = F)
kas.summary.ccm <- read.csv(file ="Results/ccm/raw_data/kas_ccm_summary.csv")
write.csv(kas.raw.ccm,file ="Results/ccm/raw_data/kas_ccm_raw.csv",row.names = F)
kas.raw.ccm <- read.csv(file ="Results/ccm/raw_data/kas_ccm_raw.csv")
write.csv(kas.lag.ccm,file ="Results/ccm/raw_data/kas_ccm_lag.csv",row.names = F)
kas.lag.ccm <- read.csv(file ="Results/ccm/raw_data/kas_ccm_lag.csv")

###########################################################################
## Save out ##
###########################################################################
summary.ccm <- rbind(kin.summary.ccm,mad.summary.ccm,
                     LZ.summary.ccm,kas.summary.ccm,
                     wind.summary.ccm)
  
raw.ccm <- rbind(kin.raw.ccm,mad.raw.ccm,
                     LZ.raw.ccm,kas.raw.ccm,
                     wind.raw.ccm)

lag.ccm <- rbind(kin.lag.ccm,mad.lag.ccm,
                     LZ.lag.ccm,kas.lag.ccm,
                     wind.lag.ccm)

write.csv(summary.ccm,file ="Results/ccm/raw_data/ccm_summary.csv",row.names = F)
save(raw.ccm,file = "Results/ccm/raw_data/ccm_raw.RData") # RData required to reduce file size compared to .csv
write.csv(lag.ccm,file ="Results/ccm/raw_data/ccm_lag.csv",row.names = F)
