### Cross Map System State with FD ###

## Preamble ##
require(tidyverse) # dplyr, ggplot etc.
require(pbmcapply) # paralled lapply 
require(data.table) # rbindlist function
require(patchwork) # plot alignment

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
                                       iter = 500,span =12*5,return.raw = T,
                                       detrend.method = "none"))
  pc <-  suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"community")],
                                    iter = 500,span =12*5,return.raw = T,
                                    detrend.method = "lm"))
  fi <-  suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"FI")],
                                   iter = 500,span =12*5,return.raw = T,
                                   detrend.method = "lm"))
  mvi <-  suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"mvi")],
                                    iter = 500,span =12*5,return.raw = T,
                                    detrend.method = "lm"))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"zp.ratio")],
                                         iter = 500,span =12*5,return.raw = T,
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
kin.phytomth.ccm.summary <- lapply(kin.phytomth.ccm, `[[`, 'summary')%>% # extract second level list elements (i.e. 'summary')
  data.table::rbindlist(idcol = "FD.metric") %>% # rbind the list with FD.metric id column
  mutate(system = "Kinneret", res = "Month",troph = "Phytoplankton") %>% # specify metadata for future plotting
  mutate(x_y.sig = ifelse(is.na(x_y.obs_value) | x_y.obs_value == "NaN","",
                          ifelse(x_y.quantile >= 0.975 | x_y.quantile <= 0.025,"*","")))%>%
  mutate(y_x.sig = ifelse(is.na(y_x.obs_value) | y_x.obs_value == "NaN","",
                          ifelse(y_x.quantile >= 0.975 | y_x.quantile <= 0.025,"*","")))# assess significance of observed cross correlation by comparing to 2.5 and 97.5 quartiles (two tailed)
kin.phytomth.ccm.raw <- lapply(kin.phytomth.ccm, `[[`, 'perm.dens')%>% # extract second level list elements (i.e. 'perm.dens')
  data.table::rbindlist(idcol = "FD.metric")%>% # rbind the list with FD.metric id column
  mutate(system = "Kinneret", res = "Month",troph = "Phytoplankton")# specify metadata for future plotting
kin.phytomth.ccm.lag <- lapply(kin.phytomth.ccm, `[[`, 'raw.obs')%>% # extract second level list elements (i.e. 'perm.dens')
  data.table::rbindlist(idcol = "FD.metric")%>% # rbind the list with FD.metric id column
  mutate(system = "Kinneret", res = "Month",troph = "Phytoplankton")# specify metadata for future plotting

kin.zoomth.ccm<- pbmcapply::pbmclapply(c("zooFDis","zooFEve","zooFRic"),function(x){
  bio <- suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"density")],
                                   detrend.method = "lm",
                                  iter = 500,span =12*5,return.raw = T))
  pc <-  suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"community")],
                                   detrend.method = "lm",
                                    iter = 500,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"FI")],
                                   detrend.method = "lm",
                                   iter = 500,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"mvi")],
                                    detrend.method = "lm",
                                    iter = 500,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"zp.ratio")],
                                         detrend.method = "lm",
                                         iter = 500,span =12*5,return.raw = T))
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
                          ifelse(x_y.quantile >= 0.975 | x_y.quantile <= 0.025,"*","")))%>%
  mutate(y_x.sig = ifelse(is.na(y_x.obs_value) | y_x.obs_value == "NaN","",
                          ifelse(y_x.quantile >= 0.975 | y_x.quantile <= 0.025,"*","")))# assess significance of observed cross correlation by comparing to 2.5 and 97.5 quartiles (two tailed)
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
                                  iter = 500,span =12*5,return.raw = T))
  pc <-  suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"community")],
                                   detrend.method = "lm",
                                    iter = 500,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"FI")],
                                   detrend.method = "lm",
                                   iter = 500,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"mvi")],
                                    detrend.method = "lm",
                                    iter = 500,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"zp.ratio")],
                                         detrend.method = "lm",
                                         iter = 500,span =12*5,return.raw = T))
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
                          ifelse(x_y.quantile >= 0.975 | x_y.quantile <= 0.025,"*","")))%>%
  mutate(y_x.sig = ifelse(is.na(y_x.obs_value) | y_x.obs_value == "NaN","",
                          ifelse(y_x.quantile >= 0.975 | y_x.quantile <= 0.025,"*","")))
mad.phytomth.ccm.raw <- lapply(mad.phytomth.ccm, `[[`, 'perm.dens')%>% 
  data.table::rbindlist(idcol = "FD.metric")%>% 
  mutate(system = "Mendota", res = "Month",troph = "Phytoplankton")
mad.phytomth.ccm.lag <- lapply(mad.phytomth.ccm, `[[`, 'raw.obs')%>% 
  data.table::rbindlist(idcol = "FD.metric")%>% 
  mutate(system = "Mendota", res = "Month",troph = "Phytoplankton")

mad.zoomth.ccm<- pbmcapply::pbmclapply(c("zooFDis","zooFEve","zooFRic"),function(x){
  bio <- suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"density")],
                                   detrend.method = "lm",
                                  iter = 500,span =12*5,return.raw = T))
  pc <-  suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"community")],
                                   detrend.method = "lm",
                                    iter = 500,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"FI")],
                                   detrend.method = "lm",
                                   iter = 500,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"mvi")],
                                    detrend.method = "lm",
                                    iter = 500,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"zp.ratio")],
                                         detrend.method = "lm",
                                         iter = 500,span =12*5,return.raw = T))
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
                          ifelse(x_y.quantile >= 0.975 | x_y.quantile <= 0.025,"*","")))%>%
  mutate(y_x.sig = ifelse(is.na(y_x.obs_value) | y_x.obs_value == "NaN","",
                          ifelse(y_x.quantile >= 0.975 | y_x.quantile <= 0.025,"*","")))
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
                                  iter = 500,span =12*5,return.raw = T))
  pc <-  suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"community")],
                                   detrend.method = "lm",
                                    iter = 500,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"FI")],
                                   detrend.method = "lm",
                                   iter = 500,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"mvi")],
                                    detrend.method = "lm",
                                    iter = 500,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"zp.ratio")],
                                         detrend.method = "lm",
                                         iter = 500,span =12*5,return.raw = T))
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
                          ifelse(x_y.quantile >= 0.975 | x_y.quantile <= 0.025,"*","")))%>%
  mutate(y_x.sig = ifelse(is.na(y_x.obs_value) | y_x.obs_value == "NaN","",
                          ifelse(y_x.quantile >= 0.975 | y_x.quantile <= 0.025,"*","")))
LZ.phytomth.ccm.raw <- lapply(LZ.phytomth.ccm, `[[`, 'perm.dens')%>% 
  data.table::rbindlist(idcol = "FD.metric")%>% 
  mutate(system = "Lower Zurich", res = "Month",troph = "Phytoplankton")
LZ.phytomth.ccm.lag <- lapply(LZ.phytomth.ccm, `[[`, 'raw.obs')%>% 
  data.table::rbindlist(idcol = "FD.metric")%>% 
  mutate(system = "Lower Zurich", res = "Month",troph = "Phytoplankton")

LZ.zoomth.ccm<- pbmcapply::pbmclapply(c("zooFDis","zooFEve","zooFRic"),function(x){
  bio <- suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"density")],
                                   detrend.method = "lm",
                                  iter = 500,span =12*5,return.raw = T))
  pc <-  suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"community")],
                                   detrend.method = "lm",
                                    iter = 500,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"FI")],
                                   detrend.method = "lm",
                                   iter = 500,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"mvi")],
                                    detrend.method = "lm",
                                    iter = 500,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"zp.ratio")],
                                         detrend.method = "lm",
                                         iter = 500,span =12*5,return.raw = T))
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
                          ifelse(x_y.quantile >= 0.975 | x_y.quantile <= 0.025,"*","")))%>%
  mutate(y_x.sig = ifelse(is.na(y_x.obs_value) | y_x.obs_value == "NaN","",
                          ifelse(y_x.quantile >= 0.975 | y_x.quantile <= 0.025,"*","")))
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
                                  iter = 500,span =12*5,return.raw = T))
  pc <-  suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"community")],
                                   detrend.method = "lm",
                                    iter = 500,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"FI")],
                                   detrend.method = "lm",
                                   iter = 500,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"mvi")],
                                    detrend.method = "lm",
                                    iter = 500,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"zp.ratio")],
                                         detrend.method = "lm",
                                         iter = 500,span =12*5,return.raw = T))
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
                          ifelse(x_y.quantile >= 0.975 | x_y.quantile <= 0.025,"*","")))%>%
  mutate(y_x.sig = ifelse(is.na(y_x.obs_value) | y_x.obs_value == "NaN","",
                          ifelse(y_x.quantile >= 0.975 | y_x.quantile <= 0.025,"*","")))
wind.phytomth.ccm.raw <- lapply(wind.phytomth.ccm, `[[`, 'perm.dens')%>% 
  data.table::rbindlist(idcol = "FD.metric")%>% 
  mutate(system = "Windermere", res = "Month",troph = "Phytoplankton")
wind.phytomth.ccm.lag <- lapply(wind.phytomth.ccm, `[[`, 'raw.obs')%>% 
  data.table::rbindlist(idcol = "FD.metric")%>% 
  mutate(system = "Windermere", res = "Month",troph = "Phytoplankton")

wind.zoomth.ccm<- pbmcapply::pbmclapply(c("zooFDis","zooFEve","zooFRic"),function(x){
  bio <- suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"density")],
                                   detrend.method = "lm",
                                  iter = 500,span =12*5,return.raw = T))
  pc <-  suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"community")],
                                   detrend.method = "lm",
                                    iter = 500,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"FI")],
                                   detrend.method = "lm",
                                   iter = 500,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"mvi")],
                                    detrend.method = "lm",
                                    iter = 500,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"zp.ratio")],
                                         detrend.method = "lm",
                                         iter = 500,span =12*5,return.raw = T))
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
                          ifelse(x_y.quantile >= 0.975 | x_y.quantile <= 0.025,"*","")))%>%
  mutate(y_x.sig = ifelse(is.na(y_x.obs_value) | y_x.obs_value == "NaN","",
                          ifelse(y_x.quantile >= 0.975 | y_x.quantile <= 0.025,"*","")))
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
                                  iter = 500,span =12*5,return.raw = T))
  pc <-  suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"community")],
                                   detrend.method = "lm",
                                    iter = 500,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"FI")],
                                   detrend.method = "lm",
                                   iter = 500,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"mvi")],
                                    detrend.method = "lm",
                                    iter = 500,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"zp.ratio")],
                                         detrend.method = "lm",
                                         iter = 500,span =12*5,return.raw = T))
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
                          ifelse(x_y.quantile >= 0.975 | x_y.quantile <= 0.025,"*","")))%>%
  mutate(y_x.sig = ifelse(is.na(y_x.obs_value) | y_x.obs_value == "NaN","",
                          ifelse(y_x.quantile >= 0.975 | y_x.quantile <= 0.025,"*","")))
kas.phytomth.ccm.raw <- lapply(kas.phytomth.ccm, `[[`, 'perm.dens')%>% 
  data.table::rbindlist(idcol = "FD.metric")%>% 
  mutate(system = "Kasumigaura", res = "Month",troph = "Phytoplankton")
kas.phytomth.ccm.lag <- lapply(kas.phytomth.ccm, `[[`, 'raw.obs')%>% 
  data.table::rbindlist(idcol = "FD.metric")%>% 
  mutate(system = "Kasumigaura", res = "Month",troph = "Phytoplankton")

kas.zoomth.ccm<- pbmcapply::pbmclapply(c("zooFDis","zooFEve","zooFRic"),function(x){
  bio <- suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"density")],
                                   detrend.method = "lm",
                                  iter = 500,span =12*5,return.raw = T))
  pc <-  suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"community")],
                                   detrend.method = "lm",
                                    iter = 500,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"FI")],
                                   detrend.method = "lm",
                                   iter = 500,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"mvi")],
                                    detrend.method = "lm",
                                    iter = 500,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"zp.ratio")],
                                         detrend.method = "lm",
                                         iter = 500,span =12*5,return.raw = T))
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
                          ifelse(x_y.quantile >= 0.975 | x_y.quantile <= 0.025,"*","")))%>%
  mutate(y_x.sig = ifelse(is.na(y_x.obs_value) | y_x.obs_value == "NaN","",
                          ifelse(y_x.quantile >= 0.975 | y_x.quantile <= 0.025,"*","")))
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

###########################################################################
## Create figures ##
###########################################################################

summary.ccm <- read.csv(file ="Results/ccm/raw_data/ccm_summary.csv")
load(file = "Results/ccm/raw_data/ccm_raw.RData") # RData required to reduce file size compared to .csv

## Cross map skill ##

pdf(file="Results/ccm/FD_perm_y_x.pdf",
    width=10, height = 8)
ggplot(raw.ccm,aes(x = state.metric, y =  y_x.skill, col = FD.metric,fill= FD.metric)) + 
  geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.05, 0.5, 0.95),scale = "width",alpha = 0.3) +
  theme_bw() + 
  geom_point(data = summary.ccm[summary.ccm$measure %in% "max.skill",],
             aes(x = state.metric, y = y_x.obs_value),position = position_dodge(width = 0.9),size=2) +
  geom_text(data = summary.ccm[summary.ccm$measure %in% "max.skill",], 
            aes(x = state.metric, y = 0.98,label = y_x.sig),col= "black",size = 4,position = position_dodge(width = 0.9))+
  geom_text(data = summary.ccm[summary.ccm$measure %in% "t.max.skill",], 
            aes(x = state.metric, y = 0.9,label = y_x.obs_value),col= "black",size = 3,position = position_dodge(width = 0.9))+
  scale_y_continuous(breaks = c(0,0.5,1.0))+
  facet_grid(system~troph,scales = "free")+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  ylab("Cross map skill") + xlab("System state proxy")+   ggtitle("Permuted cross skill of system state mapping FD")
dev.off()

pdf(file="Results/ccm/FD_perm_x_y.pdf",
    width=10, height = 8)
ggplot(raw.ccm,aes(x = state.metric, y =  x_y.skill, col = FD.metric,fill= FD.metric)) + 
  geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.05, 0.5, 0.95),scale = "width",alpha = 0.3) +
  theme_bw() + 
  geom_point(data = summary.ccm[summary.ccm$measure %in% "max.skill",],
             aes(x = state.metric, y = x_y.obs_value),position = position_dodge(width = 0.9),size=2) +
  geom_text(data = summary.ccm[summary.ccm$measure %in% "max.skill",], 
            aes(x = state.metric, y = 0.98,label = y_x.sig),col= "black",size = 4,position = position_dodge(width = 0.9))+
  geom_text(data = summary.ccm[summary.ccm$measure %in% "t.max.skill",], 
            aes(x = state.metric, y = 0.9,label = x_y.obs_value),col= "black",size = 3,position = position_dodge(width = 0.9))+
  scale_y_continuous(breaks = c(0,0.5,1.0))+
  facet_grid(system~troph,scales = "free")+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  ylab("Cross map skill") + xlab("System state proxy")+   ggtitle("Permuted cross skill of FD mapping system state")
dev.off()

## Optimal lag ##

ccm.plot.df <- summary.ccm %>%
  filter(measure != "r0.skill")%>% # keep absolute highest cross map skill
  group_by(FD.metric,system,state.metric) %>%
  select(-c(x_y.median_perm_value,y_x.median_perm_value))%>% #drop unnecessary variables
  nest(x_y = c(measure:x_y.obs_value,x_y.sig), y_x =  c(measure,y_x.quantile:y_x.obs_value,y_x.sig))%>% #nest into x_y and y_x for ease of manipulation
  mutate(x_y = map(x_y, ~.x %>% 
          pivot_longer(cols = c(x_y.quantile:x_y.obs_value), #
                       names_to = c(".value"), names_prefix = "x_y.",
                       names_repair = "unique", values_to = "value") %>% #pivot just x_y
                     mutate(lag = obs_value[2])%>% #add lag from t.absmax.skill row
                     slice(-2) %>% #drop t.absmax.skill row
                     rename(sig = x_y.sig) %>% #remove excess prefix
                     setNames(paste0('x_y.', names(.))))) %>% #add prefix for downstream wrangling
  mutate(y_x = map(y_x, ~.x %>% 
          pivot_longer(cols = c(y_x.quantile:y_x.obs_value), #repeat process for y_x
                        names_to = c(".value"),names_prefix = "y_x.",
                        names_repair = "unique", values_to = "value") %>%
                     mutate(lag = obs_value[2])%>%
                     slice(-2) %>%
                     rename(sig = y_x.sig) %>%
                     setNames(paste0('y_x.', names(.))))) %>%
  unnest(cols = c(x_y, y_x),names_repair = "unique") %>% #unnest
  pivot_longer(c(x_y.measure:y_x.lag), #pivot using prefix as reference
               names_to = c("direc",".value"),
               names_pattern = "(.*)\\.(.*)" ) %>%
  mutate(causality.direc = ifelse(grepl("^x_y",direc),"FD map State","State map FD"))%>% 
      #classify directions. Is counterintuitive as here x_y represents x map y, 
      #which if significant implies y causes x (x contains information on y and 
      #therefore is causative)
  mutate(across(quantile:lag, ~as.numeric(.))) %>% #ensure values are numeric
  mutate(lag = -1*lag) %>% #for plotting purpose convert lags from negative to positive 
  #filter(sig == "*")%>% #only keep significant relationships
  #mutate(FD.metric = ifelse(troph == "Zooplankton", paste("zoo",FD.metric,sep = ""),FD.metric))%>%
  ungroup()

count.ccmdf <- ccm.plot.df %>%
  filter(sig == "*")%>% #only keep significant relationships
  group_by(state.metric,causality.direc,FD.metric,troph)%>%
  summarise(N = n()) #significant count per group

pdf(file="Results/ccm/ccm_causality_spread.pdf",
    width=10, height = 7)
ggplot(ccm.plot.df,aes(x=lag,y= state.metric,fill= causality.direc)) + 
  geom_vline(xintercept = 0,alpha=0.6)+
  geom_boxplot(alpha = 0.8,size = 0.5)+
  geom_point(aes(x=lag,y=state.metric, group=causality.direc), alpha = 0.5, position = position_dodge(width=0.75),size = 1.3) + 
  #scale_fill_manual(values = c("#FFE7A1","#A1B4FE"),name = "Causality\ndirection",labels = c("Forward", "Reverse"))+
  scale_fill_manual(values = c("#A1B4FE","#FFE7A1"),name = "Causality\ndirection")+
  geom_text(data = count.ccmdf,
            aes(x = (max(ccm.plot.df$lag)+5),y=state.metric, label = N),
            position = position_dodge(width = 0.8))+
  facet_grid(troph~FD.metric) +
  ylab("State metric") + xlab("Optimum lag (months)")+
  theme_bw()
dev.off()

pdf(file="Results/ccm/ccm_causality_spread_alt.pdf",
    width=10, height = 7)
ggplot(ccm.plot.df,aes(x=lag,y= state.metric)) + 
  geom_vline(xintercept = 0,alpha=0.6)+
  geom_point(aes(group=causality.direc,fill=causality.direc,shape = system,alpha = as.factor(sig),col= causality.direc), 
             position = position_dodge(width=0.75),size = 3) + 
  geom_point(aes(group=causality.direc,shape = system,fill=NULL,col= causality.direc), 
             position = position_dodge(width=0.75),size = 3) +
  geom_boxplot(aes(fill = causality.direc),alpha = 0.1,size = 0.2,col="black",outlier.shape = NA)+
  geom_text(data = count.ccmdf,
            aes(x = (max(ccm.plot.df$lag)+5),y=state.metric, label = N,group =causality.direc ),
            col="black", position = position_dodge(width = 0.8))+
  scale_fill_manual(values = c("#A1B4FE","#FFE7A1"),name = "Causality\ndirection")+
  scale_color_manual(values = c("#A1B4FE","#FFE7A1"),name = "Causality\ndirection")+
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_alpha_manual(values=c(0.01,1),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+  
  facet_grid(troph~FD.metric) +
  ylab("State metric") + xlab("Optimum lag (months)")+
  #guides(fill = guide_legend(override.aes = list(col = c("#FFE7A1","#A1B4FE"))))+
  #guides(fill = guide_legend(override.aes = list(col = c("#A1B4FE","#FFE7A1"))))+
  theme_bw()
dev.off()


## Cross skill changes with lag ##

lag.ccm <- read.csv(file ="Results/ccm/raw_data/ccm_lag.csv")

ccm.lag.plot.df <- lag.ccm %>%
  group_by(FD.metric,state.metric,troph)%>%
  pivot_longer(c(x_y,y_x),
               names_to = "direc", values_to = "skill")%>%
mutate(causality.direc = ifelse(grepl("^x_y",direc),"FD map State","State map FD"))
  
pdf(file="Results/ccm/ccm_lag_changes.pdf",
    width=12, height = 10,onefile = F)
ggpubr::ggarrange(
  ggplot(filter(ccm.lag.plot.df,troph %in% "Phytoplankton"),aes(x=tp,y=skill, col =causality.direc)) +
                    geom_path() + 
                    ylab("Cross map skill") + xlab("Lag (months)")+
                    ggh4x::facet_nested(system ~ FD.metric + state.metric ) + 
                    scale_colour_manual(values = c("#A1B4FE","#DEC98C"),name = "Causality\ndirection")+
                    #scale_colour_manual(values = c("#FFE7A1","#A1B4FE","#74A180","#FF94AB","#BE86FF"),name= "Lake")+
                    scale_x_continuous(breaks = c(-30,0,30))+
                    scale_y_continuous(breaks = c(0,0.5,1.0),limits = c(0,1.0))+
                    guides(color = guide_legend(override.aes = list(alpha = 1,size=1.5) ))+
                    geom_vline(xintercept = 0,colour="black")+
                    theme_bw() + ggtitle("Phytoplankton"),
  ggplot(filter(ccm.lag.plot.df,troph %in% "Zooplankton"),aes(x=tp,y=skill, col =causality.direc)) +
    geom_path() + 
    ylab("Cross map skill") + xlab("Lag (months)")+
    ggh4x::facet_nested(system ~ FD.metric + state.metric ) + 
    scale_x_continuous(breaks = c(-30,0,30))+
    scale_y_continuous(breaks = c(0,0.5,1.0),limits = c(0,1.0))+
    scale_colour_manual(values = c("#A1B4FE","#DEC98C"),name = "Causality\ndirection")+
    #scale_colour_manual(values = c("#FFE7A1","#A1B4FE","#74A180","#FF94AB","#BE86FF"),name= "Lake")+
    guides(color = guide_legend(override.aes = list(alpha = 1,size=1.5) ))+
    geom_vline(xintercept = 0,colour="black")+
    theme_bw() + ggtitle("Zooplankton"),
                  nrow=2,common.legend=T)
dev.off()

###########################################################################
## Summary cross skills ##
###########################################################################
summary.ccm <- read.csv(file ="Results/ccm/raw_data/ccm_summary.csv")

## Proportion of forward vs reverse signficant cross maps ##
ccm.lag0.comp <- summary.ccm %>%
  filter(measure == "r0.skill")%>%  
  mutate(forward = ifelse(y_x.sig == "*" & x_y.sig != "*",TRUE,FALSE),
         reverse = ifelse(x_y.sig == "*" & y_x.sig != "*",TRUE,FALSE),
         bidirec =  ifelse(x_y.sig == "*" & y_x.sig == "*",TRUE,FALSE),
         none =  ifelse(x_y.sig != "*" & y_x.sig != "*",TRUE,FALSE),
         diff.lag = filter(summary.ccm,measure == "t.max.skill")$x_y.obs_value - filter(summary.ccm,measure == "t.max.skill")$y_x.obs_value)%>%
  group_by(troph,FD.metric,state.metric) %>%
  summarise(prop.forward=sum(forward == TRUE)/length(forward),
            prop.reverse=sum(reverse == TRUE)/length(reverse),
            prop.bidirec=sum(bidirec == TRUE)/length(bidirec),
            prop.none=sum(none == TRUE)/length(none),
            mean.lag = mean(diff.lag))
# mutate(ref.y = ifelse(FD.metric %in% "FDis",1.3, 
#                       ifelse(FD.metric %in% "FEve",1.2,1.1))) %>%
# mutate(ref.col = ifelse(FD.metric %in% "red",1.3, 
#                         ifelse(FD.metric %in% "FEve",1.2,1.1)))
write.csv(ccm.lag0.comp,file ="Results/ccm/ccm_tables/skill.comp.lag0.state.tab.csv",row.names = F)
ccm.lag0.comp <- read.csv("Results/ccm/ccm_tables/skill.comp.lag0.state.tab.csv")

ccm.lagx.comp <- summary.ccm %>%
  filter(measure == "max.skill")%>%  
  mutate(forward = ifelse(y_x.sig == "*" & x_y.sig != "*",TRUE,FALSE),
         reverse = ifelse(x_y.sig == "*" & y_x.sig != "*",TRUE,FALSE),
         bidirec =  ifelse(x_y.sig == "*" & y_x.sig == "*",TRUE,FALSE),
         none =  ifelse(x_y.sig != "*" & y_x.sig != "*",TRUE,FALSE),
         diff.lag = filter(summary.ccm,measure == "t.max.skill")$x_y.obs_value - filter(summary.ccm,measure == "t.max.skill")$y_x.obs_value)%>%
  group_by(troph,FD.metric,state.metric) %>%
  summarise(prop.forward=sum(forward == TRUE)/length(forward),
            prop.reverse=sum(reverse == TRUE)/length(reverse),
            prop.bidirec=sum(bidirec == TRUE)/length(bidirec),
            prop.none=sum(none == TRUE)/length(none),
            mean.lag = mean(diff.lag))
# mutate(ref.y = ifelse(FD.metric %in% "FDis",1.3, 
#    ifelse(FD.metric %in% "FEve",1.2,1.1))) %>%
# mutate(ref.col = ifelse(FD.metric %in% "red",1.3, 
#                       ifelse(FD.metric %in% "FEve",1.2,1.1)))
write.csv(ccm.lagx.comp,file ="Results/ccm/ccm_tables/skill.comp.lagx.state.tab.csv",row.names = F)
ccm.lagx.comp <- read.csv("Results/ccm/ccm_tables/skill.comp.lagx.state.tab.csv")

ccm.lag0.lagx.comp <- left_join(obs.ccm.y_x.lag0.state.tab,obs.ccm.y_x.lagx.state.tab,
                                by=c("troph","state.metric","FD.metric"),.groups = "rowwise",
                                suffix = c("_lag0","_lagx")) %>%
  dplyr::select(-c(mean.cor_lag0,mean.cor_lagx,cor.se_lag0,cor.se_lagx,nsig_lag0,nsig_lagx)) %>% 
  group_by()%>% rowwise()%>%
  dplyr::summarise(troph = troph,
                   FD.metric = FD.metric,
                   state.metric = state.metric,
                   cor_lag0 = median.cor_lag0,
                   cor_lagx =median.cor_lagx,
                   diff.cor = (median.cor_lagx-median.cor_lag0),
                   prop.sig_lag0 = prop.sig_lag0,
                   prop.sig_lagx =prop.sig_lagx,
                   diff.prop.sig = (prop.sig_lagx-prop.sig_lag0))

## Plots ##

pdf(file="Results/ccm/summary_ccm_r0.pdf",
    width=8, height = 5)  
plag0.fin <- ggplot(filter(summary.ccm,measure %in% "r0.skill"),
                    aes(x=state.metric,y=y_x.obs_value,col=FD.metric))+
  #geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.025, 0.5, 0.975),scale = "width",alpha = 0.3)+
  geom_hline(yintercept = 0,col="black",alpha = 0.3)+
  geom_point(position=position_dodge(width=0.75),
             aes(shape=system,alpha=y_x.sig,fill=FD.metric,group=FD.metric),size=3.5)+
  geom_point(position=position_dodge(width=0.75),
             aes(shape=system,fill=NULL,group=FD.metric),size=3.5) +
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric")+
  scale_alpha_manual(values=c(0.01,1),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  geom_boxplot(aes(fill=FD.metric),alpha=0.1,col="black",size=0.3,outlier.shape = NA)+
  ggtext::geom_richtext(data =ccm.lag0.comp,aes(x = state.metric, y =1.05,
              label = paste("<span style='color:black'>","(","</span>","<span style='color:#DEC98C'>",base::format(prop.forward,digits = 2),"</span>","<span style='color:black'>",",",base::format(prop.bidirec,digits =2),")","</span>",sep = "")),
                        alpha=0,size = 3, position = position_dodge(width = 0.75),angle = 90)+
  scale_y_continuous(breaks = seq(0,1.0,0.25),limits = c(-0.1,1.1))+
  facet_wrap(~troph)+
  ylab("Cross skill") + xlab("System state proxy")+
  guides(size = guide_legend(order = 1), 
         col = guide_legend(order = 2),
         fill = guide_legend(order = 2),
         shape = guide_legend(order = 3))+
  theme_bw()

plag0.fin
dev.off()

obs.ccm.y_x.lag0.state.tab <- summary.ccm %>%
  dplyr::select(!starts_with("x_y"))%>%
  filter(measure == "r0.skill")%>%
  group_by(troph,FD.metric,state.metric) %>%
  summarise(mean.cor = mean(y_x.obs_value),median.cor = median(y_x.obs_value),
            cor.se = sd(y_x.obs_value)/n(),
            nsig=sum(y_x.sig %in% "*"),prop.sig = sum(y_x.sig %in% "*")/length(y_x.sig)) %>%
  mutate(across(mean.cor:cor.se,~round(.x,digits=4)))
write.csv(obs.ccm.y_x.lag0.state.tab,file ="Results/ccm/ccm_tables/skill.y_x.lag0.state.tab.csv",row.names = F)
obs.ccm.y_x.lag0.state.tab <- read.csv("Results/ccm/ccm_tables/skill.y_x.lag0.state.tab.csv")

obs.ccm.x_y.lag0.state.tab <- summary.ccm %>%
  dplyr::select(!starts_with("y_x"))%>%
  filter(measure == "r0.skill")%>%
  group_by(troph,FD.metric,state.metric) %>%
  summarise(mean.cor = mean(x_y.obs_value),median.cor = median(x_y.obs_value),
            cor.se = sd(x_y.obs_value)/n(),
            nsig=sum(x_y.sig %in% "*"),prop.sig = sum(x_y.sig %in% "*")/length(x_y.sig)) %>%
  mutate(across(mean.cor:cor.se,~round(.x,digits=4)))
write.csv(obs.ccm.x_y.lag0.state.tab,file ="Results/ccm/ccm_tables/skill.x_y.lag0.state.tab.csv",row.names = F)
obs.ccm.x_y.lag0.state.tab <- read.csv("Results/ccm/ccm_tables/skill.x_y.lag0.state.tab.csv")

pdf(file="Results/ccm/summary_ccm_r0_alt.pdf",
    width=8, height = 5)  
plag0.1 <- ggplot(filter(summary.ccm,measure %in% "r0.skill"),aes(x=state.metric,y=y_x.obs_value,col=FD.metric))+
  geom_line(aes(linetype = system),position = position_dodge(width = 0.75), 
            alpha = 0, show.legend = FALSE)+
  facet_wrap(~troph)

plag0.fin <- plag0.1 + geom_segment(data = layer_data(plag0.1, 1L),
                  aes(x = xmin, xend=xmax, y = 0,yend=0, group = linetype),
                  color = "black", size = 0.5)+
  #geom_linerange(aes(xmin=state.metric,xmax=state.metric),position = position_dodge(width=0.75))+
  geom_point(position=position_dodge(width=0.75),
             aes(shape=system,alpha=y_x.sig,fill=FD.metric,group=FD.metric),size=3.5)+
  geom_point(position=position_dodge(width=0.75),
             aes(shape=system,fill=NULL,group=FD.metric),size=3.5)+
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric")+
  scale_alpha_manual(values=c(0.01,1),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  geom_point(data = obs.ccm.y_x.lag0.state.tab,aes(x=state.metric,y=mean.cor,group=FD.metric,size="mean"),position=position_dodge(width=0.75),col="black",alpha=1,shape = 21,fill="black")+
  geom_linerange(data = obs.ccm.y_x.lag0.state.tab,aes(x=state.metric,y = 0,ymin = 0,ymax=mean.cor,group=FD.metric),position=position_dodge(width=0.75),col="black",linetype = "solid")+
  # geom_text(data =count.ccf.r0.dat,aes(x = state.metric, y = ref.y,label = paste("(",NSig,"*",",",NxSig,")",sep = ""),fill=FD.metric,group = FD.metric),
  #           col= "black",size = 3, position = position_dodge(width = 0.9))+
  scale_size_manual(values = c(2.5),breaks = c("mean"),name = NULL, labels = c("Mean cross\ncorrelation"),
                    guide = guide_legend(override.aes = list(fill = c("black"),shape=c(21),size=c(3.5),alpha = c(1))))+
  ylab("Cross correlation") + xlab("System state proxy")+
  guides(size = guide_legend(order = 1), 
         col = guide_legend(order = 2),
         fill = guide_legend(order = 2),
         shape = guide_legend(order = 3))+
  theme_bw()

plag0.fin
dev.off()

obs.ccm.y_x.lagx.state.tab <- summary.ccm %>%
  dplyr::select(!starts_with("x_y"))%>%
  filter(measure %in% c("max.skill","t.max.skill"))%>%
  pivot_wider(names_from = measure,values_from = y_x.obs_value)%>%
  ungroup()%>%
  mutate(t.max.skill=ifelse(is.numeric(max.skill) & !is.na(t.max.skill),t.max.skill,
                            dplyr::lead(t.max.skill)))%>% # fill NA t.absmax with next t.absmax to associate cor with lag
  na.omit() %>%#drop duplicate rows
  group_by(troph,FD.metric,state.metric) %>%
  summarise(mean.cor = mean(max.skill),median.cor = median(max.skill),
            cor.se = sd(max.skill)/n(),
            median.lag = median(t.max.skill),lag.se = sd(t.max.skill)/n(), 
            nsig=sum(y_x.sig %in% "*"),prop.sig = sum(y_x.sig %in% "*")/length(y_x.sig))%>%
  mutate(across(mean.cor:lag.se,~round(.x,digits=4)))
write.csv(obs.ccm.y_x.lagx.state.tab,file ="Results/ccm/ccm_tables/skill.y_x.lagx.state.tab.csv",row.names = F)
obs.ccm.y_x.lagx.state.tab <- read.csv("Results/ccm/ccm_tables/skill.y_x.lagx.state.tab.csv")

obs.ccm.x_y.lagx.state.tab <- summary.ccm %>%
  dplyr::select(!starts_with("y_x"))%>%
  filter(measure %in% c("max.skill","t.max.skill"))%>%
  pivot_wider(names_from = measure,values_from = x_y.obs_value)%>%
  ungroup()%>%
  mutate(t.max.skill=ifelse(is.numeric(max.skill) & !is.na(t.max.skill),t.max.skill,
                            dplyr::lead(t.max.skill)))%>% # fill NA t.absmax with next t.absmax to associate cor with lag
  na.omit() %>%#drop duplicate rows
  group_by(troph,FD.metric,state.metric) %>%
  summarise(mean.cor = mean(max.skill),median.cor = median(max.skill),
            cor.se = sd(max.skill)/n(),
            median.lag = median(t.max.skill),lag.se = sd(t.max.skill)/n(), 
            nsig=sum(x_y.sig %in% "*"),prop.sig = sum(x_y.sig %in% "*")/length(x_y.sig))%>%
  mutate(across(mean.cor:lag.se,~round(.x,digits=4)))
write.csv(obs.ccm.x_y.lagx.state.tab,file ="Results/ccm/ccm_tables/skill.x_y.lagx.state.tab.csv",row.names = F)
obs.ccm.x_y.lagx.state.tab <- read.csv("Results/ccm/ccm_tables/skill.x_y.lagx.state.tab.csv")

pdf(file="Results/ccm/summary_ccm_lagx.pdf",
    width=10, height = 6)
pccm.lagx.1 <- ggplot(filter(summary.ccm,measure %in% "max.skill") %>% 
                        mutate(direc= ifelse(y_x.sig == "*" & x_y.sig != "*","x","y")),
                      aes(x=state.metric,y=y_x.obs_value,col=FD.metric))+
  geom_line(aes(linetype = system),position = position_dodge(width = 0.75), 
            alpha = 0, show.legend = F)+
  facet_wrap(~troph)

pccm.lagx.2 <- pccm.lagx.1 +
  # geom_segment(data = layer_data(plagx.1, 1L),
  #              aes(x = xmin, xend=xmax, y = 0,yend=0, group = linetype),
  #              color = "black", size = 0.5)+
  geom_point(position=position_dodge(width=0.75),
             aes(shape=system,alpha=y_x.sig,fill=FD.metric,group=FD.metric,col=FD.metric),size=3.5)+
  geom_point(position=position_dodge(width=0.75),
             aes(shape=system,fill=NULL,group=FD.metric,col=FD.metric),size=3.5)+
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric")+
  scale_alpha_manual(values=c(0.01,1),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  #geom_point(data = obs.ccm.y_x.lagx.state.tab,aes(x=state.metric,y=mean.cor,group=FD.metric,size="mean"),position=position_dodge(width=0.75),col="black",alpha=1,shape = 21,fill="black")+
  #geom_linerange(data = obs.ccm.y_x.lagx.state.tab,aes(x=state.metric,y = 0,ymin = 0,ymax=mean.cor,group=FD.metric),position=position_dodge(width=0.75),col="black",linetype = "solid")+
  geom_boxplot(aes(fill=FD.metric),alpha=0.1,col="black",size=0.3,outlier.shape = NA)+
  #scale_size_manual(values = c(2.5),breaks = c("mean"),name = NULL, labels = c("Mean cross\ncorrelation"),
   #                 guide = guide_legend(override.aes = list(fill = c("black"),shape=c(21),size=c(3.5),alpha = c(1))))+
  ggtext::geom_richtext(data =ccm.lagx.comp,aes(x = state.metric, y =1.05,
                                                label = paste("<span style='color:black'>","(","</span>","<span style='color:#DEC98C'>",base::format(prop.forward,digits = 2),"</span>","<span style='color:black'>",",",base::format(prop.bidirec,digits =2),")","</span>",sep = "")),
                        alpha=0,size = 3, position = position_dodge(width = 0.75),angle = 90)+
  scale_size_manual(values = c("a","b"),name = "Causality direction", breaks = c("forward","bidirec"),labels = c("Forward","Bidirectional"),
                        guide = guide_legend(override.aes = list(shape=c(21,20),size=c(3.5),alpha = c(1,1),color = c("#DEC98C","black"))))+
  scale_y_continuous(breaks = seq(0,1.0,0.25),limits = c(0,1.1))+
  facet_wrap(~troph)+
  #facet_grid(troph~FD.metric)+
  ylab("Cross skill") + xlab("System state proxy")+
  theme_bw() +
  guides(size = guide_legend(order = 1), 
         col = guide_legend(order = 2),
         fill = guide_legend(order = 2),
         shape = guide_legend(order = 3))+
  theme(axis.title.x=element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(c(2, 2, 0, 2)))

pccm.lagx.3 <- ggplot(filter(summary.ccm,measure %in% "t.max.skill") %>% 
                        mutate(sig = as.factor(filter(summary.ccm,measure %in% "max.skill")$y_x.sig)),
                      aes(x=state.metric,y=y_x.obs_value,col = FD.metric))+
  geom_hline(yintercept = 0,col="black",alpha = 0.3)+
  scale_y_binned(breaks = c(seq(60,24,-12),12,seq(-12,-60,-12)),show.limits = T)+
  geom_tile(data = expand_grid(c(seq(60,24,-12),12,seq(-12,-60,-12)),unique(summary.ccm$FD.metric),unique(summary.ccm$state.metric),unique(summary.ccm$troph)) %>%
              magrittr::set_colnames(c("y_x.obs_value","FD.metric","state.metric","troph")),
            aes(group = FD.metric),fill = "white",stat="identity",position = position_dodge(width = 0.75), col = "black", size = 0.3,width = 0.8, height = 0.9)+
  geom_point(position=position_jitterdodge(dodge.width=0.75,jitter.height = 0.2,jitter.width = 0,seed = 5),
             aes(y = y_x.obs_value,shape=system,alpha=sig,col = FD.metric,fill=FD.metric,group=FD.metric),size=1.5)+
  geom_point(position=position_jitterdodge(dodge.width=0.75,jitter.height = 0.2,jitter.width = 0,seed = 5),
             aes(y = y_x.obs_value,shape=system,col = FD.metric,fill=NULL,group=FD.metric),size=1.5)+
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake",guide = "none")+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric",guide = "none") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric",guide = "none")+
  scale_alpha_manual(values=c(0.01,1),name = "Significance", labels = c("Not significant","Significant"),
                     guide = "none")+
  facet_wrap(~troph)+
  #facet_grid(troph~FD.metric)+
  ylab("Optimal lag (months)") + xlab("System state proxy")+
  theme_bw()+
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = margin(c(2, 2, 0, 2)))

pccm.lagx.fin <- pccm.lagx.2 + pccm.lagx.3 + patchwork::plot_layout(nrow = 2,guides = "collect",heights = c(2, 1),tag_level = 'new') #tag level required for patchwork multiplots
pccm.lagx.fin
dev.off()

pdf(file="Results/ccm/summary_ccm_lagx_alt.pdf",
    width=10, height = 6)
pccm.lagx.1 <- ggplot(filter(summary.ccm,measure %in% "max.skill"),aes(x=state.metric,y=y_x.obs_value,col=FD.metric))+
  geom_line(aes(linetype = system),position = position_dodge(width = 0.75), 
            alpha = 0, show.legend = FALSE)+
  facet_wrap(~troph)

pccm.lagx.2 <- pccm.lagx.1 +
  geom_segment(data = layer_data(pccm.lagx.1, 1L),
                aes(x = xmin, xend=xmax, y = 0,yend=0, group = linetype),
                color = "black", size = 0.5)+
  geom_point(position=position_dodge(width=0.75),
             aes(shape=system,alpha=y_x.sig,fill=FD.metric,group=FD.metric),size=3.5)+
  geom_point(position=position_dodge(width=0.75),
             aes(shape=system,fill=NULL,group=FD.metric),size=3.5)+
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric")+
  scale_alpha_manual(values=c(0.01,1),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  geom_point(data = obs.ccm.y_x.lagx.state.tab,aes(x=state.metric,y=mean.cor,group=FD.metric,size="mean"),position=position_dodge(width=0.75),col="black",alpha=1,shape = 21,fill="black")+
  geom_linerange(data = obs.ccm.y_x.lagx.state.tab,aes(x=state.metric,y = 0,ymin = 0,ymax=mean.cor,group=FD.metric),position=position_dodge(width=0.75),col="black",linetype = "solid")+
  #geom_boxplot(aes(fill=FD.metric),alpha=0.1,col="black",size=0.3,outlier.shape = NA)+
  scale_size_manual(values = c(2.5),breaks = c("mean"),name = NULL, labels = c("Mean cross\ncorrelation"),
                    guide = guide_legend(override.aes = list(fill = c("black"),shape=c(21),size=c(3.5),alpha = c(1))))+
  # geom_text(data =ccm.lagx.comp,aes(x = state.metric, y =ref.y,label = paste("(",prop.forward,",",prop.bidirec,",",prop.reverse,")",sep = ""),group = FD.metric),
  #           col= "black",size = 3, position = position_dodge(width = 0.9),angle=45)+
  ggtext::geom_richtext(data =ccm.lagx.comp,aes(x = state.metric, y =ref.y,
          label = paste("<span style='color:black'>","(","</span>","<span style='color:#DEC98C'>",prop.forward,"</span>","<span style='color:black'>",",",prop.bidirec,",","</span>","<span style='color:#A1B4FE'>",prop.reverse,"</span>","<span style='color:black'>",")","</span>",sep = "")),
    alpha=0,size = 3, position = position_dodge(width = 0.9))+
  scale_y_continuous(breaks = seq(0,1.0,0.25))+
  facet_wrap(~troph)+
  #facet_grid(troph~FD.metric)+
  ylab("Cross skill") + xlab("System state proxy")+
  theme_bw() +
  guides(size = guide_legend(order = 1), 
         col = guide_legend(order = 2),
         fill = guide_legend(order = 2),
         shape = guide_legend(order = 3))+
  theme(axis.title.x=element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(c(2, 2, 0, 2)))

pccm.lagx.3 <- ggplot(filter(summary.ccm,measure %in% "t.max.skill") %>% 
                     mutate(sig = as.factor(filter(summary.ccm,measure %in% "max.skill")$y_x.sig)),
                   aes(x=state.metric,y=y_x.obs_value,col = FD.metric))+
  geom_hline(yintercept = 0,col="black",alpha = 0.3)+
  geom_point(position=position_dodge(width=0.75),
             aes(shape=system,alpha=sig,fill=FD.metric,group=FD.metric),size=2)+
  geom_point(position=position_dodge(width=0.75),
             aes(shape=system,fill=NULL,group=FD.metric),size=2) +
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake",guide = "none")+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric",guide = "none") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric",guide = "none")+
  scale_alpha_manual(values=c(0.01,1),name = "Significance", labels = c("Not significant","Significant"),
                     guide = "none")+
  #scale_x_discrete(position = "top") +
  scale_y_reverse()+
  facet_wrap(~troph)+
  #facet_grid(troph~FD.metric)+
  ylab("Optimal lag") + xlab("System state proxy")+
  theme_bw()+
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = margin(c(2, 2, 0, 2)))

pccm.lagx.fin <- pccm.lagx.2 + pccm.lagx.3 +plot_layout(nrow = 2,guides = "collect",heights = c(2, 1))
pccm.lagx.fin
dev.off()

## Combo summary figure ##
pdf(file="Results/ccm/summary_ccm_combo.pdf",
    width=15, height = 6)  
plag0.fin + pccm.lagx.fin + 
  patchwork::plot_layout(ncol = 2,nrow = 1, guides = "collect") +
  plot_annotation(tag_levels = c('a',1)) & 
  theme(plot.tag = element_text(face = "bold"))
dev.off()

