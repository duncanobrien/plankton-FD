### Cross Map System State with FD ###

## Preamble ##
require(tidyverse) # dplyr, ggplot etc.
require(pbmcapply) # paralled lapply 
require(data.table) # rbindlist function

source("Code/ccm_perm_fn.R")

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
  pc <- suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"density")],
                                       iter = 500,span =12*5,return.raw = T))
  bio <-  suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"community")],
                                    iter = 500,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"FI")],
                                   iter = 500,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"mvi")],
                                    iter = 500,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"zp.ratio")],
                                         iter = 500,span =12*5,return.raw = T))
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
  pc <- suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"density")],
                                  iter = 500,span =12*5,return.raw = T))
  bio <-  suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"community")],
                                    iter = 500,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"FI")],
                                   iter = 500,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"mvi")],
                                    iter = 500,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = kin.tot[,c("date",paste(x),"zp.ratio")],
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
write.csv(kin.raw.ccm,file ="Results/ccm/raw_data/kin_ccm_raw.csv",row.names = F)
write.csv(kin.lag.ccm,file ="Results/ccm/raw_data/kin_ccm_lag.csv",row.names = F)

## Mendota CCM ##

mad.phytomth.ccm<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"density")],
                                  iter = 500,span =12*5,return.raw = T))
  bio <-  suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"community")],
                                    iter = 500,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"FI")],
                                   iter = 500,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"mvi")],
                                    iter = 500,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"zp.ratio")],
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
  pc <- suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"density")],
                                  iter = 500,span =12*5,return.raw = T))
  bio <-  suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"community")],
                                    iter = 500,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"FI")],
                                   iter = 500,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"mvi")],
                                    iter = 500,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = mad.tot[,c("date",paste(x),"zp.ratio")],
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
  pc <- suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"density")],
                                  iter = 500,span =12*5,return.raw = T))
  bio <-  suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"community")],
                                    iter = 500,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"FI")],
                                   iter = 500,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"mvi")],
                                    iter = 500,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"zp.ratio")],
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
  pc <- suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"density")],
                                  iter = 500,span =12*5,return.raw = T))
  bio <-  suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"community")],
                                    iter = 500,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"FI")],
                                   iter = 500,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"mvi")],
                                    iter = 500,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = LZ.tot[,c("date",paste(x),"zp.ratio")],
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

write.csv(LZ.summary.ccm,file ="/Users/duncanobrien/Desktop/LZ_ccm_summary.csv",row.names = F)
write.csv(LZ.raw.ccm,file ="/Users/duncanobrien/Desktop/LZ_ccm_raw.csv",row.names = F)
write.csv(LZ.lag.ccm,file ="/Users/duncanobrien/Desktop/LZ_ccm_lag.csv",row.names = F)

write.csv(LZ.summary.ccm,file ="Results/ccm/raw_data/LZ_ccm_summary.csv",row.names = F)
write.csv(LZ.raw.ccm,file ="Results/ccm/raw_data/LZ_ccm_raw.csv",row.names = F)
write.csv(LZ.lag.ccm,file ="Results/ccm/raw_data/LZ_ccm_lag.csv",row.names = F)

## Windermere CCM ##

wind.phytomth.ccm<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"density")],
                                  iter = 500,span =12*5,return.raw = T))
  bio <-  suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"community")],
                                    iter = 500,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"FI")],
                                   iter = 500,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"mvi")],
                                    iter = 500,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"zp.ratio")],
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
  pc <- suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"density")],
                                  iter = 500,span =12*5,return.raw = T))
  bio <-  suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"community")],
                                    iter = 500,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"FI")],
                                   iter = 500,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"mvi")],
                                    iter = 500,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = wind.tot[,c("date",paste(x),"zp.ratio")],
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

ggplot(wind.raw.ccm,aes(x = state.metric, y =  y_x.skill, col = FD.metric,fill= FD.metric)) + 
  geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.05, 0.5, 0.95),scale = "width",alpha = 0.3) +
  theme_bw() + 
  geom_point(data = wind.summary.ccm[wind.summary.ccm$measure %in% "max.skill",],
             aes(x = state.metric, y = y_x.obs_value),position = position_dodge(width = 0.9),size=2) +
  geom_text(data = wind.summary.ccm[wind.summary.ccm$measure %in% "max.skill",], 
            aes(x = state.metric, y = 0.98,label = y_x.sig),col= "black",size = 4,position = position_dodge(width = 0.9))+
  geom_text(data = wind.summary.ccm[wind.summary.ccm$measure %in% "t.max.skill",], 
            aes(x = state.metric, y = 0.9,label = x_y.obs_value),col= "black",size = 3,position = position_dodge(width = 0.9))+
  scale_y_continuous(breaks = c(0,0.5,1.0))+
  facet_wrap(~troph,scales = "free")+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  ylab("Cross correlation") + xlab("System state proxy")+   ggtitle("Permuted cross skill of FD mapping system state")

## Kasumigaura CCM ##

kas.phytomth.ccm<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"density")],
                                  iter = 500,span =12*5,return.raw = T))
  bio <-  suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"community")],
                                    iter = 500,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"FI")],
                                   iter = 500,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"mvi")],
                                    iter = 500,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"zp.ratio")],
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
  pc <- suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"density")],
                                  iter = 500,span =12*5,return.raw = T))
  bio <-  suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"community")],
                                    iter = 500,span =12*5,return.raw = T))
  fi <-  suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"FI")],
                                   iter = 500,span =12*5,return.raw = T))
  mvi <-  suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"mvi")],
                                    iter = 500,span =12*5,return.raw = T))
  zp.ratio <-  suppressWarnings(ccm.perm(dat = kas.tot[,c("date",paste(x),"zp.ratio")],
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
write.csv(kas.raw.ccm,file ="Results/ccm/raw_data/kas_ccm_raw.csv",row.names = F)
write.csv(kas.lag.ccm,file ="Results/ccm/raw_data/kas_ccm_lag.csv",row.names = F)

###########################################################################
## Save out ##
###########################################################################
summary.ccm <- rbind(kin.phytomth.ccm.summary,kin.zoomth.ccm.summary,mad.phytomth.ccm.summary,mad.zoomth.ccm.summary,
                          LZ.phytomth.ccm.summary,LZ.zoomth.ccm.summary,wind.phytomth.ccm.summary,wind.zoomth.ccm.summary,
                          kas.phytomth.ccm.summary,kas.zoomth.ccm.summary)

summary.ccm <- rbind(kin.summary.ccm,mad.summary.ccm,
                     LZ.summary.ccm,kas.summary.ccm,
                     wind.summary.ccm)

raw.ccm <- rbind(kin.phytomth.ccm.raw,kin.zoomth.ccm.raw,mad.phytomth.ccm.raw,mad.zoomth.ccm.raw,
                      LZ.phytomth.ccm.raw,LZ.zoomth.ccm.raw,wind.phytomth.ccm.raw,wind.zoomth.ccm.raw,
                      kas.phytomth.ccm.raw,kas.zoomth.ccm.raw)

raw.ccm <- rbind(kin.raw.ccm,mad.raw.ccm,
                     LZ.raw.ccm,kas.raw.ccm,
                     wind.raw.ccm)

lag.ccm <- rbind(kin.phytomth.ccm.lag,kin.zoomth.ccm.lag,mad.phytomth.ccm.lag,mad.zoomth.ccm.lag,
                 LZ.phytomth.ccm.lag,LZ.zoomth.ccm.lag,wind.phytomth.ccm.lag,wind.zoomth.ccm.lag,
                 kas.phytomth.ccm.lag,kas.zoomth.ccm.lag)

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
  ylab("Cross correlation") + xlab("System state proxy")+   ggtitle("Permuted cross skill of system state mapping FD")
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
  ylab("Cross correlation") + xlab("System state proxy")+   ggtitle("Permuted cross skill of FD mapping system state")
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
  filter(sig == "*")%>% #only keep significant relationships
  #mutate(FD.metric = ifelse(troph == "Zooplankton", paste("zoo",FD.metric,sep = ""),FD.metric))%>%
  ungroup()

count.ccmdf <- ccm.plot.df %>%
  group_by(state.metric,causality.direc,FD.metric,troph)%>%
  summarise(N = n()) #significant count per group

pdf(file="Results/ccm/ccm_causality_spread.pdf",
    width=10, height = 7)
ggplot(ccm.plot.df,aes(x=lag,y= state.metric,fill= causality.direc)) + 
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
ggplot(ccm.plot.df,aes(x=lag,y= state.metric,fill= causality.direc)) + 
  geom_boxplot(alpha = 0.1,size = 0.1)+
  geom_point(aes(x=lag,y=state.metric, group=causality.direc,fill=causality.direc,shape = system),  alpha = 0.8, position = position_dodge(width=0.75),size = 3) + 
  geom_text(data = count.ccmdf,
            aes(x = (max(ccm.plot.df$lag)+5),y=state.metric, label = N),
            position = position_dodge(width = 0.8))+
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  #scale_fill_manual(values = c("#FFE7A1","#A1B4FE"),name = "Causality\ndirection",labels = c("Forward", "Reverse"))+
  scale_fill_manual(values = c("#A1B4FE","#FFE7A1"),name = "Causality\ndirection")+
  facet_grid(troph~FD.metric) +
  ylab("State metric") + xlab("Optimum lag (months)")+
  #guides(fill = guide_legend(override.aes = list(col = c("#FFE7A1","#A1B4FE"))))+
  guides(fill = guide_legend(override.aes = list(col = c("#A1B4FE","#FFE7A1"))))+
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
                    ylab("Cross map skill") + xlab("Lag")+
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
    ylab("Cross map skill") + xlab("Lag")+
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
