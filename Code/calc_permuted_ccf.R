### Cross Correlate System State with FD ###

## Preamble ##
require(tidyverse) # dplyr, ggplot etc.
require(pbmcapply) # paralled lapply 
require(data.table) # rbindlist function

source("Code/diff_perm_ccf_fn.R")

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
## Estimate cross correlation and permute ##
###########################################################################
# repeated 'diff.perm.ccf' function specified as @iter = 10000, @perm.method = "red.noise",
# @detrend.method = "lm", @span =12*5, @identical.t = F, and @lag=1.
# Results in two dataframes:
# $summary = observed absolute max correlation coef and the corresponding lag, observed correlation at lag0 and summary statistics
# $perm.dens = all permutations absolute max correlation coef and the corresponding lag plus the correlation at lag0

# Kinneret Phytoplankton and Zooplankton #

kin.phytomth.diff<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <-  suppressWarnings(diff.perm.ccf(dat = kin.tot[,c("date",paste(x),"community")],
                                        iter = 10000,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  bio <- suppressWarnings(diff.perm.ccf(dat = kin.tot[,c("date",paste(x),"density")],
                                        iter = 10000,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  fi <-  suppressWarnings(diff.perm.ccf(dat = kin.tot[,c("date",paste(x),"FI")],
                                        iter = 10000,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  mvi <-  suppressWarnings(diff.perm.ccf(dat = kin.tot[,c("date",paste(x),"mvi")],
                                         iter = 10000,span =12*5,lag=1,
                                         perm.method = "red.noise",identical.t = F,
                                         detrend.method = "lm"))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(dat = kin.tot[,c("date",paste(x),"zp.ratio")],
                                              iter = 10000,span =12*5,lag=1,
                                              perm.method = "red.noise",identical.t = F,
                                              detrend.method = "lm"))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  # extract the observed correlation coefs for FD vs each system state
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  # extract all permuted correlation coefs for FD vs each system state
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1) 
# only single core specified as 'ccm.perm' already paralled within the function. 
# Provides opportunity for further parallelisation if desired
names(kin.phytomth.diff) <- c("FDis","FEve","FRic") # name list elements
kin.phytomth.diff.summary <- lapply(kin.phytomth.diff, `[[`, 'summary')%>% # extract second level list elements (i.e. 'summary')
  data.table::rbindlist(idcol = "FD.metric") %>% # rbind the list with FD.metric id column
  mutate(system = "Kinneret", res = "Month",troph = "Phytoplankton") %>% # specify metadata for future plotting
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*",""))) # assess significance of observed cross correlation by comparing to 2.5 and 97.5 quartiles (two tailed)
kin.phytomth.diff.raw <- lapply(kin.phytomth.diff, `[[`, 'perm.dens')%>% # extract second level list elements (i.e. 'perm.dens')
  data.table::rbindlist(idcol = "FD.metric")%>% # rbind the list with FD.metric id column
  mutate(system = "Kinneret", res = "Month",troph = "Phytoplankton") # specify metadata for future plotting

kin.zoomth.diff<- pbmcapply::pbmclapply(c("zooFDis","zooFEve","zooFRic"),function(x){
  pc <-  suppressWarnings(diff.perm.ccf(dat = kin.tot[,c("date",paste(x),"community")],
                                        iter = 10000,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  bio <- suppressWarnings(diff.perm.ccf(dat = kin.tot[,c("date",paste(x),"density")],
                                        iter = 10000,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  fi <-  suppressWarnings(diff.perm.ccf(dat = kin.tot[,c("date",paste(x),"FI")],
                                        iter = 10000,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  mvi <-  suppressWarnings(diff.perm.ccf(dat = kin.tot[,c("date",paste(x),"mvi")],
                                         iter = 10000,span =12*5,lag=1,
                                         perm.method = "red.noise",identical.t = F,
                                         detrend.method = "lm"))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(dat = kin.tot[,c("date",paste(x),"zp.ratio")],
                                              iter = 10000,span =12*5,lag=1,
                                              perm.method = "red.noise",identical.t = F,
                                              detrend.method = "lm"))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(kin.zoomth.diff) <- c("FDis","FEve","FRic")
kin.zoomth.diff.summary <- lapply(kin.zoomth.diff, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Kinneret", res = "Month",troph = "Zooplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
kin.zoomth.diff.raw <- lapply(kin.zoomth.diff, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Kinneret", res = "Month",troph = "Zooplankton")

# Mendota Phytoplankton and Zooplankton #

mad.phytomth.diff<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <-  suppressWarnings(diff.perm.ccf(dat = mad.tot[,c("date",paste(x),"community")],
                                        iter = 10000,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  bio <- suppressWarnings(diff.perm.ccf(dat = mad.tot[,c("date",paste(x),"density")],
                                        iter = 10000,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  fi <-  suppressWarnings(diff.perm.ccf(dat = mad.tot[,c("date",paste(x),"FI")],
                                        iter = 10000,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  mvi <-  suppressWarnings(diff.perm.ccf(dat = mad.tot[,c("date",paste(x),"mvi")],
                                         iter = 10000,span =12*5,lag=1,
                                         perm.method = "red.noise",identical.t = F,
                                         detrend.method = "lm"))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(dat = mad.tot[,c("date",paste(x),"zp.ratio")],
                                              iter = 10000,span =12*5,lag=1,
                                              perm.method = "red.noise",identical.t = F,
                                              detrend.method = "lm"))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(mad.phytomth.diff) <- c("FDis","FEve","FRic")
mad.phytomth.diff.summary <- lapply(mad.phytomth.diff, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Mendota", res = "Month",troph = "Phytoplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
mad.phytomth.diff.raw <- lapply(mad.phytomth.diff, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Mendota", res = "Month",troph = "Phytoplankton")

mad.zoomth.diff<- pbmcapply::pbmclapply(c("zooFDis","zooFEve","zooFRic"),function(x){
  pc <-  suppressWarnings(diff.perm.ccf(dat = mad.tot[,c("date",paste(x),"community")],
                                        iter = 10000,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  bio <- suppressWarnings(diff.perm.ccf(dat = mad.tot[,c("date",paste(x),"density")],
                                        iter = 10000,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  fi <-  suppressWarnings(diff.perm.ccf(dat = mad.tot[,c("date",paste(x),"FI")],
                                        iter = 10000,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  mvi <-  suppressWarnings(diff.perm.ccf(dat = mad.tot[,c("date",paste(x),"mvi")],
                                         iter = 10000,span =12*5,lag=1,
                                         perm.method = "red.noise",identical.t = F,
                                         detrend.method = "lm"))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(dat = mad.tot[,c("date",paste(x),"zp.ratio")],
                                              iter = 10000,span =12*5,lag=1,
                                              perm.method = "red.noise",identical.t = F,
                                              detrend.method = "lm"))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(mad.zoomth.diff) <- c("FDis","FEve","FRic")
mad.zoomth.diff.summary <- lapply(mad.zoomth.diff, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Mendota", res = "Month",troph = "Zooplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
mad.zoomth.diff.raw <- lapply(mad.zoomth.diff, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Mendota", res = "Month",troph = "Zooplankton")

# Lower Zurich Phytoplankton and Zooplankton #

LZ.phytomth.diff<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <-  suppressWarnings(diff.perm.ccf(dat = LZ.tot[,c("date",paste(x),"community")],
                                        iter = 10000,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  bio <- suppressWarnings(diff.perm.ccf(dat = LZ.tot[,c("date",paste(x),"density")],
                                        iter = 10000,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  fi <-  suppressWarnings(diff.perm.ccf(dat = LZ.tot[,c("date",paste(x),"FI")],
                                        iter = 10000,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  mvi <-  suppressWarnings(diff.perm.ccf(dat = LZ.tot[,c("date",paste(x),"mvi")],
                                         iter = 10000,span =12*5,lag=1,
                                         perm.method = "red.noise",identical.t = F,
                                         detrend.method = "lm"))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(dat = LZ.tot[,c("date",paste(x),"zp.ratio")],
                                              iter = 10000,span =12*5,lag=1,
                                              perm.method = "red.noise",identical.t = F,
                                              detrend.method = "lm"))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(LZ.phytomth.diff) <- c("FDis","FEve","FRic")
LZ.phytomth.diff.summary <- lapply(LZ.phytomth.diff, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Lower Zurich", res = "Month",troph = "Phytoplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
LZ.phytomth.diff.raw <- lapply(LZ.phytomth.diff, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Lower Zurich", res = "Month",troph = "Phytoplankton")

LZ.zoomth.diff<- pbmcapply::pbmclapply(c("zooFDis","zooFEve","zooFRic"),function(x){
  pc <-  suppressWarnings(diff.perm.ccf(dat = LZ.tot[,c("date",paste(x),"community")],
                                        iter = 10000,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  bio <- suppressWarnings(diff.perm.ccf(dat = LZ.tot[,c("date",paste(x),"density")],
                                        iter = 10000,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  fi <-  suppressWarnings(diff.perm.ccf(dat = LZ.tot[,c("date",paste(x),"FI")],
                                        iter = 10000,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  mvi <-  suppressWarnings(diff.perm.ccf(dat = LZ.tot[,c("date",paste(x),"mvi")],
                                         iter = 10000,span =12*5,lag=1,
                                         perm.method = "red.noise",identical.t = F,
                                         detrend.method = "lm"))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(dat = LZ.tot[,c("date",paste(x),"zp.ratio")],
                                              iter = 10000,span =12*5,lag=1,
                                              perm.method = "red.noise",identical.t = F,
                                              detrend.method = "lm"))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(LZ.zoomth.diff) <- c("FDis","FEve","FRic")
LZ.zoomth.diff.summary <- lapply(LZ.zoomth.diff, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Lower Zurich", res = "Month",troph = "Zooplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
LZ.zoomth.diff.raw <- lapply(LZ.zoomth.diff, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Lower Zurich", res = "Month",troph = "Zooplankton")

# Windermere Phytoplankton #

wind.phytomth.diff<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <-  suppressWarnings(diff.perm.ccf(dat = wind.tot[,c("date",paste(x),"community")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  bio <- suppressWarnings(diff.perm.ccf(dat = wind.tot[,c("date",paste(x),"density")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  fi <-  suppressWarnings(diff.perm.ccf(dat = wind.tot[,c("date",paste(x),"FI")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  mvi <-  suppressWarnings(diff.perm.ccf(dat = wind.tot[,c("date",paste(x),"mvi")],
                                         iter = 999,span =12*5,lag=1,
                                         perm.method = "red.noise",identical.t = F,
                                         detrend.method = "lm"))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(dat = wind.tot[,c("date",paste(x),"zp.ratio")],
                                              iter = 999,span =12*5,lag=1,
                                              perm.method = "red.noise",identical.t = F,
                                              detrend.method = "lm"))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(wind.phytomth.diff) <- c("FDis","FEve","FRic")
wind.phytomth.diff.summary <- lapply(wind.phytomth.diff, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Windermere", res = "Month",troph = "Phytoplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
wind.phytomth.diff.raw <- lapply(wind.phytomth.diff, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Windermere", res = "Month",troph = "Phytoplankton")

# Windermere Zooplankton #

wind.zoomth.diff<- pbmcapply::pbmclapply(c("zooFDis","zooFEve","zooFRic"),function(x){
  pc <-  suppressWarnings(diff.perm.ccf(dat = wind.tot[,c("date",paste(x),"community")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  bio <- suppressWarnings(diff.perm.ccf(dat = wind.tot[,c("date",paste(x),"density")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  fi <-  suppressWarnings(diff.perm.ccf(dat = wind.tot[,c("date",paste(x),"FI")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  mvi <-  suppressWarnings(diff.perm.ccf(dat = wind.tot[,c("date",paste(x),"mvi")],
                                         iter = 999,span =12*5,lag=1,
                                         perm.method = "red.noise",identical.t = F,
                                         detrend.method = "lm"))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(dat = wind.tot[,c("date",paste(x),"zp.ratio")],
                                              iter = 999,span =12*5,lag=1,
                                              perm.method = "red.noise",identical.t = F,
                                              detrend.method = "lm"))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(wind.zoomth.diff) <- c("FDis","FEve","FRic")
wind.zoomth.diff.summary <- lapply(wind.zoomth.diff, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Windermere", res = "Month",troph = "Zooplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
wind.zoomth.diff.raw <- lapply(wind.zoomth.diff, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Windermere", res = "Month",troph = "Zooplankton")

# Kasumigaura Phytoplankton #

kas.phytomth.diff<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <-  suppressWarnings(diff.perm.ccf(dat = kas.tot[,c("date",paste(x),"community")],
                                        iter = 10000,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  bio <- suppressWarnings(diff.perm.ccf(dat = kas.tot[,c("date",paste(x),"density")],
                                        iter = 10000,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  fi <-  suppressWarnings(diff.perm.ccf(dat = kas.tot[,c("date",paste(x),"FI")],
                                        iter = 10000,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  mvi <-  suppressWarnings(diff.perm.ccf(dat = kas.tot[,c("date",paste(x),"mvi")],
                                         iter = 10000,span =12*5,lag=1,
                                         perm.method = "red.noise",identical.t = F,
                                         detrend.method = "lm"))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(dat = kas.tot[,c("date",paste(x),"zp.ratio")],
                                              iter = 10000,span =12*5,lag=1,
                                              perm.method = "red.noise",identical.t = F,
                                              detrend.method = "lm"))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(kas.phytomth.diff) <- c("FDis","FEve","FRic")
kas.phytomth.diff.summary <- lapply(kas.phytomth.diff, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Kasumigaura", res = "Month",troph = "Phytoplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
kas.phytomth.diff.raw <- lapply(kas.phytomth.diff, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Kasumigaura", res = "Month",troph = "Phytoplankton")

# Kasumigaura Zooplankton #

kas.zoomth.diff<- pbmcapply::pbmclapply(c("zooFDis","zooFEve","zooFRic"),function(x){
  pc <-  suppressWarnings(diff.perm.ccf(dat = kas.tot[,c("date",paste(x),"community")],
                                        iter = 10000,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  bio <- suppressWarnings(diff.perm.ccf(dat = kas.tot[,c("date",paste(x),"density")],
                                        iter = 10000,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  fi <-  suppressWarnings(diff.perm.ccf(dat = kas.tot[,c("date",paste(x),"FI")],
                                        iter = 10000,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "lm"))
  mvi <-  suppressWarnings(diff.perm.ccf(dat = kas.tot[,c("date",paste(x),"mvi")],
                                         iter = 10000,span =12*5,lag=1,
                                         perm.method = "red.noise",identical.t = F,
                                         detrend.method = "lm"))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(dat = kas.tot[,c("date",paste(x),"zp.ratio")],
                                              iter = 10000,span =12*5,lag=1,
                                              perm.method = "red.noise",identical.t = F,
                                              detrend.method = "lm"))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(kas.zoomth.diff) <- c("FDis","FEve","FRic")
kas.zoomth.diff.summary <- lapply(kas.zoomth.diff, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Kasumigaura", res = "Month",troph = "Zooplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
kas.zoomth.diff.raw <- lapply(kas.zoomth.diff, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Kasumigaura", res = "Month",troph = "Zooplankton")

###########################################################################
## Save out ##
###########################################################################
summary.ccf<- rbind(kin.phytomth.diff.summary,kin.zoomth.diff.summary,mad.phytomth.diff.summary,mad.zoomth.diff.summary,
                         LZ.phytomth.diff.summary,LZ.zoomth.diff.summary,wind.phytomth.diff.summary,wind.zoomth.diff.summary,
                         kas.phytomth.diff.summary,kas.zoomth.diff.summary)
raw.ccf <- rbind(kin.phytomth.diff.raw,kin.zoomth.diff.raw,mad.phytomth.diff.raw,mad.zoomth.diff.raw,
                     LZ.phytomth.diff.raw,LZ.zoomth.diff.raw,wind.phytomth.diff.raw,wind.zoomth.diff.raw,
                     kas.phytomth.diff.raw,kas.zoomth.diff.raw)

write.csv(summary.ccf,file ="Results/ccf/raw_data/summary.ccf.csv",row.names = F)
save(raw.ccf,file = "Results/ccf/raw_data/raw.ccf.RData") # RData required to reduce file size compared to .csv
load(file = "Results/ccf/raw_data/raw.ccf.RData")
summary.ccf <- read.csv("Results/ccf/raw_data/summary.ccf.csv")
