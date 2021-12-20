### Cross Correlate System State with FD ###

## Preamble ##
require(tidyverse) # dplyr, ggplot etc.
require(pbmcapply) # paralled lapply 
require(data.table) # rbindlist function
require(patchwork) # plot alignment

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
## Estimate cross correlation and permute (lm, Monthly) ##
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
## Save out (lm, Monthly) ##
###########################################################################
summary.ccf.mth1 <- rbind(kin.phytomth.diff.summary,kin.zoomth.diff.summary,mad.phytomth.diff.summary,mad.zoomth.diff.summary,
                         LZ.phytomth.diff.summary,LZ.zoomth.diff.summary,wind.phytomth.diff.summary,wind.zoomth.diff.summary,
                         kas.phytomth.diff.summary,kas.zoomth.diff.summary)
raw.ccf.mth1 <- rbind(kin.phytomth.diff.raw,kin.zoomth.diff.raw,mad.phytomth.diff.raw,mad.zoomth.diff.raw,
                     LZ.phytomth.diff.raw,LZ.zoomth.diff.raw,wind.phytomth.diff.raw,wind.zoomth.diff.raw,
                     kas.phytomth.diff.raw,kas.zoomth.diff.raw)

write.csv(summary.ccf.mth1,file ="Results/ccf/raw_data/summary.ccf.mth.lag1.csv",row.names = F)
save(raw.ccf.mth1,file = "Results/ccf/raw_data/raw.ccf.mth.lag1.RData") # RData required to reduce file size compared to .csv
load(file = "Results/ccf/raw_data/raw.ccf.mth.lag1.RData")
summary.ccf.mth1 <- read.csv("Results/ccf/raw_data/summary.ccf.mth.lag1.csv")

pdf(file="Results/ccf/FD_perm_lm_mth_absrmax.pdf",
    width=10, height = 8)
ggplot(raw.ccf.mth1,aes(x = state.metric, y =  abs.rmax, col = FD.metric,fill= FD.metric)) + 
  geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.05, 0.5, 0.95),scale = "width",alpha = 0.3) +
  theme_bw() + 
  geom_point(data = summary.ccf.mth1[summary.ccf.mth1$measure %in% "absmax.ccf",],
             aes(x = state.metric, y = obs.value),position = position_dodge(width = 0.9),size=2) +
  geom_text(data = summary.ccf.mth1[summary.ccf.mth1$measure %in% "absmax.ccf",], 
            aes(x = state.metric, y = 0.45,label = sig),col= "black",size = 4,position = position_dodge(width = 0.9))+
  geom_text(data = summary.ccf.mth1[summary.ccf.mth1$measure %in% "t.absmax.ccf",], 
            aes(x = state.metric, y = 0.4,label = obs.value),col= "black",size = 3,position = position_dodge(width = 0.9))+
  #facet_wrap(~troph,nrow=2,strip.position = "right")+
  facet_grid(system~troph)+
  scale_y_continuous(limits = c(-0.9,0.9),breaks = seq(-0.75,0.75,0.5))+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  ylab("Cross correlation") + xlab("System state proxy")+   ggtitle("Permuted CCF between detrended FD and system state")
dev.off()

pdf(file="Results/ccf/FD_perm_lm_mth_r0.pdf",
    width=10, height = 8)
ggplot(raw.ccf.mth1,aes(x = state.metric, y =  r0, col = FD.metric,fill= FD.metric)) + 
  geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.025, 0.5, 0.975),scale = "width",alpha = 0.3) +
  theme_bw() + 
  geom_point(data = summary.ccf.mth1[summary.ccf.mth1$measure %in% "r0.ccf",],
             aes(x = state.metric, y = obs.value),position = position_dodge(width = 0.9),size=2) +
  geom_text(data = summary.ccf.mth1[summary.ccf.mth1$measure %in% "r0.ccf",], 
            aes(x = state.metric, y = 0.35,label = sig),col= "black",size = 4,position = position_dodge(width = 0.9))+
  facet_grid(system~troph)+
  scale_y_continuous(limits = c(-0.9,0.9),breaks = seq(-0.75,0.75,0.5))+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  ylab("Cross correlation") + xlab("System state proxy")+   ggtitle("Permuted CCF t0 between detrended FD and system state")
dev.off()

count.ccf.r0.dat <- filter(summary.ccf.mth1,measure %in% "r0.ccf") %>%
  group_by(state.metric,FD.metric,troph)%>%
  summarise(NSig = sum(sig=="*"),
            NxSig= sum(sig!="*"))%>%
  mutate(ref.y = ifelse(FD.metric %in% "FDis",0.7, 
                        ifelse(FD.metric %in% "FEve",0.65,0.6))) #significant count per group

pdf(file="Results/ccf/summary_FD_perm_lm_mth_r0.pdf",
    width=8, height = 4)
pccf.lag0.fin <- ggplot(filter(summary.ccf.mth1,measure %in% "r0.ccf"),aes(x=state.metric,y=obs.value,col=FD.metric))+
  #geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.025, 0.5, 0.975),scale = "width",alpha = 0.3)+
  geom_hline(yintercept = 0,col="black",alpha = 0.3)+
  geom_point(position=position_dodge(width=0.75),
             aes(shape=system,alpha=sig,fill=FD.metric,group=FD.metric),size=3.5)+
  geom_point(position=position_dodge(width=0.75),
             aes(shape=system,fill=NULL,group=FD.metric),size=3.5) +
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric")+
  scale_alpha_manual(values=c(0.01,1),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  # geom_text(data =count.ccf.r0.dat,aes(x = state.metric, y = ref.y,label = paste("(",NSig,"*",",",NxSig,")",sep = ""),fill=FD.metric,group = FD.metric),
  #           col= "black",size = 3, position = position_dodge(width = 0.9))+
  geom_boxplot(aes(fill=FD.metric),alpha=0.1,col="black",size=0.3,outlier.shape = NA)+
  scale_y_continuous(limits = c(-0.9,0.9),breaks = seq(-0.75,0.75,0.5))+
  facet_wrap(~troph)+
  ylab("Cross correlation") + xlab("System state proxy")+
  guides(size = guide_legend(order = 1), 
         col = guide_legend(order = 2),
         fill = guide_legend(order = 2),
         shape = guide_legend(order = 3))+
  theme_bw()
pccf.lag0.fin
dev.off()

count.ccf.absmax.dat <- filter(summary.ccf.mth1,measure %in% "absmax.ccf") %>%
  group_by(state.metric,FD.metric,troph)%>%
  summarise(NSig = sum(sig=="*"),
            NxSig= sum(sig!="*"))%>%
  mutate(ref.y = ifelse(FD.metric %in% "FDis",0.7, 
                        ifelse(FD.metric %in% "FEve",0.65,0.6))) #significant count per group

ccf.lag1 <- ggplot(filter(summary.ccf.mth1,measure %in% "absmax.ccf"),
                   aes(x=state.metric,y=obs.value,col=FD.metric))+
  geom_hline(yintercept = 0,col="black",alpha = 0.3)+
  geom_point(position=position_dodge(width=0.75),
             aes(shape=system,alpha=sig,fill=FD.metric,group=FD.metric),size=3.5)+
  geom_point(position=position_dodge(width=0.75),
             aes(shape=system,fill=NULL,group=FD.metric),size=3.5) +
  geom_boxplot(aes(fill=FD.metric),alpha=0.1,col="black",size=0.3,outlier.shape = NA)+
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric")+
  scale_alpha_manual(values=c(0.01,1),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  # geom_text(data =count.ccf.absmax.dat,aes(x = state.metric, y = ref.y,label = paste("(",NSig,"*",",",NxSig,")",sep = ""),fill=FD.metric,group = FD.metric),
  #           col= "black",size = 3, position = position_dodge(width = 0.9))+
  scale_y_continuous(limits = c(-0.9,0.9),breaks = seq(-0.75,0.75,0.5))+
  facet_wrap(~troph)+
  ylab("Cross correlation") + xlab("System state proxy")+
  theme_bw() +
  guides(col = guide_legend(order = 1),
         fill = guide_legend(order = 1),
         shape = guide_legend(order = 2))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(), 
        panel.background = element_blank(),
        plot.margin = margin(c(2, 2, 0, 2)))

ccf.lag2 <- ggplot(filter(summary.ccf.mth1,measure %in% "t.absmax.ccf") %>% 
         mutate(sig = as.factor(filter(summary.ccf.mth1,measure %in% "absmax.ccf")$sig)),
       aes(x=state.metric,y=obs.value,col = FD.metric))+
  geom_hline(yintercept = 0,col="black",alpha = 0.3)+
  geom_point(position=position_dodge(width=0.75),
             aes(shape=system,alpha=sig,fill=FD.metric,group=FD.metric),size=2)+
  geom_point(position=position_dodge(width=0.75),
             aes(shape=system,fill=NULL,group=FD.metric),size=2) +
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake",guide = "none")+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric",guide = "none") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric",guide = "none")+
  scale_alpha_manual(values=c(0.01,1),name = "Significance", labels = c("Not significant","Significant"),
                     #guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
                     guide = "none")+
  #scale_x_discrete(position = "top") +
  scale_y_reverse()+
  facet_wrap(~troph)+
  ylab("Optimal lag") + xlab("System state proxy")+
  theme_bw()+
  #guides(shape = guide_legend(order = 1))+
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = margin(c(2, 2, 0, 2)))

ccf.lag2 <- ggplot(filter(summary.ccf.mth1,measure %in% "t.absmax.ccf") %>% 
  mutate(sig = as.factor(filter(summary.ccf.mth1,measure %in% "absmax.ccf")$sig)),
    aes(x=state.metric,y=obs.value,col = FD.metric)) + 
  geom_hline(yintercept = 0,col="black",alpha = 0.3)+
  scale_y_binned(breaks = c(seq(60,24,-12),12,seq(-12,-60,-12)),show.limits = T)+
  geom_tile(data = expand_grid(c(seq(60,24,-12),12,seq(-12,-60,-12)),unique(summary.ccf.mth1$FD.metric),unique(summary.ccf.mth1$state.metric),unique(summary.ccf.mth1$troph)) %>%
              magrittr::set_colnames(c("obs.value","FD.metric","state.metric","troph")),
            aes(group = FD.metric),fill = "white",stat="identity",position = position_dodge(width = 0.75), col = "black", size = 0.3,width = 0.8, height = 0.9)+
  geom_tile(data = expand_grid(c(12),unique(summary.ccf.mth1$FD.metric),unique(summary.ccf.mth1$state.metric),unique(summary.ccf.mth1$troph)) %>%
              magrittr::set_colnames(c("obs.value","FD.metric","state.metric","troph")),
            aes(group = FD.metric),fill = "#E5E0FF",stat="identity",position = position_dodge(width = 0.75), col = "black", size = 0.3,width = 0.8, height = 0.9)+
  geom_point(position=position_jitterdodge(dodge.width=0.75,jitter.height = 0.2,jitter.width = 0,seed = 5),
             aes(y = obs.value,shape=system,alpha=sig,col = FD.metric,fill=FD.metric,group=FD.metric),size=1.5)+
  geom_point(position=position_jitterdodge(dodge.width=0.75,jitter.height = 0.2,jitter.width = 0,seed = 5),
             aes(y = obs.value,shape=system,col = FD.metric,fill=NULL,group=FD.metric),size=1.5)+
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake",guide = "none")+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric",guide = "none") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric",guide = "none")+
  scale_alpha_manual(values=c(0.01,1),name = "Significance", labels = c("Not significant","Significant"),
                     guide = "none")+
  facet_wrap(~troph)+
  #facet_grid(troph~FD.metric)+
  ylab("Optimal lag (months)") + xlab("System state proxy")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = margin(c(2, 2, 0, 2)))


pdf(file="Results/ccf/summary_FD_perm_lm_mth_absrmax.pdf",
    width=10, height = 6)
p1 <- ccf.lag1 + ccf.lag2 + patchwork::plot_layout(nrow = 2,guides = "collect",heights = c(2, 1),tag_level = 'new')
p1
dev.off()

pdf(file="Results/ccf/summary_FD_perm_lm_mth_combo.pdf",
    width=15, height = 6)  
pccf.lag0.fin + p1 + 
  patchwork::plot_layout(ncol = 2,nrow = 1, guides = "collect") +
  patchwork::plot_annotation(tag_levels = c('a',1)) & 
  theme(plot.tag = element_text(face = "bold"))
dev.off()

###########################################################################
## Summary correlations (lm, Monthly) ##
###########################################################################
summary.ccf.mth1 <- read.csv("Results/ccf/raw_data/summary.ccf.mth.lag1.csv")

# obs.cor.lag0.lake.tab <- summary.ccf.mth1 %>%
#             filter(measure %in% "r0.ccf")%>%
#             group_by(system,troph) %>%
#             summarise(mean.cor = mean(obs.value),median.cor = median(obs.value),
#                       cor.se = sd(obs.value)/n(),
#               nsig=sum(sig %in% "*"),prop.sig = sum(sig %in% "*")/length(sig)) %>%
#             mutate(across(mean.cor:cor.se,~round(.x,digits=4)))
# write.csv(obs.cor.lag0.lake.tab,file ="Results/ccf/ccf_tables/cor.lag0.lake.tab.csv",row.names = F)
# 
# obs.cor.lag0.FD.tab <- summary.ccf.mth1 %>%
#   filter(measure %in% "r0.ccf")%>%
#   group_by(troph,FD.metric) %>% 
#   summarise(mean.cor = mean(obs.value),median.cor = median(obs.value),
#             cor.se = sd(obs.value)/n(),
#             nsig=sum(sig %in% "*"),prop.sig = sum(sig %in% "*")/length(sig)) %>%
#   mutate(across(mean.cor:cor.se,~round(.x,digits=4)))
# write.csv(obs.cor.lag0.FD.tab,file ="Results/ccf/ccf_tables/cor.lag0.FD.tab.csv",row.names = F)

obs.cor.lag0.state.tab <- summary.ccf.mth1 %>%
  filter(measure %in% "r0.ccf")%>%
  group_by(troph,FD.metric,state.metric) %>%
  summarise(mean.cor = mean(obs.value),median.cor = median(obs.value),
            cor.se = sd(obs.value)/n(),
            nsig=sum(sig %in% "*"),prop.sig = sum(sig %in% "*")/length(sig)) %>%
  mutate(across(mean.cor:cor.se,~round(.x,digits=4)))
write.csv(obs.cor.lag0.state.tab,file ="Results/ccf/ccf_tables/cor.lag0.state.tab.csv",row.names = F)
obs.cor.lag0.state.tab <- read.csv("Results/ccf/ccf_tables/cor.lag0.state.tab.csv")

pdf(file="Results/ccf/summary_FD_perm_lm_mth_r0_alt.pdf",
    width=8, height = 5)  
pccf.lag0.1 <- ggplot(filter(summary.ccf.mth1,measure %in% "r0.ccf"),aes(x=state.metric,y=obs.value,col=FD.metric))+
  #geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.025, 0.5, 0.975),scale = "width",alpha = 0.3)+
  #geom_hline(yintercept = 0,col="black",alpha = 0.3)+
  geom_line(aes(linetype = system),position = position_dodge(width = 0.75), 
            alpha = 0, show.legend = FALSE)+
  facet_wrap(~troph)

pccf.lag0.fin <- pccf.lag0.1 + geom_segment(data = layer_data(pccf.lag0.1, 1L),
                  aes(x = xmin, xend=xmax, y = 0,yend=0, group = linetype),
                  color = "black", size = 0.5)+
  #geom_linerange(aes(xmin=state.metric,xmax=state.metric),position = position_dodge(width=0.75))+
  geom_point(position=position_dodge(width=0.75),
             aes(shape=system,alpha=sig,fill=FD.metric,group=FD.metric),size=3.5)+
  geom_point(position=position_dodge(width=0.75),
             aes(shape=system,fill=NULL,group=FD.metric),size=3.5)+
  # geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width=0.75,seed=15),
  #            aes(shape=system,alpha=sig,fill=FD.metric,group=FD.metric),size=2)+
  # geom_point(position=position_jitterdodge(jitter.width = 0.1,dodge.width=0.75,seed=15),
  #            aes(shape=system,fill=NULL,group=FD.metric),size=2)+
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric")+
  scale_alpha_manual(values=c(0.01,1),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  geom_point(data = obs.cor.lag0.state.tab,aes(x=state.metric,y=mean.cor,group=FD.metric,size="mean"),position=position_dodge(width=0.75),col="black",alpha=1,shape = 21,fill="black")+
  #geom_segment(data = obs.cor.lag0.state.tab,aes(x=state.metric,y=0,xend=state.metric,yend=mean.cor,group=FD.metric),position=position_dodge(width=0.9),col="black")+
  geom_linerange(data = obs.cor.lag0.state.tab,aes(x=state.metric,y = 0,ymin = 0,ymax=mean.cor,group=FD.metric),position=position_dodge(width=0.75),col="black",linetype = "solid")+
  #geom_errorbar(data = obs.cor.lag0.state.tab,aes(x=state.metric,y = 0,ymin = 0,ymax=mean.cor,group=FD.metric), width=0.8,
  #             position=position_dodge(0.75),col="black",size=0.5) +
  # geom_text(data =count.ccf.r0.dat,aes(x = state.metric, y = ref.y,label = paste("(",NSig,"*",",",NxSig,")",sep = ""),fill=FD.metric,group = FD.metric),
  #           col= "black",size = 3, position = position_dodge(width = 0.9))+
  scale_y_continuous(limits = c(-0.9,0.9),breaks = seq(-0.75,0.75,0.5))+
  scale_size_manual(values = c(2.5),breaks = c("mean"),name = NULL, labels = c("Mean cross\ncorrelation"),
                    guide = guide_legend(override.aes = list(fill = c("black"),shape=c(21),size=c(3.5),alpha = c(1))))+
  ylab("Cross correlation") + xlab("System state proxy")+
  guides(size = guide_legend(order = 1), 
         col = guide_legend(order = 2),
         fill = guide_legend(order = 2),
         shape = guide_legend(order = 3))+
  theme_bw()

pccf.lag0.fin
dev.off()

# obs.cor.lagx.lake.tab <- summary.ccf.mth1 %>%
#   filter(measure %in% c("absmax.ccf","t.absmax.ccf"))%>%
#   pivot_wider(names_from = measure,values_from = obs.value)%>%
#   ungroup()%>%
#   mutate(t.absmax.ccf=ifelse(is.numeric(absmax.ccf) & !is.na(t.absmax.ccf),t.absmax.ccf,
#                              dplyr::lead(t.absmax.ccf)))%>% # fill NA t.absmax with next t.absmax to associate cor with lag
#   na.omit() %>%#drop duplicate rows
#   group_by(system,troph) %>%
#   summarise(mean.cor = mean(absmax.ccf),median.cor = median(absmax.ccf),
#             cor.se = sd(absmax.ccf)/n(), 
#             median.lag = median(t.absmax.ccf),lag.se = sd(t.absmax.ccf)/n(), 
#             nsig=sum(sig %in% "*"),prop.sig = sum(sig %in% "*")/length(sig)) %>%
#   mutate(across(mean.cor:lag.se,~round(.x,digits=4)))
# write.csv(obs.cor.lagx.lake.tab,file ="Results/ccf/ccf_tables/cor.lagx.lake.tab.csv",row.names = F)
# 
# obs.cor.lagx.FD.tab <- summary.ccf.mth1 %>%
#   filter(measure %in% c("absmax.ccf","t.absmax.ccf"))%>%
#   pivot_wider(names_from = measure,values_from = obs.value)%>%
#   ungroup()%>%
#   mutate(t.absmax.ccf=ifelse(is.numeric(absmax.ccf) & !is.na(t.absmax.ccf),t.absmax.ccf,
#                              dplyr::lead(t.absmax.ccf)))%>% # fill NA t.absmax with next t.absmax to associate cor with lag
#   na.omit() %>%#drop duplicate rows
#   group_by(system,troph,FD.metric) %>%
#   summarise(mean.cor = mean(absmax.ccf),median.cor = median(absmax.ccf),
#             cor.se = sd(absmax.ccf)/n(), 
#             median.lag = median(t.absmax.ccf),lag.se = sd(t.absmax.ccf)/n(), 
#             nsig=sum(sig %in% "*"),prop.sig = sum(sig %in% "*")/length(sig)) %>%
#   mutate(across(mean.cor:lag.se,~round(.x,digits=4)))
# write.csv(obs.cor.lagx.FD.tab,file ="Results/ccf/ccf_tables/cor.lagx.FD.tab.csv",row.names = F)

obs.cor.lagx.state.tab <- summary.ccf.mth1 %>%
  dplyr::select(!c(quantile,median.perm.value,obs.difference,res))%>%
  filter(measure %in% c("absmax.ccf","t.absmax.ccf"))%>%
  pivot_wider(names_from = measure,values_from = obs.value)%>%
  ungroup()%>%
  mutate(t.absmax.ccf=ifelse(is.numeric(absmax.ccf) & !is.na(t.absmax.ccf),t.absmax.ccf,
                              dplyr::lead(t.absmax.ccf)))%>% # fill NA t.absmax with next t.absmax to associate cor with lag
  na.omit() %>%#drop duplicate rows
  group_by(troph,FD.metric,state.metric) %>%
  summarise(mean.cor = mean(absmax.ccf),median.cor = median(absmax.ccf),
            cor.se = sd(absmax.ccf)/n(), 
            median.lag = median(t.absmax.ccf),lag.se = sd(t.absmax.ccf)/n(), 
            nsig=sum(sig %in% "*"),prop.sig = sum(sig %in% "*")/length(sig)) %>%
  mutate(across(mean.cor:lag.se,~round(.x,digits=4)))
write.csv(obs.cor.lagx.state.tab,file ="Results/ccf/ccf_tables/cor.lagx.state.tab.csv",row.names = F)

pdf(file="Results/ccf/summary_FD_perm_lm_mth_absrmax_alt.pdf",
    width=10, height = 6)
pccf.lagx.1 <- ggplot(filter(summary.ccf.mth1,measure %in% "absmax.ccf"),aes(x=state.metric,y=obs.value,col=FD.metric))+
  geom_line(aes(linetype = system),position = position_dodge(width = 0.75), 
            alpha = 0, show.legend = FALSE)+
  facet_wrap(~troph)

pccf.lagx.2 <- pccf.lagx.1 +
  geom_segment(data = layer_data(pccf.lagx.1, 1L),
               aes(x = xmin, xend=xmax, y = 0,yend=0, group = linetype),
               color = "black", size = 0.5)+
  geom_point(position=position_dodge(width=0.75),
             aes(shape=system,alpha=sig,fill=FD.metric,group=FD.metric),size=3.5)+
  geom_point(position=position_dodge(width=0.75),
             aes(shape=system,fill=NULL,group=FD.metric),size=3.5)+
  # geom_point(position=position_jitterdodge(jitter.width = 0.1, dodge.width=0.75,seed=15),
  #            aes(shape=system,alpha=sig,fill=FD.metric,group=FD.metric),size=2)+
  # geom_point(position=position_jitterdodge(jitter.width = 0.1,dodge.width=0.75,seed=15),
  #            aes(shape=system,fill=NULL,group=FD.metric),size=2)+
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric")+
  scale_alpha_manual(values=c(0.01,1),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  geom_point(data = obs.cor.lagx.state.tab,aes(x=state.metric,y=mean.cor,group=FD.metric,size="mean"),position=position_dodge(width=0.75),col="black",alpha=1,shape = 21,fill="black")+
  #geom_segment(data = obs.cor.lag0.state.tab,aes(x=state.metric,y=0,xend=state.metric,yend=mean.cor,group=FD.metric),position=position_dodge(width=0.75),col="black")+
  geom_linerange(data = obs.cor.lagx.state.tab,aes(x=state.metric,y = 0,ymin = 0,ymax=mean.cor,group=FD.metric),position=position_dodge(width=0.75),col="black",linetype = "solid")+
  #geom_errorbar(data = cor.lagx.state.tab,aes(x=state.metric,y = 0,ymin = 0,ymax=mean.cor,group=FD.metric), width=0.8,
  #position=position_dodge(0.75),col="black",size=0.5) +
  # geom_text(data =count.ccf.r0.dat,aes(x = state.metric, y = ref.y,label = paste("(",NSig,"*",",",NxSig,")",sep = ""),fill=FD.metric,group = FD.metric),
  #           col= "black",size = 3, position = position_dodge(width = 0.9))+
  scale_size_manual(values = c(2.5),breaks = c("mean"),name = NULL, labels = c("Mean cross\ncorrelation"),
                    guide = guide_legend(override.aes = list(fill = c("black"),shape=c(21),size=c(3.5),alpha = c(1))))+
  scale_y_continuous(limits = c(-0.9,0.9),breaks = seq(-0.75,0.75,0.5))+
  facet_wrap(~troph)+
  ylab("Cross correlation") + xlab("System state proxy")+
  theme_bw() +
  guides(size = guide_legend(order = 1), 
         col = guide_legend(order = 2),
         fill = guide_legend(order = 2),
         shape = guide_legend(order = 3))+
  theme(axis.title.x=element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(c(2, 2, 0, 2)))

pccf.lagx.fin <- pccf.lagx.2 + ccf.lag2 +plot_layout(nrow = 2,guides = "collect",heights = c(2, 1))
pccf.lagx.fin
dev.off()

lag0.lagx.comp <- left_join(obs.cor.lag0.state.tab,obs.cor.lagx.state.tab,
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

###########################################################################
## Estimate cross correlation and permute (Lag12, Monthly) ##
###########################################################################
# repeated 'diff.perm.ccf' function specified as @iter = 10000, @perm.method = "red.noise",
# @scale = T, @normalise = F, @span =12*5, @identical.t = F, @diff=T and 
# @lag=12.
# If @comp.ts = FI or mvi, these comp.ts are pre-differenced due to the rolling window 
# approach complicating the generic nature of the 'diff.perm.ccf' function
# Results in two dataframes:
# $summary = observed absolute max correlation coef and the corresponding lag, observed correlation at lag0 and summary statistics
# $perm.dens = all permutations absolute max correlation coef and the corresponding lag plus the correlation at lag0

# Kinneret Phytoplankton and Zooplankton #

kin.phytomth.diff12<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = phyto.kin.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(phyto.kin.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "red.noise", 
                                       scale = T, normalise = F,span =12*5, identical.t = F,
                                       diff=T,lag = 12,
                                       comp.ts = all.system.states$kin.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = phyto.kin.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.kin.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         diff=T,lag = 12,
                                         comp.ts = log(all.system.states$kin.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = phyto.kin.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(phyto.kin.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "red.noise",
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        diff=T,lag = 12,
                                        comp.ts = base::diff(all.system.states$kin.mth$FI[12:552],lag = 12),
                                        comp.ts.timedat = seq_along(all.system.states$kin.mth$maxt[24:552]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = phyto.kin.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.kin.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         diff=T,lag = 12,
                                         comp.ts = base::diff(log(all.system.states$kin.mth$mvi)[12:552],lag = 12),
                                         comp.ts.timedat = seq_along(all.system.states$kin.mth$maxt[24:552]),
                                         pre.diff =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = phyto.kin.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(phyto.kin.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "red.noise",
                                              scale = T, normalise = F,span =12*5,identical.t = F,
                                              diff=T,lag = 12,
                                              comp.ts = log(all.system.states$kin.mth$zp.ratio)))
  
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),"state.metric" = c(rep("Community",7),rep("Density",7),rep("FI",7),rep("MVI",7),rep("Z_P.ratio",7)))
  out.dens <- data.frame(rbind(as.data.frame(pc$perm.dens),as.data.frame(bio$perm.dens),as.data.frame(fi$perm.dens),as.data.frame(mvi$perm.dens),as.data.frame(zp.ratio$perm.dens)),
                         "state.metric" = c(rep("Community",10000),rep("Density",10000),rep("FI",10000),rep("MVI",10000),rep("Z_P.ratio",10000)))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(kin.phytomth.diff12) <- c("FDis","FEve","FRic")
kin.phytomth.diff12.summary <- lapply(kin.phytomth.diff12, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Kinneret", res = "Month",troph = "Phytoplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
kin.phytomth.diff12.raw <- lapply(kin.phytomth.diff12, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Kinneret", res = "Month",troph = "Phytoplankton")

kin.zoomth.diff12<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = zoo.kin.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(zoo.kin.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "red.noise", lag=12,
                                       scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$kin.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = zoo.kin.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.kin.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$kin.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = zoo.kin.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(zoo.kin.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "red.noise",lag=12,
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = base::diff(all.system.states$kin.mth$FI[12:552],lag=12),
                                        comp.ts.timedat = seq_along(all.system.states$kin.mth$maxt[24:552]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = zoo.kin.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.kin.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = base::diff(log(all.system.states$kin.mth$mvi)[12:552],lag=12),
                                         comp.ts.timedat = seq_along(all.system.states$kin.mth$maxt[24:552]),
                                         pre.diff =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = zoo.kin.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(zoo.kin.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "red.noise",lag=12,
                                              scale = T, normalise = F,span =12*5,identical.t = F,
                                              comp.ts = log(all.system.states$kin.mth$zp.ratio)))
  
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),"state.metric" = c(rep("Community",7),rep("Density",7),rep("FI",7),rep("MVI",7),rep("Z_P.ratio",7)))
  out.dens <- data.frame(rbind(as.data.frame(pc$perm.dens),as.data.frame(bio$perm.dens),as.data.frame(fi$perm.dens),as.data.frame(mvi$perm.dens),as.data.frame(zp.ratio$perm.dens)),
                         "state.metric" = c(rep("Community",10000),rep("Density",10000),rep("FI",10000),rep("MVI",10000),rep("Z_P.ratio",10000)))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(kin.zoomth.diff12) <- c("FDis","FEve","FRic")
kin.zoomth.diff12.summary <- lapply(kin.zoomth.diff12, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Kinneret", res = "Month",troph = "Zooplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
kin.zoomth.diff12.raw <- lapply(kin.zoomth.diff12, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Kinneret", res = "Month",troph = "Zooplankton")

# Mendota Phytoplankton and Zooplankton #

mad.phytomth.diff12<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = phyto.mad.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(phyto.mad.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "red.noise", lag=12,
                                       scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$mad.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = phyto.mad.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.mad.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$mad.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = phyto.mad.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(phyto.mad.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "red.noise",lag=12,
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = base::diff(all.system.states$mad.mth$FI[12:283],lag=12),
                                        comp.ts.timedat = seq_along(all.system.states$mad.mth$maxt[24:283]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = phyto.mad.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.mad.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = base::diff(log(all.system.states$mad.mth$mvi)[12:283],lag=12),
                                         comp.ts.timedat = seq_along(all.system.states$mad.mth$maxt[24:283]),
                                         pre.diff  =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = phyto.mad.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(phyto.mad.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "red.noise",lag=12,
                                              scale = T, normalise = F,span =12*5,identical.t = F,
                                              comp.ts = log(all.system.states$mad.mth$zp.ratio)))
  
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),"state.metric" = c(rep("Community",7),rep("Density",7),rep("FI",7),rep("MVI",7),rep("Z_P.ratio",7)))
  out.dens <- data.frame(rbind(as.data.frame(pc$perm.dens),as.data.frame(bio$perm.dens),as.data.frame(fi$perm.dens),as.data.frame(mvi$perm.dens),as.data.frame(zp.ratio$perm.dens)),
                         "state.metric" = c(rep("Community",10000),rep("Density",10000),rep("FI",10000),rep("MVI",10000),rep("Z_P.ratio",10000)))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(mad.phytomth.diff12) <- c("FDis","FEve","FRic")
mad.phytomth.diff12.summary <- lapply(mad.phytomth.diff12, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Mendota", res = "Month",troph = "Phytoplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
mad.phytomth.diff12.raw <- lapply(mad.phytomth.diff12, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Mendota", res = "Month",troph = "Phytoplankton")

mad.zoomth.diff12<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = zoo.mad.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(zoo.mad.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "red.noise", lag=12,
                                       scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$mad.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = zoo.mad.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.mad.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$mad.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = zoo.mad.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(zoo.mad.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "red.noise",lag=12,
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = base::diff(all.system.states$mad.mth$FI[12:283],lag=12),
                                        comp.ts.timedat = seq_along(all.system.states$mad.mth$maxt[24:283]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = zoo.mad.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.mad.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = base::diff(log(all.system.states$mad.mth$mvi)[12:283],lag=12),
                                         comp.ts.timedat = seq_along(all.system.states$mad.mth$maxt[24:283]),
                                         pre.diff =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = zoo.mad.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(zoo.mad.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "red.noise",lag=12,
                                              scale = T, normalise = F,span =12*5,identical.t = F,
                                              comp.ts = log(all.system.states$mad.mth$zp.ratio)))
  
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),"state.metric" = c(rep("Community",7),rep("Density",7),rep("FI",7),rep("MVI",7),rep("Z_P.ratio",7)))
  out.dens <- data.frame(rbind(as.data.frame(pc$perm.dens),as.data.frame(bio$perm.dens),as.data.frame(fi$perm.dens),as.data.frame(mvi$perm.dens),as.data.frame(zp.ratio$perm.dens)),
                         "state.metric" = c(rep("Community",10000),rep("Density",10000),rep("FI",10000),rep("MVI",10000),rep("Z_P.ratio",10000)))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(mad.zoomth.diff12) <- c("FDis","FEve","FRic")
mad.zoomth.diff12.summary <- lapply(mad.zoomth.diff12, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Mendota", res = "Month",troph = "Zooplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
mad.zoomth.diff12.raw <- lapply(mad.zoomth.diff12, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Mendota", res = "Month",troph = "Zooplankton")

# Lower Zurich Phytoplankton and Zooplankton #

LZ.phytomth.diff12<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = phyto.LZ.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(phyto.LZ.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "red.noise", lag=12,
                                       scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$LZ.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = phyto.LZ.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.LZ.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$LZ.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = phyto.LZ.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(phyto.LZ.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "red.noise",lag=12,
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = base::diff(all.system.states$LZ.mth$FI[12:391],lag=12),
                                        comp.ts.timedat = seq_along(all.system.states$LZ.mth$maxt[24:391]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = phyto.LZ.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.LZ.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = base::diff(log(all.system.states$LZ.mth$mvi)[12:391],lag=12),
                                         comp.ts.timedat = seq_along(all.system.states$LZ.mth$maxt[24:391]),
                                         pre.diff  =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = phyto.LZ.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(phyto.LZ.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "red.noise",lag=12,
                                              scale = T, normalise = F,span =12*5,identical.t = F,
                                              comp.ts = log(all.system.states$LZ.mth$zp.ratio)))
  
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),"state.metric" = c(rep("Community",7),rep("Density",7),rep("FI",7),rep("MVI",7),rep("Z_P.ratio",7)))
  out.dens <- data.frame(rbind(as.data.frame(pc$perm.dens),as.data.frame(bio$perm.dens),as.data.frame(fi$perm.dens),as.data.frame(mvi$perm.dens),as.data.frame(zp.ratio$perm.dens)),
                         "state.metric" = c(rep("Community",10000),rep("Density",10000),rep("FI",10000),rep("MVI",10000),rep("Z_P.ratio",10000)))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(LZ.phytomth.diff12) <- c("FDis","FEve","FRic")
LZ.phytomth.diff12.summary <- lapply(LZ.phytomth.diff12, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Lower Zurich", res = "Month",troph = "Phytoplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
LZ.phytomth.diff12.raw <- lapply(LZ.phytomth.diff12, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Lower Zurich", res = "Month",troph = "Phytoplankton")

LZ.zoomth.diff12<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = zoo.LZ.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(zoo.LZ.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "red.noise", lag=12,
                                       scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$LZ.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = zoo.LZ.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.LZ.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$LZ.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = zoo.LZ.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(zoo.LZ.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "red.noise",lag=12,
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = base::diff(all.system.states$LZ.mth$FI[12:391],lag=12),
                                        comp.ts.timedat = seq_along(all.system.states$LZ.mth$maxt[24:391]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = zoo.LZ.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.LZ.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = base::diff(log(all.system.states$LZ.mth$mvi)[12:391],lag=12),
                                         comp.ts.timedat = seq_along(all.system.states$LZ.mth$maxt[24:391]),
                                         pre.diff =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = zoo.LZ.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(zoo.LZ.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "red.noise",lag=12,
                                              scale = T, normalise = F,span =12*5,identical.t = F,
                                              comp.ts = log(all.system.states$LZ.mth$zp.ratio)))
  
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),"state.metric" = c(rep("Community",7),rep("Density",7),rep("FI",7),rep("MVI",7),rep("Z_P.ratio",7)))
  out.dens <- data.frame(rbind(as.data.frame(pc$perm.dens),as.data.frame(bio$perm.dens),as.data.frame(fi$perm.dens),as.data.frame(mvi$perm.dens),as.data.frame(zp.ratio$perm.dens)),
                         "state.metric" = c(rep("Community",10000),rep("Density",10000),rep("FI",10000),rep("MVI",10000),rep("Z_P.ratio",10000)))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(LZ.zoomth.diff12) <- c("FDis","FEve","FRic")
LZ.zoomth.diff12.summary <- lapply(LZ.zoomth.diff12, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Lower Zurich", res = "Month",troph = "Zooplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
LZ.zoomth.diff12.raw <- lapply(LZ.zoomth.diff12, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Lower Zurich", res = "Month",troph = "Zooplankton")

# Windermere Phytoplankton #

wind.phytomth.diff12<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = phyto.wind.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(phyto.wind.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "red.noise", lag=12,
                                       scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$wind.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = phyto.wind.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.wind.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$wind.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = phyto.wind.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(phyto.wind.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "red.noise",lag=12,
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = base::diff(all.system.states$wind.mth$FI[12:288],lag=12),
                                        comp.ts.timedat = seq_along(all.system.states$wind.mth$maxt[24:288]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = phyto.wind.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.wind.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = base::diff(log(all.system.states$wind.mth$mvi)[12:288],lag=12),
                                         comp.ts.timedat = seq_along(all.system.states$wind.mth$maxt[24:288]),
                                         pre.diff  =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = phyto.wind.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(phyto.wind.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "red.noise",lag=12,
                                              scale = T, normalise = F,span =12*5,identical.t = F,
                                              comp.ts = log(all.system.states$wind.mth$zp.ratio)))
  
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),"state.metric" = c(rep("Community",7),rep("Density",7),rep("FI",7),rep("MVI",7),rep("Z_P.ratio",7)))
  out.dens <- data.frame(rbind(as.data.frame(pc$perm.dens),as.data.frame(bio$perm.dens),as.data.frame(fi$perm.dens),as.data.frame(mvi$perm.dens),as.data.frame(zp.ratio$perm.dens)),
                         "state.metric" = c(rep("Community",10000),rep("Density",10000),rep("FI",10000),rep("MVI",10000),rep("Z_P.ratio",10000)))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(wind.phytomth.diff12) <- c("FDis","FEve","FRic")
wind.phytomth.diff12.summary <- lapply(wind.phytomth.diff12, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Windermere", res = "Month",troph = "Phytoplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
wind.phytomth.diff12.raw <- lapply(wind.phytomth.diff12, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Windermere", res = "Month",troph = "Phytoplankton")

# Windermere Zooplankton #

wind.zoomth.diff12<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = zoo.wind.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(zoo.wind.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "red.noise", lag=12,
                                       scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$wind.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = zoo.wind.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.wind.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$wind.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = zoo.wind.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(zoo.wind.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "red.noise",lag=12,
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = base::diff(all.system.states$wind.mth$FI[12:288],lag=12),
                                        comp.ts.timedat = seq_along(all.system.states$wind.mth$maxt[24:288]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = zoo.wind.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.wind.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = base::diff(log(all.system.states$wind.mth$mvi)[12:288],lag=12),
                                         comp.ts.timedat = seq_along(all.system.states$wind.mth$maxt[24:288]),
                                         pre.diff  =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = zoo.wind.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(zoo.wind.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "red.noise",lag=12,
                                              scale = T, normalise = F,span =12*5,identical.t = F,
                                              comp.ts = log(all.system.states$wind.mth$zp.ratio)))
  
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),"state.metric" = c(rep("Community",7),rep("Density",7),rep("FI",7),rep("MVI",7),rep("Z_P.ratio",7)))
  out.dens <- data.frame(rbind(as.data.frame(pc$perm.dens),as.data.frame(bio$perm.dens),as.data.frame(fi$perm.dens),as.data.frame(mvi$perm.dens),as.data.frame(zp.ratio$perm.dens)),
                         "state.metric" = c(rep("Community",10000),rep("Density",10000),rep("FI",10000),rep("MVI",10000),rep("Z_P.ratio",10000)))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(wind.zoomth.diff12) <- c("FDis","FEve","FRic")
wind.zoomth.diff12.summary <- lapply(wind.zoomth.diff12, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Windermere", res = "Month",troph = "Zooplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
wind.zoomth.diff12.raw <- lapply(wind.zoomth.diff12, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Windermere", res = "Month",troph = "Zooplankton")

# Kasumigaura Phytoplankton and Zooplankton #

kas.phytomth.diff12<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = phyto.kas.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(phyto.kas.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "red.noise", lag=12,
                                       scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$kas.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = phyto.kas.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.kas.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$kas.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = phyto.kas.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(phyto.kas.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "red.noise",lag=12,
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = base::diff(all.system.states$kas.mth$FI[12:452],lag=12),
                                        comp.ts.timedat = seq_along(all.system.states$kas.mth$maxt[24:452]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = phyto.kas.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.kas.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = base::diff(log(all.system.states$kas.mth$mvi)[12:452],lag=12),
                                         comp.ts.timedat = seq_along(all.system.states$kas.mth$maxt[24:452]),
                                         pre.diff  =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = phyto.kas.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(phyto.kas.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "red.noise",lag=12,
                                              scale = T, normalise = F,span =12*5,identical.t = F,
                                              comp.ts = log(all.system.states$kas.mth$zp.ratio)))
  
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),"state.metric" = c(rep("Community",7),rep("Density",7),rep("FI",7),rep("MVI",7),rep("Z_P.ratio",7)))
  out.dens <- data.frame(rbind(as.data.frame(pc$perm.dens),as.data.frame(bio$perm.dens),as.data.frame(fi$perm.dens),as.data.frame(mvi$perm.dens),as.data.frame(zp.ratio$perm.dens)),
                         "state.metric" = c(rep("Community",10000),rep("Density",10000),rep("FI",10000),rep("MVI",10000),rep("Z_P.ratio",10000)))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(kas.phytomth.diff12) <- c("FDis","FEve","FRic")
kas.phytomth.diff12.summary <- lapply(kas.phytomth.diff12, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Kasumigaura", res = "Month",troph = "Phytoplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
kas.phytomth.diff12.raw <- lapply(kas.phytomth.diff12, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Kasumigaura", res = "Month",troph = "Phytoplankton")

kas.zoomth.diff12<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = zoo.kas.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(zoo.kas.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "red.noise", lag=12,
                                       scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$kas.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = zoo.kas.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.kas.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$kas.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = zoo.kas.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(zoo.kas.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "red.noise",lag=12,
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = base::diff(all.system.states$kas.mth$FI[12:452],lag=12),
                                        comp.ts.timedat = seq_along(all.system.states$kas.mth$maxt[24:452]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = zoo.kas.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.kas.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = base::diff(log(all.system.states$kas.mth$mvi)[12:452],lag=12),
                                         comp.ts.timedat = seq_along(all.system.states$kas.mth$maxt[24:452]),
                                         pre.diff =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = zoo.kas.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(zoo.kas.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "red.noise",lag=12,
                                              scale = T, normalise = F,span =12*5,identical.t = F,
                                              comp.ts = log(all.system.states$kas.mth$zp.ratio)))
  
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),"state.metric" = c(rep("Community",7),rep("Density",7),rep("FI",7),rep("MVI",7),rep("Z_P.ratio",7)))
  out.dens <- data.frame(rbind(as.data.frame(pc$perm.dens),as.data.frame(bio$perm.dens),as.data.frame(fi$perm.dens),as.data.frame(mvi$perm.dens),as.data.frame(zp.ratio$perm.dens)),
                         "state.metric" = c(rep("Community",10000),rep("Density",10000),rep("FI",10000),rep("MVI",10000),rep("Z_P.ratio",10000)))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(kas.zoomth.diff12) <- c("FDis","FEve","FRic")
kas.zoomth.diff12.summary <- lapply(kas.zoomth.diff12, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Kasumigaura", res = "Month",troph = "Zooplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
kas.zoomth.diff12.raw <- lapply(kas.zoomth.diff12, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Kasumigaura", res = "Month",troph = "Zooplankton")

###########################################################################
## Save out (Lag12,Monthly) ##
###########################################################################
summary.ccf.mth12 <- rbind(kin.phytomth.diff12.summary,kin.zoomth.diff12.summary,mad.phytomth.diff12.summary,mad.zoomth.diff12.summary,
                              LZ.phytomth.diff12.summary,LZ.zoomth.diff12.summary,wind.phytomth.diff12.summary,wind.zoomth.diff12.summary,
                              kas.phytomth.diff12.summary,kas.zoomth.diff12.summary)
raw.ccf.mth12 <- rbind(kin.phytomth.diff12.raw,kin.zoomth.diff12.raw,mad.phytomth.diff12.raw,mad.zoomth.diff12.raw,
                          LZ.phytomth.diff12.raw,LZ.zoomth.diff12.raw,wind.phytomth.diff12.raw,wind.zoomth.diff12.raw,
                          kas.phytomth.diff12.raw,kas.zoomth.diff12.raw)

write.csv(summary.ccf.mth12,file ="Results/ccf/raw_data/summary.ccf.mth.lag12.csv",row.names = F)
save(raw.ccf.mth12,file = "Results/ccf/raw_data/raw.ccf.mth.lag12.RData") # RDate required to reduce file size compared to .csv

pdf(file="Results/ccf/FD_perm_lag12_diffmth_absrmax.pdf",
    width=10, height = 8)
ggplot(raw.ccf.mth12,aes(x = state.metric, y =  abs.rmax, col = FD.metric,fill= FD.metric)) + 
  geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.05, 0.5, 0.95),scale = "width",alpha = 0.3) +
  theme_bw() + 
  geom_point(data = summary.ccf.mth12[summary.ccf.mth12$measure %in% "absmax.ccf",],
             aes(x = state.metric, y = obs.value),position = position_dodge(width = 0.9),size=2) +
  geom_text(data = summary.ccf.mth12[summary.ccf.mth12$measure %in% "absmax.ccf",], 
            aes(x = state.metric, y = 0.45,label = sig),col= "black",size = 4,position = position_dodge(width = 0.9))+
  geom_text(data = summary.ccf.mth12[summary.ccf.mth12$measure %in% "t.absmax.ccf",], 
            aes(x = state.metric, y = 0.4,label = obs.value),col= "black",size = 3,position = position_dodge(width = 0.9))+
  #facet_wrap(~troph,nrow=2,strip.position = "right")+
  facet_grid(system~troph)+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  ylab("Cross correlation") + xlab("System state proxy")+   ggtitle("Permuted CCF between lag12 differenced FD and system state")
dev.off()

pdf(file="Results/ccf/FD_perm_lag12_diffmth_r0.pdf",
    width=10, height = 8)
ggplot(raw.ccf.mth12,aes(x = state.metric, y =  r0, col = FD.metric,fill= FD.metric)) + 
  geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.05, 0.5, 0.95),scale = "width",alpha = 0.3) +
  theme_bw() + 
  geom_point(data = summary.ccf.mth12[summary.ccf.mth12$measure %in% "r0.ccf",],
             aes(x = state.metric, y = obs.value),position = position_dodge(width = 0.9),size=2) +
  geom_text(data = summary.ccf.mth12[summary.ccf.mth12$measure %in% "r0.ccf",], 
            aes(x = state.metric, y = 0.45,label = sig),col= "black",size = 4,position = position_dodge(width = 0.9))+
  #facet_wrap(~troph,nrow=2,strip.position = "right")+
  facet_grid(system~troph)+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  ylab("Cross correlation") + xlab("System state proxy")+   ggtitle("Permuted CCF between lag12 differenced FD and system state")
dev.off()

###########################################################################
## Estimate cross correlation and permute (undiff, Monthly) ##
###########################################################################
# repeated 'diff.perm.ccf' function specified as @iter = 999, @perm.method = "red.noise",
# @detrend.method = "diff", @span =12*5, @identical.t = F, and @lag=1.
# Results in two dataframes:
# $summary = observed absolute max correlation coef and the corresponding lag, observed correlation at lag0 and summary statistics
# $perm.dens = all permutations absolute max correlation coef and the corresponding lag plus the correlation at lag0

# Kinneret Phytoplankton and Zooplankton #

kin.phytomth.undiff<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <-  suppressWarnings(diff.perm.ccf(dat = kin.tot[,c("date",paste(x),"community")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  bio <- suppressWarnings(diff.perm.ccf(dat = kin.tot[,c("date",paste(x),"density")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  fi <-  suppressWarnings(diff.perm.ccf(dat = kin.tot[,c("date",paste(x),"FI")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  mvi <-  suppressWarnings(diff.perm.ccf(dat = kin.tot[,c("date",paste(x),"mvi")],
                                         iter = 999,span =12*5,lag=1,
                                         perm.method = "red.noise",identical.t = F,
                                         detrend.method = "none"))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(dat = kin.tot[,c("date",paste(x),"zp.ratio")],
                                              iter = 999,span =12*5,lag=1,
                                              perm.method = "red.noise",identical.t = F,
                                              detrend.method = "none"))
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
names(kin.phytomth.undiff) <- c("FDis","FEve","FRic") # name list elements
kin.phytomth.undiff.summary <- lapply(kin.phytomth.undiff, `[[`, 'summary')%>% # extract second level list elements (i.e. 'summary')
  data.table::rbindlist(idcol = "FD.metric") %>% # rbind the list with FD.metric id column
  mutate(system = "Kinneret", res = "Month",troph = "Phytoplankton") %>% # specify metadata for future plotting
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*",""))) # assess significance of observed cross correlation by comparing to 2.5 and 97.5 quartiles (two tailed)
kin.phytomth.undiff.raw <- lapply(kin.phytomth.undiff, `[[`, 'perm.dens')%>% # extract second level list elements (i.e. 'perm.dens')
  data.table::rbindlist(idcol = "FD.metric")%>% # rbind the list with FD.metric id column
  mutate(system = "Kinneret", res = "Month",troph = "Phytoplankton") # specify metadata for future plotting

kin.zoomth.undiff<- pbmcapply::pbmclapply(c("zooFDis","zooFEve","zooFRic"),function(x){
  pc <-  suppressWarnings(diff.perm.ccf(dat = kin.tot[,c("date",paste(x),"community")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  bio <- suppressWarnings(diff.perm.ccf(dat = kin.tot[,c("date",paste(x),"density")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  fi <-  suppressWarnings(diff.perm.ccf(dat = kin.tot[,c("date",paste(x),"FI")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  mvi <-  suppressWarnings(diff.perm.ccf(dat = kin.tot[,c("date",paste(x),"mvi")],
                                         iter = 999,span =12*5,lag=1,
                                         perm.method = "red.noise",identical.t = F,
                                         detrend.method = "none"))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(dat = kin.tot[,c("date",paste(x),"zp.ratio")],
                                              iter = 999,span =12*5,lag=1,
                                              perm.method = "red.noise",identical.t = F,
                                              detrend.method = "none"))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(kin.zoomth.undiff) <- c("FDis","FEve","FRic")
kin.zoomth.undiff.summary <- lapply(kin.zoomth.undiff, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Kinneret", res = "Month",troph = "Zooplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
kin.zoomth.undiff.raw <- lapply(kin.zoomth.undiff, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Kinneret", res = "Month",troph = "Zooplankton")

# Mendota Phytoplankton and Zooplankton #

mad.phytomth.undiff<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <-  suppressWarnings(diff.perm.ccf(dat = mad.tot[,c("date",paste(x),"community")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  bio <- suppressWarnings(diff.perm.ccf(dat = mad.tot[,c("date",paste(x),"density")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  fi <-  suppressWarnings(diff.perm.ccf(dat = mad.tot[,c("date",paste(x),"FI")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  mvi <-  suppressWarnings(diff.perm.ccf(dat = mad.tot[,c("date",paste(x),"mvi")],
                                         iter = 999,span =12*5,lag=1,
                                         perm.method = "red.noise",identical.t = F,
                                         detrend.method = "none"))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(dat = mad.tot[,c("date",paste(x),"zp.ratio")],
                                              iter = 999,span =12*5,lag=1,
                                              perm.method = "red.noise",identical.t = F,
                                              detrend.method = "none"))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(mad.phytomth.undiff) <- c("FDis","FEve","FRic")
mad.phytomth.undiff.summary <- lapply(mad.phytomth.undiff, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Mendota", res = "Month",troph = "Phytoplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
mad.phytomth.undiff.raw <- lapply(mad.phytomth.undiff, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Mendota", res = "Month",troph = "Phytoplankton")

mad.zoomth.undiff<- pbmcapply::pbmclapply(c("zooFDis","zooFEve","zooFRic"),function(x){
  pc <-  suppressWarnings(diff.perm.ccf(dat = mad.tot[,c("date",paste(x),"community")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  bio <- suppressWarnings(diff.perm.ccf(dat = mad.tot[,c("date",paste(x),"density")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  fi <-  suppressWarnings(diff.perm.ccf(dat = mad.tot[,c("date",paste(x),"FI")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  mvi <-  suppressWarnings(diff.perm.ccf(dat = mad.tot[,c("date",paste(x),"mvi")],
                                         iter = 999,span =12*5,lag=1,
                                         perm.method = "red.noise",identical.t = F,
                                         detrend.method = "none"))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(dat = mad.tot[,c("date",paste(x),"zp.ratio")],
                                              iter = 999,span =12*5,lag=1,
                                              perm.method = "red.noise",identical.t = F,
                                              detrend.method = "none"))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(mad.zoomth.undiff) <- c("FDis","FEve","FRic")
mad.zoomth.undiff.summary <- lapply(mad.zoomth.undiff, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Mendota", res = "Month",troph = "Zooplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
mad.zoomth.undiff.raw <- lapply(mad.zoomth.undiff, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Mendota", res = "Month",troph = "Zooplankton")

# Lower Zurich Phytoplankton and Zooplankton #

LZ.phytomth.undiff<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <-  suppressWarnings(diff.perm.ccf(dat = LZ.tot[,c("date",paste(x),"community")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  bio <- suppressWarnings(diff.perm.ccf(dat = LZ.tot[,c("date",paste(x),"density")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  fi <-  suppressWarnings(diff.perm.ccf(dat = LZ.tot[,c("date",paste(x),"FI")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  mvi <-  suppressWarnings(diff.perm.ccf(dat = LZ.tot[,c("date",paste(x),"mvi")],
                                         iter = 999,span =12*5,lag=1,
                                         perm.method = "red.noise",identical.t = F,
                                         detrend.method = "none"))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(dat = LZ.tot[,c("date",paste(x),"zp.ratio")],
                                              iter = 999,span =12*5,lag=1,
                                              perm.method = "red.noise",identical.t = F,
                                              detrend.method = "none"))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(LZ.phytomth.undiff) <- c("FDis","FEve","FRic")
LZ.phytomth.undiff.summary <- lapply(LZ.phytomth.undiff, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Lower Zurich", res = "Month",troph = "Phytoplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
LZ.phytomth.undiff.raw <- lapply(LZ.phytomth.undiff, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Lower Zurich", res = "Month",troph = "Phytoplankton")

LZ.zoomth.undiff<- pbmcapply::pbmclapply(c("zooFDis","zooFEve","zooFRic"),function(x){
  pc <-  suppressWarnings(diff.perm.ccf(dat = LZ.tot[,c("date",paste(x),"community")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  bio <- suppressWarnings(diff.perm.ccf(dat = LZ.tot[,c("date",paste(x),"density")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  fi <-  suppressWarnings(diff.perm.ccf(dat = LZ.tot[,c("date",paste(x),"FI")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  mvi <-  suppressWarnings(diff.perm.ccf(dat = LZ.tot[,c("date",paste(x),"mvi")],
                                         iter = 999,span =12*5,lag=1,
                                         perm.method = "red.noise",identical.t = F,
                                         detrend.method = "none"))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(dat = LZ.tot[,c("date",paste(x),"zp.ratio")],
                                              iter = 999,span =12*5,lag=1,
                                              perm.method = "red.noise",identical.t = F,
                                              detrend.method = "none"))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(LZ.zoomth.undiff) <- c("FDis","FEve","FRic")
LZ.zoomth.undiff.summary <- lapply(LZ.zoomth.undiff, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Lower Zurich", res = "Month",troph = "Zooplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
LZ.zoomth.undiff.raw <- lapply(LZ.zoomth.undiff, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Lower Zurich", res = "Month",troph = "Zooplankton")

# Windermere Phytoplankton #

wind.phytomth.undiff<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <-  suppressWarnings(diff.perm.ccf(dat = wind.tot[,c("date",paste(x),"community")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  bio <- suppressWarnings(diff.perm.ccf(dat = wind.tot[,c("date",paste(x),"density")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  fi <-  suppressWarnings(diff.perm.ccf(dat = wind.tot[,c("date",paste(x),"FI")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  mvi <-  suppressWarnings(diff.perm.ccf(dat = wind.tot[,c("date",paste(x),"mvi")],
                                         iter = 999,span =12*5,lag=1,
                                         perm.method = "red.noise",identical.t = F,
                                         detrend.method = "none"))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(dat = wind.tot[,c("date",paste(x),"zp.ratio")],
                                              iter = 999,span =12*5,lag=1,
                                              perm.method = "red.noise",identical.t = F,
                                              detrend.method = "none"))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(wind.phytomth.undiff) <- c("FDis","FEve","FRic")
wind.phytomth.undiff.summary <- lapply(wind.phytomth.undiff, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Windermere", res = "Month",troph = "Phytoplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
wind.phytomth.undiff.raw <- lapply(wind.phytomth.undiff, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Windermere", res = "Month",troph = "Phytoplankton")

# Windermere Zooplankton #

wind.zoomth.undiff<- pbmcapply::pbmclapply(c("zooFDis","zooFEve","zooFRic"),function(x){
  pc <-  suppressWarnings(diff.perm.ccf(dat = wind.tot[,c("date",paste(x),"community")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  bio <- suppressWarnings(diff.perm.ccf(dat = wind.tot[,c("date",paste(x),"density")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  fi <-  suppressWarnings(diff.perm.ccf(dat = wind.tot[,c("date",paste(x),"FI")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  mvi <-  suppressWarnings(diff.perm.ccf(dat = wind.tot[,c("date",paste(x),"mvi")],
                                         iter = 999,span =12*5,lag=1,
                                         perm.method = "red.noise",identical.t = F,
                                         detrend.method = "none"))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(dat = wind.tot[,c("date",paste(x),"zp.ratio")],
                                              iter = 999,span =12*5,lag=1,
                                              perm.method = "red.noise",identical.t = F,
                                              detrend.method = "none"))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(wind.zoomth.undiff) <- c("FDis","FEve","FRic")
wind.zoomth.undiff.summary <- lapply(wind.zoomth.undiff, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Windermere", res = "Month",troph = "Zooplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
wind.zoomth.undiff.raw <- lapply(wind.zoomth.undiff, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Windermere", res = "Month",troph = "Zooplankton")

# Kasumigaura Phytoplankton #

kas.phytomth.undiff<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <-  suppressWarnings(diff.perm.ccf(dat = kas.tot[,c("date",paste(x),"community")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  bio <- suppressWarnings(diff.perm.ccf(dat = kas.tot[,c("date",paste(x),"density")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  fi <-  suppressWarnings(diff.perm.ccf(dat = kas.tot[,c("date",paste(x),"FI")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  mvi <-  suppressWarnings(diff.perm.ccf(dat = kas.tot[,c("date",paste(x),"mvi")],
                                         iter = 999,span =12*5,lag=1,
                                         perm.method = "red.noise",identical.t = F,
                                         detrend.method = "none"))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(dat = kas.tot[,c("date",paste(x),"zp.ratio")],
                                              iter = 999,span =12*5,lag=1,
                                              perm.method = "red.noise",identical.t = F,
                                              detrend.method = "none"))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(kas.phytomth.undiff) <- c("FDis","FEve","FRic")
kas.phytomth.undiff.summary <- lapply(kas.phytomth.undiff, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Kasumigaura", res = "Month",troph = "Phytoplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
kas.phytomth.undiff.raw <- lapply(kas.phytomth.undiff, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Kasumigaura", res = "Month",troph = "Phytoplankton")

# Kasumigaura Zooplankton #

kas.zoomth.undiff<- pbmcapply::pbmclapply(c("zooFDis","zooFEve","zooFRic"),function(x){
  pc <-  suppressWarnings(diff.perm.ccf(dat = kas.tot[,c("date",paste(x),"community")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  bio <- suppressWarnings(diff.perm.ccf(dat = kas.tot[,c("date",paste(x),"density")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  fi <-  suppressWarnings(diff.perm.ccf(dat = kas.tot[,c("date",paste(x),"FI")],
                                        iter = 999,span =12*5,lag=1,
                                        perm.method = "red.noise",identical.t = F,
                                        detrend.method = "none"))
  mvi <-  suppressWarnings(diff.perm.ccf(dat = kas.tot[,c("date",paste(x),"mvi")],
                                         iter = 999,span =12*5,lag=1,
                                         perm.method = "red.noise",identical.t = F,
                                         detrend.method = "none"))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(dat = kas.tot[,c("date",paste(x),"zp.ratio")],
                                              iter = 999,span =12*5,lag=1,
                                              perm.method = "red.noise",identical.t = F,
                                              detrend.method = "none"))
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),
                        "state.metric" = c(rep("Community",nrow(pc$summary)),rep("Density",nrow(bio$summary)),rep("FI",nrow(fi$summary)),rep("MVI",nrow(mvi$summary)),rep("Z_P.ratio",nrow(zp.ratio$summary))))
  out.dens <- data.frame(rbind(pc$perm.dens,bio$perm.dens,fi$perm.dens,mvi$perm.dens,zp.ratio$perm.dens),
                         "state.metric" = c(rep("Community",nrow(pc$perm.dens)),rep("Density",nrow(bio$perm.dens)),rep("FI",nrow(fi$perm.dens)),rep("MVI",nrow(mvi$perm.dens)),rep("Z_P.ratio",nrow(zp.ratio$perm.dens))))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(kas.zoomth.undiff) <- c("FDis","FEve","FRic")
kas.zoomth.undiff.summary <- lapply(kas.zoomth.undiff, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Kasumigaura", res = "Month",troph = "Zooplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
kas.zoomth.undiff.raw <- lapply(kas.zoomth.undiff, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Kasumigaura", res = "Month",troph = "Zooplankton")

###########################################################################
## Save out (undiff) ##
###########################################################################
summary.ccf.mth.undiff <- rbind(kin.phytomth.undiff.summary,kin.zoomth.undiff.summary,mad.phytomth.undiff.summary,mad.zoomth.undiff.summary,
                           LZ.phytomth.undiff.summary,LZ.zoomth.undiff.summary,wind.phytomth.undiff.summary,wind.zoomth.undiff.summary,
                           kas.phytomth.undiff.summary,kas.zoomth.undiff.summary)
raw.ccf.mth.undiff <- rbind(kin.phytomth.undiff.raw,kin.zoomth.undiff.raw,mad.phytomth.undiff.raw,mad.zoomth.undiff.raw,
                       LZ.phytomth.undiff.raw,LZ.zoomth.undiff.raw,wind.phytomth.undiff.raw,wind.zoomth.undiff.raw,
                       kas.phytomth.undiff.raw,kas.zoomth.undiff.raw)

write.csv(summary.ccf.mth.undiff,file ="Results/ccf/raw_data/summary.ccf.mth.undiff.csv",row.names = F)
save(raw.ccf.mth.undiff,file = "Results/ccf/raw_data/raw.ccf.mth.undiff.RData") # RDate required to reduce file size compared to .csv


pdf(file="Results/ccf/FD_perm_undiffmth_absrmax.pdf",
    width=10, height = 8)
ggplot(raw.ccf.mth.undiff,aes(x = state.metric, y =  abs.rmax, col = FD.metric,fill= FD.metric)) + 
  geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.05, 0.5, 0.95),scale = "width",alpha = 0.3) +
  theme_bw() + 
  geom_point(data = summary.ccf.mth.undiff[summary.ccf.mth.undiff$measure %in% "absmax.ccf",],
             aes(x = state.metric, y = obs.value),position = position_dodge(width = 0.9),size=2) +
  geom_text(data = summary.ccf.mth.undiff[summary.ccf.mth.undiff$measure %in% "absmax.ccf",], 
            aes(x = state.metric, y = 0.45,label = sig),col= "black",size = 4,position = position_dodge(width = 0.9))+
  geom_text(data = summary.ccf.mth.undiff[summary.ccf.mth.undiff$measure %in% "t.absmax.ccf",], 
            aes(x = state.metric, y = 0.4,label = obs.value),col= "black",size = 3,position = position_dodge(width = 0.9))+
  #facet_wrap(~troph,nrow=2,strip.position = "right")+
  facet_grid(system~troph)+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  ylab("Cross correlation") + xlab("System state proxy")+   ggtitle("Permuted CCF between raw FD and system state")
dev.off()

pdf(file="Results/FD_perm_un_diffmth_r0.pdf",
    width=10, height = 8)
ggplot(raw.ccf.mth.undiff,aes(x = state.metric, y =  r0, col = FD.metric,fill= FD.metric)) + 
  geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.05, 0.5, 0.95),scale = "width",alpha = 0.3) +
  theme_bw() + 
  geom_point(data = summary.ccf.mth.undiff[summary.ccf.mth.undiff$measure %in% "r0.ccf",],
             aes(x = state.metric, y = obs.value),position = position_dodge(width = 0.9),size=2) +
  geom_text(data = summary.ccf.mth.undiff[summary.ccf.mth.undiff$measure %in% "r0.ccf",], 
            aes(x = state.metric, y = 0.45,label = sig),col= "black",size = 4,position = position_dodge(width = 0.9))+
  #facet_wrap(~troph,nrow=2,strip.position = "right")+
  facet_grid(system~troph)+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  ylab("Cross correlation") + xlab("System state proxy")+   ggtitle("Permuted CCF between raw FD and system state")
dev.off()
