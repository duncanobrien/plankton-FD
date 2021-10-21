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
zoo.kas.fuzFDs.mth <- read.csv("Data/raw_FD/FD_kas_zoo_mth_raw.csv")

load("Data/all.system.states.RData")

###########################################################################
## Estimate cross correlation and permute (Lag1, Monthly) ##
###########################################################################
# repeated 'diff.perm.ccf' function specified as @iter = 10000, @perm.method = "red.noise",
# @scale = T, @normalise = F, @span =12*5, @identical.t = F, @diff=T and 
# @lag=1.
# If @comp.ts = FI or mvi, these comp.ts are pre-differenced due to the rolling window 
# approach complicating the generic nature of the 'diff.perm.ccf' function
# Results in two dataframes:
# $summary = observed absolute max correlation coef and the corresponding lag, observed correlation at lag0 and summary statistics
# $perm.dens = all permutations absolute max correlation coef and the corresponding lag plus the correlation at lag0

# Kinneret Phytoplankton and Zooplankton #

kin.phytomth.diff<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = phyto.kin.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(phyto.kin.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "red.noise", 
                                       scale = T, normalise = F,span =12*5, identical.t = F,
                                       diff=T,lag=1,
                                       comp.ts = all.system.states$kin.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = phyto.kin.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.kin.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         diff=T,lag=1,
                                         comp.ts = log(all.system.states$kin.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = phyto.kin.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(phyto.kin.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "red.noise",
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        diff=T,lag=1,
                                        comp.ts = base::diff(all.system.states$kin.mth$FI[12:552],lag=1),
                                        comp.ts.timedat = seq_along(all.system.states$kin.mth$maxt[13:552]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = phyto.kin.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.kin.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         diff=T,lag=1,
                                         comp.ts = base::diff(log(all.system.states$kin.mth$mvi)[12:552],lag=1),
                                         comp.ts.timedat = seq_along(all.system.states$kin.mth$maxt[13:552]),
                                         pre.diff =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = phyto.kin.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(phyto.kin.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "red.noise",
                                              scale = T, normalise = F,span =12*5,identical.t = F,
                                              diff=T,lag=1,
                                              comp.ts = log(all.system.states$kin.mth$zp.ratio)))
  
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),"state.metric" = c(rep("Community",7),rep("Density",7),rep("FI",7),rep("MVI",7),rep("Z_P.ratio",7)))
  # extract the observed correlation coefs for FD vs each system state
  out.dens <- data.frame(rbind(as.data.frame(pc$perm.dens),as.data.frame(bio$perm.dens),as.data.frame(fi$perm.dens),as.data.frame(mvi$perm.dens),as.data.frame(zp.ratio$perm.dens)),
                         "state.metric" = c(rep("Community",10000),rep("Density",10000),rep("FI",10000),rep("MVI",10000),rep("Z_P.ratio",10000)))
  # extract all permuted correlation coefs for FD vs each system state
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1) 
# only single core specified as 'diff.perm.ccf' already paralled within the function. 
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

kin.zoomth.diff<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = zoo.kin.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(zoo.kin.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "red.noise", lag=1,
                                       scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$kin.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = zoo.kin.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.kin.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$kin.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = zoo.kin.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(zoo.kin.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "red.noise",lag=1,
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = base::diff(all.system.states$kin.mth$FI[12:552],lag=1),
                                        comp.ts.timedat = seq_along(all.system.states$kin.mth$maxt[13:552]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = zoo.kin.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.kin.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = base::diff(log(all.system.states$kin.mth$mvi)[12:552],lag=1),
                                         comp.ts.timedat = seq_along(all.system.states$kin.mth$maxt[13:552]),
                                         pre.diff =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = zoo.kin.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(zoo.kin.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "red.noise",lag=1,
                                              scale = T, normalise = F,span =12*5,identical.t = F,
                                              comp.ts = log(all.system.states$kin.mth$zp.ratio)))
  
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),"state.metric" = c(rep("Community",7),rep("Density",7),rep("FI",7),rep("MVI",7),rep("Z_P.ratio",7)))
  out.dens <- data.frame(rbind(as.data.frame(pc$perm.dens),as.data.frame(bio$perm.dens),as.data.frame(fi$perm.dens),as.data.frame(mvi$perm.dens),as.data.frame(zp.ratio$perm.dens)),
                         "state.metric" = c(rep("Community",10000),rep("Density",10000),rep("FI",10000),rep("MVI",10000),rep("Z_P.ratio",10000)))
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
  pc <- suppressWarnings(diff.perm.ccf(ts = phyto.mad.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(phyto.mad.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "red.noise", lag=1,
                                       scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$mad.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = phyto.mad.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.mad.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$mad.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = phyto.mad.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(phyto.mad.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "red.noise",lag=1,
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = base::diff(all.system.states$mad.mth$FI,lag=1)[12:283],
                                        comp.ts.timedat = seq_along(all.system.states$mad.mth$maxt[13:283]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = phyto.mad.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.mad.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = base::diff(log(all.system.states$mad.mth$mvi),lag=1)[12:283],
                                         comp.ts.timedat = seq_along(all.system.states$mad.mth$maxt[13:283]),
                                         pre.diff  =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = phyto.mad.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(phyto.mad.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "red.noise",lag=1,
                                              scale = T, normalise = F,span =12*5,identical.t = F,
                                              comp.ts = log(all.system.states$mad.mth$zp.ratio)))
  
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),"state.metric" = c(rep("Community",7),rep("Density",7),rep("FI",7),rep("MVI",7),rep("Z_P.ratio",7)))
  out.dens <- data.frame(rbind(as.data.frame(pc$perm.dens),as.data.frame(bio$perm.dens),as.data.frame(fi$perm.dens),as.data.frame(mvi$perm.dens),as.data.frame(zp.ratio$perm.dens)),
                         "state.metric" = c(rep("Community",10000),rep("Density",10000),rep("FI",10000),rep("MVI",10000),rep("Z_P.ratio",10000)))
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

mad.zoomth.diff<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = zoo.mad.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(zoo.mad.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "red.noise", lag=1,
                                       scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$mad.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = zoo.mad.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.mad.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$mad.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = zoo.mad.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(zoo.mad.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "red.noise",lag=1,
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = base::diff(all.system.states$mad.mth$FI[12:283],lag=1),
                                        comp.ts.timedat = seq_along(all.system.states$mad.mth$maxt[13:283]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = zoo.mad.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.mad.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = base::diff(log(all.system.states$mad.mth$mvi)[12:283],lag=1),
                                         comp.ts.timedat = seq_along(all.system.states$mad.mth$maxt[13:283]),
                                         pre.diff =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = zoo.mad.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(zoo.mad.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "red.noise",lag=1,
                                              scale = T, normalise = F,span =12*5,identical.t = F,
                                              comp.ts = log(all.system.states$mad.mth$zp.ratio)))
  
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),"state.metric" = c(rep("Community",7),rep("Density",7),rep("FI",7),rep("MVI",7),rep("Z_P.ratio",7)))
  out.dens <- data.frame(rbind(as.data.frame(pc$perm.dens),as.data.frame(bio$perm.dens),as.data.frame(fi$perm.dens),as.data.frame(mvi$perm.dens),as.data.frame(zp.ratio$perm.dens)),
                         "state.metric" = c(rep("Community",10000),rep("Density",10000),rep("FI",10000),rep("MVI",10000),rep("Z_P.ratio",10000)))
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
  pc <- suppressWarnings(diff.perm.ccf(ts = phyto.LZ.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(phyto.LZ.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "red.noise", lag=1,
                                       scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$LZ.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = phyto.LZ.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.LZ.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$LZ.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = phyto.LZ.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(phyto.LZ.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "red.noise",lag=1,
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = base::diff(all.system.states$LZ.mth$FI[12:391],lag=1),
                                        comp.ts.timedat = seq_along(all.system.states$LZ.mth$maxt[13:391]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = phyto.LZ.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.LZ.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = base::diff(log(all.system.states$LZ.mth$mvi)[12:391],lag=1),
                                         comp.ts.timedat = seq_along(all.system.states$wind.mth$maxt[13:391]),
                                         pre.diff  =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = phyto.LZ.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(phyto.LZ.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "red.noise",lag=1,
                                              scale = T, normalise = F,span =12*5,identical.t = F,
                                              comp.ts = log(all.system.states$wind.mth$zp.ratio)))
  
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),"state.metric" = c(rep("Community",7),rep("Density",7),rep("FI",7),rep("MVI",7),rep("Z_P.ratio",7)))
  out.dens <- data.frame(rbind(as.data.frame(pc$perm.dens),as.data.frame(bio$perm.dens),as.data.frame(fi$perm.dens),as.data.frame(mvi$perm.dens),as.data.frame(zp.ratio$perm.dens)),
                         "state.metric" = c(rep("Community",10000),rep("Density",10000),rep("FI",10000),rep("MVI",10000),rep("Z_P.ratio",10000)))
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

LZ.zoomth.diff<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = zoo.LZ.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(zoo.LZ.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "red.noise", lag=1,
                                       scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$LZ.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = zoo.LZ.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.LZ.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$LZ.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = zoo.LZ.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(zoo.LZ.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "red.noise",lag=1,
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = base::diff(all.system.states$LZ.mth$FI[12:391],lag=1),
                                        comp.ts.timedat = seq_along(all.system.states$LZ.mth$maxt[13:391]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = zoo.LZ.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.LZ.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = base::diff(log(all.system.states$LZ.mth$mvi)[12:391],lag=1),
                                         comp.ts.timedat = seq_along(all.system.states$LZ.mth$maxt[13:391]),
                                         pre.diff =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = zoo.LZ.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(zoo.LZ.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "red.noise",lag=1,
                                              scale = T, normalise = F,span =12*5,identical.t = F,
                                              comp.ts = log(all.system.states$LZ.mth$zp.ratio)))
  
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),"state.metric" = c(rep("Community",7),rep("Density",7),rep("FI",7),rep("MVI",7),rep("Z_P.ratio",7)))
  out.dens <- data.frame(rbind(as.data.frame(pc$perm.dens),as.data.frame(bio$perm.dens),as.data.frame(fi$perm.dens),as.data.frame(mvi$perm.dens),as.data.frame(zp.ratio$perm.dens)),
                         "state.metric" = c(rep("Community",10000),rep("Density",10000),rep("FI",10000),rep("MVI",10000),rep("Z_P.ratio",10000)))
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
  pc <- suppressWarnings(diff.perm.ccf(ts = phyto.wind.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(phyto.wind.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "red.noise", lag=1,
                                       scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$wind.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = phyto.wind.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.wind.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$wind.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = phyto.wind.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(phyto.wind.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "red.noise",lag=1,
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = base::diff(all.system.states$wind.mth$FI[12:288],lag=1),
                                        comp.ts.timedat = seq_along(all.system.states$wind.mth$maxt[13:288]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = phyto.wind.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.wind.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = base::diff(log(all.system.states$wind.mth$mvi)[12:288],lag=1),
                                         comp.ts.timedat = seq_along(all.system.states$wind.mth$maxt[13:288]),
                                         pre.diff  =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = phyto.wind.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(phyto.wind.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "red.noise",lag=1,
                                              scale = T, normalise = F,span =12*5,identical.t = F,
                                              comp.ts = log(all.system.states$wind.mth$zp.ratio)))
  
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),"state.metric" = c(rep("Community",7),rep("Density",7),rep("FI",7),rep("MVI",7),rep("Z_P.ratio",7)))
  out.dens <- data.frame(rbind(as.data.frame(pc$perm.dens),as.data.frame(bio$perm.dens),as.data.frame(fi$perm.dens),as.data.frame(mvi$perm.dens),as.data.frame(zp.ratio$perm.dens)),
                         "state.metric" = c(rep("Community",10000),rep("Density",10000),rep("FI",10000),rep("MVI",10000),rep("Z_P.ratio",10000)))
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

# Kasumigaura Phytoplankton #

kas.phytomth.diff<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = phyto.kas.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(phyto.kas.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "red.noise", lag=1,
                                       scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$kas.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = phyto.kas.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.kas.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$kas.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = phyto.kas.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(phyto.kas.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "red.noise",lag=1,
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = base::diff(all.system.states$kas.mth$FI[12:452],lag=1),
                                        comp.ts.timedat = seq_along(all.system.states$kas.mth$maxt[13:452]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = phyto.kas.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.kas.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = base::diff(log(all.system.states$kas.mth$mvi)[12:452],lag=1),
                                         comp.ts.timedat = seq_along(all.system.states$kas.mth$maxt[13:452]),
                                         pre.diff  =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = phyto.kas.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(phyto.kas.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "red.noise",lag=1,
                                              scale = T, normalise = F,span =12*5,identical.t = F,
                                              comp.ts = log(all.system.states$kas.mth$zp.ratio)))
  
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),"state.metric" = c(rep("Community",7),rep("Density",7),rep("FI",7),rep("MVI",7),rep("Z_P.ratio",7)))
  out.dens <- data.frame(rbind(as.data.frame(pc$perm.dens),as.data.frame(bio$perm.dens),as.data.frame(fi$perm.dens),as.data.frame(mvi$perm.dens),as.data.frame(zp.ratio$perm.dens)),
                         "state.metric" = c(rep("Community",10000),rep("Density",10000),rep("FI",10000),rep("MVI",10000),rep("Z_P.ratio",10000)))
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

kas.zoomth.diff<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = zoo.kas.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(zoo.kas.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "red.noise", lag=1,
                                       scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$kas.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = zoo.kas.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.kas.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$kas.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = zoo.kas.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(zoo.kas.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "red.noise",lag=1,
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = base::diff(all.system.states$kas.mth$FI[12:452],lag=1),
                                        comp.ts.timedat = seq_along(all.system.states$kas.mth$maxt[13:452]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = zoo.kas.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.kas.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = base::diff(log(all.system.states$kas.mth$mvi)[12:452],lag=1),
                                         comp.ts.timedat = seq_along(all.system.states$kas.mth$maxt[13:452]),
                                         pre.diff  =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = zoo.kas.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(zoo.kas.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "red.noise",lag=1,
                                              scale = T, normalise = F,span =12*5,identical.t = F,
                                              comp.ts = log(all.system.states$kas.mth$zp.ratio)))
  
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),"state.metric" = c(rep("Community",7),rep("Density",7),rep("FI",7),rep("MVI",7),rep("Z_P.ratio",7)))
  out.dens <- data.frame(rbind(as.data.frame(pc$perm.dens),as.data.frame(bio$perm.dens),as.data.frame(fi$perm.dens),as.data.frame(mvi$perm.dens),as.data.frame(zp.ratio$perm.dens)),
                         "state.metric" = c(rep("Community",10000),rep("Density",10000),rep("FI",10000),rep("MVI",10000),rep("Z_P.ratio",10000)))
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
## Save out (Lag1, Monthly) ##
###########################################################################
summary.ccf.mth1 <- rbind(kin.phytomth.diff.summary,kin.zoomth.diff.summary,mad.phytomth.diff.summary,mad.zoomth.diff.summary,
                         LZ.phytomth.diff.summary,LZ.zoomth.diff.summary,wind.phytomth.diff.summary,kas.phytomth.diff.summary,
                         kas.zoomth.diff.summary)
raw.ccf.mth1 <- rbind(kin.phytomth.diff.raw,kin.zoomth.diff.raw,mad.phytomth.diff.raw,mad.zoomth.diff.raw,
                     LZ.phytomth.diff.raw,LZ.zoomth.diff.raw,wind.phytomth.diff.raw,kas.phytomth.diff.raw,
                     kas.zoomth.diff.raw)

write.csv(summary.ccf.mth1,file ="Results/summary.ccf.mth.lag1.csv",row.names = F)
save(raw.ccf.mth1,file = "Results/raw.ccf.mth.lag1.RData") # RData required to reduce file size compared to .csv

pdf(file="Results/FD_perm_lag1_diffmth_absrmax.pdf",
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
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  ylab("Cross correlation") + xlab("System state proxy")+   ggtitle("Permuted CCF between lag1 differenced FD and system state")
dev.off()

pdf(file="Results/FD_perm_lag1_diffmth_r0.pdf",
    width=10, height = 8)
ggplot(raw.ccf.mth1,aes(x = state.metric, y =  r0, col = FD.metric,fill= FD.metric)) + 
  geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.025, 0.5, 0.975),scale = "width",alpha = 0.3) +
  theme_bw() + 
  geom_point(data = summary.ccf.mth1[summary.ccf.mth1$measure %in% "r0.ccf",],
             aes(x = state.metric, y = obs.value),position = position_dodge(width = 0.9),size=2) +
  geom_text(data = summary.ccf.mth1[summary.ccf.mth1$measure %in% "r0.ccf",], 
            aes(x = state.metric, y = 0.35,label = sig),col= "black",size = 4,position = position_dodge(width = 0.9))+
  #facet_wrap(~troph,nrow=2,strip.position = "right")+
  facet_grid(system~troph)+
  scale_colour_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  scale_fill_manual(values=c("#969014","#22B4F5","#F07589"),name = "FD Metric") + 
  ylab("Cross correlation") + xlab("System state proxy")+   ggtitle("Permuted CCF t0 between lag1 differenced FD and system state")
dev.off()

###########################################################################
## Summary correlations (Lag1, Monthly) ##
###########################################################################
summary.ccf.mth1 <- read.csv("Results/summary.ccf.mth.lag1.csv")

obs.cor.lag0.lake.tab <- summary.ccf.mth1 %>%
            filter(measure %in% "r0.ccf")%>%
            group_by(system,troph) %>%
            summarise(mean.cor = mean(obs.value),median.cor = median(obs.value),
              se = sd(obs.value)/n(),
              nsig=sum(sig %in% "*"),prop.sig = sum(sig %in% "*")/length(sig)) %>%
            mutate(across(mean.cor:se,~signif(.x,digits=3)))
write.csv(obs.cor.lag0.lake.tab,file ="Results/ccf_tables/cor.lag0.lake.tab.csv",row.names = F)

obs.cor.lag0.FD.tab <- summary.ccf.mth1 %>%
  filter(measure %in% "r0.ccf")%>%
  group_by(troph,FD.metric) %>% 
  summarise(mean.cor = mean(obs.value),median.cor = median(obs.value),
            se = sd(obs.value)/n(),
            nsig=sum(sig %in% "*"),prop.sig = sum(sig %in% "*")/length(sig)) %>%
  mutate(across(mean.cor:se,~signif(.x,digits=3)))
write.csv(obs.cor.lag0.FD.tab,file ="Results/ccf_tables/cor.lag0.FD.tab.csv",row.names = F)

obs.cor.lag0.state.tab <- summary.ccf.mth1 %>%
  filter(measure %in% "r0.ccf")%>%
  group_by(troph,FD.metric,state.metric) %>%
  summarise(mean.cor = mean(obs.value),median.cor = median(obs.value),
            se = sd(obs.value)/n(),
            nsig=sum(sig %in% "*"),prop.sig = sum(sig %in% "*")/length(sig)) %>%
  mutate(across(mean.cor:se,~signif(.x,digits=3)))
write.csv(obs.cor.lag0.state.tab,file ="Results/ccf_tables/cor.lag0.state.tab.csv",row.names = F)

obs.cor.lagx.lake.tab <- summary.ccf.mth1 %>%
  filter(measure %in% "absmax.ccf")%>%
  group_by(system,troph) %>%
  summarise(mean.cor = mean(obs.value),median.cor = median(obs.value),
            se = sd(obs.value)/n(),
            nsig=sum(sig %in% "*"),prop.sig = sum(sig %in% "*")/length(sig)) %>%
  mutate(across(mean.cor:se,~signif(.x,digits=3)))
write.csv(obs.cor.lagx.lake.tab,file ="Results/ccf_tables/cor.lagx.lake.tab.csv",row.names = F)

obs.cor.lagx.FD.tab <- summary.ccf.mth1 %>%
  filter(measure %in% "absmax.ccf")%>%
  group_by(system,troph,FD.metric) %>%
  summarise(mean.cor = mean(obs.value),median.cor = median(obs.value),
            se = sd(obs.value)/n(),
            nsig=sum(sig %in% "*"),prop.sig = sum(sig %in% "*")/length(sig)) %>%
  mutate(across(mean.cor:se,~signif(.x,digits=3)))
write.csv(obs.cor.lagx.FD.tab,file ="Results/ccf_tables/cor.lagx.FD.tab.csv",row.names = F)

obs.cor.lagx.state.tab <- summary.ccf.mth1 %>%
  filter(measure %in% "absmax.ccf")%>%
  group_by(troph,FD.metric,state.metric) %>%
  summarise(mean.cor = mean(obs.value),median.cor = median(obs.value),
            se = sd(obs.value)/n(),
            nsig=sum(sig %in% "*"),prop.sig = sum(sig %in% "*")/length(sig)) %>%
  mutate(across(mean.cor:se,~signif(.x,digits=3)))
write.csv(obs.cor.lagx.state.tab,file ="Results/ccf_tables/cor.lagx.state.tab.csv",row.names = F)

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

kin.phytomth.diff12<- pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = phyto.kin.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(phyto.kin.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "spatial", 
                                       scale = T, normalise = F,span =12*5, identical.t = F,
                                       diff=T,lag = 12,
                                       comp.ts = all.system.states$kin.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = phyto.kin.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.kin.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "spatial",
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         diff=T,lag = 12,
                                         comp.ts = log(all.system.states$kin.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = phyto.kin.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(phyto.kin.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "spatial",
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        diff=T,lag = 12,
                                        comp.ts = base::diff(all.system.states$kin.mth$FI[12:552],lag = 12),
                                        comp.ts.timedat = seq_along(all.system.states$kin.mth$maxt[24:552]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = phyto.kin.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.kin.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "spatial",
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         diff=T,lag = 12,
                                         comp.ts = base::diff(log(all.system.states$kin.mth$mvi)[12:552],lag = 12),
                                         comp.ts.timedat = seq_along(all.system.states$kin.mth$maxt[24:552]),
                                         pre.diff =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = phyto.kin.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(phyto.kin.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "spatial",
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

kin.zoomth.diff12<- pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = zoo.kin.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(zoo.kin.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "spatial", lag=12,
                                       scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$kin.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = zoo.kin.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.kin.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "spatial",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$kin.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = zoo.kin.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(zoo.kin.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "spatial",lag=12,
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = base::diff(all.system.states$kin.mth$FI[12:552],lag=12),
                                        comp.ts.timedat = seq_along(all.system.states$kin.mth$maxt[24:552]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = zoo.kin.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.kin.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "spatial",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = base::diff(log(all.system.states$kin.mth$mvi)[12:552],lag=12),
                                         comp.ts.timedat = seq_along(all.system.states$kin.mth$maxt[24:552]),
                                         pre.diff =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = zoo.kin.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(zoo.kin.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "spatial",lag=12,
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

mad.phytomth.diff12<- pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = phyto.mad.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(phyto.mad.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "spatial", lag=12,
                                       scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$mad.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = phyto.mad.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.mad.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "spatial",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$mad.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = phyto.mad.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(phyto.mad.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "spatial",lag=12,
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = base::diff(all.system.states$mad.mth$FI[12:283],lag=12),
                                        comp.ts.timedat = seq_along(all.system.states$mad.mth$maxt[24:283]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = phyto.mad.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.mad.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "spatial",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = base::diff(log(all.system.states$mad.mth$mvi)[12:283],lag=12),
                                         comp.ts.timedat = seq_along(all.system.states$mad.mth$maxt[24:283]),
                                         pre.diff  =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = phyto.mad.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(phyto.mad.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "spatial",lag=12,
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

mad.zoomth.diff12<- pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = zoo.mad.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(zoo.mad.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "spatial", lag=12,
                                       scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$mad.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = zoo.mad.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.mad.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "spatial",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$mad.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = zoo.mad.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(zoo.mad.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "spatial",lag=12,
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = base::diff(all.system.states$mad.mth$FI[12:283],lag=12),
                                        comp.ts.timedat = seq_along(all.system.states$mad.mth$maxt[24:283]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = zoo.mad.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.mad.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "spatial",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = base::diff(log(all.system.states$mad.mth$mvi)[12:283],lag=12),
                                         comp.ts.timedat = seq_along(all.system.states$mad.mth$maxt[24:283]),
                                         pre.diff =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = zoo.mad.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(zoo.mad.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "spatial",lag=12,
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

LZ.phytomth.diff12<- pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = phyto.LZ.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(phyto.LZ.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "spatial", lag=12,
                                       scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$LZ.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = phyto.LZ.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.LZ.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "spatial",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$LZ.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = phyto.LZ.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(phyto.LZ.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "spatial",lag=12,
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = base::diff(all.system.states$LZ.mth$FI[12:391],lag=12),
                                        comp.ts.timedat = seq_along(all.system.states$LZ.mth$maxt[24:391]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = phyto.LZ.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.LZ.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "spatial",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = base::diff(log(all.system.states$LZ.mth$mvi)[12:391],lag=12),
                                         comp.ts.timedat = seq_along(all.system.states$LZ.mth$maxt[24:391]),
                                         pre.diff  =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = phyto.LZ.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(phyto.LZ.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "spatial",lag=12,
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

LZ.zoomth.diff12<- pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = zoo.LZ.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(zoo.LZ.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "spatial", lag=12,
                                       scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$LZ.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = zoo.LZ.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.LZ.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "spatial",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$LZ.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = zoo.LZ.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(zoo.LZ.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "spatial",lag=12,
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = base::diff(all.system.states$LZ.mth$FI[12:391],lag=12),
                                        comp.ts.timedat = seq_along(all.system.states$LZ.mth$maxt[24:391]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = zoo.LZ.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.LZ.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "spatial",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = base::diff(log(all.system.states$LZ.mth$mvi)[12:391],lag=12),
                                         comp.ts.timedat = seq_along(all.system.states$LZ.mth$maxt[24:391]),
                                         pre.diff =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = zoo.LZ.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(zoo.LZ.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "spatial",lag=12,
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

wind.phytomth.diff12<- pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = phyto.wind.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(phyto.wind.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "spatial", lag=12,
                                       scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$wind.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = phyto.wind.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.wind.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "spatial",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$wind.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = phyto.wind.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(phyto.wind.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "spatial",lag=12,
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = base::diff(all.system.states$wind.mth$FI[12:288],lag=12),
                                        comp.ts.timedat = seq_along(all.system.states$wind.mth$maxt[24:288]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = phyto.wind.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.wind.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "spatial",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = base::diff(log(all.system.states$wind.mth$mvi)[12:288],lag=12),
                                         comp.ts.timedat = seq_along(all.system.states$wind.mth$maxt[24:288]),
                                         pre.diff  =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = phyto.wind.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(phyto.wind.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "spatial",lag=12,
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

# Kasumigaura Phytoplankton and Zooplankton #

kas.phytomth.diff12<- pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = phyto.kas.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(phyto.kas.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "spatial", lag=12,
                                       scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$kas.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = phyto.kas.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.kas.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "spatial",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$kas.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = phyto.kas.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(phyto.kas.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "spatial",lag=12,
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = base::diff(all.system.states$kas.mth$FI[12:452],lag=12),
                                        comp.ts.timedat = seq_along(all.system.states$kas.mth$maxt[24:452]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = phyto.kas.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.kas.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "spatial",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = base::diff(log(all.system.states$kas.mth$mvi)[12:452],lag=12),
                                         comp.ts.timedat = seq_along(all.system.states$kas.mth$maxt[24:452]),
                                         pre.diff  =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = phyto.kas.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(phyto.kas.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "spatial",lag=12,
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

kas.zoomth.diff12<- pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = zoo.kas.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(zoo.kas.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "spatial", lag=12,
                                       scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$kas.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = zoo.kas.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.kas.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "spatial",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$kas.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = zoo.kas.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(zoo.kas.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "spatial",lag=12,
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = base::diff(all.system.states$kas.mth$FI[12:452],lag=12),
                                        comp.ts.timedat = seq_along(all.system.states$kas.mth$maxt[24:452]),
                                        pre.diff =T))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = zoo.kas.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.kas.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "spatial",lag=12,
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = base::diff(log(all.system.states$kas.mth$mvi)[12:452],lag=12),
                                         comp.ts.timedat = seq_along(all.system.states$kas.mth$maxt[24:452]),
                                         pre.diff =T))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = zoo.kas.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(zoo.kas.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "spatial",lag=12,
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
                              LZ.phytomth.diff12.summary,LZ.zoomth.diff12.summary,wind.phytomth.diff12.summary)
raw.ccf.mth12 <- rbind(kin.phytomth.diff12.raw,kin.zoomth.diff12.raw,mad.phytomth.diff12.raw,mad.zoomth.diff12.raw,
                          LZ.phytomth.diff12.raw,LZ.zoomth.diff12.raw,wind.phytomth.diff12.raw)

write.csv(summary.ccf.mth12,file ="Results/summary.ccf.mth.lag12.csv",row.names = F)
save(raw.ccf.mth12,file = "Results/raw.ccf.mth.lag12.RData") # RDate required to reduce file size compared to .csv

summary.ccf.mth12 <- read.csv(file ="Results/summary.ccf.mth.lag12.csv")
load(file = "Results/raw.ccf.mth.lag12.RData")
summary.ccf.mth12 <- rbind(summary.ccf.mth12,kas.phytomth.diff12.summary,kas.zoomth.diff12.summary)
raw.ccf.mth12 <- rbind(raw.ccf.mth12,kas.phytomth.diff12.raw,kas.zoomth.diff12.raw)

pdf(file="Results/FD_perm_lag12_diffmth_absrmax.pdf",
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

pdf(file="Results/FD_perm_lag12_diffmth_r0.pdf",
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
## Estimate cross correlation and permute (undiff) ##
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

kin.phytomth<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = phyto.kin.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(phyto.kin.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "red.noise", 
                                       scale = T, normalise = F,span =12*5, identical.t = F,
                                       diff=F,lag=1,
                                       comp.ts = all.system.states$kin.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = phyto.kin.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.kin.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         diff=F,lag=1,
                                         comp.ts = log(all.system.states$kin.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = phyto.kin.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(phyto.kin.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "red.noise",
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        diff=F,lag=1,
                                        comp.ts =all.system.states$kin.mth$FI[12:552],
                                        comp.ts.timedat = NULL,
                                        pre.diff = F))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = phyto.kin.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.kin.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",
                                         scale = T, normalise = F,span =12*5,identical.t = F,
                                         diff=F,lag=1,
                                         comp.ts = log(all.system.states$kin.mth$mvi)[12:552],
                                         comp.ts.timedat = NULL,
                                         pre.diff = F))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = phyto.kin.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(phyto.kin.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "red.noise",
                                              scale = T, normalise = F,span =12*5,identical.t = F,
                                              diff=F,lag=1,
                                              comp.ts = log(all.system.states$kin.mth$zp.ratio)))
  
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),"state.metric" = c(rep("Community",7),rep("Density",7),rep("FI",7),rep("MVI",7),rep("Z_P.ratio",7)))
  # extract the observed correlation coefs for FD vs each system state
  out.dens <- data.frame(rbind(as.data.frame(pc$perm.dens),as.data.frame(bio$perm.dens),as.data.frame(fi$perm.dens),as.data.frame(mvi$perm.dens),as.data.frame(zp.ratio$perm.dens)),
                         "state.metric" = c(rep("Community",10000),rep("Density",10000),rep("FI",10000),rep("MVI",10000),rep("Z_P.ratio",10000)))
  # extract all permuted correlation coefs for FD vs each system state
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1) 
# only single core specified as 'diff.perm.ccf' already paralled within the function. 
# Provides opportunity for further parallelisation if desired
names(kin.phytomth) <- c("FDis","FEve","FRic") # name list elements
kin.phytomth.summary <- lapply(kin.phytomth, `[[`, 'summary')%>% # extract second level list elements (i.e. 'summary')
  data.table::rbindlist(idcol = "FD.metric") %>% # rbind the list with FD.metric id column
  mutate(system = "Kinneret", res = "Month",troph = "Phytoplankton") %>% # specify metadata for future plotting
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*",""))) # assess significance of observed cross correlation by comparing to 2.5 and 97.5 quartiles (two tailed)
kin.phytomth.raw <- lapply(kin.phytomth, `[[`, 'perm.dens')%>% # extract second level list elements (i.e. 'perm.dens')
  data.table::rbindlist(idcol = "FD.metric")%>% # rbind the list with FD.metric id column
  mutate(system = "Kinneret", res = "Month",troph = "Phytoplankton") # specify metadata for future plotting

kin.zoomth<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = zoo.kin.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(zoo.kin.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "red.noise", lag=1,
                                       diff=F,scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$kin.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = zoo.kin.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.kin.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         diff=F,scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$kin.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = zoo.kin.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(zoo.kin.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "red.noise",lag=1,
                                        scale = T, normalise = F,span =12*5,identical.t = F,
                                        diff=F,comp.ts = all.system.states$kin.mth$FI[12:552],
                                        comp.ts.timedat = NULL,
                                        pre.diff = F))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = zoo.kin.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.kin.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         diff=F,scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$kin.mth$mvi)[12:552],
                                         comp.ts.timedat = NULL,
                                         pre.diff = F))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = zoo.kin.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(zoo.kin.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "red.noise",lag=1,
                                              scale = T, normalise = F,span =12*5,identical.t = F,
                                              comp.ts = log(all.system.states$kin.mth$zp.ratio)))
  
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),"state.metric" = c(rep("Community",7),rep("Density",7),rep("FI",7),rep("MVI",7),rep("Z_P.ratio",7)))
  out.dens <- data.frame(rbind(as.data.frame(pc$perm.dens),as.data.frame(bio$perm.dens),as.data.frame(fi$perm.dens),as.data.frame(mvi$perm.dens),as.data.frame(zp.ratio$perm.dens)),
                         "state.metric" = c(rep("Community",10000),rep("Density",10000),rep("FI",10000),rep("MVI",10000),rep("Z_P.ratio",10000)))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(kin.zoomth) <- c("FDis","FEve","FRic")
kin.zoomth.summary <- lapply(kin.zoomth, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Kinneret", res = "Month",troph = "Zooplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
kin.zoomth.raw <- lapply(kin.zoomth, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Kinneret", res = "Month",troph = "Zooplankton")


# Mendota Phytoplankton and Zooplankton #

mad.phytomth<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = phyto.mad.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(phyto.mad.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "red.noise", lag=1,
                                       diff=F,scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$mad.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = phyto.mad.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.mad.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         diff=F,scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$mad.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = phyto.mad.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(phyto.mad.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "red.noise",lag=1,
                                        diff=F, scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = all.system.states$mad.mth$FI[12:283],
                                        comp.ts.timedat = NULL,
                                        pre.diff = F))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = phyto.mad.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.mad.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         diff=F,scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$mad.mth$mvi)[12:283],
                                         comp.ts.timedat = NULL,
                                         pre.diff  = F))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = phyto.mad.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(phyto.mad.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "red.noise",lag=1,
                                              diff=F,scale = T, normalise = F,span =12*5,identical.t = F,
                                              comp.ts = log(all.system.states$mad.mth$zp.ratio)))
  
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),"state.metric" = c(rep("Community",7),rep("Density",7),rep("FI",7),rep("MVI",7),rep("Z_P.ratio",7)))
  out.dens <- data.frame(rbind(as.data.frame(pc$perm.dens),as.data.frame(bio$perm.dens),as.data.frame(fi$perm.dens),as.data.frame(mvi$perm.dens),as.data.frame(zp.ratio$perm.dens)),
                         "state.metric" = c(rep("Community",10000),rep("Density",10000),rep("FI",10000),rep("MVI",10000),rep("Z_P.ratio",10000)))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(mad.phytomth) <- c("FDis","FEve","FRic")
mad.phytomth.summary <- lapply(mad.phytomth, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Mendota", res = "Month",troph = "Phytoplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
mad.phytomth.raw <- lapply(mad.phytomth, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Mendota", res = "Month",troph = "Phytoplankton")

mad.zoomth<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = zoo.mad.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(zoo.mad.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "red.noise", lag=1,
                                       diff=F,scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$mad.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = zoo.mad.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.mad.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         diff=F,scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$mad.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = zoo.mad.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(zoo.mad.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "red.noise",lag=1,
                                        diff=F,scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = all.system.states$mad.mth$FI[12:283],
                                        comp.ts.timedat = NULL,
                                        pre.diff = F))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = zoo.mad.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.mad.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         diff=F,scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts =log(all.system.states$mad.mth$mvi)[12:283],
                                         comp.ts.timedat = NULL,
                                         pre.diff = F))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = zoo.mad.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(zoo.mad.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "red.noise",lag=1,
                                              diff=F,scale = T, normalise = F,span =12*5,identical.t = F,
                                              comp.ts = log(all.system.states$mad.mth$zp.ratio)))
  
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),"state.metric" = c(rep("Community",7),rep("Density",7),rep("FI",7),rep("MVI",7),rep("Z_P.ratio",7)))
  out.dens <- data.frame(rbind(as.data.frame(pc$perm.dens),as.data.frame(bio$perm.dens),as.data.frame(fi$perm.dens),as.data.frame(mvi$perm.dens),as.data.frame(zp.ratio$perm.dens)),
                         "state.metric" = c(rep("Community",10000),rep("Density",10000),rep("FI",10000),rep("MVI",10000),rep("Z_P.ratio",10000)))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(mad.zoomth) <- c("FDis","FEve","FRic")
mad.zoomth.summary <- lapply(mad.zoomth, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Mendota", res = "Month",troph = "Zooplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
mad.zoomth.raw <- lapply(mad.zoomth, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Mendota", res = "Month",troph = "Zooplankton")

# Lower Zurich Phytoplankton and Zooplankton #

LZ.phytomth<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = phyto.LZ.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(phyto.LZ.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "red.noise", lag=1,
                                       diff=F,scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$LZ.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = phyto.LZ.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.LZ.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         diff=F,scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$LZ.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = phyto.LZ.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(phyto.LZ.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "red.noise",lag=1,
                                        diff=F,scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = all.system.states$LZ.mth$FI[12:391],
                                        comp.ts.timedat = NULL,
                                        pre.diff = F))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = phyto.LZ.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.LZ.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         diff=F,scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$LZ.mth$mvi)[12:391],
                                         comp.ts.timedat = NULL,
                                         pre.diff  =F))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = phyto.LZ.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(phyto.LZ.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "red.noise",lag=1,
                                              diff=F,scale = T, normalise = F,span =12*5,identical.t = F,
                                              comp.ts = log(all.system.states$LZ.mth$zp.ratio)))
  
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),"state.metric" = c(rep("Community",7),rep("Density",7),rep("FI",7),rep("MVI",7),rep("Z_P.ratio",7)))
  out.dens <- data.frame(rbind(as.data.frame(pc$perm.dens),as.data.frame(bio$perm.dens),as.data.frame(fi$perm.dens),as.data.frame(mvi$perm.dens),as.data.frame(zp.ratio$perm.dens)),
                         "state.metric" = c(rep("Community",10000),rep("Density",10000),rep("FI",10000),rep("MVI",10000),rep("Z_P.ratio",10000)))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(LZ.phytomth) <- c("FDis","FEve","FRic")
LZ.phytomth.summary <- lapply(LZ.phytomth, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Lower Zurich", res = "Month",troph = "Phytoplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
LZ.phytomth.raw <- lapply(LZ.phytomth, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Lower Zurich", res = "Month",troph = "Phytoplankton")

LZ.zoomth<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = zoo.LZ.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(zoo.LZ.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "red.noise", lag=1,
                                       diff=F,scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$LZ.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = zoo.LZ.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.LZ.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         diff=F, scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$LZ.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = zoo.LZ.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(zoo.LZ.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "red.noise",lag=1,
                                        diff=F,scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = all.system.states$LZ.mth$FI[12:391],
                                        comp.ts.timedat = NULL,
                                        pre.diff = F))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = zoo.LZ.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(zoo.LZ.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         diff=F,scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = all.system.states$LZ.mth$mvi[12:391],
                                         comp.ts.timedat = NULL,
                                         pre.diff = F))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = zoo.LZ.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(zoo.LZ.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "red.noise",lag=1,
                                              diff=F,scale = T, normalise = F,span =12*5,identical.t = F,
                                              comp.ts = log(all.system.states$LZ.mth$zp.ratio)))
  
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),"state.metric" = c(rep("Community",7),rep("Density",7),rep("FI",7),rep("MVI",7),rep("Z_P.ratio",7)))
  out.dens <- data.frame(rbind(as.data.frame(pc$perm.dens),as.data.frame(bio$perm.dens),as.data.frame(fi$perm.dens),as.data.frame(mvi$perm.dens),as.data.frame(zp.ratio$perm.dens)),
                         "state.metric" = c(rep("Community",10000),rep("Density",10000),rep("FI",10000),rep("MVI",10000),rep("Z_P.ratio",10000)))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(LZ.zoomth) <- c("FDis","FEve","FRic")
LZ.zoomth.summary <- lapply(LZ.zoomth, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Lower Zurich", res = "Month",troph = "Zooplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
LZ.zoomth.raw <- lapply(LZ.zoomth, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Lower Zurich", res = "Month",troph = "Zooplankton")

# Windermere Phytoplankton #

wind.phytomth<- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  pc <- suppressWarnings(diff.perm.ccf(ts = phyto.wind.fuzFDs.mth[,paste(x)], 
                                       timedat = as.numeric(phyto.wind.fuzFDs.mth$date),
                                       iter = 10000,perm.method = "red.noise", lag=1,
                                       diff=F,scale = T, normalise = F,span =12*5,identical.t = F,
                                       comp.ts = all.system.states$wind.mth$community))
  bio <-  suppressWarnings(diff.perm.ccf(ts = phyto.wind.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.wind.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         diff=F,scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = log(all.system.states$wind.mth$density)))
  fi <-  suppressWarnings(diff.perm.ccf(ts = phyto.wind.fuzFDs.mth[,paste(x)], 
                                        timedat = as.numeric(phyto.wind.fuzFDs.mth$date),
                                        iter = 10000,perm.method = "red.noise",lag=1,
                                        diff=F,scale = T, normalise = F,span =12*5,identical.t = F,
                                        comp.ts = all.system.states$wind.mth$FI[12:288],
                                        comp.ts.timedat = seq_along(all.system.states$wind.mth$date[12:288]),                                        
                                        pre.diff = F))
  mvi <-  suppressWarnings(diff.perm.ccf(ts = phyto.wind.fuzFDs.mth[,paste(x)], 
                                         timedat = as.numeric(phyto.wind.fuzFDs.mth$date),
                                         iter = 10000,perm.method = "red.noise",lag=1,
                                         diff=F,scale = T, normalise = F,span =12*5,identical.t = F,
                                         comp.ts = all.system.states$wind.mth$mvi[12:288],
                                         comp.ts.timedat =  seq_along(all.system.states$wind.mth$date)[12:288],
                                         pre.diff  =F))
  zp.ratio <-  suppressWarnings(diff.perm.ccf(ts = phyto.wind.fuzFDs.mth[,paste(x)], 
                                              timedat = as.numeric(phyto.wind.fuzFDs.mth$date),
                                              iter = 10000,perm.method = "red.noise",lag=1,
                                              diff=F, scale = T, normalise = F,span =12*5,identical.t = F,
                                              comp.ts = log(all.system.states$wind.mth$zp.ratio)))
  
  out.val <- data.frame(rbind(pc$summary,bio$summary,fi$summary,mvi$summary,zp.ratio$summary),"state.metric" = c(rep("Community",7),rep("Density",7),rep("FI",7),rep("MVI",7),rep("Z_P.ratio",7)))
  out.dens <- data.frame(rbind(as.data.frame(pc$perm.dens),as.data.frame(bio$perm.dens),as.data.frame(fi$perm.dens),as.data.frame(mvi$perm.dens),as.data.frame(zp.ratio$perm.dens)),
                         "state.metric" = c(rep("Community",10000),rep("Density",10000),rep("FI",10000),rep("MVI",10000),rep("Z_P.ratio",10000)))
  out <- list("summary" = out.val,"perm.dens" = out.dens) 
  return(out)
},mc.cores = 1)
names(wind.phytomth) <- c("FDis","FEve","FRic")
wind.phytomth.summary <- lapply(wind.phytomth, `[[`, 'summary')%>%
  data.table::rbindlist(idcol = "FD.metric") %>% 
  mutate(system = "Windermere", res = "Month",troph = "Phytoplankton") %>%
  mutate(sig = ifelse(is.na(obs.value) | obs.value == "NaN" | obs.difference == 0,"",
                      ifelse(quantile >= 0.975 | quantile <= 0.025,"*","")))
wind.phytomth.raw <- lapply(wind.phytomth, `[[`, 'perm.dens')%>%
  data.table::rbindlist(idcol = "FD.metric")%>%
  mutate(system = "Windermere", res = "Month",troph = "Phytoplankton")

###########################################################################
## Save out (undiff) ##
###########################################################################
summary.ccf.mth.undiff <- rbind(kin.phytomth.summary,kin.zoomth.summary,mad.phytomth.summary,mad.zoomth.summary,
                           LZ.phytomth.summary,LZ.zoomth.summary,wind.phytomth.summary)
raw.ccf.mth.undiff <- rbind(kin.phytomth.raw,kin.zoomth.raw,mad.phytomth.raw,mad.zoomth.raw,
                       LZ.phytomth.raw,LZ.zoomth.raw,wind.phytomth.raw)

write.csv(summary.ccf.mth.undiff,file ="Results/summary.ccf.mth.undiff.csv",row.names = F)
save(raw.ccf.mth.undiff,file = "Results/raw.ccf.mth.undiff.RData") # RDate required to reduce file size compared to .csv


pdf(file="Results/FD_perm_undiffmth_absrmax.pdf",
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
  ylab("Cross correlation") + xlab("System state proxy")+   ggtitle("Permuted CCF between lag12 differenced FD and system state")
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
  ylab("Cross correlation") + xlab("System state proxy")+   ggtitle("Permuted CCF between lag12 differenced FD and system state")
dev.off()
