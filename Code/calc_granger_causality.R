### Granger Causality ###

## Preamble ##
require(tidyverse) # dplyr, ggplot etc.
require(vars) # VAR fitting functions 
require(lmtest) # GrangerTest function
require(pbmcapply) # paralled lapply function
require(foreach) # paralled 'for' loops
require(zoo) # na.approx function

source("Code/cross_granger_fn.R")

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
## Estimate Granger Causality ##
###########################################################################

kin.granger.TY <- pbmcapply::pbmclapply(c("FDis","FEve","FRic","zooFDis","zooFEve","zooFRic"),function(x){
  
  tmp <- foreach::foreach(i = c("density","community","FI","mvi","zp.ratio"),.combine = "rbind",.export = c("cross.granger"),.packages = c("vars","dplyr")) %do%{
    
    if(any(is.na(kin.tot[,paste(x)]))){
      kin.tot[,paste(x)] <- zoo::na.approx(kin.tot[,paste(x)],maxgap =3,na.rm=F)
    } # certain FD metrics may have NAs for single data points
    
    df <- na.omit(kin.tot)
    
    gc.df <- cross.granger.ic(ts = df[,paste(x)], comp.ts = df[,paste(i)],span = 12*5,method="TY",covariates = df[,"env"],ic = "AIC")
    gc.df$state.metric <- paste(i)
    return(gc.df)
  }
  return(tmp)
},mc.cores = 3) 
names(kin.granger.TY) <- c("FDis","FEve","FRic","zooFDis","zooFEve","zooFRic")
kin.granger.TY <- data.frame(data.table::rbindlist(kin.granger.TY,idcol = "FD.metric"))%>%
  mutate(system = "Kinneret", res = "Month") # specify metadata for future plotting

LZ.granger.TY <- pbmcapply::pbmclapply(c("FDis","FEve","FRic","zooFDis","zooFEve","zooFRic"),function(x){
  
  tmp <- foreach::foreach(i = c("density","community","FI","mvi","zp.ratio"),.combine = "rbind",.export = c("cross.granger"),.packages = c("vars","dplyr")) %do%{
    
    if(any(is.na(LZ.tot[,paste(x)]))){
      LZ.tot[,paste(x)] <- zoo::na.approx(LZ.tot[,paste(x)],maxgap =3,na.rm=F)
    } # certain FD metrics may have NAs for single data points
    
    df <- na.omit(LZ.tot)
    
    gc.df <- cross.granger.ic(ts = df[,paste(x)], comp.ts = df[,paste(i)],span = 12*5,method="TY",covariates = df[,"env"],ic = "AIC")
    gc.df$state.metric <- paste(i)
    return(gc.df)
  }
  return(tmp)
},mc.cores = 3) 
names(LZ.granger.TY) <- c("FDis","FEve","FRic","zooFDis","zooFEve","zooFRic")
LZ.granger.TY <- data.frame(data.table::rbindlist(LZ.granger.TY,idcol = "FD.metric"))%>%
  mutate(system = "Lower Zurich", res = "Month") # specify metadata for future plotting

kas.granger.TY <- pbmcapply::pbmclapply(c("FDis","FEve","FRic","zooFDis","zooFEve","zooFRic"),function(x){
  
  tmp <- foreach::foreach(i = c("density","community","FI","mvi","zp.ratio"),.combine = "rbind",.export = c("cross.granger"),.packages = c("vars","dplyr")) %do%{
    
    if(any(is.na(kas.tot[,paste(x)]))){
      kas.tot[,paste(x)] <- zoo::na.approx(kas.tot[,paste(x)],maxgap =3,na.rm=F)
    } # certain FD metrics may have NAs for single data points
    
    df <- na.omit(kas.tot)
    
    gc.df <- cross.granger.ic(ts = df[,paste(x)], comp.ts = df[,paste(i)],span = 12*5,method="TY",covariates = df[,"env"],ic = "AIC")
    gc.df$state.metric <- paste(i)
    return(gc.df)
  }
  return(tmp)
},mc.cores = 3) 
names(kas.granger.TY) <- c("FDis","FEve","FRic","zooFDis","zooFEve","zooFRic")
kas.granger.TY <- data.frame(data.table::rbindlist(kas.granger.TY,idcol = "FD.metric"))%>%
  mutate(system = "Kasumigaura", res = "Month") # specify metadata for future plotting

mad.granger.TY <- pbmcapply::pbmclapply(c("FDis","FEve","FRic","zooFDis","zooFEve","zooFRic"),function(x){
  
  tmp <- foreach::foreach(i = c("density","community","FI","mvi","zp.ratio"),.combine = "rbind",.export = c("cross.granger"),.packages = c("vars","dplyr")) %do%{
    
    if(any(is.na(mad.tot[,paste(x)]))){
      mad.tot[,paste(x)] <- zoo::na.approx(mad.tot[,paste(x)],maxgap =3,na.rm=F)
    } # certain FD metrics may have NAs for single data points
    
    df <- na.omit(mad.tot)
    
    gc.df <- cross.granger.ic(ts = df[,paste(x)], comp.ts = df[,paste(i)],span = 12*5,method="TY",covariates = df[,"env"],ic = "AIC")
    gc.df$state.metric <- paste(i)
    return(gc.df)
  }
  return(tmp)
},mc.cores = 3) 
names(mad.granger.TY) <- c("FDis","FEve","FRic","zooFDis","zooFEve","zooFRic")
mad.granger.TY <- data.frame(data.table::rbindlist(mad.granger.TY,idcol = "FD.metric"))%>%
  mutate(system = "Mendota", res = "Month") # specify metadata for future plotting

wind.granger.TY <- pbmcapply::pbmclapply(c("FDis","FEve","FRic","zooFDis","zooFEve","zooFRic"),function(x){
  
  tmp <- foreach::foreach(i = c("density","community","FI","mvi","zp.ratio"),.combine = "rbind",.export = c("cross.granger"),.packages = c("vars","dplyr")) %do%{
    
    if(any(is.na(wind.tot[,paste(x)]))){
      wind.tot[,paste(x)] <- zoo::na.approx(wind.tot[,paste(x)],maxgap =3,na.rm=F)
    } # certain FD metrics may have NAs for single data points
    
    df <- na.omit(wind.tot)
    
    gc.df <- cross.granger.ic(ts = df[,paste(x)], comp.ts = df[,paste(i)],span = 12*5,method="TY",covariates = df[,"env"],ic = "AIC")
    gc.df$state.metric <- paste(i)
    return(gc.df)
  }
  return(tmp)
},mc.cores = 3) 
names(wind.granger.TY) <- c("FDis","FEve","FRic","zooFDis","zooFEve","zooFRic")
wind.granger.TY <- data.frame(data.table::rbindlist(wind.granger.TY,idcol = "FD.metric"))%>%
  mutate(system = "Windermere", res = "Month") # specify metadata for future plotting


kin.granger.TYrestrict <- pbmcapply::pbmclapply(c("FDis","FEve","FRic","zooFDis","zooFEve","zooFRic"),function(x){
  
  tmp <- foreach::foreach(i = c("density","community","FI","mvi","zp.ratio"),.combine = "rbind",.export = c("cross.granger"),.packages = c("vars","dplyr")) %do%{
    
    if(any(is.na(kin.tot[,paste(x)]))){
      kin.tot[,paste(x)] <- zoo::na.approx(kin.tot[,paste(x)],maxgap =3,na.rm=F)
    } # certain FD metrics may have NAs for single data points
    
    df <- na.omit(kin.tot)
    
    gc.df <- cross.granger.ic(ts = df[,paste(x)], comp.ts = df[,paste(i)],span = 12*5,method="TYrestrict",covariates = df[,"env"],ic = "AIC")
    gc.df$state.metric <- paste(i)
    return(gc.df)
  }
  return(tmp)
},mc.cores = 3) 
names(kin.granger.TYrestrict) <- c("FDis","FEve","FRic","zooFDis","zooFEve","zooFRic")
kin.granger.TYrestrict <- data.frame(data.table::rbindlist(kin.granger.TYrestrict,idcol = "FD.metric"))%>%
  mutate(system = "Kinneret", res = "Month") # specify metadata for future plotting

LZ.granger.TYrestrict <- pbmcapply::pbmclapply(c("FDis","FEve","FRic","zooFDis","zooFEve","zooFRic"),function(x){
  
  tmp <- foreach::foreach(i = c("density","community","FI","mvi","zp.ratio"),.combine = "rbind",.export = c("cross.granger"),.packages = c("vars","dplyr")) %do%{
    
    if(any(is.na(LZ.tot[,paste(x)]))){
      LZ.tot[,paste(x)] <- zoo::na.approx(LZ.tot[,paste(x)],maxgap =3,na.rm=F)
    } # certain FD metrics may have NAs for single data points
    
    df <- na.omit(LZ.tot)
    
    gc.df <- cross.granger.ic(ts = df[,paste(x)], comp.ts = df[,paste(i)],span = 12*5,method="TYrestrict",covariates = df[,"env"],ic = "AIC")
    gc.df$state.metric <- paste(i)
    return(gc.df)
  }
  return(tmp)
},mc.cores = 3) 
names(LZ.granger.TYrestrict) <- c("FDis","FEve","FRic","zooFDis","zooFEve","zooFRic")
LZ.granger.TYrestrict <- data.frame(data.table::rbindlist(LZ.granger.TYrestrict,idcol = "FD.metric"))%>%
  mutate(system = "Lower Zurich", res = "Month") # specify metadata for future plotting

kas.granger.TYrestrict <- pbmcapply::pbmclapply(c("FDis","FEve","FRic","zooFDis","zooFEve","zooFRic"),function(x){
  
  tmp <- foreach::foreach(i = c("density","community","FI","mvi","zp.ratio"),.combine = "rbind",.export = c("cross.granger"),.packages = c("vars","dplyr")) %do%{
    
    if(any(is.na(kas.tot[,paste(x)]))){
      kas.tot[,paste(x)] <- zoo::na.approx(kas.tot[,paste(x)],maxgap =3,na.rm=F)
    } # certain FD metrics may have NAs for single data points
    
    df <- na.omit(kas.tot)
    
    gc.df <- cross.granger.ic(ts = df[,paste(x)], comp.ts = df[,paste(i)],span = 12*5,method="TYrestrict",covariates = df[,"env"],ic = "AIC")
    gc.df$state.metric <- paste(i)
    return(gc.df)
  }
  return(tmp)
},mc.cores = 3) 
names(kas.granger.TYrestrict) <- c("FDis","FEve","FRic","zooFDis","zooFEve","zooFRic")
kas.granger.TYrestrict <- data.frame(data.table::rbindlist(kas.granger.TYrestrict,idcol = "FD.metric"))%>%
  mutate(system = "Kasumigaura", res = "Month") # specify metadata for future plotting

mad.granger.TYrestrict <- pbmcapply::pbmclapply(c("FDis","FEve","FRic","zooFDis","zooFEve","zooFRic"),function(x){
  
  tmp <- foreach::foreach(i = c("density","community","FI","mvi","zp.ratio"),.combine = "rbind",.export = c("cross.granger"),.packages = c("vars","dplyr")) %do%{
    
    if(any(is.na(mad.tot[,paste(x)]))){
      mad.tot[,paste(x)] <- zoo::na.approx(mad.tot[,paste(x)],maxgap =3,na.rm=F)
    } # certain FD metrics may have NAs for single data points
    
    df <- na.omit(mad.tot)
    
    gc.df <- cross.granger.ic(ts = df[,paste(x)], comp.ts = df[,paste(i)],span = 12*5,method="TYrestrict",covariates = df[,"env"],ic = "AIC")
    gc.df$state.metric <- paste(i)
    return(gc.df)
  }
  return(tmp)
},mc.cores = 3) 
names(mad.granger.TYrestrict) <- c("FDis","FEve","FRic","zooFDis","zooFEve","zooFRic")
mad.granger.TYrestrict <- data.frame(data.table::rbindlist(mad.granger.TYrestrict,idcol = "FD.metric"))%>%
  mutate(system = "Mendota", res = "Month") # specify metadata for future plotting

wind.granger.TYrestrict <- pbmcapply::pbmclapply(c("FDis","FEve","FRic","zooFDis","zooFEve","zooFRic"),function(x){
  
  tmp <- foreach::foreach(i = c("density","community","FI","mvi","zp.ratio"),.combine = "rbind",.export = c("cross.granger"),.packages = c("vars","dplyr")) %do%{
    
    if(any(is.na(wind.tot[,paste(x)]))){
      wind.tot[,paste(x)] <- zoo::na.approx(wind.tot[,paste(x)],maxgap =3,na.rm=F)
    } # certain FD metrics may have NAs for single data points
    
    df <- na.omit(wind.tot)
    
    gc.df <- cross.granger.ic(ts = df[,paste(x)], comp.ts = df[,paste(i)],span = 12*5,method="TYrestrict",covariates = df[,"env"],ic = "AIC")
    gc.df$state.metric <- paste(i)
    return(gc.df)
  }
  return(tmp)
},mc.cores = 3) 
names(wind.granger.TYrestrict) <- c("FDis","FEve","FRic","zooFDis","zooFEve","zooFRic")
wind.granger.TYrestrict <- data.frame(data.table::rbindlist(wind.granger.TYrestrict,idcol = "FD.metric"))%>%
  mutate(system = "Windermere", res = "Month") # specify metadata for future plotting

###########################################################################
## Save out ##
###########################################################################
all.lakes.granger.TYRes <- rbind(kin.granger.TYrestrict,mad.granger.TYrestrict,LZ.granger.TYrestrict,wind.granger.TYrestrict,kas.granger.TYrestrict) %>%
  filter(P.value <= 0.05)

all.lakes.granger.TY <- rbind(kin.granger.TY,mad.granger.TY,LZ.granger.TY,wind.granger.TY,kas.granger.TY) %>%
  filter(P.value <= 0.05)

write.csv(all.lakes.granger.TY,file ="Results/granger_causality/all.lakes.TY.granger.csv",row.names = F)

all.lakes.granger.TY <- read.csv("Results/granger_causality/all.lakes.TY.granger.csv")

count.TYdf <- all.lakes.granger.TY %>%
  group_by(state.metric,causality.direc,FD.metric)%>%
  mutate(N = n()) # count per group

count.TYResdf <- all.lakes.granger.TYRes %>%
  group_by(state.metric,causality.direc,FD.metric)%>%
  mutate(N = n()) # count per group

pdf(file="Results/granger_causality/granger_causality_spread.pdf",
    width=10, height = 7)
ggplot(all.lakes.granger.TY,aes(x=lag,y= state.metric,fill= causality.direc)) + 
  geom_boxplot(alpha = 0.8,size = 0.5)+
  geom_point(aes(x=lag,y=state.metric, group=causality.direc), alpha = 0.5, position = position_dodge(width=0.75),size = 1.3) + 
  scale_fill_manual(values = c("#FFE7A1","#A1B4FE"),name = "Causality\ndirection",labels = c("Forward", "Reverse"))+
  geom_text(data = count.TYdf %>% distinct(FD.metric, causality.direc, state.metric, N),
           aes(x = (max(all.lakes.granger.TY$lag)+5),y=state.metric, label = N),
          position = position_dodge(width = 0.8))+
  facet_wrap(~FD.metric) +
  ylab("State metric") + xlab("Optimum lag (months)")+
  theme_bw()
dev.off()

pdf(file="Results/granger_causality/granger_causality_spread_alt.pdf",
    width=10, height = 7)
ggplot(all.lakes.granger.TY,aes(x=lag,y= state.metric,fill= causality.direc)) + 
  geom_boxplot(alpha = 0.1,size = 0.1)+
  geom_point(aes(x=lag,y=state.metric, group=causality.direc,fill=causality.direc,shape = system),  alpha = 0.8, position = position_dodge(width=0.75),size = 3) + 
  geom_text(data = count.TYdf %>% distinct(FD.metric, causality.direc, state.metric, N),
            aes(x = (max(all.lakes.granger.TY$lag)+5),y=state.metric, label = N),
            position = position_dodge(width = 0.8))+
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_fill_manual(values = c("#FFE7A1","#A1B4FE"),name = "Causality\ndirection",labels = c("Forward", "Reverse"))+
  facet_wrap(~FD.metric) +
  ylab("State metric") + xlab("Optimum lag (months)")+
  guides(fill = guide_legend(override.aes = list(col = c("#FFE7A1","#A1B4FE"))))+
  theme_bw()
dev.off()


ggplot(all.lakes.granger.TYRes,aes(x=lag,y= state.metric,fill= causality.direc)) + 
  geom_boxplot(alpha = 0.1,size = 0.1)+
  geom_point(aes(x=lag,y=state.metric, group=causality.direc,fill=causality.direc,shape = system),  alpha = 0.8, position = position_dodge(width=0.75),size = 3) + 
  geom_text(data = count.TYResdf %>% distinct(FD.metric, causality.direc, state.metric, N),
            aes(x = (max(all.lakes.granger.TYRes$lag)+5),y=state.metric, label = N),
            position = position_dodge(width = 0.8))+
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_fill_manual(values = c("#FFE7A1","#A1B4FE"),name = "Causality\ndirection",labels = c("Forward", "Reverse"))+
  facet_wrap(~FD.metric) +
  ylab("State metric") + xlab("Optimum lag (months)")+
  guides(fill = guide_legend(override.aes = list(col = c("#FFE7A1","#A1B4FE"))))+
  theme_bw()

gc.lag.changes.df <- all.lakes.granger %>%
  filter(!grepl("\\.aic$",statistic) & !grepl("\\.sig$",statistic))%>% # drop AIC and sig rows
  group_by(FD.metric,state.metric,system) %>%
  pivot_longer(c(lag_1:lag_60),names_to = "lag",values_to = "value",names_prefix = "lag_") %>% # pivot lag columns in to single column
  mutate(value = as.numeric(value),lag = as.numeric(lag))%>% # ensure values are numeric
  pivot_wider(names_from = statistic,values_from = value)%>% # pivot F/P statistic rows in to columns for ease of summarising down the line
  group_by(FD.metric,state.metric,system,lag)%>%
  pivot_longer(c(x_y.F:y_x.P),
               names_to = c("measure", ".value"),
               names_pattern = "(.*)\\.(.*)") %>%
  rename(F.value = F,P.value = P)%>% #rename problematic F and P column
  mutate(causality.direc = ifelse(grepl("^x_y",measure),"forward","reverse"))%>% # classify directions
  filter(P.value <= 0.05) # only keep significant (5% level) causality
  
pdf(file="Results/granger_causality/granger_causality_lag_changes.pdf",
    width=11, height = 5,onefile=FALSE)
ggpubr::ggarrange(ggplot(filter(gc.lag.changes.df,FD.metric %in% c("FDis","FEve","FRic")),
                         aes(x=lag,y=log(F.value), col =system)) +
  geom_point(alpha=0.4) + 
  ylab("Log F value") + xlab("Lag")+
  #ggh4x::facet_nested(state.metric ~ causality.direc + FD.metric ) + 
  ggh4x::facet_nested(causality.direc ~ FD.metric + state.metric ) + 
  #force_panelsizes(rows = c(0.5, 0.5),cols = 0.5) +
  scale_x_continuous(breaks = c(30,60))+
  scale_colour_manual(values = c("#FFE7A1","#A1B4FE","#74A180","#FF94AB","#BE86FF"),name= "Lake")+
  guides(color = guide_legend(override.aes = list(alpha = 1,size=1.5) ))+
  theme_bw() + ggtitle("Phytoplankton"),
  ggplot(filter(gc.lag.changes.df,FD.metric %in% c("zooFDis","zooFEve","zooFRic")),
         aes(x=lag,y=log(F.value), col =system)) +
    geom_point(alpha=0.4) + 
    ylab("Log F value") + xlab("Lag")+
    #ggh4x::facet_nested(state.metric ~ causality.direc + FD.metric ) + 
    ggh4x::facet_nested(causality.direc ~ FD.metric + state.metric ) + 
    #force_panelsizes(rows = c(0.5, 0.5),cols = 0.5) +
    scale_x_continuous(breaks = c(30,60))+
    scale_colour_manual(values = c("#FFE7A1","#A1B4FE","#74A180","#FF94AB","#BE86FF"),name= "Lake")+
    guides(color = guide_legend(override.aes = list(alpha = 1,size=1.5) ))+
    theme_bw()+ ggtitle("Zooplankton"),
  nrow=2,common.legend=T)
dev.off()

