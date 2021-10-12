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

zoo.kin.fuzFDs.mth <- read.csv("Data/raw_FD/FD_kin_zoo_mth_raw.csv")
zoo.LZ.fuzFDs.mth <- read.csv("Data/raw_FD/FD_LZ_zoo_mth_raw.csv")
zoo.mad.fuzFDs.mth <- read.csv("Data/raw_FD/FD_mad_zoo_mth_raw.csv")

load("Data/all.system.states.RData")

###########################################################################
## Raw FD Prep ##
###########################################################################

kin.tot <- cbind(phyto.kin.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$kin.mth[,-8])%>%
  mutate(zooFDis =  zoo.kin.fuzFDs.mth[,"FDis"],zooFEve = zoo.kin.fuzFDs.mth[,"FEve"],zooFRic = zoo.kin.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>% # log density, mvi and zo.ration to linearise
  mutate(across(-c(date,data.source,res),~scale(.x))) # center and scale to unit variance for plotting

mad.tot <- cbind(phyto.mad.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$mad.mth[,-8])%>%
  mutate(zooFDis =  zoo.mad.fuzFDs.mth[,"FDis"],zooFEve = zoo.mad.fuzFDs.mth[,"FEve"],zooFRic = zoo.mad.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

LZ.tot <- cbind(phyto.LZ.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$LZ.mth[,-8])%>%
  mutate(zooFDis =  zoo.LZ.fuzFDs.mth[,"FDis"],zooFEve = zoo.LZ.fuzFDs.mth[,"FEve"],zooFRic = zoo.LZ.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

wind.tot <- cbind(phyto.wind.fuzFDs.mth[c("FDis","FEve","FRic")],all.system.states$wind.mth[,-8])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))%>%
  mutate(zooFDis = NA, zooFEve = NA, zooFRic = NA) # add dummy zooFD variable for missing data

###########################################################################
## Estimate Granger Causality ##
###########################################################################

kin.granger <- pbmcapply::pbmclapply(c("FDis","FEve","FRic","zooFDis","zooFEve","zooFRic"),function(x){
  
  tmp <- foreach::foreach(i = c("density","community","FI","mvi","zp.ratio"),.combine = "rbind",.export = c("cross.granger"),.packages = c("vars","dplyr")) %do%{
    
    if(any(is.na(kin.tot[,paste(x)]))){
      kin.tot[,paste(x)] <- zoo::na.approx(kin.tot[,paste(x)],maxgap =3)
    } # certain FD metrics may have NAs for single data points
    
    if(i %in% c("FI","mvi")){
      df <- na.omit(kin.tot) #drop the first few rows of NAs caused by the window calculation of these metrics
    }else{
      df <- kin.tot
    }
    gc.df <- cross.granger(ts = df[,paste(x)], comp.ts =df[,paste(i)],span = 12*5,method="var")
    gc.df$state.metric <- paste(i)
    return(gc.df)
}
  return(tmp)
  },mc.cores = 1) 
names(kin.granger) <- c("FDis","FEve","FRic","zooFDis","zooFEve","zooFRic")
kin.granger <- data.frame(data.table::rbindlist(kin.granger,idcol = "FD.metric"))%>%
  mutate(system = "Kinneret", res = "Month") # specify metadata for future plotting

mad.granger <- pbmcapply::pbmclapply(c("FDis","FEve","FRic","zooFDis","zooFEve","zooFRic"),function(x){
  
  tmp <- foreach::foreach(i = c("density","community","FI","mvi","zp.ratio"),.combine = "rbind",.export = c("cross.granger"),.packages = c("vars","dplyr")) %do%{
    
    if(any(is.na(mad.tot[,paste(x)]))){
      mad.tot[,paste(x)] <- zoo::na.approx(mad.tot[,paste(x)],maxgap =3)
    } # certain FD metrics may have NAs for single data points
    
    if(i %in% c("FI","mvi")){
      df <- na.omit(mad.tot) #drop the first few rows of NAs caused by the window calculation of these metrics
    }else{
      df <- mad.tot
    }
    gc.df <- cross.granger(ts = df[,paste(x)], comp.ts =df[,paste(i)],span = 12*5,method="var")
    gc.df$state.metric <- paste(i)
    return(gc.df)
  }
  return(tmp)
},mc.cores = 1) 
names(mad.granger) <- c("FDis","FEve","FRic","zooFDis","zooFEve","zooFRic")
mad.granger <- data.frame(data.table::rbindlist(mad.granger,idcol = "FD.metric"))%>%
  mutate(system = "Mendota", res = "Month") # specify metadata for future plotting

LZ.granger <- pbmcapply::pbmclapply(c("FDis","FEve","FRic","zooFDis","zooFEve","zooFRic"),function(x){
  
  tmp <- foreach::foreach(i = c("density","community","FI","mvi","zp.ratio"),.combine = "rbind",.export = c("cross.granger"),.packages = c("vars","dplyr")) %do%{
    
    if(any(is.na(LZ.tot[,paste(x)]))){
      LZ.tot[,paste(x)] <- zoo::na.approx(LZ.tot[,paste(x)],maxgap =3)
    } # certain FD metrics may have NAs for single data points
    
    if(i %in% c("FI","mvi")){
      df <- na.omit(LZ.tot) #drop the first few rows of NAs caused by the window calculation of these metrics
    }else{
      df <- LZ.tot
    }
    gc.df <- cross.granger(ts = df[,paste(x)], comp.ts =df[,paste(i)],span = 12*5,method="var")
    gc.df$state.metric <- paste(i)
    return(gc.df)
  }
  return(tmp)
},mc.cores = 1) 
names(LZ.granger) <- c("FDis","FEve","FRic","zooFDis","zooFEve","zooFRic")
LZ.granger <- data.frame(data.table::rbindlist(LZ.granger,idcol = "FD.metric"))%>%
  mutate(system = "Lower Zurich", res = "Month") # specify metadata for future plotting

wind.granger <- pbmcapply::pbmclapply(c("FDis","FEve","FRic"),function(x){
  
  tmp <- foreach::foreach(i = c("density","community","FI","mvi","zp.ratio"),.combine = "rbind",.export = c("cross.granger"),.packages = c("vars","dplyr")) %do%{
    
    if(any(is.na(wind.tot[,paste(x)]))){
      wind.tot[,paste(x)] <- zoo::na.approx(wind.tot[,paste(x)],maxgap =3)
    } # certain FD metrics may have NAs for single data points
    
    if(i %in% c("FI","mvi")){
      df <- na.omit(wind.tot[,-c(12:14)]) #drop the first few rows of NAs caused by the window calculation of these metrics
    }else{
      df <- wind.tot
    }
    gc.df <- cross.granger(ts = df[,paste(x)], comp.ts =df[,paste(i)],span = 12*5,method="var")
    gc.df$state.metric <- paste(i)
    return(gc.df)
  }
  return(tmp)
},mc.cores = 1) 
names(wind.granger) <- c("FDis","FEve","FRic")
wind.granger <- data.frame(data.table::rbindlist(wind.granger,idcol = "FD.metric"))%>%
  mutate(system = "Windermere", res = "Month")
  # specify metadata for future plotting

###########################################################################
## Save out ##
###########################################################################

all.lakes.granger <- rbind(kin.granger,mad.granger,LZ.granger,wind.granger)
write.csv(all.lakes.granger,file ="Results/all.lakes.granger.csv",row.names = F)

all.lakes.granger <- read.csv("Results/all.lakes.granger.csv")

gc.best.lag.df <- all.lakes.granger %>%
  filter(!grepl("\\.aic$",statistic) & !grepl("\\.sig$",statistic))%>% # drop AIC and sig rows
  group_by(FD.metric,state.metric,system) %>%
  pivot_longer(c(lag_1:lag_60),names_to = "lag",values_to = "value",names_prefix = "lag_") %>% # pivot lag columns in to single column
  mutate(value = as.numeric(value),lag = as.numeric(lag))%>% # ensure values are numeric
  pivot_wider(names_from = statistic,values_from = value)%>% # pivot F/P statistic rows in to columns for ease of summarising down the line
  group_by(FD.metric,state.metric,system,lag)%>%
  pivot_longer(c(x_y.F,y_x.F),names_to = "F.measure",values_to = "F.value")%>% # pivot back down in to separate columns 
  pivot_longer(c(x_y.P,y_x.P),names_to = "P.measure",values_to = "P.value")%>%
  filter(substr(F.measure,1,3) == substr(P.measure,1,3)) %>% # match forward & reverse measures
  mutate(causality.direc = ifelse(grepl("^x_y",F.measure),"forward","reverse"))%>% # classify directions
  filter(P.value <= 0.05)%>% # only keep significant (5% level) causality
  group_by(FD.metric,state.metric,system,causality.direc)%>%
  mutate(best.value = max(F.value))%>% # find highest F value per group but keep the associated P value (was lost if not pivotted in to a separate column)
  ungroup()%>%
  filter(best.value == F.value) # filter data frame to lags with highest Granger causality

count.df <- gc.best.lag.df %>%
  group_by(state.metric,causality.direc,FD.metric)%>%
  mutate(N = n())

pdf(file="Results/granger_causality_spread.pdf",
    width=8, height = 7)
ggplot(gc.best.lag.df,aes(x= state.metric,y= lag,fill= causality.direc)) + 
  geom_boxplot(alpha = 0.8,size = 0.5)+
  geom_point(aes(y=lag, group=causality.direc), alpha = 0.5, position = position_dodge(width=0.75),size = 1.3) + 
  scale_fill_manual(values = c("#FFE7A1","#A1B4FE"),name = "Causality\ndirection",labels = c("Forward", "Reverse"))+
  geom_text(data = count.df %>% distinct(FD.metric, causality.direc, state.metric, N),
           aes(y = 65, label = N),
          position = position_dodge(width = 0.8))+
  facet_wrap(~FD.metric) +
  xlab("State metric") + ylab("Optimum lag (months)")+
  theme_bw()
dev.off()
