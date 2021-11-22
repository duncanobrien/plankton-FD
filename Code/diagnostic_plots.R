### Differenced FD Visualisation ###

## Preamble 
require(tidyverse) # dplyr, ggplot etc.
require(ggh4x) # facet_nested function

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
  #mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=1)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x))) # center and scale to unit variance for plotting

kin.diff1 <- cbind(phyto.kin.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$kin.mth[,-c(9)])%>%
  mutate(zooFDis =  zoo.kin.fuzFDs.mth[,"FDis"],zooFEve = zoo.kin.fuzFDs.mth[,"FEve"],zooFRic = zoo.kin.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=1)}))%>% # first difference (xt = xt - xt-1)
  mutate(across(-c(date,data.source,res),~scale(.x)))

kin.diff12 <- cbind(phyto.kin.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$kin.mth[,-c(9)])%>%
  mutate(zooFDis =  zoo.kin.fuzFDs.mth[,"FDis"],zooFEve = zoo.kin.fuzFDs.mth[,"FEve"],zooFRic = zoo.kin.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=12)}))%>% # first difference for monthly data (xt = xt - xt-12)
  mutate(across(-c(date,data.source,res),~scale(.x)))

kin.resid <- cbind(phyto.kin.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$kin.mth[,-c(9)])%>%
  mutate(zooFDis =  zoo.kin.fuzFDs.mth[,"FDis"],zooFEve = zoo.kin.fuzFDs.mth[,"FEve"],zooFRic = zoo.kin.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),~scale(.x))) %>%
  #mutate(across(-c(date,data.source,res),~residuals(lm(.x ~ as.numeric(date),na.action=na.exclude)))) #standardised residuals
   mutate(across(-c(date,data.source,res),function(.x){
     as.numeric(residuals(lm(.x ~ as.numeric(date),na.action=na.exclude)) - 
                  decompose(ts(.x,frequency = 12))$seasonal)})) #standardised residuals
  #mutate(across(-c(date,data.source,res),~pracma::detrend(.x))) #pracma linear detrend
  
mad.tot <- cbind(phyto.mad.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$mad.mth[,-c(9)])%>%
  mutate(zooFDis =  zoo.mad.fuzFDs.mth[,"FDis"],zooFEve = zoo.mad.fuzFDs.mth[,"FEve"],zooFRic = zoo.mad.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  #mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=1)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

mad.diff1 <- cbind(phyto.mad.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$mad.mth[,-c(9)])%>%
  mutate(zooFDis =  zoo.mad.fuzFDs.mth[,"FDis"],zooFEve = zoo.mad.fuzFDs.mth[,"FEve"],zooFRic = zoo.mad.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=1)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

mad.diff12 <- cbind(phyto.mad.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$mad.mth[,-c(9)])%>%
  mutate(zooFDis =  zoo.mad.fuzFDs.mth[,"FDis"],zooFEve = zoo.mad.fuzFDs.mth[,"FEve"],zooFRic = zoo.mad.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=12)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

mad.resid <- cbind(phyto.mad.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$mad.mth[,-c(9)])%>%
  mutate(zooFDis =  zoo.mad.fuzFDs.mth[,"FDis"],zooFEve = zoo.mad.fuzFDs.mth[,"FEve"],zooFRic = zoo.mad.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),~zoo::na.approx(.x,na.rm=F)))%>% #internal missing values
  mutate(across(-c(date,data.source,res),~scale(.x))) %>%
  #mutate(across(-c(date,data.source,res),~residuals(lm(.x ~ as.numeric(date),na.action=na.exclude)))) #standardised residuals
   mutate(across(-c(date,data.source,res),function(.x){
     as.numeric(residuals(lm(.x ~ as.numeric(date),na.action=na.exclude)) - 
                  decompose(ts(.x,frequency = 12))$seasonal)})) #standardised residuals  
  #mutate(across(-c(date,data.source,res),~pracma::detrend(.x))) #pracma linear detrend

LZ.tot <- cbind(phyto.LZ.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$LZ.mth[,-c(9)])%>%
  mutate(zooFDis =  zoo.LZ.fuzFDs.mth[,"FDis"],zooFEve = zoo.LZ.fuzFDs.mth[,"FEve"],zooFRic = zoo.LZ.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  #mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=1)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

LZ.diff1 <- cbind(phyto.LZ.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$LZ.mth[,-c(9)])%>%
  mutate(zooFDis =  zoo.LZ.fuzFDs.mth[,"FDis"],zooFEve = zoo.LZ.fuzFDs.mth[,"FEve"],zooFRic = zoo.LZ.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=1)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

LZ.diff12 <- cbind(phyto.LZ.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$LZ.mth[,-c(9)])%>%
  mutate(zooFDis =  zoo.LZ.fuzFDs.mth[,"FDis"],zooFEve = zoo.LZ.fuzFDs.mth[,"FEve"],zooFRic = zoo.LZ.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=12)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

LZ.resid <- cbind(phyto.LZ.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$LZ.mth[,-c(9)])%>%
  mutate(zooFDis =  zoo.LZ.fuzFDs.mth[,"FDis"],zooFEve = zoo.LZ.fuzFDs.mth[,"FEve"],zooFRic = zoo.LZ.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),~zoo::na.approx(.x,na.rm=F)))%>% #internal missing values
  mutate(across(-c(date,data.source,res),~scale(.x))) %>%
  #mutate(across(-c(date,data.source,res),~residuals(lm(.x ~ as.numeric(date),na.action=na.exclude)))) #standardised residuals
   mutate(across(-c(date,data.source,res),function(.x){
     as.numeric(residuals(lm(.x ~ as.numeric(date),na.action=na.exclude)) - 
                  decompose(ts(.x,frequency = 12))$seasonal)})) #standardised residuals
  #mutate(across(-c(date,data.source,res),~pracma::detrend(.x))) #pracma linear detrend

wind.tot <- cbind(phyto.wind.fuzFDs.mth[c("FDis","FEve","FRic")],all.system.states$wind.mth[,-c(9)])%>%
  mutate(zooFDis =  zoo.wind.fuzFDs.mth[,"FDis"],zooFEve = zoo.wind.fuzFDs.mth[,"FEve"],zooFRic = zoo.wind.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  #mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=1)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

wind.diff1 <- cbind(phyto.wind.fuzFDs.mth[c("FDis","FEve","FRic")],all.system.states$wind.mth[,-c(9)])%>%
  mutate(zooFDis =  zoo.wind.fuzFDs.mth[,"FDis"],zooFEve = zoo.wind.fuzFDs.mth[,"FEve"],zooFRic = zoo.wind.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=1)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

wind.diff12 <- cbind(phyto.wind.fuzFDs.mth[c("FDis","FEve","FRic")],all.system.states$wind.mth[,-c(9)])%>%
  mutate(zooFDis =  zoo.wind.fuzFDs.mth[,"FDis"],zooFEve = zoo.wind.fuzFDs.mth[,"FEve"],zooFRic = zoo.wind.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=12)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

wind.resid <- cbind(phyto.wind.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$wind.mth[,-c(9)])%>%
  mutate(zooFDis =  zoo.wind.fuzFDs.mth[,"FDis"],zooFEve = zoo.wind.fuzFDs.mth[,"FEve"],zooFRic = zoo.wind.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),~zoo::na.approx(.x,na.rm=F)))%>% #internal missing values
  mutate(across(-c(date,data.source,res),~scale(.x))) %>%
  #mutate(across(-c(date,data.source,res),~residuals(lm(.x ~ as.numeric(date),na.action=na.exclude)))) #standardised residuals
   mutate(across(-c(date,data.source,res),function(.x){
     as.numeric(residuals(lm(.x ~ as.numeric(date),na.action=na.exclude)) - 
                  decompose(ts(.x,frequency = 12))$seasonal)})) #standardised residuals
 #mutate(across(-c(date,data.source,res),~pracma::detrend(.x))) #pracma linear detrend

kas.tot <- cbind(phyto.kas.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$kas.mth[,-c(9)])%>%
  mutate(zooFDis =  zoo.kas.fuzFDs.mth[,"FDis"],zooFEve = zoo.kas.fuzFDs.mth[,"FEve"],zooFRic = zoo.kas.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  #mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=1)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

kas.diff1 <- cbind(phyto.kas.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$kas.mth[,-c(9)])%>%
  mutate(zooFDis =  zoo.kas.fuzFDs.mth[,"FDis"],zooFEve = zoo.kas.fuzFDs.mth[,"FEve"],zooFRic = zoo.kas.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=1)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

kas.diff12 <- cbind(phyto.kas.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$kas.mth[,-c(9)])%>%
  mutate(zooFDis =  zoo.kas.fuzFDs.mth[,"FDis"],zooFEve = zoo.kas.fuzFDs.mth[,"FEve"],zooFRic = zoo.kas.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=12)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

kas.resid <- cbind(phyto.kas.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$kas.mth[,-c(9)])%>%
  mutate(zooFDis =  zoo.kas.fuzFDs.mth[,"FDis"],zooFEve = zoo.kas.fuzFDs.mth[,"FEve"],zooFRic = zoo.kas.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),~zoo::na.approx(.x,na.rm=F)))%>% #internal missing values
  mutate(across(-c(date,data.source,res),~scale(.x))) %>%
  #mutate(across(-c(date,data.source,res),~residuals(lm(.x ~ as.numeric(date),na.action=na.exclude)))) #standardised residuals
   # mutate(across(-c(date,data.source,res),function(.x){
   #   as.numeric(residuals(lm(.x ~ as.numeric(date),na.action=na.exclude)) - 
   #     decompose(ts(.x,frequency = 12))$seasonal)})) #standardised residuals
  mutate(across(-c(date,data.source,res),function(.x){
               as.numeric(decompose(ts(.x,frequency = 12))$seasonal)}))
 #mutate(across(-c(date,data.source,res),~pracma::detrend(.x))) #pracma linear detrend

plot(decompose(ts(kin.tot$FDis,frequency = 12)))
# Merge lakes #
all.lakes.gam <- rbind(kin.tot,LZ.tot,mad.tot,wind.tot,kas.tot)%>%
  pivot_longer(-c(date,data.source,res),names_to = "metric",values_to = "value")%>%
  mutate(metric.type = ifelse(metric %in% c("FDis","FEve","FRic"),"Phytoplankton FD",
                              ifelse(metric %in% c("zooFDis","zooFEve","zooFRic"),"Zooplankton FD",
                                     ifelse(metric %in% c("env"),"Stress", "State Metric"))),
         metric.type = factor(metric.type,levels = c("State Metric","Phytoplankton FD","Zooplankton FD","Stress")))

all.lakes.diff1 <- rbind(kin.diff1,LZ.diff1,mad.diff1,wind.diff1,kas.diff1)%>%
  pivot_longer(-c(date,data.source,res),names_to = "metric",values_to = "value")%>%
  mutate(metric.type = ifelse(metric %in% c("FDis","FEve","FRic"),"Phytoplankton FD",
                              ifelse(metric %in% c("zooFDis","zooFEve","zooFRic"),"Zooplankton FD",
                                     ifelse(metric %in% c("env"),"Stress", "State Metric"))),
         metric.type = factor(metric.type,levels = c("State Metric","Phytoplankton FD","Zooplankton FD","Stress")))

all.lakes.diff12 <- rbind(kin.diff12,LZ.diff12,mad.diff12,wind.diff12,kas.diff12)%>%
  pivot_longer(-c(date,data.source,res),names_to = "metric",values_to = "value")%>%
  mutate(metric.type = ifelse(metric %in% c("FDis","FEve","FRic"),"Phytoplankton FD",
                              ifelse(metric %in% c("zooFDis","zooFEve","zooFRic"),"Zooplankton FD",
                                     ifelse(metric %in% c("env"),"Stress", "State Metric"))),
         metric.type = factor(metric.type,levels = c("State Metric","Phytoplankton FD","Zooplankton FD","Stress")))

all.lakes.resid <- rbind(kin.resid,mad.resid,LZ.resid,wind.resid,kas.resid)%>%
  pivot_longer(-c(date,data.source,res),names_to = "metric",values_to = "value")%>%
  mutate(metric.type = ifelse(metric %in% c("FDis","FEve","FRic"),"Phytoplankton FD",
                              ifelse(metric %in% c("zooFDis","zooFEve","zooFRic"),"Zooplankton FD",
                                     ifelse(metric %in% c("env"),"Stress", "State Metric"))),
         metric.type = factor(metric.type,levels = c("State Metric","Phytoplankton FD","Zooplankton FD","Stress")))

###########################################################################
## Diagnostic Plots ##
###########################################################################

# Distribution Skew #
hist.diff1 <- ggplot(all.lakes.diff1 %>% filter(metric.type != "State Metric" & metric.type != "Stress"),
                     aes(x=value)) +  coord_cartesian(ylim= c(0,100))+
  geom_histogram(bins = 50) + facet_grid(data.source~metric) + theme_bw() + 
  xlab("Scaled functional diversity value") + ylab("Count") + ggtitle("First differenced time series")
  
hist.diff12 <- ggplot(all.lakes.diff12 %>% filter(metric.type != "State Metric" & metric.type != "Stress"),
                      aes(x=value)) +
  geom_histogram(bins = 50) + facet_grid(data.source~metric) + theme_bw() + coord_cartesian(ylim= c(0,100))+
  xlab("Scaled functional diversity value") + ylab("Count") + ggtitle("Twelfth differenced time series")

hist.undiff <- ggplot(all.lakes.gam %>% filter(metric.type != "State Metric" & metric.type != "Stress"),
                      aes(x=value)) +  coord_cartesian(ylim= c(0,100))+
  geom_histogram(bins = 50) + facet_grid(data.source~metric) + theme_bw() + 
  xlab("Scaled functional diversity value") + ylab("Count") + ggtitle("Undifferenced time series")

hist.resid <- ggplot(all.lakes.resid %>% filter(metric.type != "State Metric" & metric.type != "Stress"),
                      aes(x=value)) +  coord_cartesian(ylim= c(0,100))+
  geom_histogram(bins = 50) + facet_grid(data.source~metric) + theme_bw() + 
  xlab("Scaled functional diversity value") + ylab("Count") + ggtitle("Time series standardised residuals")

pdf(file="Results/raw_visualisations/differenced_histograms.pdf",
    width=11, height = 7)
ggpubr::ggarrange(hist.undiff,hist.diff1,hist.diff12,hist.resid,labels = c("a","b","c","d"))
dev.off()

# Autocorrelation #

acf.undiff <-  rbind(kin.tot,LZ.tot,mad.tot,wind.tot,kas.tot) %>%
  dplyr::select(c(data.source,FDis,FEve,FRic,zooFDis,zooFEve,zooFRic))%>%
  group_by(data.source) %>%
  pivot_longer(-c(data.source), names_to = "metric",values_to ="value")%>%
  group_by(data.source,metric) %>%
  summarise(acf = acf(na.omit(value),plot = F)$acf,
            length =acf(na.omit(value),plot = F)$n.used)%>%
  mutate(lag = seq_along(1:n()))%>% 
  mutate(ciline = qnorm((1 + 0.95)/2) / sqrt(length))

acf.diff1 <- rbind(kin.diff1,LZ.diff1,mad.diff1,wind.diff1,kas.diff1) %>%
  dplyr::select(c(data.source,FDis,FEve,FRic,zooFDis,zooFEve,zooFRic))%>%
  group_by(data.source) %>%
  pivot_longer(-c(data.source), names_to = "metric",values_to ="value")%>%
  group_by(data.source,metric) %>%
  summarise(acf = acf(na.omit(value),plot = F)$acf,
            length =acf(na.omit(value),plot = F)$n.used)%>%
  mutate(lag = seq_along(1:n()))%>% 
  mutate(ciline = qnorm((1 + 0.95)/2) / sqrt(length))

acf.diff12 <- rbind(kin.diff12,LZ.diff12,mad.diff12,wind.diff12,kas.diff12) %>%
  dplyr::select(c(data.source,FDis,FEve,FRic,zooFDis,zooFEve,zooFRic))%>%
  group_by(data.source) %>%
  pivot_longer(-c(data.source), names_to = "metric",values_to ="value")%>%
  group_by(data.source,metric) %>%
  summarise(acf = acf(na.omit(value),plot = F)$acf,
           length =acf(na.omit(value),plot = F)$n.used)%>%
  mutate(lag = seq_along(1:n()))%>% 
  mutate(ciline = qnorm((1 + 0.95)/2) / sqrt(length))
  
acf.resid <- rbind(kin.resid,mad.resid,LZ.resid,wind.resid,kas.resid) %>%
  dplyr::select(c(data.source,FDis,FEve,FRic,zooFDis,zooFEve,zooFRic))%>%
  group_by(data.source) %>%
  pivot_longer(-c(data.source), names_to = "metric",values_to ="value")%>%
  group_by(data.source,metric) %>%
  summarise(acf = acf(na.omit(value),plot = F)$acf,
            length =acf(na.omit(value),plot = F)$n.used)%>%
  mutate(lag = seq_along(1:n()))%>% 
  mutate(ciline = qnorm((1 + 0.95)/2) / sqrt(length))

acf.undiff.plot <- ggplot(acf.undiff,
                         aes(x=lag,y=acf)) +  #coord_cartesian(ylim= c(0,100))+
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0))+
  geom_hline(aes(yintercept = ciline), linetype = 2, color = 'darkblue') + 
  geom_hline(aes(yintercept = -ciline), linetype = 2, color = 'darkblue') +
  facet_grid(data.source~metric) + theme_bw() + 
  xlab("Lag") + ylab("ACF") + ggtitle("Undifferenced time series ACF")

acf.diff1.plot <- ggplot(acf.diff1,
                     aes(x=lag,y=acf)) +  #coord_cartesian(ylim= c(0,100))+
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0))+
  geom_hline(aes(yintercept = ciline), linetype = 2, color = 'darkblue') + 
  geom_hline(aes(yintercept = -ciline), linetype = 2, color = 'darkblue') +
  facet_grid(data.source~metric) + theme_bw() + 
  xlab("Lag") + ylab("ACF") + ggtitle("First differenced time series ACF")

acf.diff12.plot <- ggplot(acf.diff12,
                           aes(x=lag,y=acf)) +  #coord_cartesian(ylim= c(0,100))+
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0))+
    geom_hline(aes(yintercept = ciline), linetype = 2, color = 'darkblue') + 
    geom_hline(aes(yintercept = -ciline), linetype = 2, color = 'darkblue') +
  facet_grid(data.source~metric) + theme_bw() + 
    xlab("Lag") + ylab("ACF") + ggtitle("Twelfth differenced time series ACF")

acf.resid.plot <- ggplot(acf.resid,
                          aes(x=lag,y=acf)) +  #coord_cartesian(ylim= c(0,100))+
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0))+
  geom_hline(aes(yintercept = ciline), linetype = 2, color = 'darkblue') + 
  geom_hline(aes(yintercept = -ciline), linetype = 2, color = 'darkblue') +
  facet_grid(data.source~metric) + theme_bw() + 
  xlab("Lag") + ylab("ACF") + ggtitle("Time series standardised residuals ACF")

pdf(file="Results/raw_visualisations/differenced_autocorrelation.pdf",
    width=11, height = 7)
ggpubr::ggarrange(acf.undiff.plot,acf.diff1.plot,acf.diff12.plot,acf.resid.plot,labels = c("a","b","c","d"))
dev.off()
