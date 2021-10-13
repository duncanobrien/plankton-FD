### Raw FD Visualisation ###

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

zoo.kin.fuzFDs.mth <- read.csv("Data/raw_FD/FD_kin_zoo_mth_raw.csv")
zoo.LZ.fuzFDs.mth <- read.csv("Data/raw_FD/FD_LZ_zoo_mth_raw.csv")
zoo.mad.fuzFDs.mth <- read.csv("Data/raw_FD/FD_mad_zoo_mth_raw.csv")

load("Data/all.system.states.RData")

###########################################################################
## Raw FD Prep ##
###########################################################################

kin.tot <- cbind(phyto.kin.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$kin.mth[,-c(7,9)])%>%
  mutate(zooFDis =  zoo.kin.fuzFDs.mth[,"FDis"],zooFEve = zoo.kin.fuzFDs.mth[,"FEve"],zooFRic = zoo.kin.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>% # log density, mvi and zo.ration to linearise
  #mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=1)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x))) # center and scale to unit variance for plotting

kin.diff1 <- cbind(phyto.kin.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$kin.mth[,-c(7,9)])%>%
  mutate(zooFDis =  zoo.kin.fuzFDs.mth[,"FDis"],zooFEve = zoo.kin.fuzFDs.mth[,"FEve"],zooFRic = zoo.kin.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=1)}))%>% # first difference (xt = xt - xt-1)
  mutate(across(-c(date,data.source,res),~scale(.x)))

kin.diff12 <- cbind(phyto.kin.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$kin.mth[,-c(7,9)])%>%
  mutate(zooFDis =  zoo.kin.fuzFDs.mth[,"FDis"],zooFEve = zoo.kin.fuzFDs.mth[,"FEve"],zooFRic = zoo.kin.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=12)}))%>% # first difference for monthly data (xt = xt - xt-12)
  mutate(across(-c(date,data.source,res),~scale(.x)))

mad.tot <- cbind(phyto.mad.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$mad.mth[,-c(7,9)])%>%
  mutate(zooFDis =  zoo.mad.fuzFDs.mth[,"FDis"],zooFEve = zoo.mad.fuzFDs.mth[,"FEve"],zooFRic = zoo.mad.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  #mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=1)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

mad.diff1 <- cbind(phyto.mad.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$mad.mth[,-c(7,9)])%>%
  mutate(zooFDis =  zoo.mad.fuzFDs.mth[,"FDis"],zooFEve = zoo.mad.fuzFDs.mth[,"FEve"],zooFRic = zoo.mad.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=1)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

mad.diff12 <- cbind(phyto.mad.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$mad.mth[,-c(7,9)])%>%
  mutate(zooFDis =  zoo.mad.fuzFDs.mth[,"FDis"],zooFEve = zoo.mad.fuzFDs.mth[,"FEve"],zooFRic = zoo.mad.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=12)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

LZ.tot <- cbind(phyto.LZ.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$LZ.mth[,-c(7,9)])%>%
  mutate(zooFDis =  zoo.LZ.fuzFDs.mth[,"FDis"],zooFEve = zoo.LZ.fuzFDs.mth[,"FEve"],zooFRic = zoo.LZ.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  #mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=1)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

LZ.diff1 <- cbind(phyto.LZ.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$LZ.mth[,-c(7,9)])%>%
  mutate(zooFDis =  zoo.LZ.fuzFDs.mth[,"FDis"],zooFEve = zoo.LZ.fuzFDs.mth[,"FEve"],zooFRic = zoo.LZ.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=1)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

LZ.diff12 <- cbind(phyto.LZ.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$LZ.mth[,-c(7,9)])%>%
  mutate(zooFDis =  zoo.LZ.fuzFDs.mth[,"FDis"],zooFEve = zoo.LZ.fuzFDs.mth[,"FEve"],zooFRic = zoo.LZ.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=12)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

wind.tot <- cbind(phyto.wind.fuzFDs.mth[c("FDis","FEve","FRic")],all.system.states$wind.mth[,-c(7,9)])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  # mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=1)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))%>%
  mutate(zooFDis = NA, zooFEve = NA, zooFRic = NA) # add dummy zooFD variable for missing data

wind.diff1 <- cbind(phyto.wind.fuzFDs.mth[c("FDis","FEve","FRic")],all.system.states$wind.mth[,-c(7,9)])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=1)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))%>%
  mutate(zooFDis = NA, zooFEve = NA, zooFRic = NA) # add dummy zooFD variable for missing data

wind.diff12 <- cbind(phyto.wind.fuzFDs.mth[c("FDis","FEve","FRic")],all.system.states$wind.mth[,-c(7,9)])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=12)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))%>%
  mutate(zooFDis = NA, zooFEve = NA, zooFRic = NA) # add dummy zooFD variable for missing data

# Merge lakes #

all.lakes.gam <- rbind(kin.tot,LZ.tot,mad.tot,wind.tot)%>%
  pivot_longer(-c(date,data.source,res),names_to = "metric",values_to = "value")%>%
  mutate(metric.type = ifelse(metric %in% c("FDis","FEve","FRic"),"Phytoplankton FD",
                              ifelse(metric %in% c("zooFDis","zooFEve","zooFRic"),"Zooplankton FD","State Metric")),
         metric.type = factor(metric.type,levels = c("State Metric","Phytoplankton FD","Zooplankton FD")))

all.lakes.diff1 <- rbind(kin.diff1,LZ.diff1,mad.diff1,wind.diff1)%>%
  pivot_longer(-c(date,data.source,res),names_to = "metric",values_to = "value")%>%
  mutate(metric.type = ifelse(metric %in% c("FDis","FEve","FRic"),"Phytoplankton FD",
                              ifelse(metric %in% c("zooFDis","zooFEve","zooFRic"),"Zooplankton FD","State Metric")),
         metric.type = factor(metric.type,levels = c("State Metric","Phytoplankton FD","Zooplankton FD")))

all.lakes.diff12 <- rbind(kin.diff12,LZ.diff12,mad.diff12,wind.diff12)%>%
  pivot_longer(-c(date,data.source,res),names_to = "metric",values_to = "value")%>%
  mutate(metric.type = ifelse(metric %in% c("FDis","FEve","FRic"),"Phytoplankton FD",
                              ifelse(metric %in% c("zooFDis","zooFEve","zooFRic"),"Zooplankton FD","State Metric")),
         metric.type = factor(metric.type,levels = c("State Metric","Phytoplankton FD","Zooplankton FD")))

###########################################################################
## Plot FD and State Trends ##
###########################################################################

pdf(file="Results/raw_visualisations/raw_smooth_vis.pdf",
    width=9, height = 9)
ggplot(all.lakes.gam %>% mutate(metric = ifelse(metric %in% c("zooFDis","zooFEve","zooFRic"),substr(metric,4,7),metric)),aes(x=as.numeric(date),y=value, col = metric)) + 
  geom_smooth(aes(fill = metric),method = "gam",alpha=0.3) +
  #geom_path(aes(col=metric))+
  ggh4x::facet_nested(metric.type + metric~data.source,scales = "free",
                      labeller = label_value,strip = ggh4x::strip_nested(size="constant",bleed=T),
                      space="free_x" ) +
  scale_colour_manual(values = c("#7b3294","#c2a5cf","#969014","#22B4F5","#a6dba0","#F07589","#008837","#E8E1A2"), guide = 'none')+
  scale_fill_manual(values = c("#7b3294","#c2a5cf","#969014","#22B4F5","#a6dba0","#F07589","#008837","#E8E1A2"), guide = 'none')+
  theme_bw() + xlab("Date")+ ylab("Scaled metric value")+
  geom_vline(data=all.lakes.gam %>% mutate(metric = ifelse(metric %in% c("zooFDis","zooFEve","zooFRic"),substr(metric,4,7),metric)) %>% filter(data.source == "Kinneret"), 
             aes(xintercept=1993.3), colour="black",linetype="dashed")+
  #scale_y_continuous(breaks= scales::pretty_breaks(n = 3))+
  scale_y_continuous(breaks= seq(-1,1,1))+
  #scale_x_continuous(breaks= scales::pretty_breaks(n = 3))+
  scale_x_continuous(breaks = seq(1970,2015,10))+
  theme(panel.spacing=unit(0.2,"lines"))
dev.off()

pdf(file="Results/raw_visualisations/raw_diff1_vis.pdf",
    width=9, height = 9)
ggplot(all.lakes.diff1 %>% mutate(metric = ifelse(metric %in% c("zooFDis","zooFEve","zooFRic"),substr(metric,4,7),metric)),aes(x=as.numeric(date),y=value, col = metric)) + 
  geom_path() +
  ggh4x::facet_nested(metric.type + metric~data.source,scales = "free",
                      labeller = label_value,strip = ggh4x::strip_nested(size="constant",bleed=T),
                      space="free_x") +
  scale_colour_manual(values = c("#7b3294","#c2a5cf","#969014","#22B4F5","#a6dba0","#F07589","#008837","#E8E1A2"), guide = 'none')+
  theme_bw() + xlab("Date")+ ylab("Scaled metric value")+
  geom_vline(data=all.lakes.diff1 %>% mutate(metric = ifelse(metric %in% c("zooFDis","zooFEve","zooFRic"),substr(metric,4,7),metric)) %>% filter(data.source == "Kinneret"), 
             aes(xintercept=1993.3), colour="black",linetype="dashed")+
  #scale_y_continuous(breaks= scales::pretty_breaks(n = 3))+
  scale_y_continuous(breaks= seq(-6,6,4))+
  #scale_x_continuous(breaks= scales::pretty_breaks(n = 3))+
  scale_x_continuous(breaks = seq(1970,2015,10))+
  theme(panel.spacing = unit(0.2,"line"))
dev.off()

pdf(file="Results/raw_visualisations/raw_diff12_vis.pdf",
    width=9, height = 9)
ggplot(all.lakes.diff12 %>% mutate(metric = ifelse(metric %in% c("zooFDis","zooFEve","zooFRic"),substr(metric,4,7),metric)),aes(x=as.numeric(date),y=value, col = metric)) + 
  geom_path() +
  ggh4x::facet_nested(metric.type + metric~data.source,scales = "free",
                      labeller = label_value,strip = ggh4x::strip_nested(size="constant",bleed=T),
                      space="free_x") +
  scale_colour_manual(values = c("#7b3294","#c2a5cf","#969014","#22B4F5","#a6dba0","#F07589","#008837","#E8E1A2"), guide = 'none')+
  theme_bw() + xlab("Date")+ ylab("Scaled metric value")+
  geom_vline(data=all.lakes.diff1 %>% mutate(metric = ifelse(metric %in% c("zooFDis","zooFEve","zooFRic"),substr(metric,4,7),metric)) %>% filter(data.source == "Kinneret"), 
             aes(xintercept=1993.3), colour="black",linetype="dashed")+
  #scale_y_continuous(breaks= scales::pretty_breaks(n = 3))+
  scale_y_continuous(breaks= seq(-6,6,4))+
  #scale_x_continuous(breaks= scales::pretty_breaks(n = 3))+
  scale_x_continuous(breaks = seq(1970,2015,10))+
  theme(panel.spacing = unit(0.2,"line"))
dev.off()
