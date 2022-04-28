### Figure 1 ###

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

# Merge lakes #

all.lakes.gam <- rbind(kin.tot,LZ.tot,mad.tot,wind.tot,kas.tot)%>%
  pivot_longer(-c(date,data.source,res),names_to = "metric",values_to = "value")%>%
  mutate(metric.type = ifelse(metric %in% c("FDis","FEve","FRic"),"Phytoplankton FD",
                              ifelse(metric %in% c("zooFDis","zooFEve","zooFRic"),"Zooplankton FD",
                                     ifelse(metric %in% c("env"),"Stress", "State Metric"))),
         metric.type = factor(metric.type,levels = c("State Metric","Phytoplankton FD","Zooplankton FD","Stress")))

###########################################################################
## Plot FD and State Trends ##
###########################################################################
no_env_dat <- filter(all.lakes.gam,metric != "env" )
pdf(file="Results/raw_visualisations/Figure1.pdf",
    width=9, height = 9)
ggplot(no_env_dat %>% mutate(metric = ifelse(metric %in% c("zooFDis","zooFEve","zooFRic"),substr(metric,4,7),metric)),aes(x=as.numeric(date),y=value, col = metric)) + 
  geom_point(aes(col=metric),alpha = 0.3,pch =21)+
  geom_smooth(aes(col = metric),method = "gam",formula =y ~ s(x, bs = "tp",k=15),method.args = list(method = "REML"),alpha=1,fill="grey") +
  ggh4x::facet_nested(metric.type + metric~data.source,scales = "free",
                      labeller = label_value,
                      strip = ggh4x::strip_nested(size="constant",bleed=T, 
                                                  background_y =elem_list_rect(fill = c("#D6D6D6","#EBEBEB","#EBEBEB",rep("white",12))),
                                                  background_x =elem_list_rect(fill = c(rep("white",5))),
                                                  by_layer_y = F),
                      space="free_x" ) +
  scale_colour_manual(values = c("#7b3294","#c2a5cf","#969014","#22B4F5","#a6dba0","#F07589","#008837","#E8E1A2"), guide = 'none')+
  scale_fill_manual(values = c("#7b3294","#c2a5cf","#969014","#22B4F5","#a6dba0","#F07589","#008837","#E8E1A2"), guide = 'none')+
  theme_bw() + xlab("Date")+ ylab("Scaled metric value")+
  geom_vline(data=no_env_dat %>% mutate(metric = ifelse(metric %in% c("zooFDis","zooFEve","zooFRic"),substr(metric,4,7),metric)) %>% filter(data.source == "Kasumigaura"), 
             aes(xintercept=1997.0), colour="black",linetype="dashed")+ #reported regime shifts
  geom_vline(data=no_env_dat %>% mutate(metric = ifelse(metric %in% c("zooFDis","zooFEve","zooFRic"),substr(metric,4,7),metric)) %>% filter(data.source == "Kinneret"), 
             aes(xintercept=1993.0), colour="black",linetype="dashed")+
  geom_vline(data=no_env_dat %>% mutate(metric = ifelse(metric %in% c("zooFDis","zooFEve","zooFRic"),substr(metric,4,7),metric)) %>% filter(data.source == "Mendota"), 
             aes(xintercept=2009.0), colour="black",linetype="dashed")+
  #scale_y_continuous(breaks= scales::pretty_breaks(n = 3))+
  coord_cartesian(ylim = c(-2.5,2.5))+ 
  scale_y_continuous(breaks= seq(-2,2,2))+
  #scale_x_continuous(breaks= scales::pretty_breaks(n = 3))+
  scale_x_continuous(breaks = seq(1970,2015,10))+
  theme(panel.spacing=unit(0.2,"lines"))
dev.off()