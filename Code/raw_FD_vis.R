###########################################################################
## Raw FD Visualisation ##
###########################################################################

kin.tot <- cbind(phyto.kin.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$kin.mth[,-8])%>%
  mutate(zooFDis =  zoo.kin.fuzFDs.mth[,"FDis"],zooFEve = zoo.kin.fuzFDs.mth[,"FEve"],zooFRic = zoo.kin.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  #mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=1)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

kin.diff1 <- cbind(phyto.kin.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$kin.mth[,-8])%>%
  mutate(zooFDis =  zoo.kin.fuzFDs.mth[,"FDis"],zooFEve = zoo.kin.fuzFDs.mth[,"FEve"],zooFRic = zoo.kin.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=1)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

kin.diff12 <- cbind(phyto.kin.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$kin.mth[,-8])%>%
  mutate(zooFDis =  zoo.kin.fuzFDs.mth[,"FDis"],zooFEve = zoo.kin.fuzFDs.mth[,"FEve"],zooFRic = zoo.kin.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=12)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

mad.tot <- cbind(phyto.mad.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$mad.mth[,-8])%>%
  mutate(zooFDis =  zoo.mad.fuzFDs.mth[,"FDis"],zooFEve = zoo.mad.fuzFDs.mth[,"FEve"],zooFRic = zoo.mad.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  #mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=1)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

mad.diff1 <- cbind(phyto.mad.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$mad.mth[,-8])%>%
  mutate(zooFDis =  zoo.mad.fuzFDs.mth[,"FDis"],zooFEve = zoo.mad.fuzFDs.mth[,"FEve"],zooFRic = zoo.mad.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=1)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

mad.diff12 <- cbind(phyto.mad.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$mad.mth[,-8])%>%
  mutate(zooFDis =  zoo.mad.fuzFDs.mth[,"FDis"],zooFEve = zoo.mad.fuzFDs.mth[,"FEve"],zooFRic = zoo.mad.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=12)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

LZ.tot <- cbind(phyto.LZ.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$LZ.mth[,-8])%>%
  mutate(zooFDis =  zoo.LZ.fuzFDs.mth[,"FDis"],zooFEve = zoo.LZ.fuzFDs.mth[,"FEve"],zooFRic = zoo.LZ.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  #mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=1)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

LZ.diff1 <- cbind(phyto.LZ.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$LZ.mth[,-8])%>%
  mutate(zooFDis =  zoo.LZ.fuzFDs.mth[,"FDis"],zooFEve = zoo.LZ.fuzFDs.mth[,"FEve"],zooFRic = zoo.LZ.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=1)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

LZ.diff12 <- cbind(phyto.LZ.fuzFDs.mth[,c("FDis","FEve","FRic")],all.system.states$LZ.mth[,-8])%>%
  mutate(zooFDis =  zoo.LZ.fuzFDs.mth[,"FDis"],zooFEve = zoo.LZ.fuzFDs.mth[,"FEve"],zooFRic = zoo.LZ.fuzFDs.mth[,"FRic"])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=12)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

wind.tot <- cbind(phyto.wind.fuzFDs.mth[c("FDis","FEve","FRic")],all.system.states$wind.mth[,-8])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  # mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=1)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))%>%
  mutate(zooFDis = NA, zooFEve = NA, zooFRic = NA)

wind.diff1 <- cbind(phyto.wind.fuzFDs.mth[c("FDis","FEve","FRic")],all.system.states$wind.mth[,-8])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=1)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))%>%
  mutate(zooFDis = NA, zooFEve = NA, zooFRic = NA)

wind.diff12 <- cbind(phyto.wind.fuzFDs.mth[c("FDis","FEve","FRic")],all.system.states$wind.mth[,-8])%>%
  mutate(across(c(density,mvi,zp.ratio),~log(.x)))%>%
  mutate(across(-c(date,data.source,res),function(x){x - dplyr::lag(x,n=12)}))%>%
  mutate(across(-c(date,data.source,res),~scale(.x)))

all.lakes.gam <- rbind(kin.tot,LZ.tot,mad.tot,wind.tot)%>%
  pivot_longer(-c(date,data.source,res),names_to = "metric",values_to = "value")%>%
  mutate(metric.type = ifelse(metric %in% c("FDis","FEve","FRic"),"Phytoplankton FD",
                              ifelse(metric %in% c("zooFDis","zooFEve","zooFRic"),"Zooplankton FD","State Metric")),
         metric.type = factor(metric.type,levels = c("State Metric","Phytoplankton FD","Zooplankton FD")))

png(file="Results/raw_visualisations/raw_smooth_vis.png",
    width=3000, height = 3000 ,res=300)
ggplot(all.lakes.gam %>% mutate(metric = ifelse(metric %in% c("zooFDis","zooFEve","zooFRic"),substr(metric,4,7),metric)),aes(x=as.numeric(date),y=value, col = metric)) + 
  geom_smooth(aes(fill = metric),method = "gam",alpha=0.3) +
  #facet_grid(metric.type + metric~data.source,scales = "free_x")+
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
  theme(#axis.title.y=element_blank(),
    #axis.text.y=element_blank(),
    #axis.ticks.y=element_blank(),
    panel.spacing = unit(0.2,"line"))
dev.off()

all.lakes.diff1 <- rbind(kin.diff1,LZ.diff1,mad.diff1,wind.diff1)%>%
  pivot_longer(-c(date,data.source,res),names_to = "metric",values_to = "value")%>%
  mutate(metric.type = ifelse(metric %in% c("FDis","FEve","FRic"),"Phytoplankton FD",
                              ifelse(metric %in% c("zooFDis","zooFEve","zooFRic"),"Zooplankton FD","State Metric")),
         metric.type = factor(metric.type,levels = c("State Metric","Phytoplankton FD","Zooplankton FD")))

png(file="Results/raw_visualisations/raw_diff1_vis.png",
    width=3000, height = 3000 ,res=300)
ggplot(all.lakes.diff1 %>% mutate(metric = ifelse(metric %in% c("zooFDis","zooFEve","zooFRic"),substr(metric,4,7),metric)),aes(x=as.numeric(date),y=value, col = metric)) + 
  geom_path() +
  #facet_grid(metric.type + metric~data.source,scales = "free_x")+
  ggh4x::facet_nested(metric.type + metric~data.source,scales = "free",
                      labeller = label_value,strip = strip_nested(size="constant",bleed=T),
                      space="free_x") +
  scale_colour_manual(values = c("#7b3294","#c2a5cf","#969014","#22B4F5","#a6dba0","#F07589","#008837","#E8E1A2"), guide = 'none')+
  theme_bw() + xlab("Date")+ ylab("Scaled metric value")+
  geom_vline(data=all.lakes.diff1 %>% mutate(metric = ifelse(metric %in% c("zooFDis","zooFEve","zooFRic"),substr(metric,4,7),metric)) %>% filter(data.source == "Kinneret"), 
             aes(xintercept=1993.3), colour="black",linetype="dashed")+
  #scale_y_continuous(breaks= scales::pretty_breaks(n = 3))+
  scale_y_continuous(breaks= seq(-6,6,4))+
  #scale_x_continuous(breaks= scales::pretty_breaks(n = 3))+
  scale_x_continuous(breaks = seq(1970,2015,10))+
  theme(#axis.title.y=element_blank(),
    #axis.text.y=element_blank(),
    #axis.ticks.y=element_blank(),
    panel.spacing = unit(0.2,"line"))
dev.off()

all.lakes.diff1 <- rbind(kin.diff12,LZ.diff12,mad.diff12,wind.diff12)%>%
  pivot_longer(-c(date,data.source,res),names_to = "metric",values_to = "value")%>%
  mutate(metric.type = ifelse(metric %in% c("FDis","FEve","FRic"),"Phytoplankton FD",
                              ifelse(metric %in% c("zooFDis","zooFEve","zooFRic"),"Zooplankton FD","State Metric")),
         metric.type = factor(metric.type,levels = c("State Metric","Phytoplankton FD","Zooplankton FD")))

png(file="Results/raw_visualisations/raw_diff12_vis.png",
    width=3000, height = 3000 ,res=300)
ggplot(all.lakes.diff12 %>% mutate(metric = ifelse(metric %in% c("zooFDis","zooFEve","zooFRic"),substr(metric,4,7),metric)),aes(x=as.numeric(date),y=value, col = metric)) + 
  geom_path() +
  #facet_grid(metric.type + metric~data.source,scales = "free_x")+
  ggh4x::facet_nested(metric.type + metric~data.source,scales = "free",
                      labeller = label_value,strip = strip_nested(size="constant",bleed=T),
                      space="free_x") +
  scale_colour_manual(values = c("#7b3294","#c2a5cf","#969014","#22B4F5","#a6dba0","#F07589","#008837","#E8E1A2"), guide = 'none')+
  theme_bw() + xlab("Date")+ ylab("Scaled metric value")+
  geom_vline(data=all.lakes.diff1 %>% mutate(metric = ifelse(metric %in% c("zooFDis","zooFEve","zooFRic"),substr(metric,4,7),metric)) %>% filter(data.source == "Kinneret"), 
             aes(xintercept=1993.3), colour="black",linetype="dashed")+
  #scale_y_continuous(breaks= scales::pretty_breaks(n = 3))+
  scale_y_continuous(breaks= seq(-6,6,4))+
  #scale_x_continuous(breaks= scales::pretty_breaks(n = 3))+
  scale_x_continuous(breaks = seq(1970,2015,10))+
  theme(#axis.title.y=element_blank(),
    #axis.text.y=element_blank(),
    #axis.ticks.y=element_blank(),
    panel.spacing = unit(0.2,"line"))
dev.off()
