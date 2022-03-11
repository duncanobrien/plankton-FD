### Summary Plots ###

## Preamble ##
require(tidyverse) # dplyr, ggplot etc.
require(patchwork) # plot alignment
require(ggpubr) # plot alignment

## Prepare data ##
summary.ccf.mth1 <- read.csv("Results/ccf/raw_data/summary.ccf.mth.lag1.csv")
summary.ccm <- read.csv(file ="Results/ccm/raw_data/ccm_summary.csv")

ccm.lag0.comp <- summary.ccm %>%
  filter(measure == "r0.skill")%>%  
  mutate(forward = ifelse(y_x.sig == "*" & x_y.sig != "*",TRUE,FALSE),
         reverse = ifelse(x_y.sig == "*" & y_x.sig != "*",TRUE,FALSE),
         bidirec =  ifelse(x_y.sig == "*" & y_x.sig == "*",TRUE,FALSE),
         none =  ifelse(x_y.sig != "*" & y_x.sig != "*",TRUE,FALSE))%>%
  group_by(troph,FD.metric,state.metric) %>%
  summarise(prop.forward=sum(forward == TRUE)/length(forward),
            prop.reverse=sum(reverse == TRUE)/length(reverse),
            prop.bidirec=sum(bidirec == TRUE)/length(bidirec),
            prop.none=sum(none == TRUE)/length(none))%>%
  mutate(measure = "r0.skill")

ccm.lagx.comp <- summary.ccm %>%
  filter(measure == "max.skill")%>%  
  mutate(forward = ifelse(y_x.sig == "*" & x_y.sig != "*",TRUE,FALSE),
         reverse = ifelse(x_y.sig == "*" & y_x.sig != "*",TRUE,FALSE),
         bidirec =  ifelse(x_y.sig == "*" & y_x.sig == "*",TRUE,FALSE),
         none =  ifelse(x_y.sig != "*" & y_x.sig != "*",TRUE,FALSE),
         diff.lag = filter(summary.ccm,measure == "t.max.skill")$x_y.obs_value - filter(summary.ccm,measure == "t.max.skill")$y_x.obs_value)%>%
  group_by(troph,FD.metric,state.metric) %>%
  summarise(prop.forward=sum(forward == TRUE)/length(forward),
            prop.reverse=sum(reverse == TRUE)/length(reverse),
            prop.bidirec=sum(bidirec == TRUE)/length(bidirec),
            prop.none=sum(none == TRUE)/length(none),
            mean.lag = mean(diff.lag)) %>%
            mutate(measure = "max.skill")

ccm.lag0.lagx.comp <- rbind(ccm.lag0.comp,ccm.lagx.comp)


ccm.lag.plot.df <- summary.ccm %>%
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
  #mutate(lag = -1*lag) %>% #for plotting purpose convert lags from negative to positive 
  #filter(sig == "*")%>% #only keep significant relationships
  #mutate(FD.metric = ifelse(troph == "Zooplankton", paste("zoo",FD.metric,sep = ""),FD.metric))%>%
  ungroup()

count.lag.ccmdf <- ccm.lag.plot.df %>%
  filter(sig == "*")%>% #only keep significant relationships
  group_by(state.metric,causality.direc,FD.metric,troph)%>%
  summarise(N = n()) #significant count per group

##############################################################################################################
## Cross correlations between FD and state
##############################################################################################################
FDIS.ccf <- ggplot(filter(summary.ccf.mth1,FD.metric %in% "FDis" & measure %in% c("r0.ccf","absmax.ccf")),
                        aes(x=measure,y=obs.value,col=measure))+
  # geom_rect(data = layer_data(pccf.lag0.1, 1L),
  #              aes(x = x, y = y,xmin = xmin, xmax=xmax, ymin = ymin,ymax=ymax),
  #              color = "black", size = 0.5,position = position_dodge(width=0.75))
  #geom_pointrange(aes(y=mid,ymin =lwr,ymax=upr,group=FD.metric),position= position_dodge(width=0.75),shape = 22,col="black",size=0.5)+
  #geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.025, 0.5, 0.975),scale = "width",alpha = 0.3,position = position_dodge(width=0.75))+
  geom_hline(yintercept = 0,col="black",alpha = 0.3)+
  geom_boxplot(aes(fill=measure),alpha=0.5,col="black",outlier.shape = NA,size=0.4)+
  # geom_point(position=position_dodge(width=0.75),
  #            aes(alpha=sig,group=measure,shape = system),size=2,col="black",fill="black",)+
  # geom_point(position=position_dodge(width=0.75),
  #            aes(fill=NULL,group=measure,shape = system),size=2,col="black")+
  geom_point(position=position_dodge(width=0.75),
             aes(alpha=sig,group=measure,shape = system),size=2,col="black",fill="black",)+
  geom_point(position=position_dodge(width=0.75),
             aes(fill=NULL,group=measure,shape = system),size=2,col="black")+
  #facet_grid(troph~state.metric)+
  ggh4x::facet_nested(~state.metric+troph,scales = "free",
                      labeller = label_value,space="free_x" )+
                      
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_colour_manual(values=c("#5F4B8BFF","#E69A8DFF"),name = "Lag",labels = c("LagX","Lag0")) + 
  scale_fill_manual(values=c("#5F4B8BFF","#E69A8DFF"),name = "Lag",labels = c("LagX","Lag0"))+
  scale_alpha_manual(values=c(0.01,1.0),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  # geom_text(data =count.ccf.r0.dat,aes(x = state.metric, y = ref.y,label = paste("(",NSig,"*",",",NxSig,")",sep = ""),fill=FD.metric,group = FD.metric),
  #           col= "black",size = 3, position = position_dodge(width = 0.9))+
  scale_x_discrete(labels=c("r0.ccf" = "Lag0", "absmax.ccf" = "LagX"),limits=c("r0.ccf","absmax.ccf"))+
  scale_y_continuous(limits = c(-0.9,0.9),breaks = seq(-0.75,0.75,0.5))+
  ylab("Cross correlation") + xlab("Unlagged (Lag0) vs lagged (LagX) comparison")+
  guides(size = guide_legend(order = 1), 
         col = guide_legend(order = 2),
         fill = guide_legend(order = 2),
         shape = guide_legend(order = 3))+
  ggtitle('FDis') +
  theme_bw() +  theme(plot.margin = unit(c(0,0,0.1,0), "pt"),
                      panel.spacing.x = unit(0.25, "lines"))


FEVE.ccf <- ggplot(filter(summary.ccf.mth1,FD.metric %in% "FEve" & measure %in% c("r0.ccf","absmax.ccf")),
                   aes(x=measure,y=obs.value,col=measure))+
  # geom_rect(data = layer_data(pccf.lag0.1, 1L),
  #              aes(x = x, y = y,xmin = xmin, xmax=xmax, ymin = ymin,ymax=ymax),
  #              color = "black", size = 0.5,position = position_dodge(width=0.75))
  #geom_pointrange(aes(y=mid,ymin =lwr,ymax=upr,group=FD.metric),position= position_dodge(width=0.75),shape = 22,col="black",size=0.5)+
  #geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.025, 0.5, 0.975),scale = "width",alpha = 0.3,position = position_dodge(width=0.75))+
  geom_hline(yintercept = 0,col="black",alpha = 0.3)+
  geom_boxplot(aes(fill=measure),alpha=0.5,col="black",outlier.shape = NA,size=0.4)+
  geom_point(position=position_dodge(width=0.75),
             aes(alpha=sig,group=measure,shape = system),size=2,col="black",fill="black",)+
  geom_point(position=position_dodge(width=0.75),
             aes(fill=NULL,group=measure,shape = system),size=2,col="black")+
  #facet_grid(troph~state.metric)+
  ggh4x::facet_nested(~state.metric+troph,scales = "free",
                      labeller = label_value,space="free_x")+
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_colour_manual(values=c("#5F4B8BFF","#E69A8DFF"),name = "Lag",labels = c("LagX","Lag0")) + 
  scale_fill_manual(values=c("#5F4B8BFF","#E69A8DFF"),name = "Lag",labels = c("LagX","Lag0"))+
  scale_alpha_manual(values=c(0.01,1.0),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  # geom_text(data =count.ccf.r0.dat,aes(x = state.metric, y = ref.y,label = paste("(",NSig,"*",",",NxSig,")",sep = ""),fill=FD.metric,group = FD.metric),
  #           col= "black",size = 3, position = position_dodge(width = 0.9))+
  scale_x_discrete(labels=c("r0.ccf" = "Lag0", "absmax.ccf" = "LagX"),limits=c("r0.ccf","absmax.ccf"))+
  scale_y_continuous(limits = c(-0.9,0.9),breaks = seq(-0.75,0.75,0.5))+
  ylab("Cross correlation") + xlab("Unlagged (Lag0) vs lagged (LagX) comparison")+
  guides(size = guide_legend(order = 1), 
         col = guide_legend(order = 2),
         fill = guide_legend(order = 2),
         shape = guide_legend(order = 3))+
  ggtitle('FEve') +
  theme_bw() +  theme(plot.margin = unit(c(0,0,0,0), "pt"),
                      panel.spacing.x = unit(0.25, "lines"))


FRIC.ccf <- ggplot(filter(summary.ccf.mth1,FD.metric %in% "FRic" & measure %in% c("r0.ccf","absmax.ccf")),
                   aes(x=measure,y=obs.value,col=measure))+
  # geom_rect(data = layer_data(pccf.lag0.1, 1L),
  #              aes(x = x, y = y,xmin = xmin, xmax=xmax, ymin = ymin,ymax=ymax),
  #              color = "black", size = 0.5,position = position_dodge(width=0.75))
  #geom_pointrange(aes(y=mid,ymin =lwr,ymax=upr,group=FD.metric),position= position_dodge(width=0.75),shape = 22,col="black",size=0.5)+
  #geom_violin(aes(fill = FD.metric),draw_quantiles =  c(0.025, 0.5, 0.975),scale = "width",alpha = 0.3,position = position_dodge(width=0.75))+
  geom_hline(yintercept = 0,col="black",alpha = 0.3)+
  geom_boxplot(aes(fill=measure),alpha=0.5,col="black",outlier.shape = NA,size=0.4)+
  geom_point(position=position_dodge(width=0.75),
             aes(alpha=sig,group=measure,shape = system),size=2,col="black",fill="black",)+
  geom_point(position=position_dodge(width=0.75),
             aes(fill=NULL,group=measure,shape = system),size=2,col="black")+
  #facet_grid(troph~state.metric)+
  # ggh4x::facet_nested(~state.metric+troph,scales = "free",
  #                     labeller = label_value,space="free_x",switch="x",
  #                     strip = ggh4x::strip_nested(size="constant",bleed=T, 
  #             background_x =elem_list_rect(fill = c(rep("white",5),rep("#999999",10))),
  #                                                 by_layer_x = F))+
  # 
  ggh4x::facet_nested(~state.metric+troph,scales = "free",
                      labeller = label_value,space="free_x")+
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_colour_manual(values=c("#5F4B8BFF","#E69A8DFF"),name = "Lag",labels = c("LagX","Lag0")) + 
  scale_fill_manual(values=c("#5F4B8BFF","#E69A8DFF"),name = "Lag",labels = c("LagX","Lag0"))+
  scale_alpha_manual(values=c(0.01,1.0),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  # geom_text(data =count.ccf.r0.dat,aes(x = state.metric, y = ref.y,label = paste("(",NSig,"*",",",NxSig,")",sep = ""),fill=FD.metric,group = FD.metric),
  #           col= "black",size = 3, position = position_dodge(width = 0.9))+
  scale_x_discrete(labels=c("r0.ccf" = "Lag0", "absmax.ccf" = "LagX"),limits=c("r0.ccf","absmax.ccf"))+
  scale_y_continuous(limits = c(-0.9,0.9),breaks = seq(-0.75,0.75,0.5))+
  ylab("Cross correlation") + xlab("Unlagged (Lag0) vs lagged (LagX) comparison")+
  guides(size = guide_legend(order = 1), 
         col = guide_legend(order = 2),
         fill = guide_legend(order = 2),
         shape = guide_legend(order = 3))+
  ggtitle('FRic') +
  theme_bw() +  theme(plot.margin = unit(c(0,0,0,0), "pt"),
                      panel.spacing.x = unit(0.25, "lines"))

ccf.lag.scatter.dat <- filter(summary.ccf.mth1,measure %in% "t.absmax.ccf") %>% 
  mutate(sig = as.factor(filter(summary.ccf.mth1,measure %in% "absmax.ccf")$sig),
         FD.metric = as.factor(FD.metric))

ccf.lagpanel <- ggplot(ccf.lag.scatter.dat,
      aes(x=FD.metric,y=obs.value)) + 
  geom_hline(yintercept = 0,col="black",alpha = 0.3)+
  #scale_y_binned(breaks = c(seq(60,24,-12),12,seq(-12,-60,-12)),show.limits = T)+
    geom_point(position=position_dodge(width=0.75),
           aes(alpha=sig,group=measure,shape = system),size=2,col="black",fill="black",)+
  geom_point(position=position_dodge(width=0.75),
             aes(fill=NULL,group=measure,shape = system),size=2,col="black")+
  #scale_y_binned(breaks = c(seq(60,24,-12),12,seq(-12,-60,-12)),show.limits = T)+
  geom_rect(aes(ymax = 12,ymin=-12), xmin = 0,xmax = as.numeric(ccf.lag.scatter.dat$FD.metric[[3]])+3,fill="#5F4B8BFF",alpha=0.025)+
  ggh4x::facet_nested(~state.metric+troph,scales = "free",
                      labeller = label_value,space="free_x")+
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_alpha_manual(values=c(0.01,1.0),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  # geom_text(data =count.ccf.r0.dat,aes(x = state.metric, y = ref.y,label = paste("(",NSig,"*",",",NxSig,")",sep = ""),fill=FD.metric,group = FD.metric),
  #           col= "black",size = 3, position = position_dodge(width = 0.9))+
  ylab("Strongest lag") + xlab("Functional diversity metric")+
  guides(size = guide_legend(order = 1), 
         col = guide_legend(order = 2),
         fill = guide_legend(order = 2),
         shape = guide_legend(order = 3))+
  theme_bw() +  theme(plot.margin = unit(c(0,0,0,0), "pt"),
                      panel.spacing.x = unit(0.25, "lines"))

pdf(file="Results/ccf/summary_FD_perm_lm_mth_combo_alt.pdf",
    width=14, height = 12)
FDIS.ccf + FEVE.ccf + FRIC.ccf + ccf.lagpanel +  plot_layout(nrow = 4,guides = "collect",heights = c(1,1,1,0.9))+
                    plot_annotation(tag_levels = c('a')) & 
                    theme(plot.tag = element_text(face = "bold"))&
  theme(legend.position='right')
dev.off()  

##############################################################################################################
## Cross mappings between FD and state
##############################################################################################################
FDIS.ccm <- ggplot(filter(summary.ccm,FD.metric %in% "FDis" & measure %in% c("r0.skill","max.skill")),
                   aes(x=measure,y=y_x.obs_value,col=measure))+
  geom_hline(yintercept = 0,col="black",alpha = 0.3)+
  geom_boxplot(aes(fill=measure),alpha=0.5,col="black",outlier.shape = NA,size=0.4)+
  geom_point(position=position_dodge(width=0.75),
             aes(alpha=y_x.sig,group=measure,shape = system),size=2,col="black",fill="black",)+
  geom_point(position=position_dodge(width=0.75),
             aes(fill=NULL,group=measure,shape = system),size=2,col="black")+
  ggh4x::facet_nested(~state.metric+troph,scales = "free",
                      labeller = label_value,space="free_x" )+
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_colour_manual(values=c("#5F4B8BFF","#E69A8DFF"),name = "Lag",labels = c("LagX","Lag0")) + 
  scale_fill_manual(values=c("#5F4B8BFF","#E69A8DFF"),name = "Lag",labels = c("LagX","Lag0"))+
  scale_alpha_manual(values=c(0.01,1.0),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  # ggtext::geom_richtext(data =filter(ccm.lag0.lagx.comp,FD.metric %in% "FDis"),
  #                       aes(x = measure, y =0.8,
  #                    label = paste("<span style='color:black'>","(","</span>","<span style='color:#91845C'>",base::format(prop.forward,digits = 2),"</span>","<span style='color:black'>",",",base::format(prop.bidirec,digits =2),")","</span>",sep = "")),
  #                       alpha=0,size = 3, position = position_dodge(width = 0.75),angle = 90)+
  scale_x_discrete(labels=c("r0.skill" = "Lag0", "max.skill" = "LagX"),limits=c("r0.skill","max.skill"))+
  scale_y_continuous(breaks = seq(0,1.0,0.25),limits = c(-0.17,0.8))+
  ylab("Cross map skill") + xlab("Unlagged (Lag0) vs lagged (LagX) comparison")+
  guides(size = guide_legend(order = 1), 
         col = guide_legend(order = 2),
         fill = guide_legend(order = 2),
         shape = guide_legend(order = 3))+
  ggtitle('FDis') +
  theme_bw() +  theme(plot.margin = unit(c(0,0,0.1,0), "pt"),
                      panel.spacing.x = unit(0.25, "lines"))


FEVE.ccm <- ggplot(filter(summary.ccm,FD.metric %in% "FEve" & measure %in% c("r0.skill","max.skill")),
                   aes(x=measure,y=y_x.obs_value,col=measure))+
  geom_hline(yintercept = 0,col="black",alpha = 0.3)+
  geom_boxplot(aes(fill=measure),alpha=0.5,col="black",outlier.shape = NA,size=0.4)+
  geom_point(position=position_dodge(width=0.75),
             aes(alpha=y_x.sig,group=measure,shape = system),size=2,col="black",fill="black",)+
  geom_point(position=position_dodge(width=0.75),
             aes(fill=NULL,group=measure,shape = system),size=2,col="black")+
  ggh4x::facet_nested(~state.metric+troph,scales = "free",
                      labeller = label_value,space="free_x" )+
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_colour_manual(values=c("#5F4B8BFF","#E69A8DFF"),name = "Lag",labels = c("LagX","Lag0")) + 
  scale_fill_manual(values=c("#5F4B8BFF","#E69A8DFF"),name = "Lag",labels = c("LagX","Lag0"))+
  scale_alpha_manual(values=c(0.01,1.0),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  # ggtext::geom_richtext(data =filter(ccm.lag0.lagx.comp,FD.metric %in% "FEve"),
  #                       aes(x = measure, y =0.8,
  #                             label = paste("<span style='color:black'>","(","</span>","<span style='color:#91845C'>",base::format(prop.forward,digits = 2),"</span>","<span style='color:black'>",",",base::format(prop.bidirec,digits =2),")","</span>",sep = "")),
  #                       alpha=0,size = 3, position = position_dodge(width = 0.75),angle = 90)+
  scale_x_discrete(labels=c("r0.skill" = "Lag0", "max.skill" = "LagX"),limits=c("r0.skill","max.skill"))+
  scale_y_continuous(breaks = seq(0,1.0,0.25),limits = c(-0.17,0.8))+
  ylab("Cross map skill") + xlab("Unlagged (Lag0) vs lagged (LagX) comparison")+
  guides(size = guide_legend(order = 1), 
         col = guide_legend(order = 2),
         fill = guide_legend(order = 2),
         shape = guide_legend(order = 3))+
  ggtitle('FEve') +
  theme_bw() +  theme(plot.margin = unit(c(0,0,0.1,0), "pt"),
                      panel.spacing.x = unit(0.25, "lines"))

FRIC.ccm <- ggplot(filter(summary.ccm,FD.metric %in% "FRic" & measure %in% c("r0.skill","max.skill")),
                   aes(x=measure,y=y_x.obs_value,col=measure))+
  geom_hline(yintercept = 0,col="black",alpha = 0.3)+
  geom_boxplot(aes(fill=measure),alpha=0.5,col="black",outlier.shape = NA,size=0.4)+
  geom_point(position=position_dodge(width=0.75),
             aes(alpha=y_x.sig,group=measure,shape = system),size=2,col="black",fill="black",)+
  geom_point(position=position_dodge(width=0.75),
             aes(fill=NULL,group=measure,shape = system),size=2,col="black")+
  ggh4x::facet_nested(~state.metric+troph,scales = "free",
                      labeller = label_value,space="free_x" )+
  # ggtext::geom_richtext(data =filter(ccm.lag0.lagx.comp,FD.metric %in% "FRic"),
  #                       aes(x = measure, y =0.8,
  #                         label = paste("<span style='color:black'>","(","</span>","<span style='color:#91845C'>",base::format(prop.forward,digits = 2),"</span>","<span style='color:black'>",",",base::format(prop.bidirec,digits =2),")","</span>",sep = "")),
  #                       alpha=0,size = 3, position = position_dodge(width = 0.75),angle = 90)+
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_colour_manual(values=c("#5F4B8BFF","#E69A8DFF"),name = "Lag",labels = c("LagX","Lag0")) + 
  scale_fill_manual(values=c("#5F4B8BFF","#E69A8DFF"),name = "Lag",labels = c("LagX","Lag0"))+
  scale_alpha_manual(values=c(0.01,1.0),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  scale_x_discrete(labels=c("r0.skill" = "Lag0", "max.skill" = "LagX"),limits=c("r0.skill","max.skill"))+
  scale_y_continuous(breaks = seq(0,1.0,0.25),limits = c(-0.17,0.8))+
  ylab("Cross map skill") + xlab("Unlagged (Lag0) vs lagged (LagX) comparison")+
  guides(size = guide_legend(order = 1), 
         col = guide_legend(order = 2),
         fill = guide_legend(order = 2),
         shape = guide_legend(order = 3))+
  ggtitle('FRic') +
  theme_bw() +  theme(plot.margin = unit(c(0,0,0.1,0), "pt"),
                      panel.spacing.x = unit(0.25, "lines"))


ccm.lag.scatter.dat <- filter(summary.ccm,measure %in% "t.max.skill") %>% 
  mutate(y_x.sig = as.factor(filter(summary.ccm,measure %in% "max.skill")$y_x.sig),
         FD.metric = as.factor(FD.metric))

ccm.lagpanel <- ggplot(ccm.lag.scatter.dat,
                   aes(x=FD.metric,y=y_x.obs_value)) + 
  geom_hline(yintercept = 0,col="black",alpha = 0.3)+
  geom_point(position=position_dodge(width=0.75),
             aes(alpha=y_x.sig,group=measure,shape = system),size=2,col="black",fill="black",)+
  geom_point(position=position_dodge(width=0.75),
             aes(fill=NULL,group=measure,shape = system),size=2,col="black")+
  geom_rect(aes(ymax = 12,ymin=-12), xmin = 0,xmax = as.numeric(ccm.lag.scatter.dat$FD.metric[[3]])+3,fill="#5F4B8BFF",alpha=0.025)+
  ggh4x::facet_nested(~state.metric+troph,scales = "free",
                      labeller = label_value,space="free_x")+
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_alpha_manual(values=c(0.01,1.0),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  ylab("Strongest lag") + xlab("Functional diversity metric")+
  guides(size = guide_legend(order = 1), 
         col = guide_legend(order = 2),
         fill = guide_legend(order = 2),
         shape = guide_legend(order = 3))+
  theme_bw() +  theme(plot.margin = unit(c(0,0,0.1,0), "pt"),
                      panel.spacing.x = unit(0.25, "lines"))

pdf(file="Results/ccm/summary_ccm_combo_alt.pdf",
    width=14, height = 12)
FDIS.ccm + FEVE.ccm + FRIC.ccm + ccm.lagpanel +  plot_layout(nrow = 4,guides = "collect",heights = c(1,1,1,0.9))+
  plot_annotation(tag_levels = c('a')) & 
  theme(plot.tag = element_text(face = "bold"))&
  theme(legend.position='right')
dev.off()  


##############################################################################################################
## Cross mapping lag changes
##############################################################################################################
FDIS.lagp <-  ggpubr::ggpaired(ccm.lag.plot.df %>% filter(FD.metric %in% "FDis") %>%
                             mutate(pairing = interaction(state.metric,system,FD.metric,troph)) ,
                           x = "causality.direc", y = "lag", fill = rep(c("#A1B4FE","#FFE7A1"),10), 
                           id = "pairing",alpha = 0.1,point.size =0,  line.size = 0.5,
                           ggtheme = theme_bw(),linetype = "dashed", 
                           facet.by = c("troph","state.metric")) +
  scale_x_discrete(labels =c("reverse","forward"))+
  geom_hline(yintercept = 0,alpha = 0.3)+
  geom_point(aes(group=causality.direc,shape = system, alpha = as.factor(sig)), 
             position = position_dodge(width=0.75),fill = "black",size = 2)+
  geom_point(aes(group=causality.direc,shape = system,fill=NULL),col="black",
             position = position_dodge(width=0.75),size = 2) +
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_alpha_manual(values=c(0.01,1),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  geom_text(data = count.lag.ccmdf %>% filter(FD.metric %in% "FDis"),
            aes(y = (max(ccm.lag.plot.df$lag)+10),x=causality.direc, label = N,group =causality.direc),
            col="black", position = position_dodge(width = 0.8))+
  xlab("Causality direction") + ylab("Strongest lag (months)") + 
  ggtitle('FDis') +
  scale_y_continuous(breaks = seq(-40,40,by=40),limits = c(-60,75))+
  guides(shape = guide_legend(order = 1))+
  theme(plot.margin = unit(c(1,1,0,1), "pt"),
        panel.spacing.x = unit(0.25, "lines"))


FEVE.lagp <-  ggpubr::ggpaired(ccm.lag.plot.df %>% filter(FD.metric %in% "FEve") %>%
                              mutate(pairing = interaction(state.metric,system,FD.metric,troph)) ,
                            x = "causality.direc", y = "lag", fill = rep(c("#A1B4FE","#FFE7A1"),10), 
                            id = "pairing",alpha = 0.1,point.size =0,  line.size = 0.5,
                            ggtheme = theme_bw(),linetype = "dashed", 
                            facet.by = c("troph","state.metric")) +
  scale_x_discrete(labels =c("reverse","forward"))+
  geom_hline(yintercept = 0,alpha = 0.3)+
  geom_point(aes(group=causality.direc,shape = system, alpha = as.factor(sig)), 
             position = position_dodge(width=0.75),fill = "black",size = 2)+
  geom_point(aes(group=causality.direc,shape = system,fill=NULL),col="black",
             position = position_dodge(width=0.75),size = 2) +
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_alpha_manual(values=c(0.01,1),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  geom_text(data = count.lag.ccmdf %>% filter(FD.metric %in% "FEve"),
            aes(y = (max(ccm.lag.plot.df$lag)+10),x=causality.direc, label = N,group =causality.direc),
            col="black", position = position_dodge(width = 0.8))+
  xlab("Causality direction") + ylab("Strongest lag (months)") + 
  ggtitle('FEve') +
  scale_y_continuous(breaks = seq(-40,40,by=40),limits = c(-60,75))+
  guides(shape = guide_legend(order = 1))+
  theme(plot.margin = unit(c(1,1,0,1), "pt"),
        panel.spacing.x = unit(0.25, "lines"))


FRIC.lagp <- ggpubr::ggpaired(ccm.lag.plot.df %>% filter(FD.metric %in% "FRic") %>%
                             mutate(pairing = interaction(state.metric,system,FD.metric,troph)) ,
                           x = "causality.direc", y = "lag", fill = rep(c("#A1B4FE","#FFE7A1"),10), 
                           id = "pairing",alpha = 0.1,point.size =0,  line.size = 0.5,
                           ggtheme = theme_bw(),linetype = "dashed", 
                           facet.by = c("troph","state.metric")) +
  scale_x_discrete(labels =c("reverse","forward"))+
  geom_hline(yintercept = 0,alpha = 0.3)+
  geom_point(aes(group=causality.direc,shape = system, alpha = as.factor(sig)), 
             position = position_dodge(width=0.75),fill = "black",size = 2)+
  geom_point(aes(group=causality.direc,shape = system,fill=NULL),col="black",
             position = position_dodge(width=0.75),size = 2) +
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake",
                     guide_legend(override.aes = list(fill = "white",col="white")))+
  scale_alpha_manual(values=c(0.01,1),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  geom_text(data = count.lag.ccmdf %>% filter(FD.metric %in% "FRic"),
            aes(y = (max(ccm.lag.plot.df$lag)+10),x=causality.direc, label = N,group =causality.direc),
            col="black", position = position_dodge(width = 0.8))+
  xlab("Causality direction") + ylab("Strongest lag (months)") + 
  ggtitle('FRic') +
  scale_y_continuous(breaks = seq(-40,40,by=40),limits = c(-60,75))+
  guides(shape = guide_legend(order = 1))+
  theme(plot.margin = unit(c(1,1,0,1), "pt"),
        panel.spacing.x = unit(0.25, "lines"))

pdf(file="Results/ccm/ccm_causality_spread.pdf",
    width=10, height = 10)
FDIS.lagp + FEVE.lagp + FRIC.lagp + plot_layout(nrow = 3,guides = "collect")+
  plot_annotation(tag_levels = c('a')) & 
  theme(plot.tag = element_text(face = "bold"))
dev.off()

##############################################################################################################
## Cross mapping lag changes ALT
##############################################################################################################
FDIS.lagp <-  ggpubr::ggpaired(ccm.lag.plot.df %>% filter(FD.metric %in% "FDis") %>%
                                 mutate(pairing = interaction(state.metric,system,FD.metric,troph)) ,
                               x = "causality.direc", y = "lag", fill = rep(c("#A1B4FE","#FFE7A1"),10), 
                               id = "pairing",alpha = 0.1,point.size =0,  line.size = 0.5,
                               ggtheme = theme_bw(),linetype = "dashed") +
  scale_x_discrete(labels =c("reverse","forward"))+
  geom_hline(yintercept = 0,alpha = 0.3)+
  geom_point(aes(group=causality.direc,shape = system, alpha = as.factor(sig)), 
             position = position_dodge(width=0.75),fill = "black",size = 2)+
  geom_point(aes(group=causality.direc,shape = system,fill=NULL),col="black",
             position = position_dodge(width=0.75),size = 2) +
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_alpha_manual(values=c(0.01,1),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  geom_text(data = count.lag.ccmdf %>% filter(FD.metric %in% "FDis"),
            aes(y = (max(ccm.lag.plot.df$lag)+10),x=causality.direc, label = N,group =causality.direc),
            col="black", position = position_dodge(width = 0.8))+
  xlab("Causality direction") + ylab("Strongest lag (months)") + 
  ggtitle('FDis') +
  ggh4x::facet_nested(~state.metric+troph,scales = "free",
                      labeller = label_value,space="free_x" )+
  scale_y_continuous(breaks = seq(-40,40,by=40),limits = c(-60,75))+
  guides(shape = guide_legend(order = 1))+
  theme(plot.margin = unit(c(1,1,0,1), "pt"),
        panel.spacing.x = unit(0.25, "lines")) 


FEVE.lagp <-  ggpubr::ggpaired(ccm.lag.plot.df %>% filter(FD.metric %in% "FEve") %>%
                                 mutate(pairing = interaction(state.metric,system,FD.metric,troph)) ,
                               x = "causality.direc", y = "lag", fill = rep(c("#A1B4FE","#FFE7A1"),10), 
                               id = "pairing",alpha = 0.1,point.size =0,  line.size = 0.5,
                               ggtheme = theme_bw(),linetype = "dashed") +
  scale_x_discrete(labels =c("reverse","forward"))+
  geom_hline(yintercept = 0,alpha = 0.3)+
  geom_point(aes(group=causality.direc,shape = system, alpha = as.factor(sig)), 
             position = position_dodge(width=0.75),fill = "black",size = 2)+
  geom_point(aes(group=causality.direc,shape = system,fill=NULL),col="black",
             position = position_dodge(width=0.75),size = 2) +
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_alpha_manual(values=c(0.01,1),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  geom_text(data = count.lag.ccmdf %>% filter(FD.metric %in% "FEve"),
            aes(y = (max(ccm.lag.plot.df$lag)+10),x=causality.direc, label = N,group =causality.direc),
            col="black", position = position_dodge(width = 0.8))+
  xlab("Causality direction") + ylab("Strongest lag (months)") + 
  ggtitle('FEve') +
  ggh4x::facet_nested(~state.metric+troph,scales = "free",
                      labeller = label_value,space="free_x" )+
  scale_y_continuous(breaks = seq(-40,40,by=40),limits = c(-60,75))+
  guides(shape = guide_legend(order = 1))+
  theme(plot.margin = unit(c(1,1,0,1), "pt"),
        panel.spacing.x = unit(0.25, "lines"))


FRIC.lagp <- ggpubr::ggpaired(ccm.lag.plot.df %>% filter(FD.metric %in% "FRic") %>%
                                mutate(pairing = interaction(state.metric,system,FD.metric,troph)) ,
                              x = "causality.direc", y = "lag", fill = rep(c("#A1B4FE","#FFE7A1"),10), 
                              id = "pairing",alpha = 0.1,point.size =0,  line.size = 0.5,
                              ggtheme = theme_bw(),linetype = "dashed") +
  scale_x_discrete(labels =c("reverse","forward"))+
  geom_hline(yintercept = 0,alpha = 0.3)+
  geom_point(aes(group=causality.direc,shape = system, alpha = as.factor(sig)), 
             position = position_dodge(width=0.75),fill = "black",size = 2)+
  geom_point(aes(group=causality.direc,shape = system,fill=NULL),col="black",
             position = position_dodge(width=0.75),size = 2) +
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake",
                     guide_legend(override.aes = list(fill = "white",col="white")))+
  scale_alpha_manual(values=c(0.01,1),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  geom_text(data = count.lag.ccmdf %>% filter(FD.metric %in% "FRic"),
            aes(y = (max(ccm.lag.plot.df$lag)+10),x=causality.direc, label = N,group =causality.direc),
            col="black", position = position_dodge(width = 0.8))+
  xlab("Causality direction") + ylab("Strongest lag (months)") + 
  ggtitle('FRic') +
  ggh4x::facet_nested(~state.metric+troph,scales = "free",
                      labeller = label_value,space="free_x" )+
  scale_y_continuous(breaks = seq(-40,40,by=40),limits = c(-60,75))+
  guides(shape = guide_legend(order = 1))+
  theme(plot.margin = unit(c(1,1,0,1), "pt"),
        panel.spacing.x = unit(0.25, "lines"))

pdf(file="Results/ccm/ccm_causality_spread_alt.pdf",
    width=14, height = 10)
FDIS.lagp + FEVE.lagp + FRIC.lagp + plot_layout(nrow = 3,guides = "collect")+
  plot_annotation(tag_levels = c('a')) & 
  theme(plot.tag = element_text(face = "bold"))
dev.off()
