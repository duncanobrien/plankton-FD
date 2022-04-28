### Figure 3 ###

## Preamble 
require(tidyverse) # dplyr, ggplot etc.
require(ggh4x) # facet_nested function
require(patchwork) # multi-panel plots

summary.ccm <- read.csv(file ="Results/ccm/raw_data/ccm_summary.csv")

###########################################################################
## Plot unlagged and lagged cross map skills ##
###########################################################################

ccm.boxplots <- ggplot(filter(summary.ccm, measure %in% c("r0.skill","max.skill")),
                aes(x=measure,y=y_x.obs_value,col=measure))+
  geom_hline(yintercept = 0,col="black",alpha = 0.3)+
  geom_boxplot(aes(fill=measure),alpha=0.5,col="black",outlier.shape = NA,size=0.4)+
  geom_point(position=position_dodge(width=0.75),
             aes(alpha=y_x.sig,group=measure,shape = system),size=2,col="black",fill="black",)+
  geom_point(position=position_dodge(width=0.75),
             aes(fill=NULL,group=measure,shape = system),size=2,col="black")+
  ggh4x::facet_nested(FD.metric~state.metric+troph,scales = "free",
                      labeller = label_value,space="free_x" )+
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

pdf(file="Results/ccm/Figure3.pdf",
    width=13, height = 10)
ccm.boxplots + ccf.lagpanel +  plot_layout(nrow = 2,guides = "collect",heights = c(2,0.8))+
  plot_annotation(tag_levels = c('a')) & 
  theme(plot.tag = element_text(face = "bold"))&
  theme(legend.position='right',
        strip.background = element_rect(fill="white"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))
dev.off()  