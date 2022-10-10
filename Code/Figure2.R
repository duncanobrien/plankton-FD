### Figure 2 ###

## Preamble 
require(tidyverse) # dplyr, ggplot etc.
require(ggh4x) # facet_nested function
require(patchwork) # multi-panel plots

summary.ccf <- read.csv("Results/ccf/raw_data/summary.ccf.csv")

###########################################################################
## Plot unlagged and lagged cross correlation coefficients ##
###########################################################################
ccf.boxplots <- ggplot(filter(summary.ccf, measure %in% c("r0.ccf","absmax.ccf")),
               aes(x=measure,y=obs.value,col=measure))+
  geom_hline(yintercept = 0,col="black",alpha = 0.3)+
  geom_boxplot(aes(fill=measure),alpha=0.5,col="black",outlier.shape = NA,size=0.4)+
  geom_point(position=position_dodge(width=0.75),
             aes(alpha=sig,group=measure,shape = system),size=2,col="black",fill="black",)+
  geom_point(position=position_dodge(width=0.75),
             aes(fill=NULL,group=measure,shape = system),size=2,col="black")+
  ggh4x::facet_nested(FD.metric~state.metric+troph,scales = "free",
                      labeller = label_value,space="free_x" )+
  scale_shape_manual(values = c(21,22,24,25,23),name = "Lake")+
  scale_colour_manual(values=c("#5F4B8BFF","#E69A8DFF"),name = "Lag",labels = c("LagX","Lag0")) + 
  scale_fill_manual(values=c("#5F4B8BFF","#E69A8DFF"),name = "Lag",labels = c("LagX","Lag0"))+
  scale_alpha_manual(values=c(0.01,1.0),name = "Significance", labels = c("Not significant","Significant"),
                     guide = guide_legend(override.aes = list(fill = c("white","black"),alpha = c(1,1),linetype = c("solid","solid"),shape=c(22,22))))+
  scale_x_discrete(labels=c("r0.ccf" = "Lag0", "absmax.ccf" = "LagX"),limits=c("r0.ccf","absmax.ccf"))+
  scale_y_continuous(limits = c(-0.9,0.9),breaks = seq(-0.75,0.75,0.5))+
  ylab("Cross correlation") + xlab("Unlagged (Lag0) vs lagged (LagX) comparison")+
  guides(size = guide_legend(order = 1), 
         col = guide_legend(order = 2),
         fill = guide_legend(order = 2),
         shape = guide_legend(order = 3))+
  theme_bw() +  theme(plot.margin = unit(c(0,0,1,0), "pt"),
                      panel.spacing.x = unit(0.25, "lines"))

ccf.lag.scatter.dat <- filter(summary.ccf,measure %in% "t.absmax.ccf") %>% 
  mutate(sig = as.factor(filter(summary.ccf,measure %in% "absmax.ccf")$sig),
         FD.metric = as.factor(FD.metric))

ccf.lagpanel <- ggplot(ccf.lag.scatter.dat,
                       aes(x=FD.metric,y=obs.value)) + 
  geom_hline(yintercept = 0,col="black",alpha = 0.3)+
  geom_point(position=position_dodge(width=0.75),
             aes(alpha=sig,group=measure,shape = system),size=2,col="black",fill="black",)+
  geom_point(position=position_dodge(width=0.75),
             aes(fill=NULL,group=measure,shape = system),size=2,col="black")+
  geom_rect(aes(ymax = 12,ymin=-12), xmin = 0,xmax = as.numeric(ccf.lag.scatter.dat$FD.metric[[3]])+3,fill="#5F4B8BFF",alpha=0.025)+
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

pdf(file="Results/ccf/Figure2.pdf",
    width=13, height = 10)
ccf.boxplots + ccf.lagpanel +  plot_layout(nrow = 2,guides = "collect",heights = c(2,0.8))+
  plot_annotation(tag_levels = c('a'),tag_prefix = "(",tag_suffix = ")") & 
  theme(plot.tag = element_text(face = "bold"))&
  theme(legend.position='right',
        strip.background = element_rect(fill="white"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))
dev.off()  