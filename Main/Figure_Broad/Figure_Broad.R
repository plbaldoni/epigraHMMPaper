library(ggpubr)

load('Figure_PR_H3K36me3_Encode_twocells_500bp.RData')
load('Figure_PeakCalls_H3K36me3_Encode_twocells_500bp.RData')
load('Figure_PeakCalls_H3K27me3_Encode_twocells_500bp.RData')
load('Figure_PeakCalls_EZH2_Encode_twocells_500bp.RData')
load('Figure_Memory_Time.RData')


# Fig Peak Calls
fig.calls = ggarrange(fig.h3k36me3,fig.h3k27me3,fig.ezh2,nrow = 3,ncol = 1,legend = 'none',labels = list('B','D','F'))

# Figure Metrics
fig.metrics = ggarrange(fig.ROC,fig.TPR,fig_memory_time,ncol=1,nrow=3,legend = 'none',labels = list('A','C','E'))

# Put them all together

fig.final = ggarrange(fig.metrics,fig.calls,nrow=1,ncol=2)
fig.final.label = ggarrange(get_legend(fig.ROC + theme(legend.position = 'top',legend.direction = 'horizontal',legend.box = 'vertical',legend.spacing.y = unit(-0.25, "cm")) + guides(col = guide_legend(nrow = 1,order = 1),shape = guide_legend(nrow = 1,order = 2,title = 'Nominal FDR'))),
                            fig.final,nrow=2,heights = c(0.05,0.95))

ggsave(filename = 'Figure_Broad.pdf',plot = fig.final.label,height = 11.5,width = 9,dpi = 'retina')
