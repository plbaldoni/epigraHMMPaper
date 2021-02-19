library(ggpubr)

load('Figure_Correlation_CTCF_Encode_twocells_1000bp.RData')
load('Figure_Correlation_H3K4me3_Encode_twocells_1000bp.RData')
load('Figure_PeakCalls_CTCF_Encode_twocells_1000bp.RData')
load('Figure_PeakCalls_H3K4me3_Encode_twocells_1000bp.RData')
load('Figure_PeakCalls_H3K27ac_Encode_twocells_1000bp.RData')
load('Figure_Memory_Time.RData')

# Fig Peak Calls
fig.calls = ggarrange(fig.ctcf,fig.h3k4me3,fig.h3k27ac,nrow = 3,ncol = 1,legend = 'none',labels = list('B','D','F'))

# Figure Metrics
fig.metrics = ggarrange(fig.lfc,fig.cor,fig_memory_time,ncol=1,nrow=3,legend = 'none',labels = list('A','C','E'))

# Put them all together

fig.final = ggarrange(fig.metrics,fig.calls,nrow=1,ncol=2)
fig.final.label = ggarrange(get_legend(fig.lfc + theme(legend.position = 'top',legend.direction = 'horizontal',legend.box = 'vertical',legend.spacing.y = unit(-0.25, "cm")) + guides(col = guide_legend(nrow = 1,order = 1),linetype = guide_legend(nrow = 1,order = 2,override.aes = list(col = 'black'),title = 'Observed Enrichment'))),
                            fig.final,nrow=2,heights = c(0.05,0.95))

ggsave(filename = 'Figure_Narrow.pdf',plot = fig.final.label,height = 11.5,width = 9,dpi = 'retina')
