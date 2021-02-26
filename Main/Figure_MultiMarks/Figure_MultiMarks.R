library(ggpubr)

load('Figure_Viterbi.RData')
load('Figure_Parameters.RData')
load('Figure_BIC.RData')
load('Figure_GeneExpression.RData')
load('Figure_GeneSymbols.RData')
load('Figure_Segmentation.RData')


# Figure right
fig.top <- ggarrange(fig_viterbi,fig_parameters,nrow = 1,labels = list('A','B'))
fig.mid <- ggarrange(fig_bic,fig_expression,nrow = 1,legend = 'none',labels = list('C','D'))
fig.bot <- ggarrange(fig_genesymbols,fig_segmentation[['3']],nrow = 1,legend = 'none',labels = list('E','F'))

# Put them all together

fig.final = ggarrange(fig.top,fig.mid,fig.bot,nrow=3,heights = c(0.275,0.275,0.45))
fig.final.label = ggarrange(fig.final,get_legend(fig_segmentation[['3']]),nrow=2,heights = c(0.95,0.05))

ggsave(filename = 'Figure_MultiMarks.pdf',plot = fig.final.label,height = 11.5,width = 9,dpi = 'retina')
