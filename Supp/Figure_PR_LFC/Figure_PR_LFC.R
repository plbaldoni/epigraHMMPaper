load('Figure_PR_H3K36me3_Encode_twocells_500bp_10LFC.RData')
load('Figure_PR_H3K36me3_Encode_twocells_500bp_15LFC.RData')
load('Figure_PR_H3K36me3_Encode_twocells_500bp_25LFC.RData')
load('Figure_PR_H3K36me3_Encode_twocells_500bp_30LFC.RData')

fig <- ggarrange(fig_lfc10,
                 fig_lfc15,
                 fig_lfc25,
                 fig_lfc30,
                 nrow = 2,ncol = 2,legend = 'top',common.legend = TRUE)

ggsave(fig,filename = 'Figure_PR_LFC.pdf',width = 9,height = 11.5)
