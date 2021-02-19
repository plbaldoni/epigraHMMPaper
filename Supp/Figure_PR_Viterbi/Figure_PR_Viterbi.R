load('Figure_PR_H3K36me3_Encode_twocells_250bp.RData')
load('Figure_PR_H3K36me3_Encode_twocells_500bp.RData')
load('Figure_PR_H3K36me3_Encode_twocells_750bp.RData')
load('Figure_PR_H3K36me3_Encode_twocells_1000bp.RData')

fig <- ggarrange(fig.ROC_250,fig.ROC_500,fig.ROC_750,fig.ROC_1000,
                 nrow = 2,ncol = 2,legend = 'top',common.legend = TRUE)

ggsave(fig,filename = 'Figure_PR_Viterbi.pdf',width = 9,height = 11.5)
