load('Figure_PeakCalls_H3K36me3_Encode_twocells_250bp.RData')
load('Figure_PeakCalls_H3K36me3_Encode_twocells_500bp.RData')
load('Figure_PeakCalls_H3K36me3_Encode_twocells_750bp.RData')
load('Figure_PeakCalls_H3K36me3_Encode_twocells_1000bp.RData')

fig.calls = ggarrange(fig.250,fig.500,fig.750,fig.1000,
                      nrow = 2,ncol = 2,legend = 'top',common.legend = TRUE,
                      labels = list('A','B','C','D'))

ggsave(fig.calls,filename = 'Figure_H3K36me3_Viterbi.pdf',height = 11.5,width = 9)
