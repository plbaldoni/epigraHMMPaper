# Figure A: genomic segmentation with Viterbi path

# Loading libraries
library(epigraHMM)
library(data.table)
library(magrittr)
library(plyr)
library(ggplot2)
library(ggrepel)

# Some general parameters

B <- 6

# Loading output from epigraHMM (full model)
load(
  "../../Public/epigraHMM/Helas3/Encode_H3K27me3_H3K36me3_EZH2/Output/epigraHMM_Helas3_Encode_H3K27me3_H3K36me3_EZH2_Output_500bp_6Patterns.RData"
)
epigraHMM_object <-
  epigraHMM_Helas3_Encode_H3K27me3_H3K36me3_EZH2_Output_500bp_6Patterns
rm(epigraHMM_Helas3_Encode_H3K27me3_H3K36me3_EZH2_Output_500bp_6Patterns)

# Getting parameter estimates
parhist <- metadata(epigraHMM_object)$history$parameter
parhist[, it := 1:(.N)]
parhist_delta <-
  melt(
    data = parhist,
    id.vars = 'it',
    measure.vars = paste0('delta', 1:B),
    variable.name = 'Parameter',
    value.name = 'Estimate'
  )[it == nrow(parhist),]

# Labels

label_delta <- c(
  expression(delta[1]),
  expression(delta[2]),
  expression(delta[3]),
  expression(delta[4]),
  expression(delta[5]),
  expression(delta[6])
)

names(label_delta) <- paste0('delta', 1:B)

# Mapping combinatorial patterns

labels <-
  unlist(lapply(metadata(epigraHMM_object)$control$pattern, function(x) {
    paste(unique(colData(epigraHMM_object)$condition)[x], collapse = ' &\n')
  }))

# Plotting

fig_parameters <-
  ggplot(data = parhist_delta, aes(x = Parameter, y = Estimate + 0.005)) + # Adding small number so that the rare pattern has its own bar on the plot
  theme_bw() +
  ylab('Estimate') + xlab('Mixture Model Proportions') +
  geom_bar(stat = "identity") +
  geom_text(position = 'identity', aes(label = ifelse(
    Estimate >= 0.001,
    formatC(round(Estimate, 3), format = 'f', digits = 3),
    '<0.001'
  )), vjust = -0.2) +
  annotate(
    size = 2.5,
    "text",
    x = c(1, c(1:6)),
    y = c(0.7, rep(0.625, 6)),
    hjust = c(0.25, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
    label = c('Associated Enrichment:', labels)
  ) +
  scale_x_discrete(labels = label_delta) +
  theme(panel.grid = element_blank())

save(fig_parameters, file = './Figure_Parameters.RData')