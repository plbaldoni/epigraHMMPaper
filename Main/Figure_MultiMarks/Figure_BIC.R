# Figure C: genomic segmentation with Viterbi path

# Loading libraries
library(epigraHMM)
library(data.table)
library(magrittr)
library(plyr)
library(ggplot2)
library(ggrepel)
library(ggupset)
library(S4Vectors)
library(SummarizedExperiment)

# Some general parameters

B <- 1:6

# Loading output from epigraHMM (full model)
dt_list <- list()
for (b in B) {
  name_object <- paste0('epigraHMM_Helas3_Encode_H3K27me3_H3K36me3_EZH2_Output_500bp_',b,'Patterns')
  load(file.path("../../Public/epigraHMM/Helas3/Encode_H3K27me3_H3K36me3_EZH2/Output",paste0(name_object,'.RData')))
  epigraHMM_object <- get(name_object)
  rm(list = name_object)
  
  # Mapping labels
  labels <-
    unlist(lapply(metadata(epigraHMM_object)$control$pattern, function(x) {
      paste(unique(colData(epigraHMM_object)$condition)[x], collapse = ' & ')
    }))
  
  # Organizing output
  dt_list[[b]] <- data.table(B = b,
                             BIC = tail(metadata(epigraHMM_object)$history$control$bic,1),
                             Patterns = list(labels))
}

dt_bic <- rbindlist(dt_list)

# Plotting
fig_bic <- ggplot(data = dt_bic,aes(x = Patterns,y = BIC/1e6,group = 1)) +
  geom_line() +
  geom_text(aes(label = B),vjust = -0.75,hjust = -0.75) +
  expand_limits(y = 139.25) +
  geom_point(aes(shape = (B==2),
                 color = (B==2),
                 fill = (B==2)),size = 2.5) +
  scale_shape_manual(values = c(21,23))+
  scale_color_manual(values = c('black','red')) +
  scale_fill_manual(values = c('black','red')) +
  theme_bw() + 
  theme(legend.position = "none") +
  labs(y = bquote('BIC ( x'*10^-6*')'),x = 'Combinatorial Patterns') +
  scale_x_upset(order_by = 'degree',sets = c('EZH2','H3K27me3','H3K36me3',
                                             'EZH2 & H3K27me3','EZH2 & H3K36me3','H3K27me3 & H3K36me3')) 

save(fig_bic, file = './Figure_BIC.RData')
