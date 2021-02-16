# Figure A: genomic segmentation with Viterbi path

# Loading libraries
library(epigraHMM)
library(data.table)
library(magrittr)
library(plyr)
library(ggplot2)
library(ggrepel)

# Some general parameters
mbp = 1e6

# Loading output from epigraHMM (full model)
load(
  "../../Public/epigraHMM/Helas3/Encode_H3K27me3_H3K36me3_EZH2/Output/epigraHMM_Helas3_Encode_H3K27me3_H3K36me3_EZH2_Output_500bp_6Patterns.RData"
)
epigraHMM_object <-
  epigraHMM_Helas3_Encode_H3K27me3_H3K36me3_EZH2_Output_500bp_6Patterns
rm(epigraHMM_Helas3_Encode_H3K27me3_H3K36me3_EZH2_Output_500bp_6Patterns)

# Getting viterbi sequence
viterbi <-
  rhdf5::h5read(metadata(epigraHMM_object)$output, "viterbi")[, 1]

# Loading genes
load('../../Public/Salmon/ENCODE.rnaseq.raw.RData')
ENCODE.rnaseq.raw <-
  ENCODE.rnaseq.raw[, colData(ENCODE.rnaseq.raw)$Cells %in% 'Helas3']
ENCODE.rnaseq.raw$Cells <- droplevels(ENCODE.rnaseq.raw$Cells)
ENCODE.rnaseq.raw <-
  ENCODE.rnaseq.raw[rowRanges(ENCODE.rnaseq.raw)$gene_biotype == 'protein_coding']
ENCODE.rnaseq.raw <-
  ENCODE.rnaseq.raw[overlapsAny(rowRanges(ENCODE.rnaseq.raw), rowRanges(epigraHMM_object))]

# Creating genomic ranges
gr.background <-
  rowRanges(epigraHMM_object)[viterbi == 0]
gr.background$Type = 'Background'
gr.differential <-
  rowRanges(epigraHMM_object)[viterbi == 1]
gr.differential$Type = 'Differential'
gr.enrichment <-
  rowRanges(epigraHMM_object)[viterbi == 2]
gr.enrichment$Type = 'Enrichment'

# Check whether windows overlap gene bodies (It doesn't matter now whether these genes are expressed or not, I will look into that later)
gr.background$overlapGenes <-
  overlapsAny(gr.background, rowRanges(ENCODE.rnaseq.raw))
gr.differential$overlapGenes <-
  overlapsAny(gr.differential, rowRanges(ENCODE.rnaseq.raw))
gr.enrichment$overlapGenes <-
  overlapsAny(gr.enrichment, rowRanges(ENCODE.rnaseq.raw))

# Putting the results in a data table
dt.all <-
  rbindlist(list(
    as.data.table(gr.background),
    as.data.table(gr.differential),
    as.data.table(gr.enrichment)
  ))
dt.all <-
  dt.all[, .(Windows = sum(width) / mbp), by = c('Type', 'overlapGenes')]
dt.all$overlapGenes %<>% mapvalues(from = c(T, F), to = c('Genic', 'Intergenic'))
dt.all <-
  merge(dt.all, dt.all[, .(TotalWindows = sum(Windows)), by = 'Type'], by =
          'Type')
dt.all[, Pct := Windows / TotalWindows][, TotalWindows := NULL]
setnames(dt.all, old = 'overlapGenes', new = 'Region')

# Plotting
fig_viterbi <-
  ggplot(dt.all, aes(x = Type, y = Windows, group = Region)) +
  geom_col(aes(fill = Region)) +
  geom_text_repel(data = dt.all[!Type == 'Enrichment',],aes(label = paste0(round(100 * Pct), '%')),
                  position = position_stack(vjust = 0.5),
                  direction = "y") +
  geom_text_repel(data = dt.all[Type == 'Enrichment',],aes(label = paste0(round(100 * Pct), '%')),
                  position = position_stack(vjust = 25),
                  direction = "y") +
  ylab('Mega Base Pairs') + xlab(paste0('Genomic Window Classification (HMM Viterbi Path)')) +
  theme_bw() +
  scale_fill_grey(start = 0.7, end = 0.4) +
  theme(
    legend.position = c(0.8, 0.75),
    legend.background = element_rect(fill = alpha('white', 0)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

save(fig_viterbi, file = './Figure_Viterbi.RData')
