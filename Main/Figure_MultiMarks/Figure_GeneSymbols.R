# Figure E: genomic segmentation with Viterbi path

# Loading libraries
library(epigraHMM)
library(data.table)
library(magrittr)
library(plyr)
library(ggplot2)
library(ggrepel)
library(GenomicRanges)
library(SummarizedExperiment)

# Loading output from epigraHMM (reduced)
load(
  "../../Public/epigraHMM/Helas3/Encode_H3K27me3_H3K36me3_EZH2/Output/epigraHMM_Helas3_Encode_H3K27me3_H3K36me3_EZH2_Output_500bp_2Patterns.RData"
)
epigraHMM_object <-
  epigraHMM_Helas3_Encode_H3K27me3_H3K36me3_EZH2_Output_500bp_2Patterns
rm(epigraHMM_Helas3_Encode_H3K27me3_H3K36me3_EZH2_Output_500bp_2Patterns)

# Getting viterbi sequence
viterbi <-
  rhdf5::h5read(metadata(epigraHMM_object)$output, "viterbi")[, 1]

# Getting posterior probabilities
pp <-
  as.data.table(rhdf5::h5read(metadata(epigraHMM_object)$output, 'mixtureProb'))
setnames(pp,  unlist(lapply(metadata(epigraHMM_object)$control$pattern, function(x) {
  paste(unique(colData(epigraHMM_object)$condition)[x], collapse = ' & ')
})))
pp[, Pattern := max.col(pp)]

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

# Loading genes
load('../../Public/Salmon/ENCODE.rnaseq.scaled.RData')
ENCODE.rnaseq.scaled <-
  ENCODE.rnaseq.scaled[, colData(ENCODE.rnaseq.scaled)$Cells %in% 'Helas3']
ENCODE.rnaseq.scaled$Cells <- droplevels(ENCODE.rnaseq.scaled$Cells)
ENCODE.rnaseq.scaled <-
  ENCODE.rnaseq.scaled[rowRanges(ENCODE.rnaseq.scaled)$gene_biotype == 'protein_coding']
ENCODE.rnaseq.scaled <-
  ENCODE.rnaseq.scaled[overlapsAny(rowRanges(ENCODE.rnaseq.scaled),
                                   rowRanges(epigraHMM_object))]

# Subsetting marker genes of Helas3 cell line (DOI: 10.2144/000112900)
## GNAS (ENSG00000087460) does not match range from UCSC list
## CDK4 (ENSG00000135446) does not match range from UCSC list

helas3_genes <-
  c(
    'ANAPC7',
    'ASPM',
    'CDC2L2',
    'CENPF',
    'BRD3',
    'CDK5RAP1',
    'EIF2AK1',
    'GPI',
    'GSPT2',
    'INPPL1',
    'MARK3',
    'NUP210',
    'PER1',
    'PKN3',
    'RB1',
    'SMC3',
    'SMC4',
    'STAG2',
    'ULK1',
    'YWHAE',
    'CDC25C',
    'CDKN3',
    'CHEK2',
    'EIF4A2',
    'FANCD2',
    'PCSK7',
    'SDHC',
    'SUFU',
    'TAF15',
    'YWHAZ',
    #'CDK4',#'GNAS',
    'PLK1',
    'PSMD4',
    'RAF1',
    'RAP1D',
    'BTRC',
    'MLLT3'
  )

ENCODE.rnaseq.scaled <-
  ENCODE.rnaseq.scaled[rowData(ENCODE.rnaseq.scaled)$symbol %in% helas3_genes, ]

# For each marker gene, I want to check the posterior probability of enrichment for H3K36me3 alone, given differential region

dt_overlap <-
  rbindlist(lapply(names(ENCODE.rnaseq.scaled), function(x) {
    gr.differential.overlap <-
      gr.differential[overlapsAny(gr.differential, ENCODE.rnaseq.scaled[x])]
    return(data.table(
      Gene = rowData(ENCODE.rnaseq.scaled[x])$symbol,
      Prob = pp[overlapsAny(rowRanges(epigraHMM_object), gr.differential.overlap), H3K36me3]
    ))
  }))

# Plotting

fig_genesymbols <-
  dt_overlap[, .(Mean = mean(Prob)), by = 'Gene'][order(Mean), ] %>%
  ggplot(aes(x = reorder(Gene, Mean), y = Mean)) +
  geom_segment(aes(xend = Gene, y = 0, yend = Mean), color = "black") +
  geom_point(color = "#56B4E9",
             size = 4,
             alpha = 0.6) +
  theme_light() +
  coord_flip() +
  theme(panel.grid.major.y = element_blank()) +
  labs(x = 'Marker Genes (UCSC Gene Symbols)', y = 'Differential H3K36me3 Enrichment (Mixture Prop.)')

save(fig_genesymbols, file = './Figure_GeneSymbols.RData')
