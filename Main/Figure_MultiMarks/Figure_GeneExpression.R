# Figure D: genomic segmentation with Viterbi path

# Loading libraries
library(epigraHMM)
library(data.table)
library(magrittr)
library(plyr)
library(ggplot2)
library(ggrepel)

# Some general parameters and functions

B <- 2

Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  if (length(ux[tab == max(tab)]) == 1) {
    return(ux[tab == max(tab)])
  } else{
    out = ux[tab == max(tab)]
    return(max(out))
  }
}

pat.peaks = function(x,gr.peaks,gr.genome){
  ovp = findOverlaps(gr.peaks,gr.genome,ignore.strand = T)
  out = rbindlist(lapply(unique(queryHits(ovp)), FUN = function(i){data.table(Peak = i,Mode = Modes(x[subjectHits(ovp)[queryHits(ovp) == i],]$Pattern))}))
  Pattern = rep(NA,length(gr.peaks))
  Pattern[out$Peak] = out$Mode
  return(Pattern)
}

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
pp <- as.data.table(rhdf5::h5read(metadata(epigraHMM_object)$output,'mixtureProb'))
setnames(pp,  unlist(lapply(metadata(epigraHMM_object)$control$pattern, function(x) {
  paste(unique(colData(epigraHMM_object)$condition)[x], collapse = ' & ')
})))
pp[,Pattern := max.col(pp)]

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
  ENCODE.rnaseq.scaled[overlapsAny(rowRanges(ENCODE.rnaseq.scaled), rowRanges(epigraHMM_object))]

# Check whether windows overlap gene bodies (It doesn't matter now whether these genes are expressed or not, I will look into that later)
gr.background$overlapGenes <-
  overlapsAny(gr.background, rowRanges(ENCODE.rnaseq.scaled))
gr.differential$overlapGenes <-
  overlapsAny(gr.differential, rowRanges(ENCODE.rnaseq.scaled))
gr.enrichment$overlapGenes <-
  overlapsAny(gr.enrichment, rowRanges(ENCODE.rnaseq.scaled))

# Now, I want to see what the gene expresion from genes intersecting background, and HMM Diferential states 2 and 5 are
# First, I will take the geometric mean of the RNA-seq data from the two available replicates
dt.rna <- data.table(exp(rowMeans(log(assay(ENCODE.rnaseq.scaled)))))
setnames(dt.rna,paste0('RNA.',as.character(unique(ENCODE.rnaseq.scaled$Cells))))

# Checking genes with 0 counts on RNA counts 
dt.rna[,zero_RNA := (RNA.Helas3>0)]

# Checking genes that overlap with called background, differential, and enrichment
dt.rna[,overlapBackground := overlapsAny(rowRanges(ENCODE.rnaseq.scaled),gr.background)]
dt.rna[,overlapDifferential := overlapsAny(rowRanges(ENCODE.rnaseq.scaled),gr.differential)]
dt.rna[,overlapEnrichment := overlapsAny(rowRanges(ENCODE.rnaseq.scaled),gr.enrichment)]

# How many genes extend beyond background, differential, and enrichment
dt.rna[,multiOverlap := overlapBackground+overlapDifferential+overlapEnrichment][,sum(!(multiOverlap==1))]

# Checking posterior probabilities from gene bodies
dt.rna$Pattern <- pat.peaks(x = pp,gr.peaks = rowRanges(ENCODE.rnaseq.scaled),gr.genome = rowRanges(epigraHMM_object))
dt.rna$Pattern %<>% mapvalues(from = 1:B,to = paste0('delta',1:B))
dt.rna$Pattern %<>% factor(levels = paste0('delta',B:1))

# I will subset dt.rna. I want to look only at genes that do not extend beyond background, differential, and enrichment
dt.rna.subset <- dt.rna[overlapEnrichment==F & overlapBackground==F,][(multiOverlap==1),]
dt.rna.subset[,HMM := 'Differential']
dt.rna.subset[overlapBackground==T,NewPattern := 'delta0'][overlapBackground==F,NewPattern := Pattern]

# Plotting the gene expression from genes whose predicted pattern is 2 (H3K36me3 only) or 5 (H3K27me3 & EZH2), in addition to those overlapping background

label.delta1 <- c('No\nEnrichment',gsub('&','&\n',colnames(pp)[1:2]))
names(label.delta1) = paste0('delta',c(0,1,2))

fig_expression <- ggplot(data = dt.rna.subset,
                         aes(y = NewPattern, x = log2(RNA.Helas3+1),height = ..density..,fill = ..x..))+
  geom_density_ridges_gradient(stat = "density", trim = TRUE,scale = 1)+
  scale_fill_viridis_c(option = 'B',name = 'Gene\nExpression')+
  theme_bw()+
  xlab(expression(paste("Gene Expression (log"[2],"(TPM+1))")))+
  ylab('Associated Enrichment')+
  scale_y_discrete(labels = label.delta1)+
  guides(fill = F)

save(fig_expression, file = './Figure_GeneExpression.RData')
