# Figure F: genomic segmentation with Viterbi path

# Loading libraries
library(epigraHMM)
library(data.table)
library(magrittr)
library(plyr)
library(ggplot2)
library(ggrepel)
library(GenomicRanges)
library(SummarizedExperiment)
library(tidyverse)
library(ggnewscale)

# Some general parameters and functions

bp <- 500
chromosome <- 'chr1'#'chr19'
idx.genome <- c(161626774,161943334)#c(56588768,56892457)
cellline <- 'Helas3'
size.text <- 2
size.axis <- 6
h.track = c(75,62.5,50)
h.txt = c(75,62.5,50)

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

# Loading data

load(paste0('../../Data/Encode_helas3/H3K36me3/wgEncodeBroadHistoneHelas3H3k36me3StdAlnRep1.markdup.q10.sorted.RData'))
h3k36me3.rep1 = counts[[paste0(bp)]];rm(counts)
load(paste0('../../Data/Encode_helas3/H3K36me3/wgEncodeBroadHistoneHelas3H3k36me3StdAlnRep2.markdup.q10.sorted.RData'))
h3k36me3.rep2 = counts[[paste0(bp)]];rm(counts)

load(paste0('../../Data/Encode_helas3/H3K27me3/wgEncodeBroadHistoneHelas3H3k27me3StdAlnRep1.markdup.q10.sorted.RData'))
h3k27me3.rep1 = counts[[paste0(bp)]];rm(counts)
load(paste0('../../Data/Encode_helas3/H3K27me3/wgEncodeBroadHistoneHelas3H3k27me3StdAlnRep2.markdup.q10.sorted.RData'))
h3k27me3.rep2 = counts[[paste0(bp)]];rm(counts)

load(paste0('../../Data/Encode_helas3/EZH2/wgEncodeBroadHistoneHelas3Ezh239875AlnRep1.markdup.q10.sorted.RData'))
ezh2.rep1 = counts[[paste0(bp)]];rm(counts)
load(paste0('../../Data/Encode_helas3/EZH2/wgEncodeBroadHistoneHelas3Ezh239875AlnRep2.markdup.q10.sorted.RData'))
ezh2.rep2 = counts[[paste0(bp)]]

# Loading output from epigraHMM (reduced)
load(
  "../../Public/epigraHMM/Helas3/Encode_H3K27me3_H3K36me3_EZH2/Output/epigraHMM_Helas3_Encode_H3K27me3_H3K36me3_EZH2_Output_500bp_2Patterns.RData"
)
epigraHMM_object <-
  epigraHMM_Helas3_Encode_H3K27me3_H3K36me3_EZH2_Output_500bp_2Patterns
rm(epigraHMM_Helas3_Encode_H3K27me3_H3K36me3_EZH2_Output_500bp_2Patterns)

# Getting posterior probabilities
pp <- as.data.table(rhdf5::h5read(metadata(epigraHMM_object)$output,'mixtureProb'))
setnames(pp,  unlist(lapply(metadata(epigraHMM_object)$control$pattern, function(x) {
  paste(unique(colData(epigraHMM_object)$condition)[x], collapse = ' & ')
})))
pp[,Pattern := max.col(pp)]

# Getting viterbi sequence
viterbi <-
  rhdf5::h5read(metadata(epigraHMM_object)$output, "viterbi")[, 1]

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

# Genomic ranges

countsplot = subset(counts[[paste0(bp)]],chr==chromosome,select=c('chr','start','stop'))
gr.countsplot = with(countsplot,GRanges(chr, IRanges(start=start, end=stop)))

# Matching index

idx = which.min(abs(countsplot$start - idx.genome[1])):which.min(abs(countsplot$start -
                                                                       idx.genome[2]))

# Combining data

ChIP = cbind(countsplot,Window=1:nrow(countsplot))
ChIP = cbind(ChIP,H3K27me3_1=subset(h3k27me3.rep1,chr==chromosome)$counts)
ChIP = cbind(ChIP,H3K27me3_2=subset(h3k27me3.rep2,chr==chromosome)$counts)
ChIP = cbind(ChIP,H3K36me3_1=subset(h3k36me3.rep1,chr==chromosome)$counts)
ChIP = cbind(ChIP,H3K36me3_2=subset(h3k36me3.rep2,chr==chromosome)$counts)
ChIP = cbind(ChIP,EZH2_1=subset(ezh2.rep1,chr==chromosome)$counts)
ChIP = cbind(ChIP,EZH2_2=subset(ezh2.rep2,chr==chromosome)$counts)
ChIP = as.data.table(ChIP)

ChIP.se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = as.matrix(ChIP[, c('H3K27me3_1','H3K27me3_2','H3K36me3_1','H3K36me3_2','EZH2_1','EZH2_2')])),
  rowRanges = GRanges(ChIP$chr, IRanges::IRanges(ChIP$start, ChIP$stop)),
  colData = data.frame(id = c(
    c('H3K27me3_1','H3K27me3_2','H3K36me3_1','H3K36me3_2','EZH2_1','EZH2_2')
  ))
)

ChIP.se <-
  epigraHMM::normalizeCounts(ChIP.se, epigraHMM::controlEM(),span = 1)

ChIP[,(paste0(c('H3K27me3_1','H3K27me3_2','H3K36me3_1','H3K36me3_2','EZH2_1','EZH2_2'),'.adj')) := .SD/exp(assay(ChIP.se, 'offsets')),
     .SDcols=c('H3K27me3_1','H3K27me3_2','H3K36me3_1','H3K36me3_2','EZH2_1','EZH2_2')]
ChIP[,H3K27me3 := rowSums(.SD),.SDcols = c('H3K27me3_1.adj','H3K27me3_2.adj')]
ChIP[,H3K36me3 := rowSums(.SD),.SDcols = c('H3K36me3_1.adj','H3K36me3_2.adj')]
ChIP[,EZH2 := rowSums(.SD),.SDcols = c('EZH2_1.adj','EZH2_2.adj')]

ChIP <- ChIP[,c('chr','start','stop','Window','H3K27me3','H3K36me3','EZH2')] %>%
  as_tibble() %>%
  gather(Group,Counts,H3K27me3:EZH2)
ChIP$CellLine = cellline
ChIP$Group %<>% factor(levels = c('H3K36me3','H3K27me3','EZH2'))

# Setting up peaks

gr.differential.reduced <- GenomicRanges::reduce(gr.differential)
gr.differential.reduced <- gr.differential.reduced[seqnames(gr.differential.reduced)==chromosome]
gr.differential.reduced$Pattern <- pat.peaks(x = pp,gr.differential.reduced,gr.genome = rowRanges(epigraHMM_object))
  
epigraHMM = cbind(countsplot,Window=1:nrow(countsplot),Group='H3K36me3',Method='epigraHMM')
epigraHMM$Group %<>% factor(levels = c(colnames(pp)[-ncol(pp)],'Genes'))
epigraHMM = cbind(epigraHMM,Output=1*overlapsAny(gr.countsplot,gr.differential.reduced))
epigraHMM$Counts = ifelse(epigraHMM$Output==1,h.track[1],NA)
epigraHMM <- epigraHMM %>% as_tibble()
epigraHMM$Output %<>% mapvalues(from=0:1,to=c('Non-differential','Differential'))
epigraHMM %<>% as.data.table()
epigraHMM[,Pattern := NA]
epigraHMM$Pattern[queryHits(findOverlaps(gr.countsplot,gr.differential.reduced))] <- gr.differential.reduced$Pattern[subjectHits(findOverlaps(gr.countsplot,gr.differential.reduced))]

# Loading ChromHMM

chromhmm.out <- rbindlist(lapply(3:8,function(x){
  chromhmm <- fread(paste0('../../Public/ChromHMM/Helas3/Encode_H3K27me3_H3K36me3_EZH2/Output',bp,'/ChromHMM_Helas3_Encode_H3K27me3_H3K36me3_EZH2_',bp,'bp_',x,'_segments.bed'))
  setnames(chromhmm,c('seqname','start','end','Segment'))
  gr.chromhmm <- sort(GenomicRanges::makeGRangesFromDataFrame(df = chromhmm,keep.extra.columns = TRUE))
  
  chromhmm.tmp <- as.data.table(gr.chromhmm)[,c('seqnames','start','end','Segment')]
  chromhmm.tmp[,Overlap := overlapsAny(gr.chromhmm,with(ChIP[ChIP$Window%in%idx,],GRanges(chr,IRanges(start,stop))))]
  chromhmm.tmp[Overlap==T,Counts := h.track[2]]
  chromhmm.tmp[,Group := 'H3K36me3']
  chromhmm.tmp$Group %<>% factor(levels = c(colnames(pp)[-ncol(pp)],'Genes'))
  chromhmm.tmp[,States :=  x]
  return(chromhmm.tmp)
}))

# Loading genes
load('../../Public/Salmon/ENCODE.rnaseq.scaled.RData')
ENCODE.rnaseq.scaled <-
  ENCODE.rnaseq.scaled[, colData(ENCODE.rnaseq.scaled)$Cells %in% 'Helas3']
ENCODE.rnaseq.scaled$Cells <- droplevels(ENCODE.rnaseq.scaled$Cells)
ENCODE.rnaseq.scaled <-
  ENCODE.rnaseq.scaled[rowRanges(ENCODE.rnaseq.scaled)$gene_biotype == 'protein_coding']
ENCODE.rnaseq.scaled <-
  ENCODE.rnaseq.scaled[overlapsAny(rowRanges(ENCODE.rnaseq.scaled),
                                   GRanges(chromosome, IRanges(start = idx.genome[1], end = idx.genome[2])))]
ENCODE.rnaseq.scaled <-
  ENCODE.rnaseq.scaled[end(ENCODE.rnaseq.scaled) < idx.genome[2] &
                         start(ENCODE.rnaseq.scaled) > idx.genome[1]]

# Now, I want to see what the gene expresion from genes intersecting background, and HMM Diferential states 2 and 5 are
# First, I will take the geometric mean of the RNA-seq data from the two available replicates
dt.rna <- data.table(exp(rowMeans(log(assay(ENCODE.rnaseq.scaled)))))
setnames(dt.rna,paste0('RNA.',as.character(unique(ENCODE.rnaseq.scaled$Cells))))

# Loading refseq genes

refseq.out <- as.data.table(rowRanges(ENCODE.rnaseq.scaled))[,c('seqnames','start','end','gene_name')][,RNA.Helas3 := dt.rna$RNA.Helas3]
refseq.out[,Overlap := overlapsAny(rowRanges(ENCODE.rnaseq.scaled),with(ChIP[ChIP$Window%in%idx,],GRanges(chr,IRanges(start,stop))))]
refseq.out[Overlap==T,Counts := h.track[3]]
refseq.out[Overlap==T,Lbl.Counts := 45]
refseq.out[,Group := 'H3K36me3']
refseq.out[Overlap==T,Lbl.x := start+0.5*(end-start)]
refseq.out$Group %<>% factor(levels = c(colnames(pp)[-ncol(pp)],'Genes'))

# Organizing the data

dt.segment = rbindlist(list(as.data.table(epigraHMM)))

dt.segment$Pattern %<>% as.factor() %<>% mapvalues(from = seq_len(ncol(pp)-1),to = colnames(pp)[-ncol(pp)])
dt.segment$Pattern %<>% factor(levels = colnames(pp)[-ncol(pp)])

dt.anno1 <- data.table(Group = rep('H3K36me3',2),Label = c('H3K36me3','H3K27me3 & EZH2'),
                       x = c(34100000,34400000),
                       y = rep(h.txt[1],2))
dt.anno1$Group %<>% factor(levels = colnames(pp)[-ncol(pp)])

dt.anno2 <- data.table(Group = rep('H3K36me3',3),Label = c('epigraHMM','ChromHMM','Genes'),
                       x = rep(dt.segment[Window%in%idx,min(start)-0.10*(max(start)-min(start))],3),
                       y = h.txt)
dt.anno2$Group %<>% factor(levels = colnames(pp)[-ncol(pp)])

# Plotting read counts

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

fig_segmentation <- lapply(3:8,function(x){
  ggplot(data=ChIP[ChIP$Window%in%idx,],aes(x=start,y=Counts))+
    facet_grid(rows=vars(Group))+
    geom_line()+
    geom_rect(data=chromhmm.out[Group=='H3K36me3' & States == x & Overlap == TRUE,],aes(xmin=start,xmax=end,ymin=Counts,ymax=2.5+Counts,fill = Segment))+
    geom_segment(data=refseq.out[Group=='H3K36me3' & Overlap == TRUE,],aes(x=start,xend=end,y=Counts,yend=Counts,color = log2(RNA.Helas3+1)),size = 1.5)+
    geom_text(data = dt.anno2,aes(x = x, y = y, label = Label),size = size.text) +
    geom_text_repel(data=refseq.out,aes(x = Lbl.x,y = Counts,label = gene_name),direction = 'y',size = size.text,nudge_y = -10,segment.color = 'grey') +
    theme_bw() +
    scale_fill_manual(values = cbPalette[c(7,1,4,8)],name = 'ChromHMM\nSegmentation')+
    scale_color_viridis_c(option = 'B',name = 'Gene Expression\nlog2(TPM+1)')+
    labs(x = paste0('Genomic Window (',chromosome,')'),y = 'Normalized ChIP-seq counts') +
    scale_x_continuous(labels = scales::comma,limits = c(161580000,max(ChIP[ChIP$Window%in%idx,]$start))) +
    ylim(c(0,max(dt.anno2$y)*1.1)) +
    new_scale_fill() +
    geom_rect(data=dt.segment[Window%in%idx,],aes(xmin=start,xmax=stop,ymin=Counts,ymax=2.5+Counts,fill = Pattern))+
    scale_fill_manual(values = cbPalette[c(2,3,1)],name = ' epigraHMM\n Enrichment',na.translate=FALSE) +
    theme(legend.direction = 'horizontal',legend.position = 'top',
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
})
names(fig_segmentation) <- paste0(3:8)

save(fig_segmentation, file = './Figure_Segmentation.RData')
