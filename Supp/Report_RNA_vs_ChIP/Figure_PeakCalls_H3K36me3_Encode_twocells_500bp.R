library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(reshape2)
library(tidyr)
library(magrittr)
library(plyr)
library(ggpubr)
library(data.table)
library(RColorBrewer)
library(Polychrome)
library(SummarizedExperiment)

bp = 500
mark = 'H3K36me3'
data = 'Encode_twocells'
fdr = 0.05
chromosome = paste0('chr',1:22)

# Loading epigraHMM
load(list.files(
  file.path('../../Public/epigraHMM', mark, data, 'Output'),
  paste0('Output_', bp, 'bp.RData'),
  full.names = TRUE
))

gr.epigrahmm <- epigraHMM::callPeaks(object = epigraHMM_H3K36me3_Encode_twocells_Output_500bp,
                                     hdf5 = metadata(epigraHMM_H3K36me3_Encode_twocells_Output_500bp)$output,
                                     method = 0.05)

# Loading csaw
csaw <-
  fread(list.files(
    file.path('../../Public/csaw', mark, data, paste0('Output', bp)),
    '.tsv.gz',
    full.names = TRUE
  ))
gr.csaw <-
  with(csaw, GRanges(seqnames, IRanges(start = start, end = end)))
gr.csaw$FDR <- csaw$FDR
gr.csaw <- gr.csaw[seqnames(gr.csaw) %in% chromosome]
gr.csaw <- gr.csaw[gr.csaw$FDR < fdr]

# Loading RSEG
load(list.files(
  file.path('../../Public/RSEG', mark, data, paste0('Output', bp)),
  'Output*.*.RData$',
  full.names = TRUE
))
rseg[, postprob :=  (V4 + V5)]
gr.rseg <- with(rseg, GRanges(V1, IRanges(start = V2, end = V3)))
gr.rseg$FDR = rseg$postprob
gr.rseg <-
  gr.rseg[gr.rseg$FDR >= 0 &
            gr.rseg$FDR <= 1] # There are a few > 1 cases, which I am excluding
gr.rseg <-
  gr.rseg[epigraHMM:::fdrControl(gr.rseg$FDR, fdr = fdr)] # Filtering significant windows
gr.rseg <-
  gr.rseg[seqnames(gr.rseg) %in% chromosome] # Selecting only chromosomes of interest
gr.rseg = GenomicRanges::reduce(gr.rseg)

# Loading diffReps
diffreps = fread(list.files(
  file.path('../../Public/diffReps', mark, data, paste0('Output', bp)),
  '.annotated',
  full.names = TRUE
),
header = T)
gr.diffreps = with(diffreps, GenomicRanges::GRanges(Chrom, IRanges(Start, End)))
gr.diffreps$FDR = diffreps$padj #It already gives me adjusted p-values
gr.diffreps <- gr.diffreps[seqnames(gr.diffreps) %in% chromosome]
gr.diffreps = gr.diffreps[gr.diffreps$FDR < fdr]

# Loading DiffBind
diffbind = fread(list.files(
  file.path('../../Public/DiffBind', mark, data, paste0('Output', bp)),
  '.txt',
  full.names = TRUE
), header = T)
gr.diffbind = with(diffbind, GenomicRanges::GRanges(seqnames, IRanges(start, end)))
gr.diffbind$FDR = diffbind$FDR #It already gives me adjusted p-values
gr.diffbind <- gr.diffbind[seqnames(gr.diffbind) %in% chromosome]
gr.diffbind = gr.diffbind[gr.diffbind$FDR < fdr]

# Loading peaks: ChIPComp
chipcomp = fread(list.files(
  file.path('../../Public/ChIPComp', mark, data, paste0('Output', bp)),
  '.txt',
  full.names = TRUE
), header = T)
chipcomp[, FDR := p.adjust(pvalue.wald, method = 'BH')]
gr.chipcomp = with(chipcomp, GenomicRanges::GRanges(chr, IRanges(start, end)))
gr.chipcomp$FDR = chipcomp$FDR
gr.chipcomp <- gr.chipcomp[seqnames(gr.chipcomp) %in% chromosome]
gr.chipcomp = gr.chipcomp[gr.chipcomp$FDR < fdr]

# Loading peaks: THOR
thor = fread(list.files(
  file.path('../../Public/THOR', mark, data, paste0('Output', bp)),
  paste0(bp, 'bp-diffpeaks.narrowPeak'),
  full.names = TRUE
))
thor$V8 %<>% as.numeric
gr.thor = with(thor, GenomicRanges::GRanges(V1, IRanges(V2, V3)))
gr.thor$FDR = 10 ^ (-thor$V8)  #It already gives me adjusted p-values
gr.thor <- gr.thor[seqnames(gr.thor) %in% chromosome]
gr.thor <- gr.thor[gr.thor$FDR < fdr]

# Writing BED files

writeBed = function(gr,name,fdr = 0.05){
  dt.bed <- data.frame(chrom=SummarizedExperiment::seqnames(gr),chromStart=SummarizedExperiment::start(gr),chromEnd=SummarizedExperiment::end(gr))
  
  utils::write.table(dt.bed,file="./temp.bed",row.names = FALSE,col.names = FALSE,quote = FALSE,sep="\t")
  header1 <- paste0('track name="',name,'" description="',paste0('FDR-controlled Peaks (FDR = ',fdr,')'),'"')
  header2 <- paste0('browser position ',dt.bed[1,'chrom'],':',dt.bed[1,'chromStart'],'-',dt.bed[1,'chromEnd'])
  
  system2('echo',paste0(header2,' | cat - ',"./temp.bed",' > ',"./temp1.bed"))
  system2('echo',paste0(header1,' | cat - ',"./temp1.bed",' > ',paste0(name,'.bed')))
  system2('rm',paste("./temp1.bed","./temp.bed"))
}

writeBed(gr.epigrahmm,'epigraHMM')
writeBed(gr.rseg,'RSEG')
writeBed(gr.thor,'THOR')
writeBed(gr.chipcomp,'ChIPComp-MACS2')
writeBed(gr.diffbind,'DiffBind-MACS2')
writeBed(gr.csaw,'csaw')
writeBed(gr.diffreps,'diffReps')
