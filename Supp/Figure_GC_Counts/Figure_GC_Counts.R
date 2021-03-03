library(ggplot2)
library(BSgenome.Hsapiens.UCSC.hg19)

bp = 500
ngroups = 2
mark = 'EZH2'
data = 'Encode_twocells'
chromosome = paste0('chr',c(1:22,'X'))

# Loading data
load(
  paste0(
    '../../Data/Encode_helas3/',
    mark,
    '/wgEncodeBroadHistoneHelas3Ezh239875AlnRep1.markdup.q10.sorted.RData'
  )
)
helas31 = subset(counts[[paste0(bp)]], chr %in% chromosome)
rm(counts)
load(
  paste0(
    '../../Data/Encode_helas3/',
    mark,
    '/wgEncodeBroadHistoneHelas3Ezh239875AlnRep2.markdup.q10.sorted.RData'
  )
)
helas32 = subset(counts[[paste0(bp)]], chr %in% chromosome)
rm(counts)

load(
  paste0(
    '../../Data/Encode_hepg2/',
    mark,
    '/wgEncodeBroadHistoneHepg2Ezh239875AlnRep1.markdup.q10.sorted.RData'
  )
)
hepg21 = subset(counts[[paste0(bp)]], chr %in% chromosome)
rm(counts)
load(
  paste0(
    '../../Data/Encode_hepg2/',
    mark,
    '/wgEncodeBroadHistoneHepg2Ezh239875AlnRep2.markdup.q10.sorted.RData'
  )
)
hepg22 = subset(counts[[paste0(bp)]], chr %in% chromosome)

countsraw <- as.data.table(counts[[paste0(bp)]])
countsraw <- countsraw[chr %in% chromosome, ]
gr.countsraw <- countsraw[, GRanges(chr, IRanges(start, stop))]

counts = subset(counts[[paste0(bp)]], chr %in% chromosome, select = c('chr', 'start', 'stop'))
gr.counts = with(counts, GRanges(chr, IRanges(start = start, end = stop)))

# Combining data
ChIP = cbind(counts, Window = 1:nrow(counts))
ChIP = cbind(ChIP, Helas3_1 = helas31$counts)
ChIP = cbind(ChIP, Helas3_2 = helas32$counts)
ChIP = cbind(ChIP, Hepg2_1 = hepg21$counts)
ChIP = cbind(ChIP, Hepg2_2 = hepg22$counts)
ChIP = as.data.table(ChIP)

# Getting GC content
seqs <- Views(Hsapiens,gr.counts)
gcpos <- alphabetFrequency(seqs)
gc <- rowSums(gcpos[,c('C','G')])/bp
gc[gc<0.28 | gc>0.67] <- NA

#Putting gc on data.table
ChIP[,GC := gc]

# Melting
pdf('Figure_GC_Counts.pdf',height = 6,width = 6)
par(mfrow = c(2,2))
ChIP[,{
  smoothScatter(log1p(Helas3_1)~GC,ylab = 'log1p(ChIP-seq)',main = 'Helas3 - Rep 1',ylim = c(0,6.25),xlim = c(0,1))
  smoothScatter(log1p(Helas3_2)~GC,ylab = 'log1p(ChIP-seq)',main = 'Helas3 - Rep 2',ylim = c(0,6.25),xlim = c(0,1))
  smoothScatter(log1p(Hepg2_1)~GC,ylab = 'log1p(ChIP-seq)',main = 'Hepg2 - Rep 1',ylim = c(0,6.25),xlim = c(0,1))
  smoothScatter(log1p(Hepg2_2)~GC,ylab = 'log1p(ChIP-seq)',main = 'Hepg2 - Rep 2',ylim = c(0,6.25),xlim = c(0,1))
  
}]
dev.off()
