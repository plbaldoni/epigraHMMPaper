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
library(epigraHMM)

bp = 750
chromosome = 'chr19'
ngroups = 2
mark = 'H3K36me3'
data = 'Encode_twocells'
size = 1.25
fdr = c(0.01,0.05,0.10,0.15,0.20)
methods <- c('Genes','0.01','0.05','0.10','0.15','0.20','Viterbi')
colors = c('#000000',brewer.pal(6,'Dark2'))
names(colors) <- methods 

peakpos = c(1.1,1.1+seq(0.1,by = 0.1,length.out = length(fdr)+1))
names(peakpos) <- methods 

# Loading data
load(paste0('../../Data/Encode_helas3/H3K36me3/wgEncodeBroadHistoneHelas3H3k36me3StdAlnRep1.markdup.q10.sorted.RData'))
helas31 = subset(counts[[paste0(bp)]],chr==chromosome);rm(counts)
load(paste0('../../Data/Encode_helas3/H3K36me3/wgEncodeBroadHistoneHelas3H3k36me3StdAlnRep2.markdup.q10.sorted.RData'))
helas32 = subset(counts[[paste0(bp)]],chr==chromosome);rm(counts)

load(paste0('../../Data/Encode_hepg2/H3K36me3/wgEncodeBroadHistoneHepg2H3k36me3StdAlnRep1.markdup.q10.sorted.RData'))
hepg21 = subset(counts[[paste0(bp)]],chr==chromosome);rm(counts)
load(paste0('../../Data/Encode_hepg2/H3K36me3/wgEncodeBroadHistoneHepg2H3k36me3StdAlnRep2.markdup.q10.sorted.RData'))
hepg22 = subset(counts[[paste0(bp)]],chr==chromosome)

rawcounts = as.data.table(counts[[paste0(bp)]])
rawcounts = rawcounts[chr%in%paste0('chr',1:22),]
gr.rawcounts <- rawcounts[,GRanges(chr, IRanges(start=start, end=stop))]

idx.chr = (counts[[paste0(bp)]]$chr==chromosome)

counts = subset(counts[[paste0(bp)]],chr==chromosome,select=c('chr','start','stop'))
gr.counts = with(counts,GRanges(chr, IRanges(start=start, end=stop)))

# Combining data
ChIP = cbind(counts,Window=1:nrow(counts))
ChIP = cbind(ChIP,Helas3_1=helas31$counts)
ChIP = cbind(ChIP,Helas3_2=helas32$counts)
ChIP = cbind(ChIP,Hepg2_1=hepg21$counts)
ChIP = cbind(ChIP,Hepg2_2=hepg22$counts)
ChIP <- as.data.table(ChIP)
ChIP.se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = as.matrix(ChIP[, c('Helas3_1', 'Helas3_2', 'Hepg2_1', 'Hepg2_2')])),
  rowRanges = GRanges(ChIP$chr, IRanges::IRanges(ChIP$start, ChIP$stop)),
  colData = data.frame(id = c(
    'Helas3_1', 'Helas3_2', 'Hepg2_1', 'Hepg2_2'
  ))
)
ChIP.se <-
  epigraHMM::normalizeCounts(ChIP.se, epigraHMM::controlEM(),span = 1)
ChIP[, (paste0(c(
  'Helas3_1', 'Helas3_2', 'Hepg2_1', 'Hepg2_2'
), '.adj')) := .SD / exp(assay(ChIP.se, 'offsets')),
.SDcols = c('Helas3_1', 'Helas3_2', 'Hepg2_1', 'Hepg2_2')]
ChIP[,Helas3 := rowSums(.SD),.SDcols = c('Helas3_1.adj','Helas3_2.adj')]
ChIP[,Hepg2 := rowSums(.SD),.SDcols = c('Hepg2_1.adj','Hepg2_2.adj')]

ChIP <- ChIP[,c('chr','start','stop','Window','Helas3','Hepg2')] %>%
  as.tbl() %>%
  gather(Group,Counts,Helas3:Hepg2)
ChIP$Mark = mark

# Genomic Ranges to Plot
idx.genome = c(8873645,9309280)
idx = which.min(abs(counts$start-idx.genome[1])):which.min(abs(counts$start-idx.genome[2]))

# Loading epigraHMM
load(list.files(
  file.path('../../Public/epigraHMM', mark, data, 'Output'),
  paste0('Output_', bp, 'bp.RData'),
  full.names = TRUE
))
pp.epigrahmm = epigraHMM:::getMarginalProbability(list.files(
  file.path('../../Public/epigraHMM', mark, data, 'Output'),
  paste0(bp, 'bp.h5'),
  full.names = TRUE
))
gr.epigrahmm = rowRanges(get(paste(
  'epigraHMM', mark, data, 'Output', paste0(bp, 'bp'), sep = '_'
)))

gr.epigrahmm.1 <-GenomicRanges::reduce(gr.epigrahmm[epigraHMM:::fdrControl(pp.epigrahmm[, 2], fdr = fdr[1]) & seqnames(gr.epigrahmm)==chromosome])
gr.epigrahmm.5 <-GenomicRanges::reduce(gr.epigrahmm[epigraHMM:::fdrControl(pp.epigrahmm[, 2], fdr = fdr[2]) & seqnames(gr.epigrahmm)==chromosome])
gr.epigrahmm.10 <-GenomicRanges::reduce(gr.epigrahmm[epigraHMM:::fdrControl(pp.epigrahmm[, 2], fdr = fdr[3]) & seqnames(gr.epigrahmm)==chromosome])
gr.epigrahmm.15 <-GenomicRanges::reduce(gr.epigrahmm[epigraHMM:::fdrControl(pp.epigrahmm[, 2], fdr = fdr[4]) & seqnames(gr.epigrahmm)==chromosome])
gr.epigrahmm.20 <-GenomicRanges::reduce(gr.epigrahmm[epigraHMM:::fdrControl(pp.epigrahmm[, 2], fdr = fdr[5]) & seqnames(gr.epigrahmm)==chromosome])

epigraHMM.1 = cbind(counts,Window=1:nrow(counts),Group='Helas3',Method='0.01')
epigraHMM.1 = cbind(epigraHMM.1,Output=1*overlapsAny(gr.counts,gr.epigrahmm.1))
epigraHMM.1$Counts = ifelse(epigraHMM.1$Output==1,max(ChIP[ChIP$Window%in%idx,'Counts'])*peakpos['0.01'],NA)
epigraHMM.1 <- epigraHMM.1 %>% as.tbl()
epigraHMM.1$Output %<>% mapvalues(from=0:1,to=c('Non-differential','Differential'))

epigraHMM.5 = cbind(counts,Window=1:nrow(counts),Group='Helas3',Method='0.05')
epigraHMM.5 = cbind(epigraHMM.5,Output=1*overlapsAny(gr.counts,gr.epigrahmm.5))
epigraHMM.5$Counts = ifelse(epigraHMM.5$Output==1,max(ChIP[ChIP$Window%in%idx,'Counts'])*peakpos['0.05'],NA)
epigraHMM.5 <- epigraHMM.5 %>% as.tbl()
epigraHMM.5$Output %<>% mapvalues(from=0:1,to=c('Non-differential','Differential'))

epigraHMM.10 = cbind(counts,Window=1:nrow(counts),Group='Helas3',Method='0.10')
epigraHMM.10 = cbind(epigraHMM.10,Output=1*overlapsAny(gr.counts,gr.epigrahmm.10))
epigraHMM.10$Counts = ifelse(epigraHMM.10$Output==1,max(ChIP[ChIP$Window%in%idx,'Counts'])*peakpos['0.10'],NA)
epigraHMM.10 <- epigraHMM.10 %>% as.tbl()
epigraHMM.10$Output %<>% mapvalues(from=0:1,to=c('Non-differential','Differential'))

epigraHMM.15 = cbind(counts,Window=1:nrow(counts),Group='Helas3',Method='0.15')
epigraHMM.15 = cbind(epigraHMM.15,Output=1*overlapsAny(gr.counts,gr.epigrahmm.15))
epigraHMM.15$Counts = ifelse(epigraHMM.15$Output==1,max(ChIP[ChIP$Window%in%idx,'Counts'])*peakpos['0.15'],NA)
epigraHMM.15 <- epigraHMM.15 %>% as.tbl()
epigraHMM.15$Output %<>% mapvalues(from=0:1,to=c('Non-differential','Differential'))

epigraHMM.20 = cbind(counts,Window=1:nrow(counts),Group='Helas3',Method='0.20')
epigraHMM.20 = cbind(epigraHMM.20,Output=1*overlapsAny(gr.counts,gr.epigrahmm.20))
epigraHMM.20$Counts = ifelse(epigraHMM.20$Output==1,max(ChIP[ChIP$Window%in%idx,'Counts'])*peakpos['0.20'],NA)
epigraHMM.20 <- epigraHMM.20 %>% as.tbl()
epigraHMM.20$Output %<>% mapvalues(from=0:1,to=c('Non-differential','Differential'))

# Getting viterbi sequence
viterbi <-
  rhdf5::h5read(metadata(get(paste(
    'epigraHMM', mark, data, 'Output', paste0(bp, 'bp'), sep = '_'
  )))$output, "viterbi")[, 1]
gr.viterbi <-
  rowRanges(get(paste(
    'epigraHMM', mark, data, 'Output', paste0(bp, 'bp'), sep = '_'
  )))[viterbi == 1]

Viterbi = cbind(counts,Window=1:nrow(counts),Group='Helas3',Method='Viterbi')
Viterbi = cbind(Viterbi,Output=1*overlapsAny(gr.counts,gr.viterbi))
Viterbi$Counts = ifelse(Viterbi$Output==1,max(ChIP[ChIP$Window%in%idx,'Counts'])*peakpos['Viterbi'],NA)
Viterbi <- Viterbi %>%as.tbl()
Viterbi$Output %<>% mapvalues(from=0:1,to=c('Non-differential','Differential'))

# Loading refseq genes
load('../../Public/Salmon//UCSC_RefSeq_Hg19_Genes.RData')
gr.genes <- gr.genes[overlapsAny(gr.genes, gr.counts)]

# Organizing refseq
refseq.out = cbind(counts,Window=1:nrow(counts),Group='Helas3',Method='Genes')
refseq.out = cbind(refseq.out,Output=1*overlapsAny(gr.counts,gr.genes))
refseq.out$Counts = ifelse(refseq.out$Output==1,max(ChIP[ChIP$Window%in%idx,'Counts'])*peakpos['Genes'],NA)
refseq.out <- refseq.out %>% as.tbl()
refseq.out$Output %<>% mapvalues(from=0:1,to=c('Non-differential','Differential'))

# Organizing the data
dt.segment = rbindlist(list(as.data.table(epigraHMM.1),
                            as.data.table(epigraHMM.5),
                            as.data.table(epigraHMM.10),
                            as.data.table(epigraHMM.15),
                            as.data.table(epigraHMM.20),
                            as.data.table(Viterbi),
                            as.data.table(refseq.out)))

## Plotting Peak calls

fig1.750 = ggplot(data=ChIP[ChIP$Window%in%idx,],aes(x=start,y=Counts))+
  geom_line()+
  facet_grid(rows=vars(Group))+
  annotate('rect',alpha = 0.1,xmin = 8950000,xmax = 9125000,ymin = 0,ymax = 1.05*max(dt.segment$Counts,na.rm = T)) + 
  geom_segment(data=dt.segment[Window%in%idx,],aes(x=start,xend=stop,y=Counts,yend=Counts,color = Method),size = size)+
  scale_color_manual(values = colors,breaks = names(colors))+
  scale_x_continuous(labels = scales::comma)+
  xlab(paste0('Genomic Window (',chromosome,')'))+ylab('Normalized ChIP Counts')+
  theme_bw()+
  theme(axis.title.x = element_blank(),axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank())

## Plotting Posterior Probabilities

PostProb = cbind(counts,Window=1:nrow(counts),PostProb2 = pp.epigrahmm[which(idx.chr == T), 2],Label='Prob.')
PostProb <- PostProb %>%
  as.tbl() %>%
  gather(Component,PP,PostProb2)
PostProb$Component %<>% mapvalues(from=c('PostProb2'),to=c('Differential'))

fig2.750 = ggplot(data=PostProb[PostProb$Window%in%idx,],aes(x=start,y=PP,fill=Component))+facet_grid(rows=vars(Label))+
  geom_area(position='identity',alpha=1)+
  scale_y_continuous(limits = c(0,1),breaks=c(0,1),labels = function(x) sprintf("%.1f", x))+
  scale_x_continuous(limits = range(PostProb[PostProb$Window%in%idx,'start']),labels = scales::comma)+
  xlab(paste0('Genomic Window (',chromosome,')'))+
  ylab('Prob.')+
  guides(fill=F)+
  scale_fill_manual(values = c('Differential'=as.character(colors['Genes'])))+
  theme_bw()+
  theme(legend.position="none",strip.text.y = element_text(colour = alpha('grey',0.0)))


fig.750 = ggarrange(fig1.750,fig2.750,heights = c(0.8,0.2),ncol = 1,nrow = 2,legend = 'none')
save(fig1.750,fig2.750,fig.750,file = 'Figure_PeakCalls_H3K36me3_Encode_twocells_750bp.RData')