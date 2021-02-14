library(magrittr)
library(data.table)
library(tidyverse)
library(tibble)
library(RColorBrewer)
library(GenomicRanges)
library(SummarizedExperiment)
library(vsn)
library(DESeq2)
library(ggpubr)
library(Polychrome)

# General parameters
bp = 250
chromosome = paste0('chr',1:22)
mark = 'CTCF'
data = 'Encode_twocells'
fdr = 0.05
sizetext = 7
sizetitle = 9

agg = function(x,gr.peak,gr.genome,cells = c('Helas3.adj','Hepg2.adj')){
  gr.summarized = as.data.frame(GenomicRanges::findOverlaps(gr.peak,gr.genome))
  gr.summarized = as.data.table(cbind(gr.summarized,x[gr.summarized$subjectHits,cells,with=F]));colnames(gr.summarized)=c('Peak','Genome',cells)
  gr.summarized[,Width :=  width(gr.peak[Peak])]
  gr.summarized <- gr.summarized[,.(Summary_Helas3.adj = sum(Helas3.adj,na.rm = T),
                                    Summary_Hepg2.adj = sum(Hepg2.adj,na.rm = T),
                                    Summary_Width = unique(Width)),by=Peak]
  
  return(gr.summarized)
}

getcor <- function(gr,p){
  newgr <- gr[order(-abs(log2(Summary_Hepg2.adj+1)-log2(Summary_Helas3.adj+1))),][1:round(.N*p),]
  newgr.UP = newgr[(Summary_Hepg2.adj>=Summary_Helas3.adj),]
  newgr.DOWN = newgr[(Summary_Hepg2.adj<Summary_Helas3.adj),]
  
  return(data.table(Up.LFC = ifelse(newgr.UP[,.N]>1,newgr.UP[,median(log2((Summary_Hepg2.adj+1)/(Summary_Helas3.adj+1)))],NA),
                    Down.LFC = ifelse(newgr.DOWN[,.N]>1,newgr.DOWN[,median(log2((Summary_Hepg2.adj+1)/(Summary_Helas3.adj+1)))],NA),
                    Cor = newgr[,cor(log2(Summary_Helas3.adj+1),log2(Summary_Hepg2.adj+1),method = 'spearman')]))
}

# Color of peak calls
methods <- c("ChIPComp + MACS2","csaw", "DiffBind + MACS2","diffReps","epigraHMM",'RSEG','THOR','Genes')
colors = c(kelly.colors(22)[c('red','yellow','purplishpink','yellowgreen','lightblue','buff','orange')],'black')
names(colors) <- methods 

# Loading data
load(
  paste0(
    '../../Data/Encode_helas3/',
    mark,
    '/wgEncodeBroadHistoneHelas3CtcfStdAlnRep1.markdup.q10.sorted.RData'
  )
)
helas31 = subset(counts[[paste0(bp)]], chr %in% chromosome)
rm(counts)
load(
  paste0(
    '../../Data/Encode_helas3/',
    mark,
    '/wgEncodeBroadHistoneHelas3CtcfStdAlnRep2.markdup.q10.sorted.RData'
  )
)
helas32 = subset(counts[[paste0(bp)]], chr %in% chromosome)
rm(counts)

load(
  paste0(
    '../../Data/Encode_hepg2/',
    mark,
    '/wgEncodeBroadHistoneHepg2CtcfStdAlnRep1.markdup.q10.sorted.RData'
  )
)
hepg21 = subset(counts[[paste0(bp)]], chr %in% chromosome)
rm(counts)
load(
  paste0(
    '../../Data/Encode_hepg2/',
    mark,
    '/wgEncodeBroadHistoneHepg2CtcfStdAlnRep2.markdup.q10.sorted.RData'
  )
)
hepg22 = subset(counts[[paste0(bp)]], chr %in% chromosome)

# Setting up data
counts = subset(counts[[paste0(bp)]],chr%in%chromosome,select=c('chr','start','stop'))
gr.counts = with(counts,GRanges(chr, IRanges(start=start, end=stop)))

dt <- data.table::data.table(counts,
                             Helas3_1 = helas31$counts,Helas3_2 = helas32$counts,
                             Hepg2_1 = hepg21$counts,Hepg2_2 = hepg22$counts)

# Normalizing reads and summing across replicates
dt.se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = as.matrix(dt[, c('Helas3_1', 'Helas3_2', 'Hepg2_1', 'Hepg2_2')])),
  rowRanges = GRanges(dt$chr, IRanges::IRanges(dt$start, dt$stop)),
  colData = data.frame(id = c(
    'Helas3_1', 'Helas3_2', 'Hepg2_1', 'Hepg2_2'
  ))
)
dt.se <-
  epigraHMM::normalizeCounts(dt.se, epigraHMM::controlEM(),span = 1)

dt[,(paste0(c('Helas3_1','Helas3_2','Hepg2_1','Hepg2_2'),'.adj')) := .SD/exp(assay(dt.se, 'offsets')),
   .SDcols=c('Helas3_1','Helas3_2','Hepg2_1','Hepg2_2')]

dt[,Helas3.adj := rowSums(.SD),.SDcols = paste0(c('Helas3_1','Helas3_2'),'.adj')]
dt[,Hepg2.adj := rowSums(.SD),.SDcols = paste0(c('Hepg2_1','Hepg2_2'),'.adj')]

# Loading peaks: epigraHMM
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
gr.epigrahmm <-
  GenomicRanges::reduce(gr.epigrahmm[epigraHMM:::fdrControl(pp.epigrahmm[, 2], fdr = fdr)])
dt.epigrahmm = agg(x = dt,gr.peak = gr.epigrahmm,gr.genome = gr.counts)

# Loading peaks: csaw
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
gr.csaw <- GenomicRanges::reduce(gr.csaw[gr.csaw$FDR<fdr])
dt.csaw = agg(x = dt,gr.peak = gr.csaw,gr.genome = gr.counts)

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
gr.thor <- GenomicRanges::reduce(gr.thor[gr.thor$FDR<fdr])
dt.thor = agg(x = dt,gr.peak = gr.thor,gr.genome = gr.counts)

# Loading peaks: RSEG
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
dt.rseg = agg(x = dt,gr.peak = gr.rseg,gr.genome = gr.counts)

# Loading peaks: diffReps
diffreps = fread(list.files(
  file.path('../../Public/diffReps', mark, data, paste0('Output', bp)),
  '.annotated',
  full.names = TRUE
),
header = T)
gr.diffreps = with(diffreps, GenomicRanges::GRanges(Chrom, IRanges(Start, End)))
gr.diffreps$FDR = diffreps$padj #It already gives me adjusted p-values
gr.diffreps <- gr.diffreps[seqnames(gr.diffreps) %in% chromosome]
gr.diffreps = GenomicRanges::reduce(gr.diffreps[gr.diffreps$FDR<fdr])
dt.diffreps = agg(x = dt,gr.peak = gr.diffreps,gr.genome = gr.counts)

# Loading peaks: DiffBind
diffbind = fread(list.files(
  file.path('../../Public/DiffBind', mark, data, paste0('Output', bp)),
  '.txt',
  full.names = TRUE
), header = T)
gr.diffbind = with(diffbind, GenomicRanges::GRanges(seqnames, IRanges(start, end)))
gr.diffbind$FDR = diffbind$FDR #It already gives me adjusted p-values
gr.diffbind <- gr.diffbind[seqnames(gr.diffbind) %in% chromosome]
gr.diffbind = GenomicRanges::reduce(gr.diffbind[gr.diffbind$FDR<fdr])
dt.diffbind = agg(x = dt,gr.peak = gr.diffbind,gr.genome = gr.counts)

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
gr.chipcomp = GenomicRanges::reduce(gr.chipcomp[gr.chipcomp$FDR<fdr])
dt.chipcomp = agg(x = dt,gr.peak = gr.chipcomp,gr.genome = gr.counts)

# Looping over
quant = c(seq(0.5,0.8,0.05),seq(0.825,0.9,0.025),seq(0.91,0.95,0.01),seq(0.955,0.995,0.005))
quant = c((1-quant[length(quant):1])[1:(length(quant)-1)],quant)

tpr = list()
tpr[['epigraHMM']] = list()
tpr[['thor']] = list()
tpr[['csaw']] = list()
tpr[['rseg']] = list()
tpr[['diffreps']] = list()
tpr[['diffbind']] = list()
tpr[['chipcomp']] = list()

for(i in quant){
  # cat('Quantile: ',i,'\n')
  
  ## TPR
  # epigraHMM
  tpr[['epigraHMM']] <- rbindlist(list(tpr[['epigraHMM']],data.table(Method = 'epigraHMM',Top = 100-100*i,getcor(gr = dt.epigrahmm,p = 1-i))))
  # THOR
  tpr[['thor']] <- rbindlist(list(tpr[['thor']],data.table(Method = 'THOR',Top = 100-100*i,getcor(gr = dt.thor,p = 1-i))))
  # csaw
  tpr[['csaw']] <- rbindlist(list(tpr[['csaw']],data.table(Method = 'csaw',Top = 100-100*i,getcor(gr = dt.csaw,p = 1-i))))
  # RSEG
  tpr[['rseg']] <- rbindlist(list(tpr[['rseg']],data.table(Method = 'RSEG',Top = 100-100*i,getcor(gr = dt.rseg,p = 1-i))))
  # diffReps
  tpr[['diffreps']] <- rbindlist(list(tpr[['diffreps']],data.table(Method = 'diffReps',Top = 100-100*i,getcor(gr = dt.diffreps,p = 1-i))))
  # DiffBind
  tpr[['diffbind']] <- rbindlist(list(tpr[['diffbind']],data.table(Method = 'DiffBind + MACS2',Top = 100-100*i,getcor(gr = dt.diffbind,p = 1-i))))
  # ChIPComp
  tpr[['chipcomp']] <- rbindlist(list(tpr[['chipcomp']],data.table(Method = 'ChIPComp + MACS2',Top = 100-100*i,getcor(gr = dt.chipcomp,p = 1-i))))
}

tpr.summary = rbindlist(tpr)

out.summary = tpr.summary[,c('Method','Top','Up.LFC','Down.LFC','Cor')]

out.summary.lfc <- melt(out.summary,id.vars = c('Method','Top'),measure.vars = c('Up.LFC','Down.LFC'),variable.name = 'Type',value.name = 'LFC')

out.summary.lfc <- rbindlist(list(out.summary.lfc,
                                  data.table(Method = 'Genes', Top = 99.5, Type = 'Up.LFC', LFC = 0),
                                  data.table(Method = 'Genes', Top = 99.5, Type = 'Down.LFC', LFC = 0)))

out.summary.lfc$Method %<>% factor(levels = c('epigraHMM','ChIPComp + MACS2','csaw','DiffBind + MACS2','diffReps','RSEG','THOR','Genes'))

fig.lfc = ggplot(data = out.summary.lfc,aes(y=LFC,x=Top)) +
  geom_smooth(aes(colour=Method,linetype = Type),se = F) +
  geom_line(aes(colour = Method,linetype = Type),alpha = 0)+
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'darkgrey') +
  scale_linetype_discrete(name = NULL,labels = c(expression(''%up% Hepg2 %down% Helas3),expression(''%down% Hepg2 %up% Helas3)))+
  scale_y_continuous(labels = function(x) sprintf("%.1f", x))+
  scale_color_manual(values=colors,breaks = methods)+
  theme_bw()+
  ylab('Log2 FC of CTCF Counts')+xlab('Top (%) Peaks')+
  theme(axis.text = element_text(size = sizetext),axis.title = element_text(size = sizetitle))+
  theme(legend.position = 'top',legend.direction = 'horizontal',legend.box = 'vertical',legend.spacing.y = unit(-0.25, "cm")) + guides(col = guide_legend(nrow = 1,order = 1),linetype = guide_legend(nrow = 1,order = 2,override.aes = list(col = 'black'),title = 'Observed Enrichment'))

save(fig.lfc, file = 'Figure_Correlation_CTCF_Encode_twocells_250bp.RData')
