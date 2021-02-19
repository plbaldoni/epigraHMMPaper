library(magrittr)
library(data.table)
library(tidyverse)
library(tibble)
library(RColorBrewer)
library(GenomicRanges)
library(SummarizedExperiment)
library(vsn)
library(DESeq2)
library(ggtern)

# General parameters
bp = 250
chromosome = paste0('chr',1:22)
mark = 'H3K27ac'
data = 'Encode_threecells'
top = 1
fdr = 0.05

Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  if(length(ux[tab == max(tab)])==1){
    return(ux[tab == max(tab)])
  } else{
    out = ux[tab == max(tab)]
    return(max(out))
  }
}

agg = function(x,gr.peak,gr.genome,cells=c('Helas3.adj','Hepg2.adj','Huvec.adj')){
  gr.summarized = as.data.frame(GenomicRanges::findOverlaps(gr.peak,gr.genome))
  gr.summarized = as.data.table(cbind(gr.summarized,x[gr.summarized$subjectHits,cells,with=F]));colnames(gr.summarized)=c('Peak','Genome',cells)
  gr.summarized[,Width :=  width(gr.peak[Peak])]
  gr.summarized <- gr.summarized[,.(Summary_Helas3 = sum(Helas3.adj,na.rm = T),
                                    Summary_Hepg2 = sum(Hepg2.adj,na.rm = T),
                                    Summary_Huvec = sum(Huvec.adj,na.rm = T),
                                    Summary_Width = unique(Width)),by=Peak]
  return(gr.summarized)
}

pat.peaks = function(x,gr.peaks,gr.genome){
  ovp = findOverlaps(gr.peaks,gr.genome)
  out = rbindlist(lapply(unique(queryHits(ovp)), FUN = function(i){data.table(Peak = i,Mode = Modes(x[subjectHits(ovp)[queryHits(ovp)==i],]$Pattern))}))
  Pattern = rep(NA,length(gr.peaks))
  Pattern[out$Peak] = out$Mode
  return(Pattern)
}

# Loading data
load(
  paste0(
    '../../Data/Encode_helas3/',
    mark,
    '/wgEncodeBroadHistoneHelas3H3k27acStdAlnRep1.markdup.q10.sorted.RData'
  )
)
helas31 = subset(counts[[paste0(bp)]], chr %in% chromosome)
rm(counts)
load(
  paste0(
    '../../Data/Encode_helas3/',
    mark,
    '/wgEncodeBroadHistoneHelas3H3k27acStdAlnRep2.markdup.q10.sorted.RData'
  )
)
helas32 = subset(counts[[paste0(bp)]], chr %in% chromosome)
rm(counts)

load(
  paste0(
    '../../Data/Encode_hepg2/',
    mark,
    '/wgEncodeBroadHistoneHepg2H3k27acStdAlnRep1.markdup.q10.sorted.RData'
  )
)
hepg21 = subset(counts[[paste0(bp)]], chr %in% chromosome)
rm(counts)
load(
  paste0(
    '../../Data/Encode_hepg2/',
    mark,
    '/wgEncodeBroadHistoneHepg2H3k27acStdAlnRep2.markdup.q10.sorted.RData'
  )
)
hepg22 = subset(counts[[paste0(bp)]], chr %in% chromosome)
rm(counts)

load(
  paste0(
    '../../Data/Encode_huvec/',
    mark,
    '/wgEncodeBroadHistoneHuvecH3k27acStdAlnRep1.markdup.q10.sorted.RData'
  )
)
huvec1 = subset(counts[[paste0(bp)]], chr %in% chromosome)
rm(counts)

load(
  paste0(
    '../../Data/Encode_huvec/',
    mark,
    '/wgEncodeBroadHistoneHuvecH3k27acStdAlnRep2.markdup.q10.sorted.RData'
  )
)
huvec2 = subset(counts[[paste0(bp)]], chr %in% chromosome)

# Setting up data
counts = subset(counts[[paste0(bp)]],chr%in%chromosome,select=c('chr','start','stop'))
gr.counts = with(counts,GRanges(chr, IRanges(start=start, end=stop)))

dt <- data.table::data.table(counts,
                             Helas3_1 = helas31$counts,Helas3_2 = helas32$counts,
                             Hepg2_1 = hepg21$counts,Hepg2_2 = hepg22$counts,
                             Huvec_1 = huvec1$counts,Huvec_2 = huvec2$counts)

# Normalizing reads and summing across replicates
ChIP.se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = as.matrix(dt[, c('Helas3_1', 'Helas3_2', 'Hepg2_1', 'Hepg2_2','Huvec_1','Huvec_2')])),
  rowRanges = GRanges(dt$chr, IRanges::IRanges(dt$start, dt$stop)),
  colData = data.frame(id = c(
    'Helas3_1', 'Helas3_2', 'Hepg2_1', 'Hepg2_2','Huvec_1','Huvec_2'
  ))
)
ChIP.se <-
  epigraHMM::normalizeCounts(ChIP.se, epigraHMM::controlEM(),span = 1)

dt[,(paste0(c('Helas3_1','Helas3_2','Hepg2_1','Hepg2_2','Huvec_1','Huvec_2'),'.adj')) := .SD/exp(assay(ChIP.se, 'offsets')),
   .SDcols=c('Helas3_1','Helas3_2','Hepg2_1','Hepg2_2','Huvec_1','Huvec_2'),by='chr']

dt[,Helas3.adj := rowSums(.SD),.SDcols = paste0(c('Helas3_1','Helas3_2'),'.adj')]
dt[,Hepg2.adj := rowSums(.SD),.SDcols = paste0(c('Hepg2_1','Hepg2_2'),'.adj')]
dt[,Huvec.adj := rowSums(.SD),.SDcols = paste0(c('Huvec_1','Huvec_2'),'.adj')]

# Loading peaks
load(list.files(
  file.path('../../Public/epigraHMM', mark, data, 'Output'),
  paste0('Output_', bp, 'bp.RData'),
  full.names = TRUE
))
epigraHMM_object <- get(paste(
  'epigraHMM', mark, data, 'Output', paste0(bp, 'bp'), sep = '_'
))

pp.epigrahmm = epigraHMM:::getMarginalProbability(list.files(
  file.path('../../Public/epigraHMM', mark, data, 'Output'),
  paste0(bp, 'bp.h5'),
  full.names = TRUE
))
gr.epigrahmm = rowRanges(epigraHMM_object)
gr.epigrahmm <-
  GenomicRanges::reduce(gr.epigrahmm[epigraHMM:::fdrControl(pp.epigrahmm[, 2], fdr = fdr)])
dt.epigrahmm = agg(x = dt,gr.peak = gr.epigrahmm,gr.genome = gr.counts)

# Loading posterior probabilities
pp.epigrahmm <- as.data.table(pp.epigrahmm)
pp.epigrahmm[,Pattern := max.col(pp.epigrahmm)]

# Subsetting peaks
dt.epigrahmm.top = dt.epigrahmm[order(-rowVars(as.matrix(dt.epigrahmm[,c('Summary_Helas3','Summary_Hepg2','Summary_Huvec')]))),][1:round(.N*top),]

# Getting the mode of Pattern per differntial gene

dt.epigrahmm.top$Pattern <- as.factor(pat.peaks(x = pp.epigrahmm,
                                                gr.peaks = gr.epigrahmm[dt.epigrahmm.top$Peak],
                                                gr.genome = gr.counts))

# Preparing the data to plot
group <- colData(epigraHMM_object)$condition
ref <- NULL
# Creating reference table
ref[paste0('ChIP',seq_len(length(unique(group))))] <- list(NULL)
for(i in names(ref)){ref[[i]] <- c(0,1)}
ref <- data.table::as.data.table(expand.grid(ref))
ref <- ref[order(rowSums(ref)),]
map <- data.frame(group = unique(group),col = colnames(ref))

labels <- NULL
for(i in which(!rowSums(ref)==0 & rowSums(ref)<length(unique(group)))){
  labels <- c(labels,paste(map$group[match(colnames(ref)[which(ref[i,]==1)],map$col)],collapse = ' & '))
}

dt.epigrahmm.h3k27ac <- dt.epigrahmm.top[!is.na(Pattern),]
setnames(dt.epigrahmm.h3k27ac,c('Peak','Helas3','Hepg2','Huvec','Width','Pattern'))
dt.epigrahmm.h3k27ac$Pattern %<>% plyr::mapvalues(from = 1:6,to = labels)

save(dt.epigrahmm.h3k27ac,file = 'Figure_Ternary_H3K27ac.RData')