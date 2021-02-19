library(magrittr)
library(data.table)
library(tidyverse)
library(tibble)
library(RColorBrewer)
library(GenomicRanges)
library(SummarizedExperiment)
library(vsn)
library(DESeq2)
library(ggrepel)
library(Polychrome)
library(ggrepel)

# General parameters
bp = 250
sizeanno <- 2.5
chromosome = paste0('chr', 1:22)
mark = 'H3K36me3'
data = 'Encode_twocells'
cutoff = 2
methods <-
  c(
    "ChIPComp + MACS2",
    "csaw",
    "DiffBind + MACS2",
    "diffReps",
    "epigraHMM",
    'RSEG',
    'THOR',
    'Genes'
  )
colors = c(kelly.colors(22)[c('red',
                              'yellow',
                              'purplishpink',
                              'yellowgreen',
                              'lightblue',
                              'buff',
                              'orange')], 'black')
names(colors) <- methods

getcov = function(gr.predicted,
                  gr.genome,
                  dt,
                  threshold,
                  method,
                  postprob = F) {
  out = list()
  for (x in 1:length(threshold)) {
    if (postprob) {
      gr.predicted.subset <-
        GenomicRanges::reduce(gr.predicted[epigraHMM:::fdrControl(gr.predicted$FDR, fdr = threshold[x])], ignore.strand =
                                T) # Filtering significant windows
    } else{
      gr.predicted.subset = GenomicRanges::reduce(gr.predicted[gr.predicted$FDR <=
                                                                 threshold[x]], ignore.strand = T)
    }
    
    dt[, pred := overlapsAny(gr.genome, gr.predicted.subset)]
    
    out[[x]] = data.table(
      Method = method,
      Threshold = threshold[x],
      TP = sum(dt[, DE == T & pred == T]),
      Median_Size = median(width(gr.predicted.subset)),
      Mean_Size = mean(width(gr.predicted.subset)),
      TPR = sum(dt[, DE == T &
                     pred == T]) / sum(dt[, DE == T]),
      FPR = sum(dt[, DE == F &
                     pred == T]) / sum(dt[, DE == F]),
      PPV = sum(dt[, DE == T &
                     pred == T]) / sum(dt[, pred == T])
    )
  }
  dt[, pred := NULL]
  out = rbindlist(out)
  return(out)
}

# Loading data
load(
  paste0(
    '../../Data/Encode_helas3/',
    mark,
    '/wgEncodeBroadHistoneHelas3H3k36me3StdAlnRep1.markdup.q10.sorted.RData'
  )
)
helas31 = subset(counts[[paste0(bp)]], chr %in% chromosome)
rm(counts)
load(
  paste0(
    '../../Data/Encode_helas3/',
    mark,
    '/wgEncodeBroadHistoneHelas3H3k36me3StdAlnRep2.markdup.q10.sorted.RData'
  )
)
helas32 = subset(counts[[paste0(bp)]], chr %in% chromosome)
rm(counts)

load(
  paste0(
    '../../Data/Encode_hepg2/',
    mark,
    '/wgEncodeBroadHistoneHepg2H3k36me3StdAlnRep1.markdup.q10.sorted.RData'
  )
)
hepg21 = subset(counts[[paste0(bp)]], chr %in% chromosome)
rm(counts)
load(
  paste0(
    '../../Data/Encode_hepg2/',
    mark,
    '/wgEncodeBroadHistoneHepg2H3k36me3StdAlnRep2.markdup.q10.sorted.RData'
  )
)
hepg22 = subset(counts[[paste0(bp)]], chr %in% chromosome)

# Setting up data
counts = subset(counts[[paste0(bp)]], chr %in% chromosome, select = c('chr', 'start', 'stop'))
gr.counts = with(counts, GRanges(chr, IRanges(start = start, end = stop)))

dt <- data.table::data.table(
  counts,
  Helas3_1 = helas31$counts,
  Helas3_2 = helas32$counts,
  Hepg2_1 = hepg21$counts,
  Hepg2_2 = hepg22$counts
)

# Loading genes
load('../../Public/Salmon/ENCODE.chipseq.RData')
ENCODE.chipseq <-
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = assay(ENCODE.chipseq, mark)),
    rowRanges = rowRanges(ENCODE.chipseq),
    colData = colData(ENCODE.chipseq)
  )
ENCODE.chipseq <-
  ENCODE.chipseq[, colData(ENCODE.chipseq)$Cells %in% c('Helas3', 'Hepg2')]
ENCODE.chipseq$Cells <- droplevels(ENCODE.chipseq$Cells)
ENCODE.chipseq <-
  ENCODE.chipseq[rowRanges(ENCODE.chipseq)$gene_biotype == 'protein_coding']
ENCODE.chipseq <-
  ENCODE.chipseq[overlapsAny(rowRanges(ENCODE.chipseq), gr.counts)]

# Normalizing for sequencing depth

ENCODE.chipseq <-
  epigraHMM::normalizeCounts(ENCODE.chipseq, epigraHMM::controlEM(), span = 1)
mat <-
  assay(ENCODE.chipseq, 'counts') / exp(assay(ENCODE.chipseq, 'offsets'))

# Defining DE genes based on ChIP-seq reads
dt.genes <-
  data.table(do.call(cbind, lapply(
    X = unique(ENCODE.chipseq$Cells),
    FUN = function(x) {
      apply(mat[, ENCODE.chipseq$Cells == x], 1, sum)
    }
  )))
setnames(dt.genes, as.character(unique(ENCODE.chipseq$Cells)))
dt.genes[, rowSum := ((Helas3 + Hepg2) > quantile(Helas3 + Hepg2, probs =
                                                    0.25))] #Excluding 25% of genes with low counts
dt.genes[, DE := abs(log2(Hepg2 + 1) - log2(Helas3 + 1)) > cutoff &
           rowSum][, nonDE := !DE]

# Determining differentially transcribed genes
rowRanges(ENCODE.chipseq)$DE <- dt.genes$DE

# # Setting up data
dt[, DE := overlapsAny(gr.counts, rowRanges(ENCODE.chipseq)[rowRanges(ENCODE.chipseq)$DE], ignore.strand =
                         TRUE)][, nonDE :=  !DE]

# Loading peaks: mixHMM
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
gr.epigrahmm$FDR = pp.epigrahmm[, 2]

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
  gr.rseg[seqnames(gr.rseg) %in% chromosome] # Selecting only chromosomes of interest

# # Loading peaks: diffReps
diffreps = fread(list.files(
  file.path('../../Public/diffReps', mark, data, paste0('Output', bp)),
  '.annotated',
  full.names = TRUE
),
header = T)
gr.diffreps = with(diffreps, GenomicRanges::GRanges(Chrom, IRanges(Start, End)))
gr.diffreps$FDR = diffreps$padj #It already gives me adjusted p-values
gr.diffreps <- gr.diffreps[seqnames(gr.diffreps) %in% chromosome]

# Loading peaks: DiffBind
diffbind = fread(list.files(
  file.path('../../Public/DiffBind', mark, data, paste0('Output', bp)),
  '.txt',
  full.names = TRUE
), header = T)
gr.diffbind = with(diffbind, GenomicRanges::GRanges(seqnames, IRanges(start, end)))
gr.diffbind$FDR = diffbind$FDR #It already gives me adjusted p-values
gr.diffbind <- gr.diffbind[seqnames(gr.diffbind) %in% chromosome]

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

# Calculating Metrics
fdr = list()
fdr[['epigraHMM']] = getcov(
  gr.predicted = gr.epigrahmm,
  gr.genome = gr.counts,
  dt = dt,
  threshold = c(0.01, 0.05, 0.10, 0.15, 0.20),
  method = 'epigraHMM',
  postprob = T
)
fdr[['csaw']] = getcov(
  gr.predicted = gr.csaw,
  gr.genome = gr.counts,
  dt = dt,
  threshold = c(0.01, 0.05, 0.10, 0.15, 0.20),
  method = 'csaw'
)
fdr[['thor']] = getcov(
  gr.predicted = gr.thor,
  gr.genome = gr.counts,
  dt = dt,
  threshold = c(0.01, 0.05, 0.10, 0.15, 0.20),
  method = 'THOR'
)
fdr[['rseg']] = getcov(
  gr.predicted = gr.rseg,
  gr.genome = gr.counts,
  dt = dt,
  threshold = c(0.01, 0.05, 0.10, 0.15, 0.20),
  method = 'RSEG',
  postprob = T
)
fdr[['diffreps']] = getcov(
  gr.predicted = gr.diffreps,
  gr.genome = gr.counts,
  dt = dt,
  threshold = c(0.01, 0.05, 0.10, 0.15, 0.20),
  method = 'diffReps'
)
fdr[['diffbind']] = getcov(
  gr.predicted = gr.diffbind,
  gr.genome = gr.counts,
  dt = dt,
  threshold = c(0.01, 0.05, 0.10, 0.15, 0.20),
  method = 'DiffBind + MACS2'
)
fdr[['chipcomp']] = getcov(
  gr.predicted = gr.chipcomp,
  gr.genome = gr.counts,
  dt = dt,
  threshold = c(0.01, 0.05, 0.10, 0.15, 0.20),
  method = 'ChIPComp + MACS2'
)

fdr.summary = rbindlist(fdr)
fdr.summary$Threshold %<>% as.factor()

fdr.summary <- rbindlist(list(
  fdr.summary,
  data.table(
    Method = 'Genes',
    Threshold = NA,
    TP = NA,
    Median_Size = NA,
    Mean_Size = NA,
    TPR = NA,
    FPR = NA,
    PPV = NA
  )
))

fdr.summary$Method %<>% factor(
  levels = c(
    'epigraHMM',
    'ChIPComp + MACS2',
    'csaw',
    'DiffBind + MACS2',
    'diffReps',
    'RSEG',
    'THOR',
    'Genes'
  )
)

# Plotting
fig.PR = ggplot(data = fdr.summary, aes(x = TPR, y = PPV)) +
  geom_line(aes(color = Method), size = 1.25) +
  geom_point(aes(
    x = TPR,
    y = PPV,
    fill = Method,
    shape = Threshold
  ), size = 2.25) +
  scale_shape_manual(
    values = c(
      '0.01' = 8,
      '0.05' = 21,
      '0.1' = 22,
      '0.15' = 23,
      '0.2' = 24
    ),
    na.translate = F
  ) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  guides(fill = F) +
  annotate(
    geom = 'text',
    x = 0.95,
    y = 0.95,
    label = 'Optimal',
    size = sizeanno
  ) +
  annotate(geom = 'point',
           x = 1,
           y = 1,
           size = 3) +
  labs(y = 'Observed Precision (1-FDR)', x = 'Observed Recall (TPR)', shape = 'Nominal FDR') +
  theme_bw()

fig.ROC = ggplot(data = fdr.summary, aes(x = FPR, y = TPR)) +
  geom_line(aes(color = Method), size = 1.25) +
  geom_point(aes(
    x = FPR,
    y = TPR,
    fill = Method,
    shape = Threshold
  ), size = 2.25) +
  scale_shape_manual(
    values = c(
      '0.01' = 8,
      '0.05' = 21,
      '0.1' = 22,
      '0.15' = 23,
      '0.2' = 24
    ),
    na.translate = F
  ) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  guides(fill = F) +
  annotate(
    geom = 'text',
    x = 0.05,
    y = 0.95,
    label = 'Optimal',
    size = sizeanno
  ) +
  annotate(geom = 'point',
           x = 0,
           y = 1,
           size = 3) +
  labs(y = 'Observed Sensitivity (TPR)', x = 'Observed 1 - Specificity (FPR)', shape = 'Nominal FDR') +
  theme_bw()

fig.TPR <-
  ggplot(data = fdr.summary, aes(
    x = Mean_Size / 1000,
    y = TPR,
    label = paste0(format(round(FPR, digits = 2), nsmall = 2))
  )) +
  geom_path(aes(color = Method, group = Method), size = 1.25) +
  geom_point(aes(
    x = Mean_Size / 1000,
    y = TPR,
    fill = Method,
    shape = Threshold
  ),
  size = 2.25) +
  scale_shape_manual(
    values = c(
      '0.01' = 8,
      '0.05' = 21,
      '0.1' = 22,
      '0.15' = 23,
      '0.2' = 24
    ),
    na.translate = F
  ) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  guides(fill = F) +
  labs(y = 'Observed Sensitivity (TPR)', x = 'Average Peak Size (kbp)', shape = 'Nominal FDR') +
  theme_bw() +
  coord_cartesian(xlim = c(-1, 17.5), ylim = c(0, 0.7)) +
  geom_text_repel(
    hjust = 0,
    segment.size = 0.2,
    min.segment.length = 0,
    size = 2,
    force = 15
  ) +
  annotate(
    'text',
    x = 2,
    y = 0.7,
    label = 'Observed FDR shown\nnext to each data point',
    size = sizeanno
  )

save(fig.PR, fig.ROC, fig.TPR, file = 'Figure_PR_H3K36me3_Encode_twocells_250bp.RData')
