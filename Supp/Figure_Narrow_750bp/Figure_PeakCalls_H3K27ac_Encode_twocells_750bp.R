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

bp = 750
chromosome = 'chr19'
ngroups = 2
mark = 'H3K27ac'
data = 'Encode_twocells'
size = 1.25
fdr = 0.05
methods <-
  c(
    'Genes',
    "epigraHMM",
    "ChIPComp + MACS2",
    "csaw",
    "DiffBind + MACS2",
    "diffReps",
    'RSEG',
    'THOR'
  )
colors = c('black', kelly.colors(22)[c('lightblue',
                                       'red',
                                       'yellow',
                                       'purplishpink',
                                       'yellowgreen',
                                       'buff',
                                       'orange')])
names(colors) <- methods

peakpos = c(1.1, seq(1.2, by = 0.1, length.out = 7))
names(peakpos) <- methods

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

countsraw <- as.data.table(counts[[paste0(bp)]])
countsraw <- countsraw[chr %in% paste0('chr', 1:22),]
gr.countsraw <- countsraw[, GRanges(chr, IRanges(start, stop))]

idx.chr = (counts[[paste0(bp)]]$chr == chromosome)

counts = subset(counts[[paste0(bp)]], chr == chromosome, select = c('chr', 'start', 'stop'))
gr.counts = with(counts, GRanges(chr, IRanges(start = start, end = stop)))

# Combining data
ChIP = cbind(counts, Window = 1:nrow(counts))
ChIP = cbind(ChIP, Helas3_1 = helas31$counts)
ChIP = cbind(ChIP, Helas3_2 = helas32$counts)
ChIP = cbind(ChIP, Hepg2_1 = hepg21$counts)
ChIP = cbind(ChIP, Hepg2_2 = hepg22$counts)
ChIP = as.data.table(ChIP)
ChIP.se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = as.matrix(ChIP[, c('Helas3_1', 'Helas3_2', 'Hepg2_1', 'Hepg2_2')])),
  rowRanges = GRanges(ChIP$chr, IRanges::IRanges(ChIP$start, ChIP$stop)),
  colData = data.frame(id = c(
    'Helas3_1', 'Helas3_2', 'Hepg2_1', 'Hepg2_2'
  ))
)
ChIP.se <-
  epigraHMM::normalizeCounts(ChIP.se, epigraHMM::controlEM(), span = 1)
ChIP[, (paste0(c(
  'Helas3_1', 'Helas3_2', 'Hepg2_1', 'Hepg2_2'
), '.adj')) := .SD / exp(assay(ChIP.se, 'offsets')),
.SDcols = c('Helas3_1', 'Helas3_2', 'Hepg2_1', 'Hepg2_2')]
ChIP[, Helas3 := rowSums(.SD), .SDcols = c('Helas3_1.adj', 'Helas3_2.adj')]
ChIP[, Hepg2 := rowSums(.SD), .SDcols = c('Hepg2_1.adj', 'Hepg2_2.adj')]

ChIP <-
  ChIP[, c('chr', 'start', 'stop', 'Window', 'Helas3', 'Hepg2')] %>%
  as_tibble() %>%
  gather(Group, Counts, Helas3:Hepg2)
ChIP$Mark = mark

# Genomic Ranges to Plot
idx.genome = c(37050000, 37200000)
idx = which.min(abs(counts$start - idx.genome[1])):which.min(abs(counts$start -
                                                                   idx.genome[2]))

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
gr.epigrahmm <-
  GenomicRanges::reduce(gr.epigrahmm[epigraHMM:::fdrControl(pp.epigrahmm[, 2], fdr = fdr)])
gr.epigrahmm <- gr.epigrahmm[seqnames(gr.epigrahmm) %in% chromosome]

epigraHMM = cbind(
  counts,
  Window = 1:nrow(counts),
  Group = 'Helas3',
  Method = 'epigraHMM'
)
epigraHMM = cbind(epigraHMM, Output = 1 * overlapsAny(gr.counts, gr.epigrahmm))
epigraHMM$Counts = ifelse(epigraHMM$Output == 1, max(ChIP[ChIP$Window %in%
                                                            idx, 'Counts']) * peakpos['epigraHMM'], NA)
epigraHMM <- epigraHMM %>% as_tibble()
epigraHMM$Output %<>% mapvalues(from = 0:1,
                                to = c('Non-differential', 'Differential'))

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

Csaw = cbind(
  counts,
  Window = 1:nrow(counts),
  Group = 'Helas3',
  Method = 'csaw'
)
Csaw = cbind(Csaw, Output = 1 * overlapsAny(gr.counts, gr.csaw))
Csaw$Counts = ifelse(Csaw$Output == 1, max(ChIP[ChIP$Window %in% idx, 'Counts']) *
                       peakpos['csaw'], NA)
Csaw <- Csaw %>% as.tbl()
Csaw$Output %<>% mapvalues(from = 0:1,
                           to = c('Non-differential', 'Differential'))

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

RSEG = cbind(
  counts,
  Window = 1:nrow(counts),
  Group = 'Helas3',
  Method = 'RSEG'
)
RSEG = cbind(RSEG, Output = 1 * overlapsAny(gr.counts, gr.rseg))
RSEG$Counts = ifelse(RSEG$Output == 1, max(ChIP[ChIP$Window %in% idx, 'Counts']) *
                       peakpos['RSEG'], NA)
RSEG <- RSEG %>% as.tbl()
RSEG$Output %<>% mapvalues(from = 0:1,
                           to = c('Non-differential', 'Differential'))

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

diffReps = cbind(
  counts,
  Window = 1:nrow(counts),
  Group = 'Helas3',
  Method = 'diffReps'
)
diffReps = cbind(diffReps, Output = 1 * overlapsAny(gr.counts, gr.diffreps))
diffReps$Counts = ifelse(diffReps$Output == 1, max(ChIP[ChIP$Window %in% idx, 'Counts']) *
                           peakpos['diffReps'], NA)
diffReps <- diffReps %>% as.tbl()
diffReps$Output %<>% mapvalues(from = 0:1,
                               to = c('Non-differential', 'Differential'))

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

DiffBind = cbind(
  counts,
  Window = 1:nrow(counts),
  Group = 'Helas3',
  Method = 'DiffBind + MACS2'
)
DiffBind = cbind(DiffBind, Output = 1 * overlapsAny(gr.counts, gr.diffbind))
DiffBind$Counts = ifelse(DiffBind$Output == 1, max(ChIP[ChIP$Window %in% idx, 'Counts']) *
                           peakpos['DiffBind + MACS2'], NA)
DiffBind <- DiffBind %>% as.tbl()
DiffBind$Output %<>% mapvalues(from = 0:1,
                               to = c('Non-differential', 'Differential'))

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

ChIPComp = cbind(
  counts,
  Window = 1:nrow(counts),
  Group = 'Helas3',
  Method = 'ChIPComp + MACS2'
)
ChIPComp = cbind(ChIPComp, Output = 1 * overlapsAny(gr.counts, gr.chipcomp))
ChIPComp$Counts = ifelse(ChIPComp$Output == 1, max(ChIP[ChIP$Window %in% idx, 'Counts']) *
                           peakpos['ChIPComp + MACS2'], NA)
ChIPComp <- ChIPComp %>% as.tbl()
ChIPComp$Output %<>% mapvalues(from = 0:1,
                               to = c('Non-differential', 'Differential'))

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

THOR = cbind(
  counts,
  Window = 1:nrow(counts),
  Group = 'Helas3',
  Method = 'THOR'
)
THOR = cbind(THOR, Output = 1 * overlapsAny(gr.counts, gr.thor))
THOR$Counts = ifelse(THOR$Output == 1, max(ChIP[ChIP$Window %in% idx, 'Counts']) *
                       peakpos['THOR'], NA)
THOR <- THOR %>% as.tbl()
THOR$Output %<>% mapvalues(from = 0:1,
                           to = c('Non-differential', 'Differential'))

# Loading refseq genes
load('../../Public/Salmon//UCSC_RefSeq_Hg19_Genes.RData')
gr.genes <- gr.genes[overlapsAny(gr.genes, gr.counts)]

# Organizing refseq
refseq.out = cbind(
  counts,
  Window = 1:nrow(counts),
  Group = 'Helas3',
  Method = 'Genes'
)
refseq.out = cbind(refseq.out, Output = 1 * overlapsAny(gr.counts, gr.genes))
refseq.out$Counts = ifelse(refseq.out$Output == 1, max(ChIP[ChIP$Window %in%
                                                              idx, 'Counts']) * peakpos['Genes'], NA)
refseq.out <- refseq.out %>% as.tbl()
refseq.out$Output %<>% mapvalues(from = 0:1,
                                 to = c('Non-differential', 'Differential'))

# Organizing the data
dt.segment = rbindlist(
  list(
    as.data.table(epigraHMM),
    as.data.table(Csaw),
    as.data.table(RSEG),
    as.data.table(diffReps),
    as.data.table(DiffBind),
    as.data.table(ChIPComp),
    as.data.table(THOR),
    as.data.table(refseq.out)
  )
)

dt.segment$Method %<>% factor(
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

## Plotting Peak calls

fig1.h3k27ac = ggplot(data = ChIP[ChIP$Window %in% idx,], aes(x = start, y =
                                                                Counts)) +
  geom_line() +
  facet_grid(rows = vars(Group)) +
  annotate(
    'rect',
    alpha = 0.1,
    xmin = 37095129 - 10 * bp,
    xmax = 37097128 + 10 * bp,
    ymin = 0,
    ymax = 1.05 * max(dt.segment$Counts, na.rm = T)
  ) +
  geom_segment(
    data = dt.segment[Window %in% idx,],
    aes(
      x = start,
      xend = stop,
      y = Counts,
      yend = Counts,
      color = Method
    ),
    size = size
  ) +
  scale_color_manual(values = colors) +
  scale_x_continuous(labels = scales::comma) +
  xlab(paste0('Genomic Window (', chromosome, ')')) + ylab('Normalized ChIP Counts') +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  )

## Plotting Posterior Probabilities

PostProb = cbind(
  counts,
  Window = 1:nrow(counts),
  PostProb2 = pp.epigrahmm[which(idx.chr == T), 2],
  Label = 'Prob.'
)
PostProb <- PostProb %>%
  as_tibble() %>%
  gather(Component, PP, PostProb2)
PostProb$Component %<>% mapvalues(from = c('PostProb2'), to = c('Differential'))

fig2.h3k27ac = ggplot(data = PostProb[PostProb$Window %in% idx,], aes(x =
                                                                        start, y = PP, fill = Component)) + facet_grid(rows = vars(Label)) +
  geom_area(position = 'identity', alpha = 1) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 1),
    labels = function(x)
      sprintf("%.1f", x)
  ) +
  scale_x_continuous(limits = range(PostProb[PostProb$Window %in% idx, 'start']), labels = scales::comma) +
  xlab(paste0('Genomic Window (', chromosome, ')')) +
  ylab('Prob.') +
  guides(fill = F) +
  scale_fill_manual(values = c('Differential' = as.character(colors['epigraHMM']))) +
  theme_bw() +
  theme(legend.position = "none",
        strip.text.y = element_text(colour = alpha('grey', 0.0)),panel.grid = element_blank())

fig.h3k27ac = ggarrange(
  fig1.h3k27ac,
  fig2.h3k27ac,
  heights = c(0.8, 0.2),
  ncol = 1,
  nrow = 2,
  legend = 'none'
)

save(fig1.h3k27ac, fig2.h3k27ac, fig.h3k27ac, file = 'Figure_PeakCalls_H3K27ac_Encode_twocells_750bp.RData')
