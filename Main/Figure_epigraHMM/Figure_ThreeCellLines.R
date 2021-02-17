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
library(SummarizedExperiment)

bp = 500
chromosome = 'chr7'
idx.genome = c(130720310, 133063299)
ngroups = 3
mark = 'H3K27me3'
data = 'Encode_threecells'
size = 2
fdr = 0.05

# Loading data
load(
  paste0(
    '../../Data/Encode_helas3/H3K27me3/wgEncodeBroadHistoneHelas3H3k27me3StdAlnRep1.markdup.q10.sorted.RData'
  )
)
helas31 = subset(counts[[paste0(bp)]])
rm(counts)
load(
  paste0(
    '../../Data/Encode_helas3/H3K27me3/wgEncodeBroadHistoneHelas3H3k27me3StdAlnRep2.markdup.q10.sorted.RData'
  )
)
helas32 = subset(counts[[paste0(bp)]])
rm(counts)

load(
  paste0(
    '../../Data/Encode_hepg2/H3K27me3/wgEncodeBroadHistoneHepg2H3k27me3StdAlnRep1.markdup.q10.sorted.RData'
  )
)
hepg21 = subset(counts[[paste0(bp)]])
rm(counts)
load(
  paste0(
    '../../Data/Encode_hepg2/H3K27me3/wgEncodeBroadHistoneHepg2H3k27me3StdAlnRep2.markdup.q10.sorted.RData'
  )
)
hepg22 = subset(counts[[paste0(bp)]])
rm(counts)

load(
  paste0(
    '../../Data/Encode_huvec/H3K27me3/wgEncodeBroadHistoneHuvecH3k27me3StdAlnRep1.markdup.q10.sorted.RData'
  )
)
huvec1 = subset(counts[[paste0(bp)]])
rm(counts)
load(
  paste0(
    '../../Data/Encode_huvec/H3K27me3/wgEncodeBroadHistoneHuvecH3k27me3StdAlnRep2.markdup.q10.sorted.RData'
  )
)
huvec2 = subset(counts[[paste0(bp)]])

# Genomic ranges
counts = subset(counts[[paste0(bp)]], chr == chromosome, select = c('chr', 'start', 'stop'))
gr.counts = with(counts, GRanges(chr, IRanges(start = start, end = stop)))

# Combining data
ChIP = cbind(counts, Window = 1:nrow(counts))
ChIP = cbind(ChIP, Helas3_1 = subset(helas31, chr == chromosome)$counts)
ChIP = cbind(ChIP, Helas3_2 = subset(helas32, chr == chromosome)$counts)
ChIP = cbind(ChIP, Hepg2_1 = subset(hepg21, chr == chromosome)$counts)
ChIP = cbind(ChIP, Hepg2_2 = subset(hepg22, chr == chromosome)$counts)
ChIP = cbind(ChIP, Huvec_1 = subset(huvec1, chr == chromosome)$counts)
ChIP = cbind(ChIP, Huvec_2 = subset(huvec2, chr == chromosome)$counts)
ChIP = as.data.table(ChIP)

ChIP.se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = as.matrix(ChIP[, c('Helas3_1',
                                            'Helas3_2',
                                            'Hepg2_1',
                                            'Hepg2_2',
                                            'Huvec_1',
                                            'Huvec_2')])),
  rowRanges = GRanges(ChIP$chr, IRanges::IRanges(ChIP$start, ChIP$stop)),
  colData = data.frame(id = c(
    c(
      'Helas3_1',
      'Helas3_2',
      'Hepg2_1',
      'Hepg2_2',
      'Huvec_1',
      'Huvec_2'
    )
  ))
)

ChIP.se <-
  epigraHMM::normalizeCounts(ChIP.se, epigraHMM::controlEM(), span = 1)

ChIP[, (paste0(
  c(
    'Helas3_1',
    'Helas3_2',
    'Hepg2_1',
    'Hepg2_2',
    'Huvec_1',
    'Huvec_2'
  ),
  '.adj'
)) := .SD / exp(assay(ChIP.se, 'offsets')),
.SDcols = c('Helas3_1',
            'Helas3_2',
            'Hepg2_1',
            'Hepg2_2',
            'Huvec_1',
            'Huvec_2')]
ChIP[, Helas3 := rowSums(.SD), .SDcols = c('Helas3_1.adj', 'Helas3_2.adj')]
ChIP[, Hepg2 := rowSums(.SD), .SDcols = c('Hepg2_1.adj', 'Hepg2_2.adj')]
ChIP[, Huvec := rowSums(.SD), .SDcols = c('Huvec_1.adj', 'Huvec_2.adj')]

ChIP <-
  ChIP[, c('chr', 'start', 'stop', 'Window', 'Helas3', 'Hepg2', 'Huvec')] %>%
  as_tibble() %>%
  gather(Group, Counts, Helas3:Huvec)
ChIP$Mark = mark

# Genomic Ranges to Plot
idx = which.min(abs(counts$start - idx.genome[1])):which.min(abs(counts$start -
                                                                   idx.genome[2]))

# Color of peak calls
methods <-
  c(
    'Genes',
    "ChIPComp + MACS2",
    "csaw",
    "DiffBind + MACS2",
    "diffReps",
    "epigraHMM",
    'RSEG',
    'THOR'
  )
colors = c('#000000', Polychrome::kelly.colors(22)[c('red',
                                                     'yellow',
                                                     'purplishpink',
                                                     'yellowgreen',
                                                     'lightblue',
                                                     'buff',
                                                     'orange')])
names(colors) <- methods

# Position of peak calls (e.g, 1.2 means that the peak call will be placed 1.2 times the maximum observed read count of the plotted region)
peakpos = c(1.05, 1.15, 1.25, 1.35, NA, NA, NA, NA)
names(peakpos) = methods

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

csaw = cbind(
  counts,
  Window = 1:nrow(counts),
  Group = 'Helas3',
  Method = 'csaw'
)
csaw = cbind(csaw, Output = 1 * overlapsAny(gr.counts, gr.csaw))
csaw$Counts = ifelse(csaw$Output == 1, max(ChIP[ChIP$Window %in% idx, 'Counts']) *
                       peakpos['csaw'], NA)
csaw <- csaw %>% as_tibble()
csaw$Output %<>% mapvalues(from = 0:1,
                           to = c('Non-differential', 'Differential'))

# Loading ChIPComp
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

# Loading refseq genes
load('../../Public/Salmon/ENCODE.rnaseq.scaled.RData')
ENCODE.rnaseq.scaled <-
  ENCODE.rnaseq.scaled[, colData(ENCODE.rnaseq.scaled)$Cells %in% c('Helas3', 'Hepg2', 'Huvec')]
ENCODE.rnaseq.scaled$Cells <- droplevels(ENCODE.rnaseq.scaled$Cells)
ENCODE.rnaseq.scaled <-
  ENCODE.rnaseq.scaled[rowRanges(ENCODE.rnaseq.scaled)$gene_biotype == 'protein_coding']
ENCODE.rnaseq.scaled <-
  ENCODE.rnaseq.scaled[overlapsAny(rowRanges(ENCODE.rnaseq.scaled),
                                   GRanges(chromosome, IRanges(start = idx.genome[1], end = idx.genome[2])))]
ENCODE.rnaseq.scaled <-
  ENCODE.rnaseq.scaled[end(ENCODE.rnaseq.scaled) < idx.genome[2] &
                         start(ENCODE.rnaseq.scaled) > idx.genome[1]]

refseq.out = cbind(
  counts,
  Window = 1:nrow(counts),
  Group = 'Helas3',
  Method = 'Genes'
)
refseq.out = cbind(refseq.out,
                   Output = 1 * overlapsAny(gr.counts, ENCODE.rnaseq.scaled))
refseq.out$Counts = ifelse(refseq.out$Output == 1, max(ChIP[ChIP$Window %in%
                                                              idx, 'Counts']) * peakpos['Genes'], NA)
refseq.out <- refseq.out %>% as.tbl()
refseq.out$Output %<>% mapvalues(from = 0:1,
                                 to = c('Non-differential', 'Differential'))

# Organizing the data
dt.segment = rbindlist(list(
  as.data.table(ChIPComp),
  as.data.table(csaw),
  as.data.table(DiffBind),
  as.data.table(refseq.out)
))
dt.segment$Method %<>% factor(levels = c('ChIPComp + MACS2', 'csaw', 'DiffBind + MACS2', 'Genes'))

## Plotting Peak calls

fig_example <-
  ggplot(data = ChIP[ChIP$Window %in% idx, ], aes(x = start, y = Counts)) +
  facet_grid(rows = vars(Group)) +
  geom_segment(
    x = 130785000,
    xend = 130785000,
    y = 0,
    yend = 42.5,
    color = 'grey',
    linetype = 2
  ) +
  geom_segment(
    x = 131015000,
    xend = 131015000,
    y = 0,
    yend = 42.5,
    color = 'grey',
    linetype = 2
  ) +
  geom_segment(
    x = 131200000,
    xend = 131200000,
    y = 0,
    yend = 42.5,
    color = 'grey',
    linetype = 3
  ) +
  geom_segment(
    x = 131350000,
    xend = 131350000,
    y = 0,
    yend = 42.5,
    color = 'grey',
    linetype = 3
  ) +
  geom_segment(
    x = 132250000,
    xend = 132250000,
    y = 0,
    yend = 42.5,
    color = 'grey',
    linetype = 4
  ) +
  geom_segment(
    x = 131795000,
    xend = 131795000,
    y = 0,
    yend = 42.5,
    color = 'grey',
    linetype = 4
  ) +
  annotate(
    'rect',
    alpha = 0.10,
    xmin = 130785000,
    xmax = 131015000,
    ymin = 0,
    ymax = 42.5
  ) +
  annotate(
    'rect',
    alpha = 0.10,
    xmin = 131200000,
    xmax = 131350000,
    ymin = 0,
    ymax = 42.5
  ) +
  annotate(
    'rect',
    alpha = 0.10,
    xmin = 132250000,
    xmax = 131795000,
    ymin = 0,
    ymax = 42.5
  ) +
  geom_line() +
  geom_segment(
    data = dt.segment[Window %in% idx, ],
    aes(
      x = start,
      xend = stop,
      y = Counts,
      yend = Counts,
      color = Method
    ),
    size = size
  ) +
  scale_x_continuous(
    limits = range(ChIP[ChIP$Window %in% idx, 'start']),
    labels = scales::comma,
    position = 'top'
  ) +
  scale_color_manual(values = colors) +
  theme_bw() + xlab(paste0('Genomic Window (', chromosome, ')')) + ylab('Normalized ChIP-seq Counts') +
  guides(col = guide_legend(nrow = 1)) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = 'top',
    legend.direction = 'horizontal'
  )

save(fig_example,file = 'Figure_ThreeCellLines.RData')
ggsave(fig_example,filename = 'Figure_ThreeCellLines.pdf',height = 4.5,width = 9,dpi = 'retina')
