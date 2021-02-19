library(data.table)
library(ggplot2)
library(Polychrome)
library(magrittr)
library(forcats)
library(plyr)
library(GenomicRanges)
library(ggpubr)

# This plots the specificity and sensitivity values relevant to a given *.tsv file.

prefix = 'hist'
incoming = paste0('../../Sim/SequencingReads/', prefix, "_result.tsv")
methods <-
  c("epigraHMM",
    "ChIPComp + HOMER",
    "csaw",
    "DiffBind + HOMER",
    "diffReps",
    'RSEG',
    'THOR')
colors = c(kelly.colors(22)[c('lightblue',
                              'red',
                              'yellow',
                              'purplishpink',
                              'yellowgreen',
                              'buff',
                              'orange')])
names(colors) <- methods

cutoff <- 0.20
sizeanno <- 2.5
units = 1e3
heights = 1 + seq(0.15, by = 0.15, length.out = 7)
names(heights) = c(methods)
size = 1.25

# Loading in the data.
ided <- fread(incoming)

# How many times RSEG failed out of 100 simulated data sets?
ided[Method == 'RSEG', .(Error = sum(Error == 1)), by = c('It')][, sum(Error >
                                                                         0)] #26

# Out of the number of times that RSEG did not fail, how many times it called the entire genome as differential?
ided[Method == 'RSEG' &
       Error == 0, ][, sum(CallSize > 1e6), by = 'Cutoff'] #3

# Plotting Precision/Recall curve

summarized.roc.ided <- ided[, .(
  TPR = mean(TPR, na.rm = T),
  FPR = mean(FPR, na.rm = T),
  PPV = mean(PPV, na.rm = T)
), by = c('Method', 'Cutoff')]
summarized.roc.ided$Cutoff %<>% as.character()

figA = ggplot(data = summarized.roc.ided, aes(x = (1 - PPV), y = TPR)) +
  geom_line(aes(color = Method), size = 1.25) +
  geom_point(aes(
    x = (1 - PPV),
    y = TPR,
    fill = Method,
    shape = Cutoff
  ), size = 2.25) +
  scale_shape_manual(values = c(
    '0.01' = 8,
    '0.05' = 21,
    '0.1' = 22,
    '0.15' = 23,
    '0.2' = 24
  )) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  guides(fill = F) +
  annotate(geom = 'point',
           x = 0,
           y = 1,
           size = 3) +
  coord_cartesian(ylim = c(0, 1)) +
  annotate(
    geom = 'text',
    x = 0.1,
    y = 1,
    label = 'Optimal',
    size = sizeanno
  ) +
  labs(x = 'Observed False Discovery Rate', y = 'Observed True Positive Rate', shape = 'Nominal\nFDR') +
  theme_bw()

summarized.size <-
  ided[, .(
    HitsPerPeak = mean(HitsPerPeak, na.rm = T),
    NCallsNPeaks = mean(NCalls / NPeaks, na.rm = T)
  ), by = c('Method', 'Cutoff')]

figB <-
  ggplot(data = summarized.size, aes(x = HitsPerPeak, y = NCallsNPeaks)) +
  geom_line(aes(color = Method), size = 1.25) +
  geom_point(aes(
    x = HitsPerPeak,
    y = NCallsNPeaks,
    fill = Method,
    shape = as.factor(Cutoff)
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
  coord_cartesian(xlim = c(0, 11), ylim = c(-1, 31)) +
  geom_vline(xintercept = 1,
             linetype = 2,
             color = 'grey') +
  geom_hline(yintercept = 1,
             linetype = 2,
             color = 'grey') +
  annotate(
    geom = 'text',
    x = 0.5,
    y = -0.5,
    label = 'Optimal',
    size = sizeanno
  ) +
  labs(y = 'Called/Simulated Peaks', x = 'Hits per peak', shape = 'Nominal FDR') +
  theme_bw()


# Plotting peak calls

gr.true = with(
  read.table('../../Sim/SequencingReads/autosim1/hist1_log.txt', header = T),
  GRanges(chr, IRanges(start, end))
)

gr.tile = with(
  fread(
    '../../Sim/SequencingReads/autosim1/results/epigraHMM_ranges.txt'
  ),
  GRanges(seqnames, IRanges(start, end))
)

dt.tile = data.table(
  chr = 'chrA',
  start = start(gr.tile),
  stop = end(gr.tile),
  window = 1:length(start(gr.tile)),
  counts1 = bamsignals::bamCount(
    bampath = '../../Sim/SequencingReads/autosim1/hist1_out_1.bam',
    gr = gr.tile,
    verbose = T,
    shift = 0
  ),
  counts2 = bamsignals::bamCount(
    bampath = '../../Sim/SequencingReads/autosim1/hist1_out_2.bam',
    gr = gr.tile,
    verbose = T,
    shift = 0
  ),
  counts3 = bamsignals::bamCount(
    bampath = '../../Sim/SequencingReads/autosim1/hist1_out_3.bam',
    gr = gr.tile,
    verbose = T,
    shift = 0
  ),
  counts4 = bamsignals::bamCount(
    bampath = '../../Sim/SequencingReads/autosim1/hist1_out_4.bam',
    gr = gr.tile,
    verbose = T,
    shift = 0
  ),
  control1 = bamsignals::bamCount(
    bampath = '../../Sim/SequencingReads/autosim1/control1_out_1.bam',
    gr = gr.tile,
    verbose = T,
    shift = 0
  ),
  control2 = bamsignals::bamCount(
    bampath = '../../Sim/SequencingReads/autosim1/control1_out_2.bam',
    gr = gr.tile,
    verbose = T,
    shift = 0
  ),
  control3 = bamsignals::bamCount(
    bampath = '../../Sim/SequencingReads/autosim1/control1_out_3.bam',
    gr = gr.tile,
    verbose = T,
    shift = 0
  ),
  control4 = bamsignals::bamCount(
    bampath = '../../Sim/SequencingReads/autosim1/control1_out_4.bam',
    gr = gr.tile,
    verbose = T,
    shift = 0
  )
)

dt.tile[, Simulated := overlapsAny(gr.tile, gr.true)]

dt.tile[, countsA := counts1 + counts2]
dt.tile[, countsB := counts3 + counts4]

dt.melted.tile <-
  data.table::melt(
    dt.tile,
    id.vars = c('start', 'stop', 'window'),
    measure.vars = c('countsA', 'countsB'),
    variable.name = 'condition',
    value.name = 'counts'
  )
dt.melted.tile$condition %<>% mapvalues(from = c('countsA', 'countsB'),
                                        to = c('Condition A', 'Condition B'))

idx = (dt.tile$window[which.min(abs(124800000 - dt.tile$start))] - 0):(dt.tile$window[which.min(abs(124950000 -
                                                                                                      dt.tile$stop))] + 0)

# Loading peaks
gr.ChIPComp <-
  GenomicRanges::reduce(with(
    fread(
      '../../Sim/SequencingReads/autosim1/results/ChIPComp_table.txt'
    )[FDR <= cutoff, ],
    GRanges(chr, IRanges(start, end))
  ))
ChIPComp = cbind(
  dt.tile[, c('chr', 'start', 'stop', 'window')],
  Group = 'Condition A',
  Method = 'ChIPComp + HOMER',
  Output = 1 * overlapsAny(gr.tile, gr.ChIPComp)
)
ChIPComp[, Counts := ifelse(ChIPComp$Output == 1, max(dt.tile[window %in% idx, countsA]) *
                              heights['ChIPComp + HOMER'], NA)]
ChIPComp$Output %<>% mapvalues(from = 0:1,
                               to = c('Non-differential', 'Differential'))

gr.csaw <-
  GenomicRanges::reduce(with(
    fread('../../Sim/SequencingReads/autosim1/results/csaw_ranges.txt')[fread('../../Sim/SequencingReads/autosim1/results/csaw_table.txt')$FDR <=
                                                                       cutoff, ],
    GRanges(seqnames, IRanges(start, end))
  ))
csaw = cbind(
  dt.tile[, c('chr', 'start', 'stop', 'window')],
  Group = 'Condition A',
  Method = 'csaw',
  Output = 1 * overlapsAny(gr.tile, gr.csaw)
)
csaw[, Counts := ifelse(csaw$Output == 1, max(dt.tile[window %in% idx, countsA]) *
                          heights['csaw'], NA)]
csaw$Output %<>% mapvalues(from = 0:1,
                           to = c('Non-differential', 'Differential'))

gr.DiffBind <-
  GenomicRanges::reduce(with(
    fread(
      '../../Sim/SequencingReads/autosim1/results/DiffBind_ranges.txt'
    )[FDR <= cutoff, ],
    GRanges(seqnames, IRanges(start, end))
  ))
DiffBind = cbind(
  dt.tile[, c('chr', 'start', 'stop', 'window')],
  Group = 'Condition A',
  Method = 'DiffBind + HOMER',
  Output = 1 * overlapsAny(gr.tile, gr.DiffBind)
)
DiffBind[, Counts := ifelse(DiffBind$Output == 1, max(dt.tile[window %in% idx, countsA]) *
                              heights['DiffBind + HOMER'], NA)]
DiffBind$Output %<>% mapvalues(from = 0:1,
                               to = c('Non-differential', 'Differential'))

gr.diffReps <-
  GenomicRanges::reduce(with(
    fread(
      '../../Sim/SequencingReads/autosim1/results/diffReps_table.txt'
    )[FDR <= cutoff, ],
    GRanges(Chrom, IRanges(Start, End))
  ))
diffReps = cbind(
  dt.tile[, c('chr', 'start', 'stop', 'window')],
  Group = 'Condition A',
  Method = 'diffReps',
  Output = 1 * overlapsAny(gr.tile, gr.diffReps)
)
diffReps[, Counts := ifelse(diffReps$Output == 1, max(dt.tile[window %in% idx, countsA]) *
                              heights['diffReps'], NA)]
diffReps$Output %<>% mapvalues(from = 0:1,
                               to = c('Non-differential', 'Differential'))

gr.epigraHMM <-
  GenomicRanges::reduce(with(
    fread(
      '../../Sim/SequencingReads/autosim1/results/epigraHMM_table.txt'
    )[epigraHMM:::fdrControl(prob = Postprob, fdr = cutoff), ],
    GRanges(seqnames, IRanges(start, end))
  ))
epigraHMM = cbind(
  dt.tile[, c('chr', 'start', 'stop', 'window')],
  Group = 'Condition A',
  Method = 'epigraHMM',
  Output = 1 * overlapsAny(gr.tile, gr.epigraHMM)
)
epigraHMM[, Counts := ifelse(epigraHMM$Output == 1, max(dt.tile[window %in%
                                                                  idx, countsA]) * heights['epigraHMM'], NA)]
epigraHMM$Output %<>% mapvalues(from = 0:1,
                                to = c('Non-differential', 'Differential'))

gr.RSEG <-
  GenomicRanges::reduce(with(
    fread('../../Sim/SequencingReads/autosim1/results/RSEG_table.txt')[epigraHMM:::fdrControl(prob = Postprob, fdr = cutoff), ],
    GRanges(V1, IRanges(V2, V3))
  ))
RSEG = cbind(
  dt.tile[, c('chr', 'start', 'stop', 'window')],
  Group = 'Condition A',
  Method = 'RSEG',
  Output = 1 * overlapsAny(gr.tile, gr.RSEG)
)
RSEG[, Counts := ifelse(RSEG$Output == 1, max(dt.tile[window %in% idx, countsA]) *
                          heights['RSEG'], NA)]
RSEG$Output %<>% mapvalues(from = 0:1,
                           to = c('Non-differential', 'Differential'))

gr.THOR <-
  GenomicRanges::reduce(with(
    fread('../../Sim/SequencingReads/autosim1/results/THOR_table.txt')[FDR <= cutoff, ],
    GRanges(V1, IRanges(V2, V3))
  ))
THOR = cbind(
  dt.tile[, c('chr', 'start', 'stop', 'window')],
  Group = 'Condition A',
  Method = 'THOR',
  Output = 1 * overlapsAny(gr.tile, gr.THOR)
)
THOR[, Counts := ifelse(THOR$Output == 1, max(dt.tile[window %in% idx, countsA]) *
                          heights['THOR'], NA)]
THOR$Output %<>% mapvalues(from = 0:1,
                           to = c('Non-differential', 'Differential'))

dt.peaks = rbindlist(list(ChIPComp[window %in% idx, ],
                          csaw[window %in% idx, ],
                          DiffBind[window %in% idx, ],
                          diffReps[window %in% idx, ],
                          epigraHMM[window %in% idx, ],
                          RSEG[window %in% idx, ],
                          THOR[window %in% idx, ]))

dt.peaks$Track = dt.peaks$Counts

dt.melted.tile$Group = dt.melted.tile$condition
dt.melted.tile$ChIP = dt.melted.tile$counts

dt.all = merge(
  dt.melted.tile[window %in% idx, c('start', 'stop', 'Group', 'ChIP')],
  dt.peaks[, c('start', 'stop', 'Group', 'Method', 'Track')],
  by = c('start', 'stop', 'Group'),
  all.x = T
)

figD = ggplot(data = dt.all, aes(x = start, y = ChIP)) +
  facet_wrap( ~ Group, nrow = 2, ncol = 1) +
  annotate(
    'rect',
    alpha = 0.10,
    xmin = 124852109 - 0.1 * abs(diff(c(
      124862108 , 124852109
    ))),
    xmax = 124862108 + 0.1 * abs(diff(c(
      124862108 , 124852109
    ))),
    ymin = 0,
    ymax = 1.05 * max(dt.all$Track, na.rm = T)
  ) +
  annotate(
    'rect',
    alpha = 0.10,
    xmin = 124905109 - 0.05 * abs(diff(c(
      124905109, 124914708
    ))),
    xmax = 124914708 + 0.05 * abs(diff(c(
      124905109, 124914708
    ))),
    ymin = 0,
    ymax = 1.05 * max(dt.all$Track, na.rm = T)
  ) +
  geom_line() +
  geom_segment(aes(
    x = start,
    xend = stop,
    y = Track,
    yend = Track,
    color = Method
  ),
  size = size) +
  scale_color_manual(values = colors, breaks = names(colors)) +
  scale_x_continuous(breaks = c(124825000,124925000), labels = scales::comma) +
  labs(y = 'Total ChIP counts', x = 'Genomic Position') +
  theme_bw()

# Plotting computing time

ided.time <-
  ided[Cutoff == cutoff, ] #It does not matter which cutoff we pick, we ran each method only once per simulated data sets anyway
ided.time$Method %<>% mapvalues(
  from = c('DiffBind + HOMER', 'ChIPComp + HOMER'),
  to = c('DiffBind +\nHOMER', 'ChIPComp +\nHOMER')
)

newcolors <- colors
names(newcolors) %<>% mapvalues(
  from = c('DiffBind + HOMER', 'ChIPComp + HOMER'),
  to = c('DiffBind +\nHOMER', 'ChIPComp +\nHOMER')
)

figC <- ggplot(data = ided.time, aes(x = Method, y = Time, fill = Method)) +
  theme_bw() +
  labs(y = 'Time (minutes)', x = 'Method') +
  geom_jitter(width = 0.2, alpha = 0.25) +
  geom_boxplot(outlier.shape = NA,alpha = 0.75) + 
  guides(x = guide_axis(angle = 30)) +
  scale_fill_manual(values = newcolors)

# Creating figure

fig <-
  ggarrange(
    figA + theme(
      legend.position = 'top',
      legend.direction = 'horizontal',
      legend.box = 'vertical',
      legend.spacing.y = unit(-0.25, "cm")
    ) + guides(
      col = guide_legend(nrow = 1, order = 1),
      shape = guide_legend(
        nrow = 1,
        order = 2,
        title = 'Nominal FDR'
      )
    ),
    figB + theme(
      legend.position = 'top',
      legend.direction = 'horizontal',
      legend.box = 'vertical',
      legend.spacing.y = unit(-0.25, "cm")
    ) + guides(
      col = guide_legend(nrow = 1, order = 1),
      shape = guide_legend(
        nrow = 1,
        order = 2,
        title = 'Nominal FDR'
      )
    ),
    figC + theme(
      legend.position = 'top',
      legend.direction = 'horizontal',
      legend.box = 'vertical',
      legend.spacing.y = unit(-0.25, "cm")
    ) + guides(
      col = guide_legend(nrow = 1, order = 1),
      shape = guide_legend(
        nrow = 1,
        order = 2,
        title = 'Nominal FDR'
      )
    ),
    figD + theme(
      legend.position = 'top',
      legend.direction = 'horizontal',
      legend.box = 'vertical',
      legend.spacing.y = unit(-0.25, "cm")
    ) + guides(
      col = guide_legend(nrow = 1, order = 1),
      shape = guide_legend(
        nrow = 1,
        order = 2,
        title = 'Nominal FDR'
      )
    ),
    nrow = 2,
    ncol = 2,
    legend = 'top',
    common.legend = T,
    labels = list('A', 'B', 'C', 'D')
  )

# Render
ggsave(
  plot = fig,
  filename = paste0('Figure_Simulation_SequencingReads_',100*cutoff,'.pdf'),
  height = 8,
  width = 8,
  dpi = 'retina'
)
