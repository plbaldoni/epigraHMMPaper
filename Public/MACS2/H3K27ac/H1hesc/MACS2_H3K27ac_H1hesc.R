library(GenomicRanges)
library(microbenchmark)

method = 'MACS2'
mark = 'H3K27ac'
cell = 'H1hesc'

chip = c(
  'wgEncodeBroadHistoneH1hescH3k27acStdAlnRep1.markdup.q10.sorted.bam',
  'wgEncodeBroadHistoneH1hescH3k27acStdAlnRep1.markdup.q10.sorted.bam'
)

control = c(
  'wgEncodeBroadHistoneH1hescControlStdAlnRep1.markdup.q10.sorted.bam',
  'wgEncodeBroadHistoneH1hescControlStdAlnRep2.markdup.q10.sorted.bam'
)

dirdata = file.path('../../../../Data', paste0('Encode_', tolower(cell)))
cptime = list()

### Organizing other parameters

opt = ifelse(
  mark %in% c('H3K36me3', 'H3K27me3', 'EZH2'),
  '--broad -f BAM -g 2.80e+09 -B --broad-cutoff 0.1',
  ifelse(
    mark %in% c('H3K4me3', 'H3K27ac', 'CTCF'),
    '-f BAM -g 2.80e+09 -B -q 0.01',
    NA
  )
)

system(paste(
  'mkdir -p',
  paste0('./temporary_data/chip'),
  paste0('./temporary_data/control')
))
system(paste('mkdir', paste0('./Output/')))

### Copying files

cmd = paste('cp',
            paste0(file.path(dirdata, mark, chip), collapse = ' '),
            paste0('./temporary_data/chip'))
cat('Command: ', cmd, '\n')
system(cmd)
chipfiles = list.files(
  path = paste0('./temporary_data/chip'),
  pattern = '*.bam$',
  full.names = T
)

cmd = paste('cp',
            paste0(file.path(dirdata, 'Control', control), collapse = ' '),
            paste0('./temporary_data/control'))
cat('Command: ', cmd, '\n')
system(cmd)
controlfiles = list.files(
  path = paste0('./temporary_data/control'),
  pattern = '*.bam$',
  full.names = T
)

### Running MACS2
cmd = paste(
  'macs2 callpeak',
  opt,
  '-t',
  paste(chipfiles, collapse = ' '),
  '-c',
  paste(controlfiles, collapse = ' '),
  '--outdir',
  paste0('./Output/'),
  '-n',
  paste(method, mark, cell, 'Output', sep = '_')
)
cat('Command: ', cmd, '\n')

cptime[[paste(method, mark, cell, 'Output', sep = '_')]] = microbenchmark::microbenchmark(system(cmd), times = 1)

### Organizing output
if (mark %in% c('H3K36me3', 'H3K27me3', 'EZH2')) {
  dat = read.table(paste0(
    paste0('./Output/'),
    paste(method, mark, cell, 'Output', sep = '_'),
    '_peaks.broadPeak'
  ))[, c('V1', 'V2', 'V3')]
  colnames(dat) = c('chr', 'start', 'stop')
  write.table(
    dat,
    file = paste0(
      paste0('./Output/'),
      paste(method, mark, cell, 'Output', sep = '_')
    ),
    quote = F,
    row.names = F,
    col.names = T
  )
}
if (mark %in% c('H3K4me3', 'H3K27ac', 'CTCF')) {
  dat = read.table(paste0(
    paste0('./Output/'),
    paste(method, mark, cell, 'Output', sep = '_'),
    '_peaks.narrowPeak'
  ))[, c('V1', 'V2', 'V3')]
  colnames(dat) = c('chr', 'start', 'stop')
  write.table(
    dat,
    file = paste0(
      paste0('./Output/'),
      paste(method, mark, cell, 'Output', sep = '_')
    ),
    quote = F,
    row.names = F,
    col.names = T
  )
}

### Saving computing time
save(cptime, file = paste0(
  paste0('./Output/'),
  paste(method, mark, cell, 'Time.RData', sep = '_')
))

### Removing tmp
system(paste('rm -r',
             paste0('./temporary_data')))