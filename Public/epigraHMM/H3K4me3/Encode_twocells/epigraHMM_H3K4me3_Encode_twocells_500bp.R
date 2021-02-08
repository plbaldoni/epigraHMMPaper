rm(list = ls())
library(epigraHMM)
library(data.table)
library(GenomicRanges)
library(rtracklayer)

bin <- 500
condition <- rep(c('helas3', 'hepg2'), each = 2)
replicate <- rep(c(1, 2), times = 2)
rdata <-
  c(
    'wgEncodeBroadHistoneHelas3H3k4me3StdAlnRep1.markdup.q10.sorted.RData',
    'wgEncodeBroadHistoneHelas3H3k4me3StdAlnRep2.markdup.q10.sorted.RData',
    'wgEncodeBroadHistoneHepg2H3k4me3StdAlnRep1.markdup.q10.sorted.RData',
    'wgEncodeBroadHistoneHepg2H3k4me3StdAlnRep2.markdup.q10.sorted.RData'
  )
mark <- 'H3K4me3'
method <- 'epigraHMM'


# Fixed parameters
type <- ifelse(
  length(unique(condition)) == 2,
  'Encode_twocells',
  ifelse(length(unique(condition)) == 3, 'Encode_threecells', 'Encode_fourcells')
)
chromosome <- paste0('chr', 1:22)

# Organizing input data

for (idx in seq_len(length(rdata))) {
  load(file.path(
    '../../../../Data',
    paste0('Encode_', condition[idx]),
    mark,
    rdata[idx]
  ))
  dat <- data.table(counts[[paste0(bin)]])
  setnames(dat, c(
    'chr',
    'start',
    'stop',
    paste(condition[idx], replicate[idx], sep = '.')
  ))
  rm(counts)
  if (idx == 1) {
    colData <- data.table(dat)
  } else{
    colData <-
      merge(
        colData,
        data.table(dat),
        sort = FALSE,
        all.x = TRUE,
        by = c('chr', 'start', 'stop')
      )
  }
  rm(dat)
}

colData <- colData[chr %in% chromosome,]
rowRanges <-
  makeGRangesFromDataFrame(
    colData,
    keep.extra.columns = FALSE,
    seqnames.field = 'chr',
    start.field = 'start',
    end.field = 'stop'
  )

# Output directory

diroutput = './Output/'
dir.create(diroutput)

# Error function for output
saveerror = function(name, text) {
  fileConn <- file(paste0(diroutput, name))
  writeLines(c(text), fileConn)
  close(fileConn)
}

# epigraHMM
cptime <- list()
cptime[[paste(method, mark, type, 'Time', paste0(bin, 'bp'), sep = '_')]] <-
  microbenchmark::microbenchmark({
    # Creating object
    out.epigraHMM <- epigraHMM::epigraHMMDataSetFromMatrix(
      countData = as.matrix(colData[, paste(condition, replicate, sep = '.'), with = FALSE]),
      colData = data.frame(condition = condition, replicate = replicate),
      rowRanges = rowRanges
    )
    
    # Removing unecessary objects
    rm(colData, rowRanges)
    
    # Normalization
    out.epigraHMM <-
      epigraHMM::normalizeCounts(out.epigraHMM, control = controlEM(), span = 1)
    
    # Initialization
    out.epigraHMM <-
      epigraHMM::initializer(object = out.epigraHMM,
                             control = epigraHMM::controlEM(maxIterEM = 15))
    
    assign(
      paste(method, mark, type, 'Output', paste0(bin, 'bp'), sep = '_'),
      epigraHMM::epigraHMM(
        object = out.epigraHMM,
        type = 'differential',
        dist = 'nb',
        control = epigraHMM::controlEM(
          tempDir = normalizePath(diroutput),
          fileName = paste(method, mark, type, 'Output', paste0(bin, 'bp'), sep = '_')
        )
      )
    )
  }, times = 1)

### Saving computing time
save(cptime,file = paste0(diroutput, paste(
  method, mark, type, 'Time', paste0(bin, 'bp', '.RData'), sep = '_'
)))

### Saving output
do.call(save, list(
  paste(method, mark, type, 'Output', paste0(bin, 'bp'), sep = '_'),
  file = paste0(diroutput, paste(
    method, mark, type, 'Output', paste0(bin, 'bp', '.RData'), sep = '_'
  ))
))

### Saving bed file

tobed <-
  get(paste(method, mark, type, 'Output', paste0(bin, 'bp'), sep = '_'))
tobed <- callPeaks(
  object = tobed,
  hdf5 = file.path(diroutput, paste(
    method, mark, type, 'Output', paste0(bin, 'bp.h5'), sep = '_'
  )),
  method = 'viterbi',
  saveToFile = FALSE
)

rtracklayer::export.bed(object = tobed, con = file.path(diroutput, paste(
  method, mark, type, 'Bed', paste0(bin, 'bp.bed'), sep = '_'
)))
