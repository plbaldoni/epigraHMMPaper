rm(list = ls())
library(epigraHMM)
library(data.table)
library(GenomicRanges)
library(rtracklayer)

pattern <- list(c(1,2),3,2) #k27-ezh2 k36 k27
bin <- 500
cellline <- 'Helas3'
type <- tail(strsplit(getwd(), '/')[[1]], 1)
condition <- rep(c('H3K27me3', 'H3K36me3', 'EZH2'), each = 2)
replicate <- rep(c(1, 2), times = 3)
rdata <-
  c(
    'wgEncodeBroadHistoneHelas3H3k27me3StdAlnRep1.markdup.q10.sorted.RData',
    'wgEncodeBroadHistoneHelas3H3k27me3StdAlnRep2.markdup.q10.sorted.RData',
    'wgEncodeBroadHistoneHelas3H3k36me3StdAlnRep1.markdup.q10.sorted.RData',
    'wgEncodeBroadHistoneHelas3H3k36me3StdAlnRep2.markdup.q10.sorted.RData',
    'wgEncodeBroadHistoneHelas3Ezh239875AlnRep1.markdup.q10.sorted.RData',
    'wgEncodeBroadHistoneHelas3Ezh239875AlnRep2.markdup.q10.sorted.RData'
  )
method <- 'epigraHMM'

# Fixed parameters
chromosome <- paste0('chr', 1:22)

# Organizing input data

for (idx in seq_len(length(rdata))) {
  load(file.path(
    '../../../../Data',
    paste0('Encode_', tolower(cellline)),
    condition[idx],
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

colData <- colData[chr %in% chromosome, ]
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
cptime[[paste(method, cellline, type, 'Time', paste0(bin, 'bp'),paste0(length(pattern),'Patterns'), sep = '_')]] <- microbenchmark::microbenchmark({
  # Creating object
  out.epigraHMM <- epigraHMM::epigraHMMDataSetFromMatrix(
    countData = as.matrix(colData[, paste(condition, replicate, sep = '.'), with = FALSE]),
    colData = data.frame(condition = condition, replicate = replicate),
    rowRanges = rowRanges
  )
  
  # Removing unnecessary objects
  rm(colData, rowRanges)
  
  # Normalization
  out.epigraHMM <-
    epigraHMM::normalizeCounts(out.epigraHMM, control = controlEM(), span = 1)
  
  # Initialization
  out.epigraHMM <-
    epigraHMM::initializer(object = out.epigraHMM,
                           control = epigraHMM::controlEM(maxIterEM = 15))
  
  assign(
    paste(method, cellline, type, 'Output', paste0(bin, 'bp'),paste0(length(pattern),'Patterns'), sep = '_'),
    epigraHMM::epigraHMM(
      object = out.epigraHMM,
      type = 'differential',
      dist = 'nb',
      control = epigraHMM::controlEM(
        pattern = pattern,
        tempDir = normalizePath(diroutput),
        fileName = paste(method, cellline, type, 'Output', paste0(bin, 'bp'),paste0(length(pattern),'Patterns'), sep = '_')
      )
    )
  )
},times = 1)

### Saving computing time
save(cptime,file = paste0(diroutput, paste(
  method, cellline, type, 'Time', paste0(bin, 'bp'),paste0(length(pattern),'Patterns','.RData'), sep = '_'
)))

### Saving output
do.call(save, list(
  paste(method, cellline, type, 'Output', paste0(bin, 'bp'),paste0(length(pattern),'Patterns'), sep = '_'),
  file = paste0(diroutput, paste(
    method, cellline, type, 'Output', paste0(bin, 'bp'),paste0(length(pattern),'Patterns','.RData'), sep = '_'
  ))
))

### Saving bed file
tobed <-
  get(paste(method, cellline, type, 'Output', paste0(bin, 'bp'),paste0(length(pattern),'Patterns'), sep = '_'))
tobed <- callPeaks(
  object = tobed,
  hdf5 = file.path(
    diroutput,
    paste(method, cellline, type, 'Output', paste0(bin, 'bp'),paste0(length(pattern),'Patterns','.h5'), sep = '_')
  ),
  method = 'viterbi',
  saveToFile = FALSE
)

rtracklayer::export.bed(object = tobed, con = file.path(
  diroutput,
  paste(method, cellline, type, 'Bed', paste0(bin, 'bp'),paste0(length(pattern),'Patterns','.bed'), sep = '_')
))
