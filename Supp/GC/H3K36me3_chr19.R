library(data.table)
library(SummarizedExperiment)
library(ggplot2)
library(magrittr)
library(epigraHMM)

foo <- function(mark,windowSize,chrlist = c(paste0('chr',1:22),'chrX'),gcRange = c(0.28,0.67),
                dataDir = "/pine/scr/b/a/baldoni/Rashid/epigraHMMPaper/Data"){
  
  # Setting up paths
  allDir <- normalizePath(dataDir)
  allDir <- list.dirs(path = allDir,full.names = TRUE,recursive = TRUE)
  allDir <- allDir[grep(pattern = mark,x = allDir)]
  
  # Getting bam files and Filtering relevant cell lines
  bamDir <- list.files(path = allDir,pattern = 'wgEncode.*.bam$',full.names = TRUE,recursive = TRUE)
  bamDir <- bamDir[grep("helas3|hepg2",bamDir)]
  
  # Setting up input for epigraHMM
  colData <- data.table(condition = unlist(lapply(strsplit(bamDir,'/'),function(x){x[10]})))
  colData[,replicate := sequence(.N),by = 'condition']
  colData[,bam := bamDir]
  colData <- as.data.frame(colData)
  
  # Tiling up the specified genome
  object <- epigraHMM::epigraHMMDataSetFromBam(colData = colData,
                                               genome = 'hg19',
                                               windowSize = windowSize,
                                               gapTrack = TRUE,
                                               blackList = TRUE)
  
  # Subsetting
  object <- object[GenomeInfoDb::seqnames(object)%in%chrlist]
  
  # Getting GC
  gc <- rowSums(Biostrings::alphabetFrequency(IRanges::Views(BSgenome.Hsapiens.UCSC.hg19::Hsapiens,object@rowRanges))[,c('C','G')])/width(object@rowRanges)
  gc.trimmed <- gc
  gc.trimmed[(gc.trimmed < gcRange[1]) | (gc.trimmed > gcRange[2])] <- NA
  
  # Plotting Counts vs. GC
  figRawCounts <- data.table(Counts = c(assay(object)),
                             GC = rep(gc.trimmed,ncol(object)),
                             Condition = rep(object@colData$condition,each = length(object)),
                             Replicate = rep(object@colData$replicate,each = length(object))) %>%
    ggplot(aes(y = log1p(Counts),x = GC))+
    facet_grid(rows = vars(Replicate),cols = vars(Condition))+
    geom_point(shape = ".",alpha = 0.25)+
    theme_bw()
  
  # Estimating GC
  cts.GCadj <- gcapc::refineSites(counts = assay(object),
                                  sites = object@rowRanges,
                                  gcrange = c(0,1),
                                  flank = 0,
                                  plot = FALSE,
                                  genome = 'hg19')
  
  figAdjCounts <- data.table(`Adjusted Counts` = c(cts.GCadj),
                             GC = rep(gc.trimmed,ncol(object)),
                             Condition = rep(object@colData$condition,each = length(object)),
                             Replicate = rep(object@colData$replicate,each = length(object))) %>%
    ggplot(aes(y = log1p(`Adjusted Counts`),x = GC))+
    facet_grid(rows = vars(Replicate),cols = vars(Condition))+
    geom_point(shape = ".",alpha = 0.25)+
    theme_bw()
  
  # Calling peaks WITHOUT GC
  
  ## Creating control
  
  control <- epigraHMM::controlEM()
  
  ## Initializing EM
  
  outputInitial <- epigraHMM::initializer(object = object,
                                          control = control)
  
  ## Sample-specific peaks
  
  samplePeaks <- lapply(seq_len(ncol(outputInitial)),function(x){
    output <- epigraHMM::epigraHMM(object = outputInitial[,x],
                                   type = 'consensus',
                                   control = control,
                                   dist = 'zinb',
                                   random = FALSE)
    return(metadata(output)$prob$Enrichment)
  })
  
  ## Single-condition consensus peaks
  
  consensusPeaks <- lapply(unique(outputInitial@colData$condition),function(x){
    
    cols <- (outputInitial@colData$condition == x)
    message('Cols: ',cols)
    
    output <- outputInitial[,cols]
    output <- epigraHMM::normalizeCounts(object = output,
                                         type = 'gam',
                                         by = NULL)
    output <- epigraHMM::epigraHMM(object = output,
                                   control = control,
                                   type = 'consensus',
                                   dist = 'zinb',
                                   random = FALSE)
    return(metadata(output)$prob$Enrichment)
  })
  
  ## Multiple-condition differential peaks
  
  differentialPeaks <- lapply('donothing',function(x){
    
    output <- outputInitial
    output <- epigraHMM::normalizeCounts(object = output,
                                         type = 'gam',
                                         by = NULL)
    output <- epigraHMM::epigraHMM(object = output,
                                   control = control,
                                   type = 'differential',
                                   dist = 'zinb',
                                   random = FALSE)
    return(metadata(output)$prob$Differential)
  })
  
  # Calling peaks WITH GC
  
  ## Preparing offsets
  
  offset.GCadj <- log((assay(outputInitial)+1)/(cts.GCadj+1))
  outputInitial.GCadj <- outputInitial
  
  dimnames(offset.GCadj) <- dimnames(assay(outputInitial.GCadj))
  assay(outputInitial.GCadj,'offset',withDimnames=TRUE) <- offset.GCadj
  
  ## Sample-specific peaks 
  
  samplePeaks.GC <- lapply(seq_len(ncol(outputInitial.GCadj)),function(x){
    
    output <- epigraHMM::epigraHMM(object = outputInitial.GCadj[,x],
                                   type = 'consensus',
                                   control = control,
                                   dist = 'zinb',
                                   random = FALSE)
    return(metadata(output)$prob$Enrichment)
  })
  
  ## Single-condition consensus peaks WITH GC
  
  consensusPeaks.GC <- lapply(unique(outputInitial.GCadj@colData$condition),function(x){
    
    cols <- (outputInitial.GCadj@colData$condition == x)
    message('Cols: ',cols)
    
    output <- outputInitial.GCadj[,cols]
    output <- epigraHMM::normalizeCounts(object = output,
                                         type = 'gam',
                                         by = NULL)
    output <- epigraHMM::epigraHMM(object = output,
                                   control = control,
                                   type = 'consensus',
                                   dist = 'zinb',
                                   random = FALSE)
    return(metadata(output)$prob$Enrichment)
  })
  
  ## Multiple-condition differential peaks WITH GC
  
  differentialPeaks.GC <- lapply('donothing',function(x){
    
    output <- outputInitial.GCadj
    output <- epigraHMM::normalizeCounts(object = output,
                                         type = 'gam',
                                         by = NULL)
    output <- epigraHMM::epigraHMM(object = output,
                                   control = control,
                                   type = 'differential',
                                   dist = 'zinb',
                                   random = FALSE)
    return(metadata(output)$prob$Differential)
  })
  
  # Returning
  return(list('figRawCounts' = figRawCounts,
              'figAdjCounts' = figAdjCounts,
              'samplePeaks' = samplePeaks,
              'consensusPeaks' = consensusPeaks,
              'differentialPeaks' = differentialPeaks,
              'outputInitial.GCadj' = outputInitial.GCadj,
              'samplePeaks.GC' = samplePeaks.GC,
              'consensusPeaks.GC' = consensusPeaks.GC,
              'differentialPeaks.GC' = differentialPeaks.GC))
  
}

# H3K36me3

foo_H3K36me3 <- lapply(c(250,500,1000,10000),function(x){
  foo(mark = 'H3K36me3',windowSize = x,chrlist = 'chr19')
})

save(foo_H3K36me3,file = 'H3K36me3_chr19.RData',compress = 'xz')
rm(foo_H3K36me3)
