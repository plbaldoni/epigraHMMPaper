library(csaw)
library(edgeR)
library(rtracklayer)
library(statmod)
library(xscss)
library(microbenchmark)

mark = 'H3K27me3'
bpvec = 250 # Window size
method = 'csaw'
cell = 'Encode_threecells'

chip1 = c(
  'wgEncodeBroadHistoneHelas3H3k27me3StdAlnRep1.markdup.q10.sorted.bam',
  'wgEncodeBroadHistoneHelas3H3k27me3StdAlnRep2.markdup.q10.sorted.bam'
)
chip2 = c(
  'wgEncodeBroadHistoneHepg2H3k27me3StdAlnRep1.markdup.q10.sorted.bam',
  'wgEncodeBroadHistoneHepg2H3k27me3StdAlnRep2.markdup.q10.sorted.bam'
)

chip3 = c(
  'wgEncodeBroadHistoneHuvecH3k27me3StdAlnRep1.markdup.q10.sorted.bam',
  'wgEncodeBroadHistoneHuvecH3k27me3StdAlnRep2.markdup.q10.sorted.bam'
)

cell.type = c('Helas3', 'Helas3', 'Hepg2', 'Hepg2', 'Huvec', 'Huvec') # Cell type
contrast = 2:3 # Contrast for test differential

for (bp in bpvec) {
  cptime = list()
  tmpdir = paste0('temporary_data', bp)
  outdir = paste0('Output', bp)
  
  system(paste('mkdir', tmpdir))
  system(paste('mkdir', paste0(tmpdir, '/chip'), paste0(tmpdir, '/control')))
  system(paste('mkdir', outdir))
  
  ### Loading data
  cmd = paste('cp',
              paste(
                list.files(
                  file.path('../../../../Data/Encode_helas3/', mark),
                  'wg*.*.bam',
                  full.names = TRUE
                ),
                collapse = ' '
              ),
              paste0(tmpdir, '/chip'))
  cat('Command: ', cmd, '\n')
  system(cmd)
  
  cmd = paste('cp',
              paste(
                list.files(
                  file.path('../../../../Data/Encode_hepg2/', mark),
                  'wg*.*.bam',
                  full.names = TRUE
                ),
                collapse = ' '
              ),
              paste0(tmpdir, '/chip'))
  cat('Command: ', cmd, '\n')
  system(cmd)
  
  cmd = paste('cp',
              paste(
                list.files(
                  file.path('../../../../Data/Encode_huvec/', mark),
                  'wg*.*.bam',
                  full.names = TRUE
                ),
                collapse = ' '
              ),
              paste0(tmpdir, '/chip'))
  cat('Command: ', cmd, '\n')
  system(cmd)
  
  cptime[[paste(method, mark, cell, 'Output', bp, sep = '_')]] = microbenchmark({
    widebp = 10 * bp #Size of large bins to compute normalization
    max.delay = 250 #Max. shift to attempt in cross-correlation analysis
    common.length = 200 #Rescale fragment length to this value
    fold.change = log2(2) #Fold change for filter
    filterct = 20 #An integer scalar for the minimum count sum across libraries for each window.
    tol = 100 #Parameters for mergeWindows
    max.width = 5000 #Parameters for mergeWindows
    
    ### csaw starts now
    # List of bam files
    bam.files = list.files(
      path = paste0(tmpdir, '/chip'),
      pattern = '*.bam$',
      full.names = T
    )
    bam.files
    
    # Design matrix
    design <- model.matrix(~ factor(cell.type))
    colnames(design) <- c("intercept", unique(cell.type)[-1])
    design
    
    # Parameters (PCR duplicates already removed and quality score filtered)
    param <- readParam(dedup = F)
    param
    
    # Estimating the average fragment length (rescaling all to 200bp)
    x = lapply(bam.files,
               correlateReads,
               param = param,
               max.dist = max.delay)
    multi.frag.lens = list(unlist(lapply(x, maximizeCcf)), common.length)
    multi.frag.lens
    
    # Counting reads
    data <-
      windowCounts(
        bam.files,
        width = bp,
        ext = multi.frag.lens,
        param = param,
        filter = filterct
      )
    data
    
    # Filtering data
    data.large <-
      windowCounts(bam.files,
                   width = widebp,
                   bin = T,
                   param = param)
    
    bin.ab <-
      scaledAverage(data.large, scale = median(getWidths(data.large)) / median(getWidths(data)))
    
    threshold <- median(bin.ab) + fold.change
    
    keep.global <- aveLogCPM(asDGEList(data)) >  threshold
    
    sum(keep.global)
    
    # Creating filtered data
    filtered.data <- data[keep.global,]
    
    # Testing for DB
    
    y <-
      DGEList(assay(filtered.data), lib.size = filtered.data$totals)
    y$samples$norm.factors <- 1
    y$offset <- NULL
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust = TRUE)
    out <- glmQLFTest(fit,coef  = contrast)
    tabres <- topTags(out, nrow(out))$table
    tabres <- tabres[order(as.integer(rownames(tabres))),]
    
    merged <-
      mergeWindows(rowRanges(filtered.data),
                   tol = tol,
                   max.width = max.width)
    tabneg <- combineTests(merged$id, tabres)
    
    # Organizing output
    
    output = merged$region
    output$PValue = tabneg$PValue
    output$FDR = tabneg$FDR
    
    # Saving file in tsv format with all details preserved
    ofile <-
      gzfile(paste0(
        outdir,
        "/",
        paste(method, mark, cell, paste0(bp, 'bp'), sep = '_'),
        ".tsv.gz"
      ), open = "w")
    write.table(
      as.data.frame(output),
      file = ofile,
      row.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )
    close(ofile)
    
    # Saving file in bed format for Genome Browser visualization
    test <- output
    names(test) <- paste0("region", 1:length(output))
    export(test, paste0(outdir, "/temp.bed"))
    system(paste(
      "echo '",
      paste0(
        'track name="csaw (',
        bp,
        'bp)" description="',
        mark,
        '" color=150,150,0'
      ),
      "' | cat - ",
      paste0(outdir, "/temp.bed"),
      " > ",
      paste0(
        outdir,
        "/",
        paste(method, mark, cell, paste0(bp, 'bp'), sep = '_'),
        ".bed"
      )
    ))
    system(paste0('rm ', paste0(outdir, "/temp.bed")))
    
    ### Removing files
    system(paste('rm -r', tmpdir))
  }, times = 1)
  
  ### Saving computing time
  save(cptime, file = paste0(outdir, paste(
    method, mark, cell, 'Time', paste0(bp, 'bp.RData'), sep = '_'
  )))
}