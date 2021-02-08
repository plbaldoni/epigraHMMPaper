library(data.table)

bpsize = 750

chrom = '../../../../Data/human-hg19-size.bed'
deadzone = '../../../../Data/deadzones-k36-hg19.bed'

method = 'RSEG'
mark = 'H3K27me3'
cell = 'Encode_twocells'

chip1.pool = paste0('Encode_helas3_', mark, '.sorted.bed')
chip2.pool = paste0('Encode_hepg2_', mark, '.sorted.bed')

datadir = '../../../../Data/'

for (bp in bpsize) {
  cptime = list()
  
  system(paste('mkdir', paste0('./temporary_data', bp)))
  system(paste(
    'mkdir',
    paste0('./temporary_data', bp, '/chip1'),
    paste0('./temporary_data', bp, '/chip2')
  ))
  
  #######################
  ### Pooled analysis ###
  #######################
  cmd = paste(
    'cp',
    paste0(datadir, 'Encode_helas3/', mark, '/', chip1.pool, collapse = ' '),
    paste0('./temporary_data', bp, '/chip1')
  )
  cat('Command: ', cmd, '\n')
  system(cmd)
  
  chipfiles = list.files(
    path = paste0('./temporary_data', bp, '/chip1'),
    pattern = '*.bed$',
    full.names = T
  )
  system('export LC_ALL=C')
  system(paste(
    'sort -k1,1 -k3,3n -k2,2n -k6,6r',
    chipfiles,
    paste('-o', paste0(
      './temporary_data', bp, '/chip1/chip1.bed'
    ))
  ))
  
  cmd = paste(
    'cp',
    paste0(datadir, 'Encode_hepg2/', mark, '/', chip2.pool, collapse = ' '),
    paste0('./temporary_data', bp, '/chip2')
  )
  cat('Command: ', cmd, '\n')
  system(cmd)
  
  chipfiles = list.files(
    path = paste0('./temporary_data', bp, '/chip2'),
    pattern = '*.bed$',
    full.names = T
  )
  system('export LC_ALL=C')
  system(paste(
    'sort -k1,1 -k3,3n -k2,2n -k6,6r',
    chipfiles,
    paste('-o', paste0(
      './temporary_data', bp, '/chip2/chip2.bed'
    ))
  ))
  
  outdir = paste0('./Output', bp, '/')
  system(paste('mkdir', outdir))
  
  cptime[[paste(method, mark, cell, 'Output', 'Pooled', paste0(bp, 'bp'), sep =
                  '_')]] =
    microbenchmark::microbenchmark({
      cmd = paste(
        'rseg-diff -verbose -mode 3',
        '-out',
        paste0(
          outdir,
          paste(
            method,
            mark,
            cell,
            'Output',
            'Pooled',
            paste0(bp, 'bp.bed'),
            sep = '_'
          )
        ),
        '-score',
        paste0(
          outdir,
          paste(
            method,
            mark,
            cell,
            'Output',
            'Pooled',
            paste0(bp, 'bp.wig'),
            sep = '_'
          )
        ),
        '-chrom',
        chrom,
        '-bin-size',
        bp,
        #'-deadzones', # With -deadzones, RSEG throws an error : 'rseg-diff: ../common/ReadCounts.cpp:262: void AdjustBinSize(std::vector<SimpleGenomicRegion>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<long unsigne d int>&, size_t, size_t): Assertion `*std::min_element(nondead_scales.begin(), nondead_scales.end())>=0.0' failed.
        #deadzone,
        paste0('./temporary_data', bp, '/chip1/chip1.bed'),
        paste0('./temporary_data', bp, '/chip2/chip2.bed')
      )
      cat('Command: ', cmd)
      system(cmd)
    },times = 1)
  
  ### Saving computing time
  save(cptime, file = paste0(outdir, paste(
    method, mark, cell, 'Time', paste0(bp, 'bp.RData'), sep = '_'
  )))
  
  ### Saving output in RData format
  rseg <-
    fread(paste0(
      outdir,
      paste(
        method,
        mark,
        cell,
        'Output',
        'Pooled',
        paste0(bp, 'bp.wig'),
        sep = '_'
      )
    ))
  save(rseg, file = paste0(
    outdir,
    paste(
      method,
      mark,
      cell,
      'Output',
      'Pooled',
      paste0(bp, 'bp.wig.RData'),
      sep = '_'
    )
  ))
  
  ### Deleting temporary folder
  system(paste('rm -r', paste0('./temporary_data', bp)))
  
  cat('Done!\n')
}
