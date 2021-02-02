# (genome and chromosome files come from RGT installation http://www.regulatory-genomics.org/rgt/rgt-data-folder/)

library(microbenchmark)

method <- 'THOR'
mark <- 'EZH2'
cell <- 'Encode_twocells'
bp <- 750

outdir = paste0('Output', bp, '/')
system(paste('mkdir', outdir))

cmd = paste(
  'rgt-THOR',
  paste0('THOR_', mark, '_', cell, '.config'),
  '--name',
  paste0('THOR_', mark, '_', cell, '_', bp, 'bp'),
  '-b',
  bp,
  '--pvalue 1.0',
  '--output-dir',
  outdir
)
cat('Command: ', cmd)

cptime <- microbenchmark({
  system(cmd)
}, times = 1)

save(cptime, file = paste0(outdir, paste(
  method, mark, cell, 'Time', paste0(bp, 'bp.RData'), sep = '_'
)))