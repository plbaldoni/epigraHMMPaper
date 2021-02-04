method = 'diffReps'
mark = 'H3K27ac'
cell = 'Encode_twocells'
bp = 750

chip1 = c(
  'wgEncodeBroadHistoneHelas3H3k27acStdAlnRep1.markdup.q10.sorted.bed',
  'wgEncodeBroadHistoneHelas3H3k27acStdAlnRep2.markdup.q10.sorted.bed'
)

control1 = c(
  'wgEncodeBroadHistoneHelas3ControlStdAlnRep1.markdup.q10.sorted.bed',
  'wgEncodeBroadHistoneHelas3ControlStdAlnRep2.markdup.q10.sorted.bed'
)

chip2 = c(
  'wgEncodeBroadHistoneHepg2H3k27acStdAlnRep1.markdup.q10.sorted.bed',
  'wgEncodeBroadHistoneHepg2H3k27acStdAlnRep2.markdup.q10.sorted.bed'
)

control2 = c(
  'wgEncodeBroadHistoneHepg2ControlStdAlnRep1.markdup.q10.sorted.bed',
  'wgEncodeBroadHistoneHepg2ControlStdAlnRep2.markdup.q10.sorted.bed'
)

dirchip1 = '../../../../Data/Encode_helas3/'
dirchip2 = '../../../../Data/Encode_hepg2/'

dircontrol1 = '../../../../Data/Encode_helas3/'
dircontrol2 = '../../../../Data/Encode_hepg2/'

for (j in bp) {
  outdir = paste0('./Output', j, '/')
  system(paste('mkdir', outdir))
  
  cptime = list()
  
  cptime[[paste(method, mark, cell, 'Output', paste0(j, 'bp'), sep = '_')]] =
    system.time({
      cmd = paste(
        'diffReps.pl --gname hg19 --report',
        paste0(
          outdir,
          paste(method, mark, cell, 'Output', paste0(j, 'bp.txt'), sep = '_')
        ),
        '--treatment',
        paste0(dirchip1, mark, '/', chip1, collapse = ' '),
        '--control',
        paste0(dirchip2, mark, '/', chip2, collapse = ' '),
        '--btr',
        paste0(dircontrol1, 'Control/', control1, collapse = ' '),
        '--bco',
        paste0(dircontrol2, 'Control/', control2, collapse = ' '),
        '--window',
        j,
        '--pval 1',
        '--nsd',
        ifelse(
          mark %in% c('H3K4me3', 'H3K27ac', 'CTCF'),
          'sharp',
          ifelse(mark %in% c('H3K27me3', 'H3K36me3', 'EZH2'), 'broad', NA)
        ),
        '--meth nb'
      )
      cat('Command: ', cmd)
      system(cmd)
    })
  
  ### Saving computing time
  save(cptime, file = paste0(outdir, paste(
    method, mark, cell, 'Time', paste0(j, 'bp.RData'), sep = '_'
  )))
}

cat('Done!\n')