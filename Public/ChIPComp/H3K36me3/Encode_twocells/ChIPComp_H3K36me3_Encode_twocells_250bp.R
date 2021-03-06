library(ChIPComp)
library(microbenchmark)

bin = 250
method = 'ChIPComp'
mark = 'H3K36me3'
cell = 'Encode_twocells'

chip1 = c('wgEncodeBroadHistoneHelas3H3k36me3StdAlnRep1.markdup.q10.sorted.bam',
          'wgEncodeBroadHistoneHelas3H3k36me3StdAlnRep2.markdup.q10.sorted.bam')

chip2 = c('wgEncodeBroadHistoneHepg2H3k36me3StdAlnRep1.markdup.q10.sorted.bam',
          'wgEncodeBroadHistoneHepg2H3k36me3StdAlnRep2.markdup.q10.sorted.bam')

control1 = c('wgEncodeBroadHistoneHelas3ControlStdAlnRep1.markdup.q10.sorted.bam',
             'wgEncodeBroadHistoneHelas3ControlStdAlnRep2.markdup.q10.sorted.bam')

control2 = c('wgEncodeBroadHistoneHepg2ControlStdAlnRep1.markdup.q10.sorted.bam',
             'wgEncodeBroadHistoneHepg2ControlStdAlnRep2.markdup.q10.sorted.bam')

dirchip1 = '../../../../Data/Encode_helas3/'
dirchip2 = '../../../../Data/Encode_hepg2/'

### Organizing other parameters

cptime = list()

if(mark%in%c('H3K36me3','H3K27me3','EZH2')){
  dirpeak1 = paste0('../../../MACS2/',mark,'/Helas3/Output/MACS2_',mark,'_Helas3_Output_peaks.broadPeak')
  dirpeak2 = paste0('../../../MACS2/',mark,'/Hepg2/Output/MACS2_',mark,'_Hepg2_Output_peaks.broadPeak') 
}
if(mark%in%c('H3K4me3','H3K27ac','CTCF')){
  dirpeak1 = paste0('../../../MACS2/',mark,'/Helas3/Output/MACS2_',mark,'_Helas3_Output_peaks.narrowPeak')
  dirpeak2 = paste0('../../../MACS2/',mark,'/Hepg2/Output/MACS2_',mark,'_Hepg2_Output_peaks.narrowPeak') 
}

outdir = paste0('./Output', bin, '/')
system(paste('mkdir', outdir))

### Input for ChIPComp

conf <- data.frame(SampleID = 1:4,
                   condition = c('Helas3','Helas3','Hepg2','Hepg2'),
                   factor = rep(mark,4),
                   ipReads = c(paste0(dirchip1,mark,'/',chip1[1]),paste0(dirchip1,mark,'/',chip1[2]),paste0(dirchip2,mark,'/',chip2[1]),paste0(dirchip2,mark,'/',chip2[2])),
                   ctReads = c(paste0(dirchip1,'Control','/',control1[1]),paste0(dirchip1,'Control','/',control1[2]),paste0(dirchip2,'Control','/',control2[1]),paste0(dirchip2,'Control','/',control2[2])),
                   peaks = c(dirpeak1,dirpeak1,dirpeak2,dirpeak2))

design = as.data.frame(lapply(conf[, c("condition", "factor")], function(x){as.numeric(as.factor(x))})) - 1
design = as.data.frame(model.matrix( ~ condition, design))

### Running ChIPComp

for(bp in bin){
  cptime[[paste(method,mark,cell,'Output',bp,sep='_')]] = microbenchmark({
    countSet = makeCountSet(conf,design,filetype="bam",species="hg19",binsize=bp)
    countSet = ChIPComp(countSet)
    write.table(countSet$db,file=paste0(outdir,paste(method,mark,cell,'Output',paste0(bp,'bp.txt'),sep='_')),quote=F,row.names=F)
  },times=1)
}

### Saving computing time
save(cptime, file = paste0(outdir, paste(
  method, mark, cell, 'Time', paste0(bin, 'bp.RData'), sep = '_'
)))