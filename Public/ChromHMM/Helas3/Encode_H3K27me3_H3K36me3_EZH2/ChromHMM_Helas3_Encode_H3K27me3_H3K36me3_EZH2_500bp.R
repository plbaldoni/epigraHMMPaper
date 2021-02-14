### ChromHMM - Helas3
### General parameters
bp = 500

statenum <- c(3,4,5,6,7,8)

cell1 = c('wgEncodeBroadHistoneHelas3H3k27me3StdAlnRep1.markdup.q10.sorted.bam',
          'wgEncodeBroadHistoneHelas3H3k27me3StdAlnRep2.markdup.q10.sorted.bam')

cell2 = c('wgEncodeBroadHistoneHelas3H3k36me3StdAlnRep1.markdup.q10.sorted.bam',
          'wgEncodeBroadHistoneHelas3H3k36me3StdAlnRep2.markdup.q10.sorted.bam')

cell3 = c('wgEncodeBroadHistoneHelas3Ezh239875AlnRep1.markdup.q10.sorted.bam',
          'wgEncodeBroadHistoneHelas3Ezh239875AlnRep2.markdup.q10.sorted.bam')

control = c('wgEncodeBroadHistoneHelas3ControlStdAlnRep1.markdup.q10.sorted.bam',
            'wgEncodeBroadHistoneHelas3ControlStdAlnRep2.markdup.q10.sorted.bam')

### Creating temporary input/output file directory
tmpdir = paste0('./Input',bp)
outdir = paste0('./Output',bp)
bindir = paste0('./BinarizeBed',bp)

system(paste('mkdir',tmpdir))
system(paste('mkdir',bindir))
system(paste('mkdir',outdir))

### Moving input files

dirdt = '../../../../Data/'

method = 'ChromHMM'
cellline = 'Helas3'
marks = c('H3K27me3','H3K36me3','EZH2')
type <- paste(c('Encode',marks),collapse = '_')

### Saving the data.frame with the inpuit information

dt <- data.frame(Cell=rep(paste(method,cellline,type,paste0(bp,'bp'),sep='_'),6),
                 Mark=c(rep('H3K27me3',length(cell1)),rep('H3K36me3',length(cell2)),rep('EZH2',length(cell3))),
                 ChIP=c(cell1,cell2,cell3),
                 Control=c(control,control,control))

write.table(dt,file=paste0('Input',bp,'.txt'),row.names=F,col.names=F,quote=F,sep="\t")

# Mark 1

cmd = paste('cp',paste0(dirdt,'Encode_helas3/H3K27me3/',cell1,collapse = ' '),
            paste0(dirdt,'Encode_helas3/H3K27me3/',gsub('.bam','.bam.bai',cell1),collapse = ' '),
            paste0(tmpdir))
cat('Command: ',cmd,'\n')
system(cmd)

# Mark 2

cmd = paste('cp',paste0(dirdt,'Encode_helas3/H3K36me3/',cell2,collapse = ' '),
            paste0(dirdt,'Encode_helas3/H3K36me3/',gsub('.bam','.bam.bai',cell2),collapse = ' '),
            paste0(tmpdir))
cat('Command: ',cmd,'\n')
system(cmd)

# Mark 3

cmd = paste('cp',paste0(dirdt,'Encode_helas3/EZH2/',cell3,collapse = ' '),
            paste0(dirdt,'Encode_helas3/EZH2/',gsub('.bam','.bam.bai',cell3),collapse = ' '),
            paste0(tmpdir))
cat('Command: ',cmd,'\n')
system(cmd)

# Controls

cmd = paste('cp',paste0(dirdt,'Encode_helas3/Control/',control,collapse = ' '),
            paste0(dirdt,'Encode_helas3/Control/',gsub('.bam','.bam.bai',control),collapse = ' '),
            paste0(tmpdir))
cat('Command: ',cmd,'\n')
system(cmd)

### Runnign ChromHMM-BinarizeBed

cmd = paste('java -jar ~/ChromHMM/ChromHMM.jar BinarizeBam -b',bp,'-c',tmpdir,
            '~/ChromHMM/CHROMSIZES/hg19.txt',tmpdir,paste0('Input',bp,'.txt'),bindir)
cat('Command: ',cmd,'\n')
system(cmd)

### Runnign ChromHMM-LearnModel

for(num in statenum){
  cmd = paste('unset DISPLAY && java -jar ~/ChromHMM/ChromHMM.jar LearnModel -p 0 -b',bp,
              bindir,outdir,num,'hg19')
  cat('Command: ',cmd,'\n')
  system(cmd)   
}

### Deleting extra files
system(paste('rm -r',paste0('Input',bp,'.txt'),tmpdir,bindir))