library(csaw)
library(bamsignals)

load('./hg19/human.hg19.ranges.mappable.RData')
load('./hg19/human.hg19.ranges.blacklist.RData')
size = c(100,250,500,750,1000)

paths <- list.dirs('.',recursive = FALSE)
paths <- paths[grep('Encode',paths)]

for(folder in paths){
  
  dirs = list.dirs(path = folder, full.names = TRUE, recursive = TRUE)
  
  # Sorting and removing duplicates
  
  for(i in dirs){
    files = list.files(path = i,pattern = ".*.bam$",full.names = TRUE)
    if(length(files)>0){
      for(j in files){
        # Sorting
        cmd = paste('samtools sort -n -o',gsub('.bam','.sort.bam',j),j)
        cat(paste('Command:',cmd,'\n'))
        system(cmd)
        
        cmd = paste('rm',j)
        cat(paste('Command:',cmd,'\n'))
        system(cmd)
        
        # Add ms and MC tags for markdup to use later
        cmd = paste('samtools fixmate -m',gsub('.bam','.sort.bam',j),gsub('.bam','.sort.fixmate.bam',j))
        cat(paste('Command:',cmd,'\n'))
        system(cmd)
        
        cmd = paste('rm',gsub('.bam','.sort.bam',j))
        cat(paste('Command:',cmd,'\n'))
        system(cmd)
        
        # Markdup needs position order
        cmd = paste('samtools sort -o',gsub('.bam','.sort.fixmate.positionsort.bam',j),gsub('.bam','.sort.fixmate.bam',j))
        cat(paste('Command:',cmd,'\n'))
        system(cmd)
        
        cmd = paste('rm',gsub('.bam','.sort.fixmate.bam',j))
        cat(paste('Command:',cmd,'\n'))
        system(cmd)
        
        # Finally mark duplicates
        cmd = paste('samtools markdup -r',gsub('.bam','.sort.fixmate.positionsort.bam',j),gsub('.bam','.markdup.bam',j))
        cat(paste('Command:',cmd,'\n'))
        system(cmd)
        
        cmd = paste('rm',gsub('.bam','.sort.fixmate.positionsort.bam',j))
        cat(paste('Command:',cmd,'\n'))
        system(cmd)
        
        # Remove low quality reads
        cmd = paste('samtools view -q 10 -b',gsub('.bam','.markdup.bam',j),'>',gsub('.bam','.markdup.q10.bam',j))
        cat(paste('Command:',cmd,'\n'))
        system(cmd)
        
        cmd = paste('rm',gsub('.bam','.markdup.bam',j))
        cat(paste('Command:',cmd,'\n'))
        system(cmd)
      }
    }
  }
  
  # Sorting again and Indexing
  
  for(i in dirs){
    files = list.files(path = i,pattern = ".*.markdup.q10.bam$",full.names = T)
    if(length(files)>0){
      for(j in files){
        # Sorting by leftmost coordinates
        cmd = paste('samtools sort -o',gsub('.bam','.sorted.bam',j),j)
        cat(paste('Command:',cmd,'\n'))
        system(cmd)
        
        cmd = paste('rm',j)
        cat(paste('Command:',cmd,'\n'))
        system(cmd)
        
        # Indexing
        cmd = paste('samtools index -b',gsub('.bam','.sorted.bam',j),gsub('.bam','.sorted.bam.bai',j))
        cat(paste('Command:',cmd,'\n'))
        system(cmd)
      }
    }
  }
  
  # Transforming .bam to .bed (to be used by JAMM)
  
  for(i in dirs){
    files = list.files(path = i,pattern = ".*.bam$",full.names = T)
    if(length(files)>0){
      for(j in files){
        cmd = paste('bedtools bamtobed -i',j,'>',gsub('.bam','.bed',j))
        cat(paste('Command:',cmd,'\n'))
        system(cmd)
      }
    }
  }
  
  # Estimating fragment length and tabulating read counts
  
  for(i in dirs){
    files = list.files(path = i,pattern = ".*.bam$",full.names = T)
    if(length(files)>0){
      for(j in files){
        cat(paste('File:',j,'\n'))
        
        # Creating list
        counts = list()
        
        # Fragment length
        fleng = maximizeCcf(correlateReads(j,param=readParam(discard=hg19.discard)))
        cat(paste('Frag. length:',fleng,'\n'))
        sft = ifelse(strsplit(j,'/')[[1]][3]=='Control',0,fleng/2)
        cat(paste('Shift:',sft,'\n'))
        
        # Looping through the window sizes
        for(k in size){
          cat(paste('Window:',k,'\n'))
          gr.tile = with(data.frame(GenomicRanges::tile(hg19,width=k))[,c('seqnames','start','end')],
                         GenomicRanges::GRanges(seqnames,IRanges::IRanges(start,end)))
          
          # Checking for valid chromosomes
          valid.tile = as.logical(GenomeInfoDb::seqnames(gr.tile) %in% GenomeInfoDb::seqlevels(Rsamtools::BamFile(j)))
          
          counts[[paste0(k)]] = data.frame(chr=as.character(seqnames(gr.tile[valid.tile])),start=start(gr.tile[valid.tile]),stop=end(gr.tile[valid.tile]),
                                           counts=bamCount(bampath=j,gr=gr.tile[valid.tile],verbose=T,shift=sft))
        }
        # Writing .RData
        save(counts,file=gsub('.bam','.RData',j),compress = "xz")
      }
    }
  }
  
  # Pooling .bam files, sorting and transforming .bam to .bed (to be used by JAMM)
  
  for(i in dirs){
    files = list.files(path = i,pattern = ".*.bam$",full.names = T)
    if(length(files)>0){
      # Setting name
      namebam = paste0(strsplit(i,'/')[[1]][2:3],collapse='_')
      
      # Pooling
      cmd = paste('samtools merge',paste0(i,'/',namebam,'.bam'),paste0(files,collapse=' '))
      cat(paste('Command:',cmd,'\n'))
      system(cmd)
      
      # Sorting by leftmost coordinates
      cmd = paste('samtools sort -o',paste0(i,'/',namebam,'.sorted.bam'),paste0(i,'/',namebam,'.bam'))
      cat(paste('Command:',cmd,'\n'))
      system(cmd)
      
      cmd = paste('rm',paste0(i,'/',namebam,'.bam'))
      cat(paste('Command:',cmd,'\n'))
      system(cmd)
      
      # Indexing
      cmd = paste('samtools index -b',paste0(i,'/',namebam,'.sorted.bam'),paste0(i,'/',namebam,'.sorted.bam.bai'))
      cat(paste('Command:',cmd,'\n'))
      system(cmd)
      
      # Transforming .bam to .bed
      cmd = paste('bedtools bamtobed -i',paste0(i,'/',namebam,'.sorted.bam'),'>',paste0(i,'/',namebam,'.sorted.bed'))
      cat(paste('Command:',cmd,'\n'))
      system(cmd)
    }
  }
  
  # Tabulating read counts for the pooled data
  
  for(i in dirs){
    files = list.files(path = i,pattern = ".*.bam$",full.names = T)
    if(length(files)>0){
      # Setting name
      namebam = paste0(strsplit(i,'/')[[1]][2:3],collapse='_')
      
      cat(paste('File:',paste0(i,'/',namebam,'.sorted.bam'),'\n'))
      
      # Creating list
      counts = list()
      
      # Fragment length
      fleng = maximizeCcf(correlateReads(paste0(i,'/',namebam,'.sorted.bam'),param=readParam(discard=hg19.discard)))
      cat(paste('Frag. length:',fleng,'\n'))
      sft = ifelse(strsplit(i,'/')[[1]][3]=='Control',0,fleng/2)
      cat(paste('Shift:',sft,'\n'))
      
      # Looping through the window sizes
      for(k in size){
        cat(paste('Window:',k,'\n'))
        gr.tile = with(data.frame(GenomicRanges::tile(hg19,width=k))[,c('seqnames','start','end')],
                       GenomicRanges::GRanges(seqnames,IRanges::IRanges(start,end)))
        
        # Checking for valid chromosomes
        valid.tile = as.logical(GenomeInfoDb::seqnames(gr.tile) %in% GenomeInfoDb::seqlevels(Rsamtools::BamFile(paste0(i,'/',namebam,'.sorted.bam'))))
        
        counts[[paste0(k)]] = data.frame(chr=as.character(seqnames(gr.tile[valid.tile])),start=start(gr.tile[valid.tile]),stop=end(gr.tile[valid.tile]),
                                         counts=bamCount(bampath=j,gr=gr.tile[valid.tile],verbose=T,shift=sft))
      }
      # Writing .RData
      save(counts,file=gsub('.bam','.RData',paste0(i,'/',namebam,'.sorted.bam')),compress = "xz")
    }
  } 
}
