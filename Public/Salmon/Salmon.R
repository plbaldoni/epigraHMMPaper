# Creating transcriptome index (from https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/04_quasi_alignment_salmon.html)
if(!file.exists(file.path('Homo_sapiens.GRCh37.cdna.all.fa.salmon.index','versionInfo.json'))){
    system('mkdir Homo_sapiens.GRCh37.cdna.all.fa.salmon.index')
    system('salmon index -t ../../Data/transcriptome/Homo_sapiens.GRCh37.75.cdna.all.fa -i Homo_sapiens.GRCh37.cdna.all.fa.salmon.index')
}

# Run Salmon

cellline <- c('Helas3','Hepg2','Huvec','H1hesc')
names(cellline) <- c('Encode_helas3','Encode_hepg2','Encode_huvec','Encode_h1hesc')


for(i in seq_len(length(cellline))){
    system(paste('mkdir',paste0(cellline[i])))
    
    cmd = paste('salmon quant --gcBias -i Homo_sapiens.GRCh37.cdna.all.fa.salmon.index -l A -1',
                paste0(paste0('../../Data/',names(cellline[i]),'/RNAseq/'),paste0('wgEncodeCshlLongRnaSeq',cellline[i],'CellPapFastqRd1Rep1.fastq.gz')),'-2',
                paste0(paste0('../../Data/',names(cellline[i]),'/RNAseq/'),paste0('wgEncodeCshlLongRnaSeq',cellline[i],'CellPapFastqRd2Rep1.fastq.gz')),
                '--validateMappings -p 1 -o',paste0(cellline[i],'/Output1'))
    
    if(!file.exists(file.path(paste0(cellline[i],'/Output1'),'quant.sf'))){system(cmd)}
    
    cmd = paste('salmon quant --gcBias -i Homo_sapiens.GRCh37.cdna.all.fa.salmon.index -l A -1',
                paste0(paste0('../../Data/',names(cellline[i]),'/RNAseq/'),paste0('wgEncodeCshlLongRnaSeq',cellline[i],'CellPapFastqRd1Rep2.fastq.gz')),'-2',
                paste0(paste0('../../Data/',names(cellline[i]),'/RNAseq/'),paste0('wgEncodeCshlLongRnaSeq',cellline[i],'CellPapFastqRd2Rep2.fastq.gz')),
                '--validateMappings -p 1 -o',paste0(cellline[i],'/Output2'))
    
    if(!file.exists(file.path(paste0(cellline[i],'/Output2'),'quant.sf'))){system(cmd)}
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Now, summarizing salmon results with tximport
rm(list=ls())

library(EnsDb.Hsapiens.v75)
library(tximport)
library(SummarizedExperiment)
library(DESeq2)

edb = EnsDb.Hsapiens.v75
Tx <- transcripts(edb,return.type='data.frame')
tx2gene = subset(Tx,select=c('tx_id','gene_id'))

cells = c('H1hesc','Helas3','Hepg2','Huvec')

files = NULL
for(i in cells){
    cat(i)
    files = c(files,list.files(path = i,pattern = '^quant.sf$',recursive = T,full.names = T))
}

for(i in files){
    cellname = strsplit(i,'/')[[1]][length(strsplit(i,'/')[[1]])-2]
    repname = strsplit(i,'/')[[1]][length(strsplit(i,'/')[[1]])-1]; repname = substr(repname,nchar(repname),nchar(repname))
    names(files)[which(files==i)] = paste0(cellname,'.Rep',repname)
}

txi <- tximport(files, type = "salmon", tx2gene = tx2gene,countsFromAbundance='no')
txi.scaled <- tximport(files, type = "salmon", tx2gene = tx2gene,countsFromAbundance='lengthScaledTPM')

### Creating DESeqDataSet
#### Raw
eset.DESeq <- DESeqDataSetFromTximport(txi,colData = data.frame(Cells=unlist(lapply(strsplit(colnames(txi$counts),"\\."), FUN = function(x){x[1]})),
                                                                Replicates = unlist(lapply(strsplit(colnames(txi$counts),"\\."), FUN = function(x){substr(x[2],nchar(x[2]),nchar(x[2]))}))),design=~Cells)
rowRanges(eset.DESeq) <- genes(edb)[rownames(eset.DESeq)]

ENCODE.rnaseq.raw = eset.DESeq
seqlevelsStyle(ENCODE.rnaseq.raw) <- 'UCSC'

save(ENCODE.rnaseq.raw,file='ENCODE.rnaseq.raw.RData',compress = 'xz')

#### Scaled
eset.DESeq.scaled <- DESeqDataSetFromTximport(txi.scaled,
                                              colData = data.frame(Cells=unlist(lapply(strsplit(colnames(txi.scaled$counts),"\\."), FUN = function(x){x[1]})),
                                                                   Replicates = unlist(lapply(strsplit(colnames(txi.scaled$counts),"\\."), FUN = function(x){substr(x[2],nchar(x[2]),nchar(x[2]))}))),design=~Cells)
rowRanges(eset.DESeq.scaled) <- genes(edb)[rownames(eset.DESeq.scaled)]

ENCODE.rnaseq.scaled = eset.DESeq.scaled
seqlevelsStyle(ENCODE.rnaseq.scaled) <- 'UCSC'

save(ENCODE.rnaseq.scaled,file='ENCODE.rnaseq.scaled.RData',compress = 'xz')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Now, I want to check whether gene bodies intersect with ChIP-seq peaks

rm(list=ls())

library(SummarizedExperiment)
library(AnnotationHub)
library(GenomicRanges)
library(rtracklayer)
library(bamsignals)
library(csaw)
library(EnsDb.Hsapiens.v75)
library(tximport)
library(DESeq2)

dirchip = '../../Data/'
chromosome = c(paste0('chr',1:22),'chrX')

# ChIP-seq data

chip.ezh2 <- c(paste0(dirchip,'Encode_h1hesc/EZH2/wgEncodeBroadHistoneH1hescEzh239875AlnRep1.markdup.q10.sorted.bam'),
               paste0(dirchip,'Encode_h1hesc/EZH2/wgEncodeBroadHistoneH1hescEzh239875AlnRep2.markdup.q10.sorted.bam'),
               paste0(dirchip,'Encode_helas3/EZH2/wgEncodeBroadHistoneHelas3Ezh239875AlnRep1.markdup.q10.sorted.bam'),
               paste0(dirchip,'Encode_helas3/EZH2/wgEncodeBroadHistoneHelas3Ezh239875AlnRep2.markdup.q10.sorted.bam'),
               paste0(dirchip,'Encode_hepg2/EZH2/wgEncodeBroadHistoneHepg2Ezh239875AlnRep1.markdup.q10.sorted.bam'),
               paste0(dirchip,'Encode_hepg2/EZH2/wgEncodeBroadHistoneHepg2Ezh239875AlnRep2.markdup.q10.sorted.bam'),
               paste0(dirchip,'Encode_huvec/EZH2/wgEncodeBroadHistoneHuvecEzh239875AlnRep1.markdup.q10.sorted.bam'),
               paste0(dirchip,'Encode_huvec/EZH2/wgEncodeBroadHistoneHuvecEzh239875AlnRep2.markdup.q10.sorted.bam'))

chip.h3k4me3 <- c(paste0(dirchip,'Encode_h1hesc/H3K4me3/wgEncodeBroadHistoneH1hescH3k4me3StdAlnRep1.markdup.q10.sorted.bam'),
                  paste0(dirchip,'Encode_h1hesc/H3K4me3/wgEncodeBroadHistoneH1hescH3k4me3StdAlnRep2.markdup.q10.sorted.bam'),
                  paste0(dirchip,'Encode_helas3/H3K4me3/wgEncodeBroadHistoneHelas3H3k4me3StdAlnRep1.markdup.q10.sorted.bam'),
                  paste0(dirchip,'Encode_helas3/H3K4me3/wgEncodeBroadHistoneHelas3H3k4me3StdAlnRep2.markdup.q10.sorted.bam'),
                  paste0(dirchip,'Encode_hepg2/H3K4me3/wgEncodeBroadHistoneHepg2H3k4me3StdAlnRep1.markdup.q10.sorted.bam'),
                  paste0(dirchip,'Encode_hepg2/H3K4me3/wgEncodeBroadHistoneHepg2H3k4me3StdAlnRep2.markdup.q10.sorted.bam'),
                  paste0(dirchip,'Encode_huvec/H3K4me3/wgEncodeBroadHistoneHuvecH3k4me3StdAlnRep1.markdup.q10.sorted.bam'),
                  paste0(dirchip,'Encode_huvec/H3K4me3/wgEncodeBroadHistoneHuvecH3k4me3StdAlnRep2.markdup.q10.sorted.bam'))

chip.h3k27ac <- c(paste0(dirchip,'Encode_h1hesc/H3K27ac/wgEncodeBroadHistoneH1hescH3k27acStdAlnRep1.markdup.q10.sorted.bam'),
                  paste0(dirchip,'Encode_h1hesc/H3K27ac/wgEncodeBroadHistoneH1hescH3k27acStdAlnRep2.markdup.q10.sorted.bam'),
                  paste0(dirchip,'Encode_helas3/H3K27ac/wgEncodeBroadHistoneHelas3H3k27acStdAlnRep1.markdup.q10.sorted.bam'),
                  paste0(dirchip,'Encode_helas3/H3K27ac/wgEncodeBroadHistoneHelas3H3k27acStdAlnRep2.markdup.q10.sorted.bam'),
                  paste0(dirchip,'Encode_hepg2/H3K27ac/wgEncodeBroadHistoneHepg2H3k27acStdAlnRep1.markdup.q10.sorted.bam'),
                  paste0(dirchip,'Encode_hepg2/H3K27ac/wgEncodeBroadHistoneHepg2H3k27acStdAlnRep2.markdup.q10.sorted.bam'),
                  paste0(dirchip,'Encode_huvec/H3K27ac/wgEncodeBroadHistoneHuvecH3k27acStdAlnRep1.markdup.q10.sorted.bam'),
                  paste0(dirchip,'Encode_huvec/H3K27ac/wgEncodeBroadHistoneHuvecH3k27acStdAlnRep2.markdup.q10.sorted.bam'))

chip.h3k36me3 <- c(paste0(dirchip,'Encode_h1hesc/H3K36me3/wgEncodeBroadHistoneH1hescH3k36me3StdAlnRep1.markdup.q10.sorted.bam'),
                   paste0(dirchip,'Encode_h1hesc/H3K36me3/wgEncodeBroadHistoneH1hescH3k36me3StdAlnRep2.markdup.q10.sorted.bam'),
                   paste0(dirchip,'Encode_helas3/H3K36me3/wgEncodeBroadHistoneHelas3H3k36me3StdAlnRep1.markdup.q10.sorted.bam'),
                   paste0(dirchip,'Encode_helas3/H3K36me3/wgEncodeBroadHistoneHelas3H3k36me3StdAlnRep2.markdup.q10.sorted.bam'),
                   paste0(dirchip,'Encode_hepg2/H3K36me3/wgEncodeBroadHistoneHepg2H3k36me3StdAlnRep1.markdup.q10.sorted.bam'),
                   paste0(dirchip,'Encode_hepg2/H3K36me3/wgEncodeBroadHistoneHepg2H3k36me3StdAlnRep2.markdup.q10.sorted.bam'),
                   paste0(dirchip,'Encode_huvec/H3K36me3/wgEncodeBroadHistoneHuvecH3k36me3StdAlnRep1.markdup.q10.sorted.bam'),
                   paste0(dirchip,'Encode_huvec/H3K36me3/wgEncodeBroadHistoneHuvecH3k36me3StdAlnRep2.markdup.q10.sorted.bam'))

chip.h3k27me3 <- c(paste0(dirchip,'Encode_h1hesc/H3K27me3/wgEncodeBroadHistoneH1hescH3k27me3StdAlnRep1.markdup.q10.sorted.bam'),
                   paste0(dirchip,'Encode_h1hesc/H3K27me3/wgEncodeBroadHistoneH1hescH3k27me3StdAlnRep2.markdup.q10.sorted.bam'),
                   paste0(dirchip,'Encode_helas3/H3K27me3/wgEncodeBroadHistoneHelas3H3k27me3StdAlnRep1.markdup.q10.sorted.bam'),
                   paste0(dirchip,'Encode_helas3/H3K27me3/wgEncodeBroadHistoneHelas3H3k27me3StdAlnRep2.markdup.q10.sorted.bam'),
                   paste0(dirchip,'Encode_hepg2/H3K27me3/wgEncodeBroadHistoneHepg2H3k27me3StdAlnRep1.markdup.q10.sorted.bam'),
                   paste0(dirchip,'Encode_hepg2/H3K27me3/wgEncodeBroadHistoneHepg2H3k27me3StdAlnRep2.markdup.q10.sorted.bam'),
                   paste0(dirchip,'Encode_huvec/H3K27me3/wgEncodeBroadHistoneHuvecH3k27me3StdAlnRep1.markdup.q10.sorted.bam'),
                   paste0(dirchip,'Encode_huvec/H3K27me3/wgEncodeBroadHistoneHuvecH3k27me3StdAlnRep2.markdup.q10.sorted.bam'))

# Loading blacklisted regions
load('../../Data/hg19/human.hg19.ranges.blacklist.RData')

# Loading genes
load('ENCODE.rnaseq.raw.RData')

ENCODE.rnaseq.raw.subset <- ENCODE.rnaseq.raw[seqnames(ENCODE.rnaseq.raw)%in%chromosome]

# Mapping ChIP-seq reads
chipEzh2 = list()
chipH3K4me3 = list()
chipH3K27ac = list()
chipH3K36me3 = list()
chipH3K27me3 = list()

for(i in 1:length(chip.ezh2)){
    if(!is.na(chip.ezh2[i])){
        cat('ChIP-seq data: ',chip.ezh2[i],'\n')
        fleng <- maximizeCcf(correlateReads(chip.ezh2[i],param=readParam(discard=hg19.discard)))
        cat('Fragment length: ',fleng,'\n')
        chipEzh2[[i]] <- bamCount(bampath = chip.ezh2[i],gr = rowRanges(ENCODE.rnaseq.raw.subset),shift = fleng/2) # Mapping onto gene bodies
    } else{
        chipEzh2[[i]] <- rep(NA,length(rowRanges(ENCODE.rnaseq.raw.subset)))
    }
}

for(i in 1:length(chip.h3k4me3)){
    if(!is.na(chip.h3k4me3[i])){
        cat('ChIP-seq data: ',chip.h3k4me3[i],'\n')
        fleng <- maximizeCcf(correlateReads(chip.h3k4me3[i],param=readParam(discard=hg19.discard)))
        cat('Fragment length: ',fleng,'\n')
        chipH3K4me3[[i]] <- bamCount(bampath = chip.h3k4me3[i],gr = promoters(rowRanges(ENCODE.rnaseq.raw.subset)),shift = fleng/2) # Mapping onto gene promoters
    } else{
        chipH3K4me3[[i]] <- rep(NA,length(rowRanges(ENCODE.rnaseq.raw.subset)))
    }
}

for(i in 1:length(chip.h3k27ac)){
    if(!is.na(chip.h3k27ac[i])){
        cat('ChIP-seq data: ',chip.h3k27ac[i],'\n')
        fleng <- maximizeCcf(correlateReads(chip.h3k27ac[i],param=readParam(discard=hg19.discard)))
        cat('Fragment length: ',fleng,'\n')
        chipH3K27ac[[i]] <- bamCount(bampath = chip.h3k27ac[i],gr = promoters(rowRanges(ENCODE.rnaseq.raw.subset)),shift = fleng/2) # Mapping onto gene promoters
    } else{
        chipH3K27ac[[i]] <- rep(NA,length(rowRanges(ENCODE.rnaseq.raw.subset)))
    }
}

for(i in 1:length(chip.h3k36me3)){
    if(!is.na(chip.h3k36me3[i])){
        cat('ChIP-seq data: ',chip.h3k36me3[i],'\n')
        fleng <- maximizeCcf(correlateReads(chip.h3k36me3[i],param=readParam(discard=hg19.discard)))
        cat('Fragment length: ',fleng,'\n')
        chipH3K36me3[[i]] <- bamCount(bampath = chip.h3k36me3[i],gr = rowRanges(ENCODE.rnaseq.raw.subset),shift = fleng/2) # Mapping onto gene bodies
    } else{
        chipH3K36me3[[i]] <- rep(NA,length(rowRanges(ENCODE.rnaseq.raw.subset)))
    }
}

for(i in 1:length(chip.h3k27me3)){
    if(!is.na(chip.h3k27me3[i])){
        cat('ChIP-seq data: ',chip.h3k27me3[i],'\n')
        fleng <- maximizeCcf(correlateReads(chip.h3k27me3[i],param=readParam(discard=hg19.discard)))
        cat('Fragment length: ',fleng,'\n')
        chipH3K27me3[[i]] <- bamCount(bampath = chip.h3k27me3[i],gr = rowRanges(ENCODE.rnaseq.raw.subset),shift = fleng/2) # Mapping onto gene bodies
    } else{
        chipH3K27me3[[i]] <- rep(NA,length(rowRanges(ENCODE.rnaseq.raw.subset)))
    }
}


ENCODE.chipseq <- ENCODE.rnaseq.raw.subset

assay(ENCODE.chipseq,'EZH2',withDimnames = FALSE) <- do.call(cbind,chipEzh2)
assay(ENCODE.chipseq,'H3K4me3',withDimnames = FALSE) <- do.call(cbind,chipH3K4me3)
assay(ENCODE.chipseq,'H3K27ac',withDimnames = FALSE) <- do.call(cbind,chipH3K27ac)
assay(ENCODE.chipseq,'H3K36me3',withDimnames = FALSE) <- do.call(cbind,chipH3K36me3)
assay(ENCODE.chipseq,'H3K27me3',withDimnames = FALSE) <- do.call(cbind,chipH3K27me3)

save(ENCODE.chipseq,file = 'ENCODE.chipseq.RData',compress = 'xz')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Now, creating RData element with RefSeq gene bodies to be used when plotting
rm(list=ls())

library("AnnotationHub")

# Loading genome
load('../../Data/Encode_helas3/H3K36me3/wgEncodeBroadHistoneHelas3H3k36me3StdAlnRep1.markdup.q10.sorted.RData')

# Setting up data
counts = subset(counts[[paste0(100)]],chr%in%c(paste0('chr',1:22),'chrX','chrY'),select=c('chr','start','stop'))
gr.counts = with(counts,GRanges(chr, IRanges(start=start, end=stop)))

# Loading genes
ah <- AnnotationHub()
qhs <- query(ah, c("RefSeq", "Homo sapiens", "hg19"))
genes <- qhs[['AH5040']]
gr.genes <- genes[overlapsAny(genes,gr.counts)]
gr.genes <- gr.genes[seqnames(gr.genes)%in%c(paste0('chr',1:22),'chrX','chrY')]

# Saving genes
save(gr.genes,file='UCSC_RefSeq_Hg19_Genes.RData',compress = 'xz')


