library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)

#Chromosome sizes
chr = paste0("chr", c(1:22, "X", "Y"))
hg19.full = GRanges(chr,IRanges(0,seqlengths(Hsapiens)[chr]))

#Defining blacklist and gap track
hg19.blacklist.dac = import('ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz')
hg19.blacklist.duke = import('ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityRegionsExcludable.bed.gz')
hg19.gap = with(read.table(textConnection(readLines(gzcon(url("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gap.txt.gz"))))),GRanges(V2,IRanges(V3,V4)))

hg19.blacklist.dac = hg19.blacklist.dac[seqnames(hg19.blacklist.dac)%in%chr]
hg19.blacklist.duke = hg19.blacklist.duke[seqnames(hg19.blacklist.duke)%in%chr]
hg19.gap = hg19.gap[seqnames(hg19.gap)%in%chr]

hg19.discard = suppressWarnings(union(hg19.gap,union(hg19.blacklist.dac,hg19.blacklist.duke)))
start(hg19.discard)[start(hg19.discard)==0] <- 1

#Defining the mappable hg19
hg19 = suppressWarnings(setdiff(hg19.full,hg19.discard))

#Saving .RData
system('mkdir hg19')
save(hg19,file='./hg19/human.hg19.ranges.mappable.RData')
save(hg19.discard,file='./hg19/human.hg19.ranges.blacklist.RData')
