#This creates a figure that shows MA plots and adjusted read counts
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(reshape2)
library(tidyr)
library(magrittr)
library(plyr)
library(ggpubr)
library(edgeR)
library(epigraHMM)

bp = 500
chromosome = 'chr19'
ngroups=3
chipnames = c('h1hesc_rep1','h1hesc_rep2','helas3_rep1','helas3_rep2','hepg2_rep1','hepg2_rep2')

# Loading data

## Group 1: H1hesc
load('../../Data/Encode_h1hesc/H3K36me3/wgEncodeBroadHistoneH1hescH3k36me3StdAlnRep1.markdup.q10.sorted.RData')
h1hesc_1 = subset(counts[[paste0(bp)]],chr==chromosome)$counts;rm(counts)
load('../../Data/Encode_h1hesc/H3K36me3/wgEncodeBroadHistoneH1hescH3k36me3StdAlnRep1.markdup.q10.sorted.RData')
h1hesc_2 = subset(counts[[paste0(bp)]],chr==chromosome)$counts;rm(counts)

## Group 2: Helas3
load('../../Data/Encode_helas3/H3K36me3/wgEncodeBroadHistoneHelas3H3k36me3StdAlnRep1.markdup.q10.sorted.RData')
helas3_1 = subset(counts[[paste0(bp)]],chr==chromosome)$counts;rm(counts)
load('../../Data/Encode_helas3/H3K36me3/wgEncodeBroadHistoneHelas3H3k36me3StdAlnRep2.markdup.q10.sorted.RData')
helas3_2 = subset(counts[[paste0(bp)]],chr==chromosome)$counts;rm(counts)

## Group 2: Hepg2
load('../../Data/Encode_hepg2/H3K36me3/wgEncodeBroadHistoneHepg2H3k36me3StdAlnRep1.markdup.q10.sorted.RData')
hepg2_1 = subset(counts[[paste0(bp)]],chr==chromosome)$counts;rm(counts)
load('../../Data/Encode_hepg2/H3K36me3/wgEncodeBroadHistoneHepg2H3k36me3StdAlnRep2.markdup.q10.sorted.RData')
hepg2_2 = subset(counts[[paste0(bp)]],chr==chromosome)$counts

counts = subset(counts[[paste0(bp)]],chr==chromosome,select=c('chr','start','stop'))
gr.counts = with(counts,GRanges(chr, IRanges(start=start, end=stop)))

# Combining data
# Y = cbind(counts,Window=1:nrow(counts))
# Y = cbind(Y,H1hesc_1=h1hesc_1,H1hesc_2=h1hesc_2,Helas3_1=helas3_1,Helas3_2=helas3_2,Hepg2_1=hepg2_1,Hepg2_2=hepg2_2)
Y = cbind(H1hesc_1=h1hesc_1,H1hesc_2=h1hesc_2,Helas3_1=helas3_1,Helas3_2=helas3_2,Hepg2_1=hepg2_1,Hepg2_2=hepg2_2)
# ChIP <- ChIP %>%
#     as.tbl() %>%
#     gather(Group,Counts,H1hesc:Hepg2)
N = ncol(Y)
M = nrow(Y)

### Creating offset ###
log.Y = log(Y+1)
avg.Y = exp(rowMeans(log.Y))
log.avg.Y = log(avg.Y)

log.Y.M = log.Y - matrix(log.avg.Y,nrow=M,ncol=N,byrow = F) #M, from MA plot
log.Y.A = 0.5*(log.Y + matrix(log.avg.Y,nrow=M,ncol=N,byrow = F)) #A, from MA plot

# Calculating Loess
offset = sapply(1:N,function(i){loessFit(y=log.Y.M[,i],x=log.Y.A[,i],span=1)$fitted})

# Setting up names
colnames(log.Y.M) <- paste0('M_',colnames(log.Y.M))
colnames(log.Y.A) <- paste0('A_',colnames(log.Y.A))
colnames(offset) <- gsub('A_','O_',colnames(log.Y.A))

# Plotting MA plots
tb <- as_tibble(log.Y.M) %>%
  bind_cols(as_tibble(log.Y.A)) %>%
  bind_cols(as_tibble(offset))

h1hesc = ggplot(data = tb,aes(y=M_H1hesc_1,x=A_H1hesc_1))+
  stat_density2d(aes(fill = ..density..^0.25),geom = "tile", contour = FALSE,n=200)+
  scale_fill_continuous(low = "white", high = "dodgerblue4")+
  theme_bw() + theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  ylab('Minus')+xlab('Average')+
  scale_x_continuous(limits=c(-0.35,4.25))+
  scale_y_continuous(limits=c(-2.75,2.75))+
  geom_line(aes(x=A_H1hesc_1,y=O_H1hesc_1),color='blue')+
  geom_hline(yintercept=0,color='red')+
  ggtitle('H1hesc')+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size=rel(0.75)),axis.title.x = element_text(size=rel(0.75)),
        axis.text.y = element_text(size=rel(0.75)),axis.title.y = element_text(size=rel(0.75)),
        plot.title = element_text(size=rel(0.75)))

helas3 = ggplot(data = tb,aes(y=M_Helas3_1,x=A_Helas3_1))+
  stat_density2d(aes(fill = ..density..^0.25),geom = "tile", contour = FALSE,n=200)+
  scale_fill_continuous(low = "white", high = "dodgerblue4")+
  theme_bw() + theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab('Average')+
  scale_x_continuous(limits=c(-0.35,4.25))+
  scale_y_continuous(limits=c(-2.75,2.75))+
  geom_line(aes(x=A_Helas3_1,y=O_Helas3_1),color='blue')+
  geom_hline(yintercept=0,color='red')+
  theme(axis.title.y = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  ggtitle('Helas3')+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size=rel(0.75)),axis.title.x = element_text(size=rel(0.75)),
        plot.title = element_text(size=rel(0.75)))

hepg2 = ggplot(data = tb,aes(y=M_Hepg2_1,x=A_Hepg2_1))+
  stat_density2d(aes(fill = ..density..^0.25),geom = "tile", contour = FALSE,n=200)+
  scale_fill_continuous(low = "white", high = "dodgerblue4")+
  theme_bw() + theme(legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  xlab('Average')+
  scale_x_continuous(limits=c(-0.35,4.25))+
  scale_y_continuous(limits=c(-2.75,2.75))+
  geom_line(aes(x=A_Hepg2_1,y=O_Hepg2_1),color='blue')+
  geom_hline(yintercept=0,color='red')+
  theme(axis.title.y = element_blank(),axis.line.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())+
  ggtitle('Hepg2')+theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x = element_text(size=rel(0.75)),axis.title.x = element_text(size=rel(0.75)),
        plot.title = element_text(size=rel(0.75)))

maplot = ggarrange(h1hesc,helas3,hepg2,ncol = 3,nrow=1)

### Now, plotting transformed read counts
idx = which.min(abs(counts$start-56642201)):which.min(abs(counts$start-56947948)) #103550:106200

ChIP = cbind(counts,Window=1:nrow(counts))
ChIP = cbind(ChIP,H1hesc=h1hesc_1)
ChIP = cbind(ChIP,Helas3=helas3_1)
ChIP = cbind(ChIP,Hepg2=hepg2_1)
ChIP <- ChIP %>%
  as_tibble() %>%
  gather(Group,Counts,H1hesc:Hepg2) %>%
  bind_cols(Offset=c(offset[,c(1,3,5)]))

h3k36me3 = ggplot(data=ChIP[ChIP$Window%in%idx,],aes(x=start,y=Counts/exp(Offset)))+geom_line()+facet_grid(cols=vars(Group))+
  theme_bw()+xlab('Genomic Window (chr19)')+ylab('Offset-adjusted ChIP Counts')+
  scale_x_continuous(limits = range(ChIP[ChIP$Window%in%idx,'start']))+
  scale_y_continuous(limits = with(ChIP[ChIP$Window%in%idx,],range(Counts)))+
  theme(strip.background = element_blank(),strip.text.x = element_blank())+
  theme(axis.text.x = element_text(size=rel(0.75)),axis.title.x = element_text(size=rel(0.75)),
        axis.text.y = element_text(size=rel(0.75)),axis.title.y = element_text(size=rel(0.75)))

h3k36me3.unadj = ggplot(data=ChIP[ChIP$Window%in%idx,],aes(x=start,y=Counts))+geom_line()+facet_grid(cols=vars(Group))+
  theme_bw()+xlab('Genomic Window (chr19)')+ylab('ChIP Counts')+
  scale_x_continuous(limits = range(ChIP[ChIP$Window%in%idx,'start']))+
  scale_y_continuous(limits = with(ChIP[ChIP$Window%in%idx,],range(Counts)))+theme(axis.text.x = element_text(size=rel(0.75)),axis.title.x = element_text(size=rel(0.75)),
                                                                                   axis.text.y = element_text(size=rel(0.75)),axis.title.y = element_text(size=rel(0.75)),
                                                                                   strip.text.x = element_text(size=rel(0.75)))+
  theme(axis.title.x = element_blank(),axis.line.x = element_blank(),axis.ticks.x = element_blank(),axis.text.x = element_blank())

### Putting everything together

fig = ggarrange(maplot,h3k36me3.unadj,h3k36me3,ncol=1,nrow=3)

ggsave(fig,filename = 'Figure_Normalization.pdf',height = 11*0.4,width = 8.5*0.75)
