library(data.table)
library(ggplot2)
library(Polychrome)
library(magrittr)
library(forcats)
library(plyr)
library(GenomicRanges)
library(ggpubr)

size = 1.25
cutoff <- 0.05
methods <- c("epigraHMM","ChIPComp + HOMER","csaw", "DiffBind + HOMER","diffReps",'RSEG','THOR')
heights = 1+seq(0.15,by = 0.15,length.out = 7)
names(heights) = c(methods)
colors = c(kelly.colors(22)[c('lightblue','red','yellow','purplishpink','yellowgreen','buff','orange')])
names(colors) <- methods 

gr.true = with(read.table('../../Sim/SequencingReads/autosim1//hist1_log.txt',header=T),GRanges(chr,IRanges(start,end)))

gr.tile = with(fread('../../Sim/SequencingReads/autosim1/results/epigraHMM_ranges.txt'),GRanges(seqnames,IRanges(start,end)))

dt.tile = data.table(chr='chrA',start=start(gr.tile),
                     stop=end(gr.tile),
                     window = 1:length(start(gr.tile)),
                     counts1=bamsignals::bamCount(bampath='../../Sim/SequencingReads/autosim1/hist1_out_1.bam',gr=gr.tile,verbose=T,shift=0),
                     counts2=bamsignals::bamCount(bampath='../../Sim/SequencingReads/autosim1/hist1_out_2.bam',gr=gr.tile,verbose=T,shift=0),
                     counts3=bamsignals::bamCount(bampath='../../Sim/SequencingReads/autosim1/hist1_out_3.bam',gr=gr.tile,verbose=T,shift=0),
                     counts4=bamsignals::bamCount(bampath='../../Sim/SequencingReads/autosim1/hist1_out_4.bam',gr=gr.tile,verbose=T,shift=0),
                     control1=bamsignals::bamCount(bampath='../../Sim/SequencingReads/autosim1/control1_out_1.bam',gr=gr.tile,verbose=T,shift=0),
                     control2=bamsignals::bamCount(bampath='../../Sim/SequencingReads/autosim1/control1_out_2.bam',gr=gr.tile,verbose=T,shift=0),
                     control3=bamsignals::bamCount(bampath='../../Sim/SequencingReads/autosim1/control1_out_3.bam',gr=gr.tile,verbose=T,shift=0),
                     control4=bamsignals::bamCount(bampath='../../Sim/SequencingReads/autosim1/control1_out_4.bam',gr=gr.tile,verbose=T,shift=0))

dt.tile[,Simulated := overlapsAny(gr.tile,gr.true)]

dt.tile[,countsA := counts1+counts2]
dt.tile[,countsB := counts3+counts4]

dt.melted.tile <- data.table::melt(dt.tile,id.vars = c('start','stop','window'),measure.vars = c('countsA','countsB'),variable.name = 'condition',value.name = 'counts')
dt.melted.tile$condition %<>% mapvalues(from = c('countsA','countsB'),to = c('Condition A','Condition B'))


plotfunc = function(j){
  idx = (dt.tile$window[which.min(abs(124800000+j*150000-dt.tile$start))]-0):(dt.tile$window[which.min(abs(124950000+j*150000-dt.tile$stop))]+0) #j = 0 gives the paper results
  
  # Loading peaks
  gr.ChIPComp <- GenomicRanges::reduce(with(fread('../../Sim/SequencingReads/autosim1/results/ChIPComp_table.txt')[FDR<=cutoff,],GRanges(chr,IRanges(start,end))))
  ChIPComp = cbind(dt.tile[,c('chr','start','stop','window')],Group='Condition A',Method='ChIPComp + HOMER',Output=1*overlapsAny(gr.tile,gr.ChIPComp))
  ChIPComp[,Counts := ifelse(ChIPComp$Output==1,max(dt.tile[window%in%idx,countsA])*heights['ChIPComp + HOMER'],NA)]
  ChIPComp$Output %<>% mapvalues(from=0:1,to=c('Non-differential','Differential'))
  
  gr.csaw<- GenomicRanges::reduce(with(fread('../../Sim/SequencingReads/autosim1/results/csaw_ranges.txt')[fread('../../Sim/SequencingReads/autosim1/results/csaw_table.txt')$FDR<=cutoff,],GRanges(seqnames,IRanges(start,end))))
  csaw = cbind(dt.tile[,c('chr','start','stop','window')],Group='Condition A',Method='csaw',Output=1*overlapsAny(gr.tile,gr.csaw))
  csaw[,Counts := ifelse(csaw$Output==1,max(dt.tile[window%in%idx,countsA])*heights['csaw'],NA)]
  csaw$Output %<>% mapvalues(from=0:1,to=c('Non-differential','Differential'))
  
  gr.DiffBind<- GenomicRanges::reduce(with(fread('../../Sim/SequencingReads/autosim1/results/DiffBind_ranges.txt')[FDR<=cutoff,],GRanges(seqnames,IRanges(start,end))))
  DiffBind = cbind(dt.tile[,c('chr','start','stop','window')],Group='Condition A',Method='DiffBind + HOMER',Output=1*overlapsAny(gr.tile,gr.DiffBind))
  DiffBind[,Counts := ifelse(DiffBind$Output==1,max(dt.tile[window%in%idx,countsA])*heights['DiffBind + HOMER'],NA)]
  DiffBind$Output %<>% mapvalues(from=0:1,to=c('Non-differential','Differential'))
  
  gr.diffReps <- GenomicRanges::reduce(with(fread('../../Sim/SequencingReads/autosim1/results/diffReps_table.txt')[FDR<=cutoff,],GRanges(Chrom,IRanges(Start,End))))
  diffReps = cbind(dt.tile[,c('chr','start','stop','window')],Group='Condition A',Method='diffReps',Output=1*overlapsAny(gr.tile,gr.diffReps))
  diffReps[,Counts := ifelse(diffReps$Output==1,max(dt.tile[window%in%idx,countsA])*heights['diffReps'],NA)]
  diffReps$Output %<>% mapvalues(from=0:1,to=c('Non-differential','Differential'))
  
  gr.epigraHMM <- GenomicRanges::reduce(with(fread('../../Sim/SequencingReads/autosim1/results/epigraHMM_table.txt')[epigraHMM:::fdrControl(prob = Postprob,fdr = cutoff),],GRanges(seqnames,IRanges(start,end))))
  epigraHMM = cbind(dt.tile[,c('chr','start','stop','window')],Group='Condition A',Method='epigraHMM',Output=1*overlapsAny(gr.tile,gr.epigraHMM))
  epigraHMM[,Counts := ifelse(epigraHMM$Output==1,max(dt.tile[window%in%idx,countsA])*heights['epigraHMM'],NA)]
  epigraHMM$Output %<>% mapvalues(from=0:1,to=c('Non-differential','Differential'))
  
  gr.RSEG <- GenomicRanges::reduce(with(fread('../../Sim/SequencingReads/autosim1/results/RSEG_table.txt')[epigraHMM:::fdrControl(prob = Postprob,fdr = cutoff),],GRanges(V1,IRanges(V2,V3))))
  RSEG = cbind(dt.tile[,c('chr','start','stop','window')],Group='Condition A',Method='RSEG',Output=1*overlapsAny(gr.tile,gr.RSEG))
  RSEG[,Counts := ifelse(RSEG$Output==1,max(dt.tile[window%in%idx,countsA])*heights['RSEG'],NA)]
  RSEG$Output %<>% mapvalues(from=0:1,to=c('Non-differential','Differential'))
  
  gr.THOR <- GenomicRanges::reduce(with(fread('../../Sim/SequencingReads/autosim1/results/THOR_table.txt')[FDR<=cutoff,],GRanges(V1,IRanges(V2,V3))))
  THOR = cbind(dt.tile[,c('chr','start','stop','window')],Group='Condition A',Method='THOR',Output=1*overlapsAny(gr.tile,gr.THOR))
  THOR[,Counts := ifelse(THOR$Output==1,max(dt.tile[window%in%idx,countsA])*heights['THOR'],NA)]
  THOR$Output %<>% mapvalues(from=0:1,to=c('Non-differential','Differential'))
  
  dt.peaks = rbindlist(list(ChIPComp[window%in%idx,],
                            csaw[window%in%idx,],
                            DiffBind[window%in%idx,],
                            diffReps[window%in%idx,],
                            epigraHMM[window%in%idx,],
                            RSEG[window%in%idx,],
                            THOR[window%in%idx,]))
  
  dt.peaks$Track = dt.peaks$Counts
  
  dt.melted.tile$Group = dt.melted.tile$condition
  dt.melted.tile$ChIP = dt.melted.tile$counts
  
  dt.all = merge(dt.melted.tile[window%in%idx,c('start','stop','Group','ChIP')],
                 dt.peaks[,c('start','stop','Group','Method','Track')],by=c('start','stop','Group'),all.x=T)
  
  figD = ggplot(data = dt.all,aes(x = start,y = ChIP))+
    facet_wrap(~Group,nrow=2,ncol=1)+
    geom_line()+
    geom_segment(aes(x=start,xend=stop,y=Track,yend=Track,color=Method),size=size)+
    scale_color_manual(values = colors,breaks = names(colors))+
    scale_x_continuous(labels = scales::comma)+
    labs(y='Total ChIP counts',x='Genomic Position')+
    theme_bw()
  
  return(figD)
}


fig1 <- plotfunc(j = 5)
fig2 <- plotfunc(j = 10)
fig3 <- plotfunc(j = -5)
fig4 <- plotfunc(j = -60)

fig <- ggarrange(fig1,fig2,fig3,fig4,nrow = 2,ncol = 2,common.legend = TRUE,legend = 'top')

ggsave(fig,filename = 'Figure_SequencingReads_Example.pdf',height = 11.5,width = 9)
