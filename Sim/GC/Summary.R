rm(list=ls())
workdir <- getwd()

library(data.table)
library(ggplot2)
library(magrittr)
library(plyr)
library(epigraHMM)

dirdt <- '/pine/scr/b/a/baldoni/Rashid/epigraHMMPaper/Sim/GC/Analysis'
datdir <- '/pine/scr/b/a/baldoni/Rashid/epigraHMMPaper/Sim/GC/Data'
dirout <- 'Summary'
idx = 7500:9500
fdr <- c(0.01,0.05,0.10,0.20)
dt.plot <- list()
system(paste('mkdir',dirout))

dirs <- list.dirs(path = dirdt,full.names = FALSE,recursive = FALSE)

res <- list()
res.fdr <- list()

for(i in dirs){
  subres <- list()
  subres.fdr <- list()
  if(length(list.files(file.path(dirdt,i),'*.RData')) == 100){
    for(sim in 1:100){
      message(i,', simulation ',sim)
      
      tryCatch({
        load(file.path(dirdt,i,paste0('Simulation',sim,'.RData')))
        
        dt[,Window := 1:.N]
        dt.viterbi <- melt(dt,id.vars = c('Window','GoldStandard'),
                           measure.vars = c('Viterbi.naive','Viterbi.gc','Viterbi.norm','Viterbi.gc.norm'),
                           variable.name = 'Method',value.name = 'Call')
        subres[[sim]] <- dt.viterbi[,.(Sensitivity = sum(Call=='D' & GoldStandard)/sum(GoldStandard),
                                       Specificity = sum(!Call=='D' & !GoldStandard)/sum(!GoldStandard)),by = 'Method'][,Scenario := i][,Simulation := sim]
        
        for(cutoff in fdr){
          dt.fdr <- dt[,c('Window','GoldStandard','PP.naive','PP.gc','PP.norm','PP.gc.norm')]
          dt.fdr[,PP.naive.fdr := 1*epigraHMM:::fdrControl(PP.naive,fdr = cutoff)]
          dt.fdr[,PP.gc.fdr := 1*epigraHMM:::fdrControl(PP.gc,fdr = cutoff)]
          dt.fdr[,PP.norm.fdr := 1*epigraHMM:::fdrControl(PP.norm,fdr = cutoff)]
          dt.fdr[,PP.gc.norm.fdr := 1*epigraHMM:::fdrControl(PP.gc.norm,fdr = cutoff)]
          
          dt.fdr <- melt(dt.fdr,id.vars = c('Window','GoldStandard'),measure.vars = c('PP.naive.fdr','PP.gc.fdr','PP.norm.fdr','PP.gc.norm.fdr'),
                         variable.name = 'Method',value.name = 'Call')
          subres.fdr[[paste0(sim,'-',cutoff)]] <- dt.fdr[,.(Sensitivity = sum(Call==1 & GoldStandard)/sum(GoldStandard),
                                                            Obs_FDR = sum(Call==1 & !GoldStandard)/sum(Call==1)),by = 'Method'][,Scenario := i][,Simulation := sim][,Nominal_FDR := cutoff]
        }
        
        # if(sim==1){
        #   load(file.path(datdir,i,paste0(i,"_",sim,".RData")))
        #   dt.plot[[paste0(i,"_",sim)]] <- rbindlist(list(as.data.table(as.data.frame(dat$ChIP))[,Scenario := i][,Window := 1:.N]))
        #   setnames(dt.plot[[paste0(i,"_",sim)]],c('A.1','A.2','B.1','B.2','Scenario','Window'))
        #   rm(dat)
        # 
        # }
        
      },error = function(cond){message(cond)})
    }
  }
  res[[i]] <- rbindlist(subres)
  res.fdr[[i]] <- rbindlist(subres.fdr)
}

res <- rbindlist(res)
res.fdr <- rbindlist(res.fdr)
dt.plot <- rbindlist(dt.plot)

# Figure 1

res$Method %<>% plyr::mapvalues(from = c('Viterbi.naive','Viterbi.gc','Viterbi.norm','Viterbi.gc.norm'),
                                to = c('Naive','GC-adj.','Normalization','GC-adj. + Normalization'))
res.fdr$Method %<>% plyr::mapvalues(from = c('PP.naive.fdr','PP.gc.fdr','PP.norm.fdr','PP.gc.norm.fdr'),
                                    to = c('Naive','GC-adj.','Normalization','GC-adj. + Normalization'))

for(scenario in unique(res$Scenario)){
  fig <- ggplot(data = res[Scenario == scenario,],aes(x = 1-Specificity,y = Sensitivity,color = Method,group = Simulation))+
    geom_line(col = 'gray',linetype = 'dashed',alpha = 0.5)+
    geom_point(alpha = 0.5)+
    theme_bw() +
    scale_color_brewer(palette = 'Set1')+
    theme(legend.position = 'right',legend.direction = 'vertical')+
    guides(color = guide_legend(title="Offset Method"))+
    labs(title = paste0('Scenario: ',scenario))
  ggsave(fig,filename = file.path(dirout,paste0('Summary_ROC_',scenario,'.png')),height = 7,width = 9)
}

for(scenario in unique(res.fdr$Scenario)){
  summarized.roc.ided <- res.fdr[Scenario == scenario,][,.(Sensitivity = mean(Sensitivity,na.rm = T),Obs_FDR = mean(Obs_FDR,na.rm = T)),by = c('Method','Nominal_FDR')]
  summarized.roc.ided$Nominal_FDR %<>% as.character()
  
  fig.fdr <- ggplot(data = summarized.roc.ided,aes(x = Obs_FDR,y = Sensitivity))+
    geom_line(aes(color = Method),size=1.25)+
    geom_point(aes(x = Obs_FDR,y = Sensitivity,fill = Method,shape = Nominal_FDR),size = 2.25)+
    scale_shape_manual(values = c('0.01'=8,'0.05'=21,'0.1'=22,'0.15'=23,'0.2'=24))+
    guides(fill = F)+
    labs(x='Observed False Discovery Rate',y = 'Observed True Positive Rate',shape = 'Nominal\nFDR',title = paste0('Scenario: ',scenario))+
    theme_bw()
  
  ggsave(fig.fdr,filename = file.path(dirout,paste0('Summary_FDR_',scenario,'.png')),height = 7,width = 9)
}

# Figure 2

for(scenario in dirs){
  
  load(file.path(datdir,scenario,paste0(scenario,"_",1,".RData")))
  dt.plot <- rbindlist(list(as.data.table(as.data.frame(dat$ChIP))[,Scenario := scenario][,Window := 1:.N][,Z := as.numeric(dat$z)][,GC := as.numeric(dat$gc)]))
  setnames(dt.plot,c('A.1','A.2','B.1','B.2','Scenario','Window','Z','GC'))
  rm(dat)
  
  dt.plot.melted <- melt(dt.plot,id.vars = c('Scenario','Window','Z','GC'),measure.vars = c('A.1','A.2','B.1','B.2'),variable.name = 'Sample',value.name = 'Counts')
  dt.plot.melted[,Replicate := substr(x = Sample,3,3)][,Condition := substr(x = Sample,1,1)]
  
  dt.plot.melted[,Zdeconv := ifelse((Condition == 'A' & (Z %in% c(2,4))) | (Condition == 'B' & (Z %in% c(3,4))),'Enrichment','Background')]

  fig.counts <- ggplot(data = dt.plot.melted[Window %in% idx,],aes(x = Window,y = Counts))+
    facet_grid(rows = vars(Replicate),cols = vars(Condition))+
    geom_line()+
    ylab('ChIP Counts')+xlab('Genomic Window')+
    labs(title = paste0('Scenario: ',scenario))+
    theme_bw()
  
  fig.lm <- ggplot(data = dt.plot.melted,aes(x = GC,y = log1p(Counts),color = Zdeconv))+
    facet_grid(rows = vars(Replicate),cols = vars(Condition))+
    geom_point(alpha = 0.25,size = 1) + 
    geom_smooth(se = FALSE, method = 'gam')+
    ylab('log(ChIP + 1)')+xlab('GC')+
    theme_bw()+
    guides(color = guide_legend(title="Latent HMM state"))+
    scale_color_brewer(palette = 'Set1') +
    labs(title = paste0('Scenario: ',scenario))+
    theme(legend.position = 'right',legend.direction = 'vertical')
  
  figAB <- ggpubr::ggarrange(fig.counts,fig.lm,nrow = 2,ncol = 1,common.legend = FALSE,labels = list('A','B'))
  ggsave(figAB,filename = file.path(dirout,paste0('Summary_Scatter_',scenario,'.png')),height = 9,width = 7)
}
