rm(list=ls())
workdir <- getwd()

library(data.table)
library(ggplot2)
library(magrittr)
library(plyr)

dirdt <- '/pine/scr/b/a/baldoni/Rashid/Project2/Simulation_Control/Analysis'
datdir <- '/pine/scr/b/a/baldoni/Rashid/Project2/Simulation_Control/Data'
dirout <- 'Summary'
idx = 7500:9500
dt.plot <- list()
system(paste('mkdir',dirout))

dirs <- list.dirs(path = dirdt,full.names = FALSE,recursive = FALSE)

res <- list()

for(i in dirs){
    subres <- list()
    for(sim in 1:100){
        message(i,', simulation ',sim)
        
        tryCatch({
            load(file.path(dirdt,i,paste0('Simulation',sim,'.RData')))
            
            dt <- as.data.table(dt)[,Window := 1:.N]
            dt <- melt(dt,id.vars = c('Window','GoldStandard'),measure.vars = c('WithoutControl','WithControl','WithGroupedControl'),variable.name = 'Method',value.name = 'Call')
            subres[[sim]] <- dt[,.(Sensitivity = sum(Call=='D' & GoldStandard==1)/sum(GoldStandard==1),
                                   Specificity = sum(!Call=='D' & !GoldStandard==1)/sum(!GoldStandard==1)),by = 'Method'][,Scenario := i][,Simulation := sim]
            
            if(sim==1){
                load(file.path(datdir,i,paste0(i,"_",sim,".RData")))
                dat1 <- dat
                rm(dat)
                load(file.path(datdir,i,paste0(i,"_",sim+100,".RData")))
                dat2 <- dat
                rm(dat)
                
                dt.plot[[paste0(i,"_",sim)]] <- rbindlist(list(as.data.table(dat1[,c('y.1','y.2','z','x.1','x.2')])[,Scenario := i][,Condition := 'A'][,Window := 1:nrow(dat1)],
                                                               as.data.table(dat2[,c('y.1','y.2','z','x.1','x.2')])[,Scenario := i][,Condition := 'B'][,Window := 1:nrow(dat2)]))
            }
            
        },error = function(cond){message(cond)})
    }
    res[[i]] <- rbindlist(subres)
}

res <- rbindlist(res)
dt.plot <- rbindlist(dt.plot)

# Figure 1

res$Method %<>% plyr::mapvalues(from = c('WithoutControl','WithControl','WithGroupedControl'),
                                to = c('Without Control','With Control','With Grouped Control\n(Grouped on Latent HMM State)'))

fig <- ggplot(data = res,aes(x = 1-Specificity,y = Sensitivity,color = Method,group = Simulation))+
    geom_line(col = 'gray',linetype = 'dashed',alpha = 0.5)+
    geom_point(alpha = 0.5)+
    theme_bw() +
    scale_color_brewer(palette = 'Set1')+
    theme(legend.position = 'right',legend.direction = 'vertical')+
    guides(color = guide_legend(title="Offset Method"))

# Figure 2

dt.plot.melted <- melt(as.data.table(dt.plot)[,Window :=  rep(1:nrow(dat1),2)],
                       id.vars = c('Condition','Window'),measure.vars = list(c('y.1','y.2'),c('x.1','x.2'),c('z','z')),
                       variable.name = 'Replicate',value.name = c('ChIP.Counts','Control.Counts','Z'))

dt.plot.melted$Replicate %<>% plyr::mapvalues(from = c(1,2),to = paste0('Replicate ',1:2))
dt.plot.melted$Condition %<>% plyr::mapvalues(from = c('A','B'),to = c('Condition A','Condition B'))
dt.plot.melted$Z %<>% plyr::mapvalues(from = c(0,1),to = c('Background','Enrichment'))


fig.counts <- ggplot(data = dt.plot.melted[Window %in% idx,],aes(x = Window,y = ChIP.Counts))+
    facet_grid(rows = vars(Replicate),cols = vars(Condition))+
    geom_line()+
    ylab('ChIP Counts')+xlab('Genomic Window')+
    theme_bw()

fig.lm <- ggplot(data = dt.plot.melted,aes(x = log(Control.Counts+1),y = log(ChIP.Counts+1),color = as.factor(Z)))+
    facet_grid(rows = vars(Replicate),cols = vars(Condition))+
    geom_point(alpha = 0.25,size = 1) + 
    geom_smooth(se = FALSE, method = lm)+
    ylab('log(ChIP + 1)')+xlab('log(Control + 1)')+
    theme_bw()+
    guides(color = guide_legend(title="Latent HMM state"))+
    scale_color_brewer(palette = 'Set1') +
    theme(legend.position = 'right',legend.direction = 'vertical')

### Puting everything together

figABC <- ggpubr::ggarrange(fig.counts,fig.lm,fig,nrow = 3,ncol = 1,common.legend = FALSE,labels = list('A','B','C'))

ggsave(figABC,filename = file.path(dirout,'Summary.pdf'),height = 9,width = 7)
ggsave(figABC,filename = file.path(dirout,'Summary.jpg'),height = 9,width = 7)

