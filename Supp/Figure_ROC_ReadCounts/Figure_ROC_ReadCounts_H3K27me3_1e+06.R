library(data.table)
library(ggplot2)
library(magrittr)
library(plyr)
library(ggpubr)

mark = 'H3K27me3'
size = 1e6

### Simulation scnearios
# Marker:   H3K27me3    H3K36me3
# Groups:   2, 3, 4
# Replicates: 1, 2, 4
# Nwindow: 1e5, 5e5, 1e6
# Nsim: 100

scenario = expand.grid(Marker=c('H3K27me3','H3K36me3'),Groups=c(2,3,4),Replicates=c(1,2,4),Nwindow=c(1e5,5e5,1e6),Nsim=100,Pct=round(100*c(0.7,0.8,0.9,1)))
scenario$Label = with(scenario,paste(Marker,paste0(Groups,'G'),paste0(Replicates,'R'),paste0(Nwindow,'W'),paste0(Pct,'SNR'),sep="_"))
scenario <- scenario[order(scenario$Groups,scenario$Replicates,scenario$Nwindow,-scenario$Pct),]
# print(scenario)

### Loading FDRpp

pct = c(70,80,90,100)

fdrpp.lst = list()
for(i in pct){
  outRData <-
    paste0('../../Sim/ReadCounts/Metrics_',
           i,
           'SNR/',
           'Metrics_',
           i,
           'SNR_fdrpp.RData')
  if (file.exists(outRData)) {
    load(outRData)
    fdrpp.lst[[paste0(i)]] = fdrpp.list
  }
}
fdrpp.lst <- unlist(fdrpp.lst,recursive = F)

fdrpp <- rbindlist(lapply(fdrpp.lst,FUN = function(x){
  y = rbindlist(x,idcol = 'rn') #as.data.table(do.call(rbind,x),keep.rownames = T)
  y[,Simulation := as.numeric(gsub("Sim","",rn))][,rn := NULL]
  return(y)
}),idcol = 'rn')

label <- as.data.table(do.call(rbind,strsplit(fdrpp$rn,"_")))
setnames(label,c('Marker','Groups','Replicates','Nwindow','Pct'))
label[,Groups := as.numeric(gsub("G","",Groups))]
label[,Replicates := as.numeric(gsub("R","",Replicates))]
label[,Nwindow := as.numeric(gsub("W","",Nwindow))]
label[,Pct := as.numeric(gsub("S.*","",Pct))]
label[,Marker := (gsub(".*\\.","",Marker))]

dt <- data.table(fdrpp[,rn := NULL],label)

### Merging the results

dt$Groups %<>% mapvalues(from = 2:4,to = paste(2:4,'Conditions'))
dt$Pct %<>% mapvalues(from = c(70,80,90,100),to = c('70% SNR','80% SNR','90% SNR','Obs. SNR'))

dt.subset <- dt[Marker==mark & Nwindow==size & FDR>=0.01 & FDR<=0.2,]
dt.subset <- dt.subset[,.(FPR = mean(FPR,na.rm=T),Precision = mean(Precision,na.rm=T),Recall = mean(Recall,na.rm=T)),by = c('Marker','Groups','Replicates','Nwindow','Pct','FDR')]
dt.subset <- dt.subset[as.character(FDR)%in%paste0(c(0.01,0.05,0.10,0.15,0.20)),]

dt.subset$Pct %<>% mapvalues(from = c("70% SNR","80% SNR","90% SNR","Obs. SNR"),to = c("70% of Observed SNR","80% of Observed SNR","90% of Observed SNR","Observed SNR")) %<>% factor(levels = c("70% of Observed SNR","80% of Observed SNR","90% of Observed SNR","Observed SNR"))
dt.subset$FDR %<>% factor(levels = paste0(c(0.01,0.05,0.10,0.15,0.20))) %<>% mapvalues(from = paste0(c(0.01,0.05,0.10,0.15,0.20)),to = sprintf("%.2f",c(0.01,0.05,0.10,0.15,0.20)))
dt.subset$Replicates %<>% factor(levels = c(1,2,4)) %<>% mapvalues(from = c(1,2,4), to = c('1','2','4'))

### Plotting the results

fig_roc = ggplot(data = dt.subset,aes(x = (1-Precision),y = Recall))+
  facet_grid(Groups~Pct) +
  geom_line(aes(color = Replicates),size=1.25)+
  geom_point(aes(x = (1-Precision),y = Recall,fill = Replicates,shape = FDR),size = 2.25)+
  scale_shape_manual(values = c('0.01'=8,'0.05'=21,'0.10'=22,'0.15'=23,'0.20'=24))+
  scale_color_brewer(palette = 'Dark2') +
  scale_fill_brewer(palette = 'Dark2') +
  guides(fill = F)+
  labs(x='Observed False Discovery Rate',y = 'Observed True Positive Rate',shape = 'Nominal\nFDR',color = 'Replicates',fill = 'Replicates')+
  theme_bw()

rm(list=setdiff(ls(),c('fig_roc','mark','size')))

### Simulation scnearios
# Marker:   H3K27me3    H3K36me3
# Groups:   2, 3, 4
# Replicates: 1, 2, 4
# Nwindow: 1e5, 5e5, 1e6
# Nsim: 100

scenario = expand.grid(Marker=c('H3K27me3','H3K36me3'),Groups=c(2,3,4),Replicates=c(1,2,4),Nwindow=c(1e5,5e5,1e6),Nsim=100,Pct=round(100*c(0.7,0.8,0.9,1)))
scenario$Label = with(scenario,paste(Marker,paste0(Groups,'G'),paste0(Replicates,'R'),paste0(Nwindow,'W'),paste0(Pct,'SNR'),sep="_"))
scenario <- scenario[order(scenario$Groups,scenario$Replicates,scenario$Nwindow,-scenario$Pct),]
# print(scenario)

### Loading Predicted Pattern

pct = c(70,80,90,100)

preddiff.lst = list()
for(i in pct){
  outRData <-
    paste0('../../Sim/ReadCounts/Metrics_',
           i,
           'SNR/',
           'Metrics_',
           i,
           'SNR_preddiff.RData')
  if (file.exists(outRData)) {
    load(outRData)
    preddiff.lst[[paste0(i)]] = preddiff.list
  }
}
preddiff.lst <- unlist(preddiff.lst,recursive = F)
names(preddiff.lst) <- gsub('.*\\.','',names(preddiff.lst))

preddiff <- lapply(preddiff.lst, FUN = function(x){
  ### Adjusting cases where the model predicted fewer differential states
  newx = lapply(x,FUN = function(y){
    if(nrow(y)<ncol(y)){
      z = matrix(0,nrow = ncol(y),ncol = ncol(y))
      colnames(z) = paste0(1:ncol(y));rownames(z) = paste0(1:ncol(y))
      for(idx in rownames(y)){
        z[as.numeric(idx),] = y[idx,]
      }
      return(as.matrix(z))
    } else{
      return(unclass(y))
    }
  })
  
  newx.prob = lapply(newx,FUN = function(y){
    newy = as.matrix(y)
    newy = newy[-c(1,nrow(newy)),-c(1,ncol(newy))]
    #newy = newy/sum(newy)
    return(newy)
  })
  
  newx.prob.norm = Reduce("+",newx.prob)/length(newx.prob)
  dt = as.data.table(melt(newx.prob.norm))
  setnames(dt,c('Predicted','Simulated','Prob'))
  return(dt)
})

preddiff = rbindlist(preddiff,idcol = 'Label')

label <- as.data.table(do.call(rbind,strsplit(preddiff$Label,"_")))
setnames(label,c('Marker','Groups','Replicates','Nwindow','Pct'))
label[,Groups := as.numeric(gsub("G","",Groups))]
label[,Replicates := as.numeric(gsub("R","",Replicates))]
label[,Nwindow := as.numeric(gsub("W","",Nwindow))]
label[,Pct := as.numeric(gsub("S.*","",Pct))]
label[,Marker := (gsub(".*\\.","",Marker))]

preddiff <- data.table(preddiff[,Label := NULL],label)

### Cleaning up ls()

# rm(list=setdiff(ls(),c('preddiff')))

### Merging the results

preddiff$Groups %<>% mapvalues(from = 2:4,to = paste(2:4,'Conditions'))
preddiff$Pct %<>% mapvalues(from = c(70,80,90,100),to = c('70% SNR','80% SNR','90% SNR','Obs. SNR'))
preddiff$Predicted %<>% as.factor %<>% factor(levels = paste0(1:16))
preddiff$Simulated %<>% as.factor %<>% factor(levels = paste0(16:1))
preddiff$Replicates %<>% mapvalues(from = c(1,2,4),to = c('1 Replicate','2 Replicates','4 Replicates'))

### Plotting the results

preddiff.subset = preddiff[Marker==mark & Nwindow==size & Groups=='3 Conditions',]
preddiff.subset$PredictedFactor <- plyr::mapvalues(preddiff.subset$Predicted,from = 2:7,to = c('100','010','001','110','101','011'))
preddiff.subset$SimulatedFactor <- plyr::mapvalues(preddiff.subset$Simulated,from = 2:7,to = c('100','010','001','110','101','011'))
preddiff.subset$Pct %<>% 
  mapvalues(from = c("70% SNR","80% SNR","90% SNR","Obs. SNR"),to = c("70% of Observed SNR","80% of Observed SNR","90% of Observed SNR","Observed SNR")) %<>% 
  factor(levels = c("70% of Observed SNR","80% of Observed SNR","90% of Observed SNR","Observed SNR"))

fig_heatmap <- ggplot(data = preddiff.subset,
                      aes(x = PredictedFactor,y = SimulatedFactor,fill = Prob))+
  facet_grid(Replicates~Pct) +
  geom_tile(colour="grey",size=0.25) + 
  scale_fill_gradient(high = "black",low = "white",breaks = c(10,45000,90000),labels = c(0,45000,90000))+
  ylab('Simulated')+xlab('Classified')+
  theme_bw() + 
  labs(fill = 'Genomic Windows')+
  theme(axis.text.x=element_text(angle = 90),legend.position = 'bottom',legend.direction = 'horizontal')+
  scale_x_discrete(labels = c('EBB','BEB','BBE','EEB','EBE','BEE'))+
  scale_y_discrete(labels = c('BEE','EBE','EEB','BBE','BEB','EBB'))

fig <- ggarrange(fig_roc+theme(legend.direction = 'vertical'),
                 fig_heatmap+theme(legend.direction = 'vertical')+labs(fill = 'Genomic\nWindows'),
                 ncol=1,nrow = 2,labels = list('A','B'),common.legend = F,legend = 'right')

ggsave(fig,filename = paste0('Figure_ROC_ReadCounts_',mark,'_',size,'.pdf'),height = 11.5,width = 9)
