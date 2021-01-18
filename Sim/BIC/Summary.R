rm(list=ls())

library(data.table)
library(epigraHMM)
library(ggplot2)
library(magrittr)
library(grid)
library(viridis)

# Scenarios data frame
scenario = expand.grid(Marker=c('H3K36me3','H3K27me3'),Groups=c(2,3,4),Replicates=c(1,2,4),Nwindow=c(1e5),Nsim=100,Pct=round(100*c(0.7,0.8,0.9,1)))
scenario$Label = with(scenario,paste(Marker,paste0(Groups,'G'),paste0(Replicates,'R'),paste0(Nwindow,'W'),paste0(Pct,'SNR'),sep="_"))
scenario <- scenario[order(scenario$Groups,scenario$Replicates,scenario$Nwindow,-scenario$Pct),]
print(scenario)

# Checking if simulations completed
scenario$complete <- NA
for(i in seq_len(nrow(scenario))){
  scenario$numcomplete[i] <- sum(file.exists(paste0('./Sim/BIC/Output_',scenario$Pct[i],'SNR/',scenario$Label[i],'/',scenario$Label[i],'_',1:scenario$Nsim[i],'/',scenario$Label[i],'_',1:scenario$Nsim[i],'_epigraHMM.RData')))
  scenario$complete[i] <- all(file.exists(paste0('./Sim/BIC/Output_',scenario$Pct[i],'SNR/',scenario$Label[i],'/',scenario$Label[i],'_',1:scenario$Nsim[i],'/',scenario$Label[i],'_',1:scenario$Nsim[i],'_epigraHMM.RData')))
}
scenario$index <- seq(1,nrow(scenario))

# Function to extract BIC
getBIC <- function(object){
  return(rbindlist(lapply(seq_len(length(object)),function(x){
    patname <- names(object[x])
    k <- length(strsplit(patname,'_')[[1]])
    return(data.table(Pattern = patname,K = k,BIC = object[x][[1]][['BIC']]))
  })))
}

for(i in seq_len(nrow(scenario))){
  output = list()
  
  # Creating directory
  summaryDir <- paste0('./Sim/BIC/Summary_',scenario$Pct[i],'SNR')
  if(!dir.exists(summaryDir)){
    dir.create(summaryDir,recursive = TRUE)
  }
  
  # Loading output data
  if(scenario$complete[i] & !file.exists(file.path(summaryDir,paste0(scenario$Label[i],'.RData')))){
    
    # Preparing lists
    output[['BIC']] <- list()
    output[['Metrics']] <- list()
    
    outfiles = paste0('./Sim/BIC/Output_',scenario$Pct[i],'SNR/',scenario$Label[i],'/',scenario$Label[i],'_',1:scenario$Nsim[i],'/',scenario$Label[i],'_',1:scenario$Nsim[i],'_epigraHMM.RData')
    for(j in 1:length(outfiles)){
      message('Output: ',outfiles[j],'\n')
      
      # Loading dataset
      load(outfiles[j])
      
      # Output
      
      ## BIC
      
      dtBIC <- rbindlist(list(data.table(scenario[i,c('Marker','Groups','Replicates','Nwindow','Pct')],Sim = j,Method = 'Fast',getBIC(bic_list[['fast']])),
                              data.table(scenario[i,c('Marker','Groups','Replicates','Nwindow','Pct')],Sim = j,Method = 'Full',getBIC(bic_list[['full']]))))
      
      setkey(dtBIC,'K')
      dtBIC$sortedPattern <- unlist(lapply(seq_len(nrow(dtBIC)),function(x){paste(sort(strsplit(dtBIC$Pattern[x],'_')[[1]]),collapse = '_')}))
      
      tbPatterns <- rbindlist(lapply(unique(dtBIC$K),function(x){subtest <- unique(dtBIC[K==x,'sortedPattern'])[order(sortedPattern),][,patternOrder := seq_len(.N)][,K := x]}))
      
      output[['BIC']][[j]] <- merge(dtBIC,tbPatterns,by = c('K','sortedPattern'),all.x = TRUE)
      
      # ## Metrics
      # output[['Metrics']][[j]] <- data.table(scenario[i,c('Marker','Groups','Replicates','Nwindow','Pct')],Sim = j,data.table::rbindlist(lapply(seq_len(length(bic_list)),function(x){
      #   DT <- data.table(Sim = !(Z %in% range(Z)),Pred = NA)
      #   data.table(Pattern = names(bic_list)[x],data.table::rbindlist(lapply(c(0.01,0.05,0.10,0.15,0.20),function(y){
      #     DT[,Pred := epigraHMM:::fdrControl(prob = bic_list[[x]]$prob,fdr = y)]
      #     return(data.table(FDR = y,
      #                       FPR = sum(DT[,Sim==F & Pred==T])/sum(DT[,Sim==F]),
      #                       Precision = sum(DT[,Sim==T & Pred==T])/sum(DT[,Pred==T]),
      #                       Recall = sum(DT[,Sim==T & Pred==T])/sum(DT[,Sim==T])))
      #   })))
      # })))
    }
    # Combining output
    output[['BIC']] <- data.table::rbindlist(output[['BIC']])
    # output[['Metrics']] <- data.table::rbindlist(output[['Metrics']])
    
    # Saving output
    save(output,file = file.path(summaryDir,paste0(scenario$Label[i],'.RData')))
  } else{
    if(scenario$complete[i] & file.exists(file.path(summaryDir,paste0(scenario$Label[i],'.RData')))){
      
      load(file.path(summaryDir,paste0(scenario$Label[i],'.RData'))) 
    }
  }
  
  if(scenario$complete[i] & file.exists(file.path(summaryDir,paste0(scenario$Label[i],'.RData')))){
    # # Creating plots w/ metrics
    # metrics <- output[['Metrics']][,list('Sensitivity' = mean(Recall),
    #                                      '1-Specificity' = mean(FPR),
    #                                      'Observed FDR' = mean(1-Precision)),by = c('Marker','Groups','Replicates','Nwindow','Pct','Pattern','FDR')]
    # metrics$`Nominal FDR` <- factor(format(round(metrics$FDR, digits=2), nsmall = 2),levels = c('0.01','0.05','0.10','0.15','0.20'))
    # metrics[Pattern == ifelse(Groups == 2,'2',ifelse(Groups == 3,'1_2/3',NA)),Pattern := paste0(Pattern,' (Optimal)')]
    # 
    # fig.metrics <- ggplot(data = metrics,aes(x = `Observed FDR`,y = `Sensitivity`,color = Pattern)) + 
    #   geom_point(aes(shape = `Nominal FDR`)) +
    #   geom_line() + 
    #   theme_bw() + theme(panel.grid = element_blank())+
    #   ylim(0,1) + scale_x_continuous(breaks = c(0.01,0.05,0.10,0.15,0.20)) +
    #   geom_vline(xintercept = c(0.01,0.05,0.10,0.15,0.20),linetype='dashed',size = 0.5)
    # 
    # ggsave(filename = file.path(summaryDir,paste0(scenario$Label[i],'_metrics.png')),
    #        plot = fig.metrics,
    #        width = 7,
    #        height = ifelse(scenario$Groups[i]==4,7,5),
    #        dpi = "retina")
    # 
    # # ROC
    # 
    # fig.roc <- ggplot(data = metrics,aes(x = `1-Specificity`,y = `Sensitivity`,color = Pattern)) + 
    #   geom_point(aes(shape = `Nominal FDR`)) +
    #   geom_line() + 
    #   theme_bw() + theme(panel.grid = element_blank())+
    #   ylim(0,1) + xlim(0,1)
    # 
    # ggsave(filename = file.path(summaryDir,paste0(scenario$Label[i],'_roc.png')),
    #        plot = fig.roc,
    #        width = 7,
    #        height = ifelse(scenario$Groups[i]==4,7,5),
    #        dpi = "retina")
    
    # Creating plots with BIC
    Test <- output[['BIC']][Pattern == '1_2/3' & Method == 'Fast',][,.(optimalBIC = min(BIC)),by = c('Sim')]
    Test <- merge(output[['BIC']],Test,by = 'Sim',all.x = TRUE)
    avgTest <- Test[,.(BIC = mean(BIC),diffBIC = mean(BIC - optimalBIC)),by = c('K','patternOrder','sortedPattern')]
    lineTest <- Test[Method == 'Fast' & nchar(sortedPattern) >= nchar('1'),]
    lineTest <- unique(merge(lineTest,unique(lineTest[,c('K','patternOrder')])[order(-K),patternRank := 1:.N],
                             by = c('K','patternOrder'),all.x = TRUE)[,c('K','patternOrder','patternRank')])
    setkey(lineTest,'K')
    
    fig.metrics.bic <- ggplot(avgTest) + 
      theme_classic()+theme(axis.line.y=element_blank(),
                            axis.text.y=element_blank(),
                            axis.ticks.y=element_blank())+
      geom_tile(aes(x = K, y = patternOrder, fill = BIC),color = "black", size = 0.25)+
      geom_text(aes(x = K, y = patternOrder,label = gsub('_',' ',sortedPattern)),
                size = 3,color = 'white',alpha = 0.5,fontface = "bold")+
      ylab('Combinatorial Patterns')+xlab('Number of HMM differential mixture components (L)')+
      scale_x_continuous(breaks = seq_len(max(output[['BIC']]$K)))+
      scale_fill_viridis(option = 'viridis')+
      geom_segment(inherit.aes = FALSE,data = lineTest,aes(xend = K,x=c(tail(K, n=-1), NA), yend = patternOrder,y=c(tail(patternOrder, n=-1), NA),color = 'Model Selection\nSteps'),linetype = 'dashed')+
      geom_text(inherit.aes = FALSE,data = lineTest,aes(x = K, y = patternOrder,label = patternRank),
                alpha = ifelse(lineTest$patternRank==5,1,0.9),color = 'white',size =5,fontface = "bold")+
      scale_color_manual(values = 'gray')+
      labs(fill = 'Average BIC\n(100 datasets)',color = NULL)
    
    ggsave(filename = file.path(summaryDir,paste0(scenario$Label[i],'_bic.png')),
           plot = fig.metrics.bic,
           width = 9,
           height = 5,
           dpi = "retina")
    
    fig.metrics.diffBic <- ggplot(avgTest) + 
      theme_classic()+theme(axis.line.y=element_blank(),
                            axis.text.y=element_blank(),
                            axis.ticks.y=element_blank())+
      geom_tile(aes(x = K, y = patternOrder, fill = diffBIC),color = "black", size = 0.25)+
      geom_text(aes(x = K, y = patternOrder,label = gsub('_',' ',sortedPattern)),
                size = 3,color = 'black',alpha = 0.5,fontface = "bold")+
      ylab('Combinatorial Patterns')+xlab('Number of HMM differential mixture components (L)')+
      scale_x_continuous(breaks = seq_len(max(output[['BIC']]$K)))+
      # scale_fill_distiller(palette = "RdBu")+
      scale_fill_gradient2()+
      geom_segment(inherit.aes = FALSE,data = lineTest,aes(xend = K,x=c(tail(K, n=-1), NA), yend = patternOrder,y=c(tail(patternOrder, n=-1), NA),color = 'Model Selection\nSteps'),linetype = 'dashed')+
      geom_text(inherit.aes = FALSE,data = lineTest,aes(x = K, y = patternOrder,label = patternRank),
                alpha = ifelse(lineTest$patternRank==5,1,0.9),color = 'black',size =5,fontface = "bold")+
      scale_color_manual(values = 'gray')+
      labs(fill = 'Average BIC\nDifference From\nOptimal Model\n(100 datasets)',color = NULL)
    
    ggsave(filename = file.path(summaryDir,paste0(scenario$Label[i],'_diffbic.png')),
           plot = fig.metrics.diffBic,
           width = 9,
           height = 5,
           dpi = "retina")
  }
}

