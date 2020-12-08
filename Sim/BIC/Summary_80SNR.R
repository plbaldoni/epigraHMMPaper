rm(list=ls())

library(data.table)
library(epigraHMM)
library(ggplot2)
library(magrittr)

pct = 0.8

# Scenarios data frame
scenario = expand.grid(Marker=c('H3K36me3','H3K27me3'),Groups=c(2,3,4),Replicates=c(1,2,4),Nwindow=c(1e5),Nsim=100,Pct=round(100*pct))
scenario$Label = with(scenario,paste(Marker,paste0(Groups,'G'),paste0(Replicates,'R'),paste0(Nwindow,'W'),paste0(Pct,'SNR'),sep="_"))
scenario <- scenario[order(scenario$Groups,scenario$Replicates,scenario$Nwindow,-scenario$Pct),]
print(scenario)

# Creating output directory
summaryDir <- paste0('./Summary_',round(100*pct),'SNR')
if(!dir.exists(summaryDir)){
    dir.create(summaryDir)
}

# Checking if simulations completed
scenario$complete <- NA
for(i in seq_len(nrow(scenario))){
    scenario$numcomplete[i] <- sum(file.exists(paste0('./Output_',scenario$Pct[i],'SNR/',scenario$Label[i],'/',scenario$Label[i],'_',1:scenario$Nsim[i],'/',scenario$Label[i],'_',1:scenario$Nsim[i],'_epigraHMM.RData')))
    scenario$complete[i] <- all(file.exists(paste0('./Output_',scenario$Pct[i],'SNR/',scenario$Label[i],'/',scenario$Label[i],'_',1:scenario$Nsim[i],'/',scenario$Label[i],'_',1:scenario$Nsim[i],'_epigraHMM.RData')))
}
scenario$index <- seq(1,nrow(scenario))

for(i in seq_len(nrow(scenario))){
    output = list()
    
    # Loading output data
    if(scenario$complete[i] & !file.exists(file.path(summaryDir,paste0(scenario$Label[i],'.RData')))){
        
        # Preparing lists
        output[['BIC']] <- list()
        output[['Metrics']] <- list()
        
        outfiles = paste0('./Output_',scenario$Pct[i],'SNR/',scenario$Label[i],'/',scenario$Label[i],'_',1:scenario$Nsim[i],'/',scenario$Label[i],'_',1:scenario$Nsim[i],'_epigraHMM.RData')
        for(j in 1:length(outfiles)){
            message('Output: ',outfiles[j],'\n')
            
            # Loading dataset
            load(outfiles[j])
            
            # Output
            output[['BIC']][[j]] <- data.table::melt(data.table(scenario[i,c('Marker','Groups','Replicates','Nwindow','Pct')],Sim = j,as.data.table(lapply(bic_list,function(x){x$BIC}))),
                                                     id.vars = c('Marker','Groups','Replicates','Nwindow','Pct','Sim'),variable.name = "Pattern",value.name = "BIC")[,Optimal := ifelse(Groups == 2 & Pattern == '2',TRUE,ifelse(Groups == 3 & Pattern == '1_2/3',TRUE,FALSE))][,Seq := 1:.N]
            
            output[['Metrics']][[j]] <- data.table(scenario[i,c('Marker','Groups','Replicates','Nwindow','Pct')],Sim = j,data.table::rbindlist(lapply(seq_len(length(bic_list)),function(x){
                DT <- data.table(Sim = !(Z %in% range(Z)),Pred = NA)
                data.table(Pattern = names(bic_list)[x],data.table::rbindlist(lapply(c(0.01,0.05,0.10,0.15,0.20),function(y){
                    DT[,Pred := epigraHMM:::fdrControl(prob = bic_list[[x]]$prob,fdr = y)]
                    return(data.table(FDR = y,
                                      FPR = sum(DT[,Sim==F & Pred==T])/sum(DT[,Sim==F]),
                                      Precision = sum(DT[,Sim==T & Pred==T])/sum(DT[,Pred==T]),
                                      Recall = sum(DT[,Sim==T & Pred==T])/sum(DT[,Sim==T])))
                })))
            })))
        }
        
        # Combining output
        output[['BIC']] <- data.table::rbindlist(output[['BIC']])
        output[['Metrics']] <- data.table::rbindlist(output[['Metrics']])
        
        # Saving output
        save(output,file = file.path(summaryDir,paste0(scenario$Label[i],'.RData')))
        
        # Creating plots w/ metrics
        metrics <- output[['Metrics']][,list('Sensitivity' = mean(Recall),
                                             '1-Specificity' = mean(FPR),
                                             'Observed FDR' = mean(1-Precision)),by = c('Marker','Groups','Replicates','Nwindow','Pct','Pattern','FDR')]
        metrics$`Nominal FDR` <- factor(format(round(metrics$FDR, digits=2), nsmall = 2),levels = c('0.01','0.05','0.10','0.15','0.20'))
        metrics[Pattern == ifelse(Groups == 2,'2',ifelse(Groups == 3,'1_2/3',NA)),Pattern := paste0(Pattern,' (Optimal)')]
        
        fig.metrics <- ggplot(data = metrics,aes(x = `Observed FDR`,y = `Sensitivity`,color = Pattern)) + 
            geom_point(aes(shape = `Nominal FDR`)) +
            geom_line() + 
            theme_bw() + theme(panel.grid = element_blank())+
            ylim(0,1) + scale_x_continuous(breaks = c(0.01,0.05,0.10,0.15,0.20)) +
            geom_vline(xintercept = c(0.01,0.05,0.10,0.15,0.20),linetype='dashed',size = 0.5)
        
        ggsave(filename = file.path(summaryDir,paste0(scenario$Label[i],'_metrics.png')),
               plot = fig.metrics,
               width = 7,
               height = ifelse(scenario$Groups[i]==4,7,5),
               dpi = "retina")
        
        # ROC
        
        fig.roc <- ggplot(data = metrics,aes(x = `1-Specificity`,y = `Sensitivity`,color = Pattern)) + 
            geom_point(aes(shape = `Nominal FDR`)) +
            geom_line() + 
            theme_bw() + theme(panel.grid = element_blank())+
            ylim(0,1) + xlim(0,1)
        
        ggsave(filename = file.path(summaryDir,paste0(scenario$Label[i],'_roc.png')),
               plot = fig.roc,
               width = 7,
               height = ifelse(scenario$Groups[i]==4,7,5),
               dpi = "retina")
        
        # Creating plots with BIC
        metrics.bic <- output[['BIC']][,list(BIC_mean = mean(BIC),BIC_sd = sd(BIC)),by = c('Marker','Groups','Replicates','Nwindow','Pct','Pattern')]
        metrics.bic[,BIC_minusSD := BIC_mean -BIC_sd][,BIC_plusSD := BIC_mean + BIC_sd]
        metrics.bic[Pattern == ifelse(Groups == 2,'2',ifelse(Groups == 3,'1_2/3',NA)),Pattern := paste0(Pattern,' (Optimal)')]
        metrics.bic$Pattern %<>% factor(levels = unique(metrics.bic$Pattern))
        
        fig.metrics.bic <- ggplot(data = metrics.bic,aes(x = Pattern,y = BIC_mean))+
            xlab('Combinatorial Patterns in Mixture Model') + ylab('Average ± SD BIC (100 Simulated Datasets)')+
            geom_line(aes(group = 1)) +
            geom_vline(aes(xintercept = 2,color = '1_2/3 (Optimal)'),linetype='dotted')+
            geom_pointrange(aes(ymin = BIC_minusSD, ymax = BIC_plusSD))+
            scale_color_manual(name = "Optimal Pattern", values = c(`1_2/3 (Optimal)` = "red"))+
            theme_bw() + theme(panel.grid = element_blank(),legend.position = c(0.8,0.8),axis.text.x = element_text(angle = 30, hjust = 1))
        
        ggsave(filename = file.path(summaryDir,paste0(scenario$Label[i],'_bic.png')),
               plot = fig.metrics.bic,
               width = 7,
               height = ifelse(scenario$Groups[i]==4,7,5),
               dpi = "retina")
    }
}

