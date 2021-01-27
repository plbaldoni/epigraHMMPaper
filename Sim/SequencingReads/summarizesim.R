library(data.table)
iters <- 1:100
results <- list()

# Loading results
for (it in iters) {
    cat(paste('Iteration:',it,'\n'))
    if (file.exists(paste0('./autosim',it,'/autosim',it,'.R'))) {
        results[[it]] <- data.table(It = it,fread(paste0('./autosim',it,'/hist',it,'_result.tsv'),na.strings = "NA"))
    }
}

# Combining results
results <- rbindlist(results)
setnames(results,c('It','Method','Cutoff','OneMinusPrecision','Recall','Time','CallSize',
                   'PeakSize','HitsPerPeak','NCalls','NPeaks','TPR','FPR','PPV','Error'))
write.table(results,file = "hist_result.tsv",quote = F,row.names = F,col.names = T,sep = '\t')
