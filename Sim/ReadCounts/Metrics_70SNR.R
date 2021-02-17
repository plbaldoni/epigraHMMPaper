rm(list=ls())

library(data.table)

pct = 0.7

### Simulation scnearios
# Marker:   H3K27me3    H3K36me3
# Groups:   2, 3, 4
# Replicates: 1, 2, 4
# Nwindow: 1e5, 5e5, 1e6
# Nsim: 100

scenario = expand.grid(Marker=c('H3K27me3','H3K36me3'),Groups=c(2,3,4),Replicates=c(1,2,4),Nwindow=c(1e5,5e5,1e6),Nsim=100,Pct=round(100*pct))
scenario$Label = with(scenario,paste(Marker,paste0(Groups,'G'),paste0(Replicates,'R'),paste0(Nwindow,'W'),paste0(Pct,'SNR'),sep="_"))
scenario <- scenario[order(scenario$Groups,scenario$Replicates,scenario$Nwindow,-scenario$Pct),]
print(scenario)

acc.list = list()
prec.list = list()
rec.list = list()
f1.list = list()
kappa.list = list()
par.list = list()
preddiff.list = list()
bias.list = list()
fdrpp.list = list()

genpar = function(Marker,Groups,pct){
    ## Initial probability (always begin with background) and Transition Matrix
    pi = c(1,0)
    ## Number of states
    N.states = 2^Groups
    if(Marker=='H3K36me3'){
        ## See 01_H3K36me3_Parameters.R
        ## Transition Probability Matrix
        d1 = 1-1/(131.92967+1)  
        d2 = 1-1/(43.75663+1) 
        gamma = matrix(0,N.states,N.states)
        diag(gamma) = c(d1,rep(d2,N.states-1))
        gamma[1,2:N.states] = (1-d1)/(N.states-1)
        gamma[gamma==0] = (1-d2)/(N.states-1)
        
        ## ChIP parameters
        par.background = c(0.7560612,0,3.230474)
        par.enrichment = c(2.4169519,0,2.891908)
        
        ## Offsets
        #I am using the same offset for replicates from the same group
        offsets = NULL
        offsets[paste0('Group',1:Groups)] = list(NULL)
        offsets[paste0('Group',1:Groups)] = rep(0,Groups) 
    }
    if(Marker=='H3K27me3'){
        ## See 02_H3K27me3_Parameters.R
        ## Transition Probability Matrix
        d1 = 1-1/(143.57818+1) 
        d2 = 1-1/(49.47422+1)
        gamma = matrix(0,N.states,N.states)
        diag(gamma) = c(d1,rep(d2,N.states-1))
        gamma[1,2:N.states] = (1-d1)/(N.states-1)
        gamma[gamma==0] = (1-d2)/(N.states-1)
        
        ## ChIP parameters
        par.background = c(1.11594274,0,3.599259)
        par.enrichment = c(2.2810865,0,4.076325)
        
        ## Offsets
        #I am using the same offset for replicates from the same group
        offsets = NULL
        offsets[paste0('Group',1:Groups)] = list(NULL)
        offsets[paste0('Group',1:Groups)] = rep(0,Groups) 
    }
    # Adjusting for SNR pct
    snr = exp(par.enrichment[1])/exp(par.background[1])
    var.enrichment = exp(par.enrichment[1])*(1+exp(par.enrichment[1])/par.enrichment[3])
    meanvar.enrichment = exp(par.enrichment[1])/(var.enrichment)
    par.enrichment.new = c(log(pct*snr*exp(par.background[1])),0,(pct*snr*exp(par.background[1]))*meanvar.enrichment/(1-meanvar.enrichment))
    par.enrichment = par.enrichment.new
    
    return(list('pi'=pi,'gamma'=gamma,'par.background'=par.background,'par.enrichment'=par.enrichment,'offsets'=offsets))
}

for(i in seq_len(nrow(scenario))){
    output = list()
    output[['z']] = list()
    output[['viterbi']] = list()
    output[['psi']] = list()
    output[['preddiff']] = list()
    output[['fdrpp']] = list()
    
    # Loading simulated data
    datafiles = paste0('./Data_',scenario$Pct[i],'SNR/',scenario$Label[i],'/',scenario$Label[i],'_',1:scenario$Nsim[i],'.RData')
    for(j in 1:length(datafiles)){
        load(datafiles[j])
        output[['z']][[j]] = dat$z
    }
    
    # Loading output data
    outfiles = paste0('./Output_',scenario$Pct[i],'SNR/',scenario$Label[i],'/',scenario$Label[i],'_',1:scenario$Nsim[i],'/',scenario$Label[i],'_',1:scenario$Nsim[i],'_mixHMMConstr.RData')
    for(j in 1:length(outfiles)){
        if(file.exists(outfiles[j])){
            cat('Output: ',outfiles[j],'\n')
            load(outfiles[j])
            output[['viterbi']][[j]] = peaks$Viterbi
            output[['psi']][[j]] = peaks$Psi
            output[['preddiff']][[j]] = apply(peaks$Mix.Prob,1,which.max) + 1
            output[['fdrpp']][[j]] = peaks$Prob
        }
    }
    
    # Calculating measures now
    acc.list[[scenario$Label[i]]] = list()
    prec.list[[scenario$Label[i]]] = list()
    rec.list[[scenario$Label[i]]] = list()
    f1.list[[scenario$Label[i]]] = list()
    kappa.list[[scenario$Label[i]]] = list()
    par.list[[scenario$Label[i]]] = list()
    preddiff.list[[scenario$Label[i]]] = list()
    bias.list[[scenario$Label[i]]] = list()
    fdrpp.list[[scenario$Label[i]]] = list()
    
    # Model-specifc parameters
    par = genpar(scenario$Marker[i],scenario$Groups[i],pct=scenario$Pct[i]/100)
    
    # Looping through simulated datasets
    for(j in 1:scenario$Nsim[i]){
        if(file.exists(outfiles[j])){
            z = output[['z']][[j]]
            viterbi = output[['viterbi']][[j]]
            psi = output[['psi']][[j]]
            preddiff = output[['preddiff']][[j]]
            fdrpp = output[['fdrpp']][[j]]
            
            viterbi.preddiff = ifelse(viterbi==0,1,ifelse(viterbi==2,2^scenario$Groups[i],preddiff))
            z.agg = ifelse(z==1,0,ifelse(z==2^scenario$Groups[i],2,1))
            
            cm = table(viterbi,z.agg)
            n = sum(cm) # number of windows
            nc = nrow(cm) # number of states
            diag = diag(cm) # number of correctly classified windows per state
            rowsums = apply(cm, 1, sum) # sum of predicted condition
            colsums = apply(cm, 2, sum) # sum of true condition
            p = rowsums / n # predicted prevalence
            q = colsums / n # true prevalence
            
            accuracy = sum(diag) / n
            precision = diag / rowsums
            recall = diag /  colsums
            f1 = 2 * precision * recall / (precision + recall)
            kappa = (accuracy - sum(p*q)) / (1 - sum(p*q))
            
            acc.list[[scenario$Label[i]]][[paste0('Sim',j)]] = accuracy
            prec.list[[scenario$Label[i]]][[paste0('Sim',j)]] = precision
            rec.list[[scenario$Label[i]]][[paste0('Sim',j)]] = recall
            f1.list[[scenario$Label[i]]][[paste0('Sim',j)]] = f1
            kappa.list[[scenario$Label[i]]][[paste0('Sim',j)]] = kappa
            par.list[[scenario$Label[i]]][[paste0('Sim',j)]] = psi[c('HMM.Mean.Int','HMM.Mean.Diff','HMM.Disp.Int','HMM.Disp.Diff')]
            preddiff.list[[scenario$Label[i]]][[paste0('Sim',j)]] = as.matrix(table(viterbi.preddiff,z)) #Accuracy of differential pattern prediction
            bias.list[[scenario$Label[i]]][[paste0('Sim',j)]] = c(psi['HMM.Mean.Int']-par$par.background[1],
                                                                  psi['HMM.Mean.Diff']-(par$par.enrichment[1]-par$par.background[1]),
                                                                  psi['HMM.Disp.Int']-log(par$par.background[3]),
                                                                  psi['HMM.Disp.Diff']-log(par$par.enrichment[3]/par$par.background[3]))
            
            #### FDR thresholding for Precision vs. Recall
            fdr = seq(0.005,0.250,0.005)
            DT = data.table(Sim = z.agg,Postprob = fdrpp$PostProb2)
            fdrpp.list[[scenario$Label[i]]][[paste0('Sim',j)]] = rbindlist(lapply(fdr,FUN = function(x){
                
                threshold = tryCatch(DT[,mixHMM::cutoff(postprob = Postprob,alpha = x)], error=function(err) NA)
                
                if(is.na(threshold)){
                    return(data.table(FDR = x,FPR = NA,Precision = NA,Recall = NA))
                } else{
                    DT[,Pred := 1*(Postprob>threshold)]
                    return(data.table(FDR = x,
                                      FPR = sum(DT[,Sim==F & Pred==T])/sum(DT[,Sim==F]),
                                      Precision = sum(DT[,Sim==T & Pred==T])/sum(DT[,Pred==T]),
                                      Recall = sum(DT[,Sim==T & Pred==T])/sum(DT[,Sim==T])))   
                }
            }))
        }
    }
}

cmd = paste0('mkdir Metrics_',round(100*pct),'SNR')
system(cmd)

save(acc.list,file=paste0('Metrics_',round(100*pct),'SNR','/Metrics_',round(100*pct),'SNR_acc.RData'))
save(prec.list,file=paste0('Metrics_',round(100*pct),'SNR','/Metrics_',round(100*pct),'SNR_prec.RData'))
save(rec.list,file=paste0('Metrics_',round(100*pct),'SNR','/Metrics_',round(100*pct),'SNR_rec.RData'))
save(f1.list,file=paste0('Metrics_',round(100*pct),'SNR','/Metrics_',round(100*pct),'SNR_f1.RData'))
save(kappa.list,file=paste0('Metrics_',round(100*pct),'SNR','/Metrics_',round(100*pct),'SNR_kappa.RData'))
save(par.list,file=paste0('Metrics_',round(100*pct),'SNR','/Metrics_',round(100*pct),'SNR_par.RData'))
save(preddiff.list,file=paste0('Metrics_',round(100*pct),'SNR','/Metrics_',round(100*pct),'SNR_preddiff.RData'))
save(bias.list,file=paste0('Metrics_',round(100*pct),'SNR','/Metrics_',round(100*pct),'SNR_bias.RData'))
save(fdrpp.list,file=paste0('Metrics_',round(100*pct),'SNR','/Metrics_',round(100*pct),'SNR_fdrpp.RData'))
