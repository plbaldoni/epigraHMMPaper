#### mixHMM Simulation with reduced SNR####
library(Rcpp)
library(ZIMHMM)
library(plyr)
library(mixHMM)

### Simulation scnearios
# Marker:   H3K27me3    H3K36me3
# Groups:   2, 3, 4
# Replicates: 1, 2, 4
# Nwindow: 1e5, 5e5, 1e6
# Nsim: 100

scenario = expand.grid(Marker=c('H3K27me3','H3K36me3'),Groups=c(2,3,4),Replicates=c(1,2,4),Nwindow=c(1e5,5e5,1e6),Nsim=100,Pct=round(100*c(0.7,0.8,0.9,1)))
scenario$Label = with(scenario,paste(Marker,paste0(Groups,'G'),paste0(Replicates,'R'),paste0(Nwindow,'W'),paste0(Pct,'SNR'),sep="_"))
print(scenario)

## Creating folders fot datasets
for(i in unique(scenario$Pct)){
    if(!dir.exists(paste0('Data_',i,'SNR/'))){
        system(paste('mkdir',paste0('Data_',i,'SNR/')))
    }
}

## Control Parameters: mean (mux), dispersion (phix)
# Current simulation does not account for control, so it is useless here
mux = 9
phix = 2.5

## ChIP Parameters:
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

for(i in 1:nrow(scenario)){
    scenario.i = scenario[i,]
    
    if(!dir.exists(paste0(paste0('Data_',scenario.i$Pct,'SNR/'),scenario.i$Label))){
        cmd = paste('mkdir',paste0(paste0('Data_',scenario.i$Pct,'SNR/'),scenario.i$Label))
        cat('Command: ',cmd)
        system(cmd)   
    }
    
    ### Simulation begins
    # Sample size and number of windows
    ngroups = scenario.i$Groups
    nreps = scenario.i$Replicates
    nwindow = scenario.i$Nwindow
    nstates = 2^ngroups
    
    # Model-specifc parameters
    par = genpar(scenario.i$Marker,ngroups,pct=scenario.i$Pct/100)
    pi = par$pi
    gamma = par$gamma
    par.background = par$par.background
    par.enrichment = par$par.enrichment
    offsets = par$offsets
    
    #Data frame with differential paths
    difflist = NULL
    difflist[paste0('Group',1:ngroups)] = list(NULL)
    for(i in names(difflist)){difflist[[i]] = c(NA,i)}
    difflist = as.data.frame(expand.grid(difflist))
    difflist = difflist[order(rowSums(is.na(difflist)),decreasing=T),]
    difflist
    
    for(sim in 1:scenario.i$Nsim)
    {
        if(!file.exists(paste0('./',paste0('Data_',scenario.i$Pct,'SNR/'),scenario.i$Label,'/',scenario.i$Label,'_',sim,'.RData'))){
            cat("\r",paste('Simulation: '),sim)
            
            # Simulate hidden states
            z = mixHMM::simulateHMM(gamma,nwindow)
            
            # Now simulating control. Again, control is irrelevant here but I will keep it
            # to avoid rewriting the entire simulation.
            controllist = NULL
            controllist[paste0('Group',1:ngroups)] = list(NULL)
            for(i in names(controllist)){
                for(j in 1:nreps){
                    controllist[[i]][[paste0('Replicate',j)]] = rnbinom(nwindow,mu=mux,size=phix)
                }
            }
            str(controllist)
            
            # Now simulating ChIP
            chiplist = NULL
            chiplist[paste0('Group',1:ngroups)] = list(NULL)
            for(i in names(chiplist)){
                for(k in 1:nreps){
                    chiplist[[i]][[paste0('Replicate',k)]] = rep(0,nwindow)
                    for(j in 1:nstates){
                        if(is.na(difflist[[i]][[j]])){
                            chiplist[[i]][[paste0('Replicate',k)]][which(z==j)] = rnbinom(sum(z==j),mu=exp(par.background[1]+par.background[2]*log(controllist[[i]][[paste0('Replicate',k)]][which(z==j)]+1)+offsets[[i]]),size=par.background[3])
                        } else{
                            chiplist[[i]][[paste0('Replicate',k)]][which(z==j)] = rnbinom(sum(z==j),mu=exp(par.enrichment[1]+par.enrichment[2]*log(controllist[[i]][[paste0('Replicate',k)]][which(z==j)]+1)+offsets[[i]]),size=par.enrichment[3])
                        }
                    }
                }
            }
            str(chiplist)
            
            dat = list('sim'=sim,'z'=z,'ChIP'=chiplist)
            save(dat,file=paste0('./',paste0('Data_',scenario.i$Pct,'SNR/'),scenario.i$Label,'/',scenario.i$Label,'_',sim,'.RData'))
        }
    }
}
