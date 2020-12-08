rm(list=ls())
wd = getwd()
### Simulation scnearios
# Marker:   H3K27me3    H3K36me3
# Groups:   2, 3, 4
# Replicates: 1, 2, 4
# Nwindow: 1e5, 5e5, 1e6
# Nsim: 100

scenario = expand.grid(Marker=c('H3K36me3','H3K27me3'),Groups=c(2,3,4),Replicates=c(1,2,4),Nwindow=c(1e5),Nsim=100,Pct=round(100*c(0.7,0.8,0.9,1)))
scenario$Label = with(scenario,paste(Marker,paste0(Groups,'G'),paste0(Replicates,'R'),paste0(Nwindow,'W'),paste0(Pct,'SNR'),sep="_"))
scenario <- scenario[order(scenario$Groups,scenario$Replicates,scenario$Nwindow,-scenario$Pct),]
print(scenario)

for(i in unique(scenario$Pct)){
    if(!dir.exists(paste0('Analysis_',i,'SNR/'))){
        system(paste('mkdir',paste0('Analysis_',i,'SNR/')))
        system(paste('mkdir',paste0('Output_',i,'SNR/')))   
    }
}

for(i in 1:nrow(scenario)){
    scenario.i = scenario[i,]
    
    if(!dir.exists(paste0(paste0('Analysis_',scenario.i$Pct,'SNR/'),scenario.i$Label))){
        # Creating code directory
        dir.create(paste0(paste0('Analysis_',scenario.i$Pct,'SNR/'),scenario.i$Label),recursive = TRUE)
    }
    
    if(!dir.exists(paste0(paste0('Output_',scenario.i$Pct,'SNR/'),scenario.i$Label))){
        # Creating output directory
        dir.create(paste0(paste0('Output_',scenario.i$Pct,'SNR/'),scenario.i$Label),recursive = TRUE)
    }
    
    for(sim in 1:scenario.i$Nsim){
        if(!file.exists(paste0(wd,'/',paste0('Analysis_',scenario.i$Pct,'SNR/'),scenario.i$Label,'/',scenario.i$Label,'_',sim,'/',paste0(scenario.i$Label,"_",sim,".R")))){
            # Copying Base.R to tmp
            cmd = paste0('cp Base.R tmp.R')
            cat('Command: ',cmd)
            system(cmd)
            
            # Creating code directory
            cmd = paste('mkdir',paste0(paste0('Analysis_',scenario.i$Pct,'SNR/'),scenario.i$Label,'/',scenario.i$Label,'_',sim))
            cat('Command: ',cmd)
            system(cmd)
            
            # Writing the specification to the top of the file
            cmd = paste0("echo '",paste0('sim = "_',sim,'"'),paste0(';label = "',scenario.i$Label,'"'),paste0(';pct = "',scenario.i$Pct,'"'),"' | cat - tmp.R > ",paste0('./',paste0('Analysis_',scenario.i$Pct,'SNR/'),scenario.i$Label,'/',scenario.i$Label,'_',sim,'/',scenario.i$Label,"_",sim,".R"))
            cat('Command: ',cmd)
            system(cmd)
            
            # Remove tmp
            cmd = paste0('rm tmp.R')
            cat('Command: ',cmd)
            system(cmd)
        }
        if(!file.exists(paste0(paste0('Output_',scenario.i$Pct,'SNR/'),scenario.i$Label,'/',scenario.i$Label,'_',sim,'/',scenario.i$Label,'_',sim,'_epigraHMM.RData')) &
           file.exists(paste0(paste0('../ReadCounts/Data_',scenario.i$Pct,'SNR/',scenario.i$Label,'/'),scenario.i$Label,'_',sim,'.RData'))){

            # Creating output directory
            dir.create(paste0(paste0('Output_',scenario.i$Pct,'SNR/'),scenario.i$Label,'/',scenario.i$Label,'_',sim),recursive = TRUE)
            
            # Submitting jobs
            setwd(paste0(wd,'/',paste0('Analysis_',scenario.i$Pct,'SNR/'),scenario.i$Label,'/',scenario.i$Label,'_',sim,'/'))
            cmd = paste0('sbatch -t 2:00:00 --mem=24g R CMD BATCH ./',scenario.i$Label,"_",sim,".R")
            cat('Command: ',cmd)
            system(cmd)
            setwd(wd)   
        }
    }
}

