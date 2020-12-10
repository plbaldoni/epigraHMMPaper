
iters <- 1:100
is.tf <- F
wd <- getwd()

for(tf in is.tf){
    for(it in iters){
        # Creating code directory
        cmd = paste0('mkdir autosim',it)
        cat('Command: ',cmd)
        system(cmd)
        
        # Copying autosim.R to tmp.R
        cmd = paste0('cp autosim.R tmp.R')
        cat('Command: ',cmd)
        system(cmd)
        
        # Writing the specification to the top of the file autosim.R
        cmd = paste0("echo '",paste0('is.tf = ',tf),paste0(';iters = ',it),"' | cat - tmp.R > ",paste0('./autosim',it,".R"))
        cat('Command: ',cmd)
        system(cmd)
        
        # Moving new autosim
        cmd = paste0('mv autosim',it,'.R ./autosim',it,'/')
        cat('Command: ',cmd)
        system(cmd)
        
        # Remove tmp
        cmd = paste0('rm tmp.R')
        cat('Command: ',cmd)
        system(cmd)
        
        # Submitting jobs, change it if not using SLURM
        setwd(paste0('./autosim',it,'/'))
        cmd = paste0('sbatch -t 12:00:00 --mem=16g R CMD BATCH ./autosim',it,".R")
        cat('Command: ',cmd)
        system(cmd)
        setwd(wd)
    }
}


