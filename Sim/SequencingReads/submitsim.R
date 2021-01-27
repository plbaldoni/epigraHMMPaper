
iters <- 1:100
is.tf <- F
wd <- normalizePath('.')

for (tf in is.tf) {
    for (it in iters) {
        # Creating code directory
        dir.create(file.path(wd,paste0('autosim',it)))
        
        # Copying autosim.R to tmp.R
        cmd = paste('cp',file.path(wd,'autosim.R'),file.path(wd,'tmp.R'))
        cat('Command: ',cmd)
        system(cmd)
        
        # Writing the specification to the top of the file autosim.R
        cmd = paste0("echo '",paste0('is.tf = ',tf),paste0(';iters = ',it),"' | cat - ",file.path(wd,'tmp.R')," > ",
                     file.path(wd,paste0('autosim',it,'.R')))
        cat('Command: ',cmd)
        system(cmd)
        
        # Moving new autosim
        cmd = paste('mv',file.path(wd,paste0('autosim',it,'.R')),file.path(wd,paste0('autosim',it)))
        cat('Command: ',cmd)
        system(cmd)
        
        # Remove tmp
        cmd = paste('rm',file.path(wd,'tmp.R'))
        cat('Command: ',cmd)
        system(cmd)
        
        # Submitting jobs, change it if not using SLURM
        setwd(file.path(wd,paste0('autosim',it)))
        
        sink(file = paste0('autosim',it,'.sh'))
        cat(
'#!/bin/bash
 
#SBATCH -t 12:00:00
#SBATCH --mem=12g
 
source ~/anaconda3/etc/profile.d/conda.sh
conda activate py36
R CMD BATCH ',paste0('autosim',it,'.R'))
        sink()
        
        cmd = paste('sbatch',paste0('autosim',it,'.sh'))
        cat('Command: ',cmd)
        system(cmd)
        setwd(wd)
    }
}


