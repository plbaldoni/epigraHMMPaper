rm(list=ls())
workdir <- getwd()

system('mkdir ./Analysis/')

dirdt <- '/pine/scr/b/a/baldoni/Rashid/Project2/Simulation_Control/Data'

dirs <- list.dirs(path = dirdt,full.names = FALSE,recursive = FALSE)

for(i in dirs){
    
    diranalysis <- paste0('./Analysis/',i)
    
    system(paste('mkdir',diranalysis))
    
    for(sim in 1:100){
        sink(file = file.path(diranalysis,paste0('Simulation',sim,'.R')))
        cat("

library(epigraHMM)

load('",file.path(dirdt,i,paste0(i,'_',sim,'.RData')),"')
dat1 <- dat
load('",file.path(dirdt,i,paste0(i,'_',sim+100,'.RData')),"')
dat2 <- dat
rm(dat)

M <- nrow(dat1)
controlEM <- epigraHMM::controlEM(quiet = F)

chip <- cbind(as.matrix(dat1[,c('y.1','y.2')]),as.matrix(dat2[,c('y.1','y.2')]))
control <- cbind(as.matrix(dat1[,c('x.1','x.2')]),as.matrix(dat2[,c('x.1','x.2')]))
z <- dat1$z+dat2$z

# Model without control

object.nocontrol <- epigraHMM::epigraHMMDataSetFromMatrix(countData = chip,
                                       colData = data.frame(condition = c('A','A','B','B'),replicate = c(1,2,1,2)),
                                       rowRanges = GenomicRanges::GRanges('chr1',IRanges::IRanges(1+500*0:(M-1),500+500*0:(M-1))))

object.nocontrol <- createOffset(object = object.nocontrol,type = 'loess',span = 1)

object.nocontrol <- epigraHMM(object.nocontrol,control = controlEM,type = 'differential',dist = 'nb',random = F)


# Model with control

object.withcontrol <- epigraHMM::epigraHMMDataSetFromMatrix(countData = chip,
                                       colData = data.frame(condition = c('A','A','B','B'),replicate = c(1,2,1,2)),
                                       rowRanges = GenomicRanges::GRanges('chr1',IRanges::IRanges(1+500*0:(M-1),500+500*0:(M-1))))

metadata(object.withcontrol)$adjustment <- list(control)
metadata(object.withcontrol)$grouping <- rep('A',nrow(control))
metadata(object.withcontrol)$normalize <- T

object.withcontrol <- createOffset(object = object.withcontrol,type = 'loess')

object.withcontrol <- epigraHMM(object.withcontrol,control = controlEM,type = 'differential',dist = 'nb',random = F)

# Model with grouped control

object.withgroupedcontrol <- epigraHMM::epigraHMMDataSetFromMatrix(countData = chip,
                                       colData = data.frame(condition = c('A','A','B','B'),replicate = c(1,2,1,2)),
                                       rowRanges = GenomicRanges::GRanges('chr1',IRanges::IRanges(1+500*0:(M-1),500+500*0:(M-1))))

metadata(object.withgroupedcontrol)$adjustment <- list(control)
metadata(object.withgroupedcontrol)$grouping <- as.character(z)
metadata(object.withgroupedcontrol)$normalize <- T

object.withgroupedcontrol <- createOffset(object = object.withgroupedcontrol,type = 'loess')

object.withgroupedcontrol <- epigraHMM(object.withgroupedcontrol,control = controlEM,type = 'differential',dist = 'nb',random = F)

# Summarizing the results

dt <- data.table(WithoutControl = metadata(object.nocontrol)$viterbi,
                 WithControl = metadata(object.withcontrol)$viterbi,
                 WithGroupedControl = metadata(object.withgroupedcontrol)$viterbi,
                 GoldStandard = z)
                 
save(dt,file = './Simulation",sim,".RData',compress = 'xz')

",sep = "")
        sink()
        # Changing directories and submitting job
        setwd(diranalysis)
        system(paste0('sbatch -t 1:00:00 --mem=4g R CMD BATCH ',paste0('./Simulation',sim,'.R')))
        setwd(workdir)
    }
}

