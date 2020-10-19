rm(list=ls())
workdir <- getwd()

system('mkdir ./Analysis/')

dirdt <- paste0(workdir,'/Data')

dirs <- list.dirs(path = dirdt,full.names = FALSE,recursive = FALSE)

for(i in dirs){
  
  diranalysis <- paste0('./Analysis/',i)
  
  system(paste('mkdir',diranalysis))
  
  for(sim in 1:100){
    sink(file = file.path(diranalysis,paste0('Simulation',sim,'.R')))
    cat("

library(epigraHMM)
library(data.table)

my_line <- function(x,y,...){
    smoothScatter(x,y, nrpoints = 0, add = TRUE)
    abline(a = 0,b = 1,col = 'red',...)
    abline(h = mean(y),col = 'purple')
}

load('",file.path(dirdt,i,paste0(i,'_',sim,'.RData')),"')
chip <- do.call(cbind,lapply(dat$ChIP,function(x){matrix(unlist(x),byrow = FALSE,ncol = length(x))}))
z <- dat$z
gc <- dat$gc
rm(dat)

M <- nrow(chip)
controlEM <- epigraHMM::controlEM(quiet = F,maxit.em = 100,epsilon.em = c(1e-4,1e-4,1e-4,1e-4))

# Initializing

object <- epigraHMMDataSetFromMatrix(countData = chip,
                                     colData = data.frame(condition = c('A','A','B','B'),replicate = c(1,2,1,2)),
                                     rowRanges = GenomicRanges::GRanges('chr1',IRanges::IRanges(1+500*0:(M-1),500+500*0:(M-1))))

object <- initializer(object = object,control = controlEM)

# Model naive (without gc/normalization)

object.naive <- epigraHMM(object,control = controlEM,type = 'differential',dist = 'nb',random = F)

# Model with gc

object.gc <- cleanGC(object = object,gc = gc)
object.gc <- epigraHMM(object.gc,control = controlEM,type = 'differential',dist = 'nb',random = F)

# Model with normalization

object.norm <- normalizeCounts(object = object,type = 'gam',by = NULL)
object.norm <- epigraHMM(object.norm,control = controlEM,type = 'differential',dist = 'nb',random = F)

# Model with gc + normalization

object.gc.norm <- cleanGC(object = object,gc = gc)
object.gc.norm <- normalizeCounts(object = object.gc.norm,type = 'gam',by = NULL)
object.gc.norm <- epigraHMM(object.gc.norm,control = controlEM,type = 'differential',dist = 'nb',random = F)

# Summarizing the results

dt <- data.table(GoldStandard = (z==2 | z==3),
                 Viterbi.naive = metadata(object.naive)$viterbi,
                 Viterbi.gc = metadata(object.gc)$viterbi,
                 Viterbi.norm = metadata(object.norm)$viterbi,
                 Viterbi.gc.norm = metadata(object.gc.norm)$viterbi,
                 PP.naive = as.numeric(metadata(object.naive)$prob$Differential),
                 PP.gc = as.numeric(metadata(object.gc)$prob$Differential),
                 PP.norm = as.numeric(metadata(object.norm)$prob$Differential),
                 PP.gc.norm = as.numeric(metadata(object.gc.norm)$prob$Differential))

save(dt,file = './Simulation",sim,".RData',compress = 'xz')

# Plotting data

cts <- log1p(assay(object,'counts'))
cts.gc <- log1p(assay(object.gc,'counts')/exp(assay(object.gc,'offset')))
cts.norm <- log1p(assay(object.norm,'counts')/exp(assay(object.norm,'offset')))
cts.gc.norm <- log1p(assay(object.gc.norm,'counts')/exp(assay(object.gc.norm,'offset')))
idx <- 10000:11000

pdf(file = './Simulation",sim,".pdf',height = 7,width = 7)
pairs(cbind(log1p(assay(object.naive,'counts')),gc), lower.panel = NULL, upper.panel = my_line,main = 'Naive model - Counts ~ GC')
pairs(cbind(log1p(assay(object.naive,'counts')),ref=rowMeans(log1p(assay(object.naive,'counts')))), lower.panel = NULL, upper.panel = my_line,main = 'Naive model - Counts ~ Ref.')

pairs(cbind(cts.gc,gc), lower.panel = NULL, upper.panel = my_line,main = 'GC-adj. model - Counts ~ GC')
pairs(cbind(cts.gc,ref=rowMeans(cts.gc)), lower.panel = NULL, upper.panel = my_line,main = 'GC-adj. model - Counts ~ Ref.')

pairs(cbind(cts.norm,gc), lower.panel = NULL, upper.panel = my_line,main = 'Ref.-adj. model - Counts ~ GC')
pairs(cbind(cts.norm,ref=rowMeans(cts.norm)), lower.panel = NULL, upper.panel = my_line,main = 'Ref.-adj. model - Counts ~ Ref.')

pairs(cbind(cts.gc.norm,gc), lower.panel = NULL, upper.panel = my_line,main = '(GC + Ref.)-adj. model - Counts ~ GC')
pairs(cbind(cts.gc.norm,ref=rowMeans(cts.gc.norm)), lower.panel = NULL, upper.panel = my_line,main = '(GC + Ref.)-adj. model - Counts ~ Ref.')

par(mfrow = c(2,2),oma = c(0, 0, 2, 0))
plot(exp(cts)[idx,1],ylim=c(0,max(exp(cts)[idx,])),type='l',main = 'Sample A1',ylab = 'Counts')
plot(exp(cts)[idx,2],ylim=c(0,max(exp(cts)[idx,])),type='l',main = 'Sample A2',ylab = 'Counts')
plot(exp(cts)[idx,3],ylim=c(0,max(exp(cts)[idx,])),type='l',main = 'Sample B1',ylab = 'Counts')
plot(exp(cts)[idx,4],ylim=c(0,max(exp(cts)[idx,])),type='l',main = 'Sample B1',ylab = 'Counts')
mtext('Naive model', outer = TRUE, cex = 1.5)

par(mfrow = c(2,2),oma = c(0, 0, 2, 0))
plot(exp(cts.gc)[idx,1],ylim=c(0,max(exp(cts.gc)[idx,])),type='l',main = 'Sample A1',ylab = 'Counts')
plot(exp(cts.gc)[idx,2],ylim=c(0,max(exp(cts.gc)[idx,])),type='l',main = 'Sample A2',ylab = 'Counts')
plot(exp(cts.gc)[idx,3],ylim=c(0,max(exp(cts.gc)[idx,])),type='l',main = 'Sample B1',ylab = 'Counts')
plot(exp(cts.gc)[idx,4],ylim=c(0,max(exp(cts.gc)[idx,])),type='l',main = 'Sample B1',ylab = 'Counts')
mtext('GC-adj. model', outer = TRUE, cex = 1.5)

par(mfrow = c(2,2),oma = c(0, 0, 2, 0))
plot(exp(cts.norm)[idx,1],ylim=c(0,max(exp(cts.norm)[idx,])),type='l',main = 'Sample A1',ylab = 'Counts')
plot(exp(cts.norm)[idx,2],ylim=c(0,max(exp(cts.norm)[idx,])),type='l',main = 'Sample A2',ylab = 'Counts')
plot(exp(cts.norm)[idx,3],ylim=c(0,max(exp(cts.norm)[idx,])),type='l',main = 'Sample B1',ylab = 'Counts')
plot(exp(cts.norm)[idx,4],ylim=c(0,max(exp(cts.norm)[idx,])),type='l',main = 'Sample B1',ylab = 'Counts')
mtext('Ref.-adj. model', outer = TRUE, cex = 1.5)

par(mfrow = c(2,2),oma = c(0, 0, 2, 0))
plot(exp(cts.gc.norm)[idx,1],ylim=c(0,max(exp(cts.gc.norm)[idx,])),type='l',main = 'Sample A1',ylab = 'Counts')
plot(exp(cts.gc.norm)[idx,2],ylim=c(0,max(exp(cts.gc.norm)[idx,])),type='l',main = 'Sample A2',ylab = 'Counts')
plot(exp(cts.gc.norm)[idx,3],ylim=c(0,max(exp(cts.gc.norm)[idx,])),type='l',main = 'Sample B1',ylab = 'Counts')
plot(exp(cts.gc.norm)[idx,4],ylim=c(0,max(exp(cts.gc.norm)[idx,])),type='l',main = 'Sample B1',ylab = 'Counts')
mtext('(GC+Ref.)-adj. model', outer = TRUE, cex = 1.5)

dev.off()


",sep = "")
    sink()
    # Changing directories and submitting job
    setwd(diranalysis)
    system(paste0('sbatch -t 2:00:00 --mem=8g R CMD BATCH ',paste0('./Simulation',sim,'.R')))
    setwd(workdir)
  }
}

