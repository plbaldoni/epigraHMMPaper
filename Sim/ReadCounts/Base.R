rm(list=setdiff(ls(),c('sim','label','pct')))
dirdata <- file.path(normalizePath('../..'),paste0('Data_',pct,'SNR'),label)
diroutput <- file.path(normalizePath('../..'),paste0('Output_',pct,'SNR'),label,paste0(label,sim))

if(!dir.exists(dirdata)){dir.create(dirdata,recursive = TRUE)}
if(!dir.exists(diroutput)){dir.create(diroutput,recursive = TRUE)}

library(microbenchmark)
library(epigraHMM)

# Error function for output
saveerror = function(name,text){
    fileConn<-file(paste0(diroutput,name))
    writeLines(c(text), fileConn)
    close(fileConn)
}

# Loading data
load(file.path(dirdata,paste0(label,sim,'.RData')))

# Creating input elements
countData <- as.matrix(data.frame(dat$ChIP))
colnames(countData) <- paste0(names(lapply(dat$ChIP,function(x){names(x)})),'.',unlist(lapply(dat$ChIP,function(x){names(x)})))
colData <- data.frame(condition = unlist(lapply(strsplit(x=colnames(countData),split='\\.'),function(x){x[1]})),
                      replicate = unlist(lapply(strsplit(x=colnames(countData),split='\\.'),function(x){x[2]})))
rowRanges <- GenomicRanges::GRanges('chrA',IRanges::IRanges(start = seq(1,100*nrow(countData),by = 100),width = 100))
control <- epigraHMM::controlEM(epsilon.em = c(1e-3,1e-3,1e-6,1e-3),maxit.innerem = 3,quiet = FALSE)

# Gold standard
Z = dat$z

# mixHMMConstr
tryCatch({
    time <- microbenchmark({
        assign('object',epigraHMM::epigraHMMDataSetFromMatrix(countData = countData,colData = colData,rowRanges = rowRanges))
        assign('object',epigraHMM::initializer(object = get('object'), control = control))
        assign('object',epigraHMM::epigraHMM(object = get('object'), control = control, type = 'differential', dist = 'nb', random = FALSE))
    },times=1)
    save(time,object,Z,file=file.path(diroutput,paste0(label,sim,'_epigraHMM.RData')))
    rm(list=ls())
},warning=function(e){cat('WARNING: ',conditionMessage(e), "\n");saveerror(paste0(label,sim,'_epigraHMM.txt'),'Error')},error=function(e){cat('ERROR: ',conditionMessage(e), "\n");saveerror(paste0(label,sim,'_epigraHMM.txt'),'Error')})
