rm(list=setdiff(ls(),c('sim','label','pct')))
dirdata <- file.path(normalizePath('../../../../ReadCounts'),paste0('Data_',pct,'SNR'),label)
diroutput <- file.path(normalizePath('../../../'),paste0('Output_',pct,'SNR'),label,paste0(label,sim))

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
countData <- as.matrix(data.frame(dat$ChIPBIC))
colnames(countData) <- unlist(lapply(seq_len(length(dat$ChIPBIC)),function(x){
    paste0(rep(names(dat$ChIPBIC)[x],length(dat$ChIPBIC[[1]])),".",names(dat$ChIPBIC[[1]]))
}))
colData <- data.frame(condition = unlist(lapply(strsplit(x=colnames(countData),split='\\.'),function(x){x[1]})),
                      replicate = unlist(lapply(strsplit(x=colnames(countData),split='\\.'),function(x){x[2]})))
rowRanges <- GenomicRanges::GRanges('chrA',IRanges::IRanges(start = seq(1,100*nrow(countData),by = 100),width = 100))

# Gold standard
Z <- dat$z
ZBic <- dat$zBic

# Defining differential patterns now
ngroups = length(unique(colData$condition))
difflist = NULL
difflist[paste0('Group',1:ngroups)] = list(NULL)
for(i in names(difflist)){difflist[[i]] = c(NA,i)}
difflist = as.data.frame(expand.grid(difflist))
difflist = difflist[order(rowSums(is.na(difflist)),decreasing=T),]

# Pattern list
if(ngroups == 2){
    patternlist <- list(list(2), # Optimal
                        list(1,2)) 
} else{
    if(ngroups == 3){
        patternlist <- list(list(1), 
                            list(1,c(2,3)), # Optimal
                            list(1,c(2,3),2),
                            list(1,c(2,3),2,3),
                            list(1,c(2,3),2,3,c(1,2)),
                            list(1,c(2,3),2,3,c(1,2),c(1,3)))
    } else{
        patternlist <- list(list(1),
                            list(1,c(2,3,4)), # Optimal
                            list(1,c(2,3,4),2),
                            list(1,c(2,3,4),2,3),
                            list(1,c(2,3,4),2,3,4),
                            list(1,c(2,3,4),2,3,4,c(1,2)),
                            list(1,c(2,3,4),2,3,4,c(1,2),c(1,3)),
                            list(1,c(2,3,4),2,3,4,c(1,2),c(1,3),c(1,4)),
                            list(1,c(2,3,4),2,3,4,c(1,2),c(1,3),c(1,4),c(2,3)),
                            list(1,c(2,3,4),2,3,4,c(1,2),c(1,3),c(1,4),c(2,3),c(2,4)),
                            list(1,c(2,3,4),2,3,4,c(1,2),c(1,3),c(1,4),c(2,3),c(2,4),c(3,4)),
                            list(1,c(2,3,4),2,3,4,c(1,2),c(1,3),c(1,4),c(2,3),c(2,4),c(3,4),c(1,2,3)),
                            list(1,c(2,3,4),2,3,4,c(1,2),c(1,3),c(1,4),c(2,3),c(2,4),c(3,4),c(1,2,3),c(1,2,4)),
                            list(1,c(2,3,4),2,3,4,c(1,2),c(1,3),c(1,4),c(2,3),c(2,4),c(3,4),c(1,2,3),c(1,2,4),c(1,3,4)))
    }
}

# mixHMMConstr
tryCatch({
    
    initalObject <- epigraHMM::epigraHMMDataSetFromMatrix(countData = countData,colData = colData,rowRanges = rowRanges)
    initalObject <- epigraHMM::initializer(object = initalObject,control = epigraHMM::controlEM(epsilon.em = c(1e-3,1e-3,1e-6,1e-3),maxit.innerem = 3,quiet = FALSE))
    
    bic_list <- list()
    
    for(i in seq_len(length(patternlist))){
        # Fitting model with reduced patterns
        object <- epigraHMM::epigraHMM(object = initalObject,
                                       control = epigraHMM::controlEM(epsilon.em = c(1e-3,1e-3,1e-6,1e-3),maxit.innerem = 3,quiet = FALSE,pattern = patternlist[[i]]),
                                       type = 'differential', dist = 'nb', random = FALSE)
        
        # Creating output
        namelist <- paste(lapply(patternlist[[i]],function(x){paste(x,collapse = '/')}),collapse = '_')
        bic_list[[namelist]] <- list()
        bic_list[[namelist]][['prob']] <- metadata(object)$prob$Differential
        bic_list[[namelist]][['viterbi']] <- metadata(object)$viterbi
        bic_list[[namelist]][['BIC']] <- tail(metadata(object)$parHist$BIC,1)
    }
    
    save(object,bic_list,Z,ZBic,file=file.path(diroutput,paste0(label,sim,'_epigraHMM.RData')))
    rm(list=ls())
},warning=function(e){cat('WARNING: ',conditionMessage(e), "\n");saveerror(paste0(label,sim,'_epigraHMM.txt'),'Error')},error=function(e){cat('ERROR: ',conditionMessage(e), "\n");saveerror(paste0(label,sim,'_epigraHMM.txt'),'Error')})
