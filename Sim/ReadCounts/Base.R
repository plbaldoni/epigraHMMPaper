rm(list=setdiff(ls(),c('sim','label','pct')))
dirdata=paste0('/pine/scr/b/a/baldoni/Rashid/Project2/Simulation/Data_',pct,'SNR/',label,'/')
diroutput=paste0('/pine/scr/b/a/baldoni/Rashid/Project2/Simulation/Output_',pct,'SNR/',label,'/',label,sim,'/')

library(Rcpp)
library(microbenchmark)
library(Biobase)
library(mixHMM)

# Error function for output
saveerror = function(name,text){
    fileConn<-file(paste0(diroutput,name))
    writeLines(c(text), fileConn)
    close(fileConn)
}

# Loading data
load(paste0(dirdata,label,sim,'.RData'))

# Creating input elements
ChIP = as.matrix(data.frame(dat$ChIP))
offset = matrix(0,ncol=ncol(ChIP),nrow=nrow(ChIP))
group = rep(1:length(dat$ChIP),each=length(dat$ChIP$Group1))
B = 2^length(dat$ChIP)-2 
control = mixHMM::controlPeaks(epsilon.em = c(1e-3,1e-3,1e-6,1e-3),maxit.innerem = 3)

Z = dat$z

# mixHMMConstr
tryCatch({
    time = microbenchmark(assign('tmp',mixHMM::mixHMMConstr(ChIP = ChIP,Control = NULL,offset = offset,group = group,B = B,control = control)),times=1)
    peaks = tmp
    save(time,peaks,Z,file=paste0(diroutput,paste0(label,sim,'_mixHMMConstr.RData')))
    rm(time);rm(peaks);rm(tmp)
},warning=function(e){cat('WARNING: ',conditionMessage(e), "\n");saveerror(paste0(label,sim,'_mixHMMConstr.txt'),'Error')},error=function(e){cat('ERROR: ',conditionMessage(e), "\n");saveerror(paste0(label,sim,'_mixHMMConstr.txt'),'Error')})
