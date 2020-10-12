library(Rcpp)
library(plyr)
library(data.table)
library(magrittr)
library(ggplot2)
### My function

my_line <- function(x,y,...){
  smoothScatter(x,y, nrpoints = 0, add = TRUE)
  abline(a = 0,b = 1,col = 'red',...)
  abline(h = mean(y),col = 'purple')
}

eq <- function(xy0,xy1,phi0 = 10,phi1 = 2.5,prob = c(0.025,0.975),alpha = 11,plot = FALSE){
  a0 <- (xy0[1,1]*(xy0[3,2]-xy0[2,2])+xy0[2,1]*(xy0[1,2]-xy0[3,2])+xy0[3,1]*(xy0[2,2]-xy0[1,2]))/((xy0[1,1]-xy0[2,1])*(xy0[1,1]-xy0[3,1])*(xy0[2,1]-xy0[3,1]))
  b0 <- (xy0[2,2]-xy0[1,2])/(xy0[2,1]-xy0[1,1]) - a0*(xy0[1,1]+xy0[2,1])
  c0 <- xy0[1,2] - a0*xy0[1,1]^2-b0*xy0[1,1]
  
  a1 <- (xy1[1,1]*(xy1[3,2]-xy1[2,2])+xy1[2,1]*(xy1[1,2]-xy1[3,2])+xy1[3,1]*(xy1[2,2]-xy1[1,2]))/((xy1[1,1]-xy1[2,1])*(xy1[1,1]-xy1[3,1])*(xy1[2,1]-xy1[3,1]))
  b1 <- (xy1[2,2]-xy1[1,2])/(xy1[2,1]-xy1[1,1]) - a1*(xy1[1,1]+xy1[2,1])
  c1 <- xy1[1,2] - a1*xy1[1,1]^2-b1*xy1[1,1]
  
  if(plot){
    ## Means
    
    xseq <- seq(0,1,by = 0.001)
    y0seq <- log(exp(a0*xseq^2 + b0*xseq + c0)+1)
    y1seq <- log(exp(a1*xseq^2 + b1*xseq + c1)+1)
    
    ## Plotting
    par(mfrow = c(2,2))
    plot(y1seq~xseq,type='l',xlim = c(0,1.05),ylim = c(0,5.1),xlab = 'GC',ylab = 'log(mean+1)')
    points(x = xy1[,1],y = log(exp(xy1[,2])+1))
    segments(1,0,1,5,lty = "dashed",col = 'grey')
    segments(0,5,1,5,lty = "dashed",col = 'grey')
    segments(0,0,1,0,lty = "dashed",col = 'grey')
    segments(0,0,0,5,lty = "dashed",col = 'grey')
    segments(0.3,0,0.3,5,lty = "dashed",col = 'grey')
    segments(0.7,0,0.7,5,lty = "dashed",col = 'grey')
    lines(x = xseq,y = log(qnbinom(p = prob[1],size = phi1,mu = exp(y1seq))+1),col = 'blue', type = 'l')
    lines(x = xseq,y = log(qnbinom(p = prob[2],size = phi1,mu = exp(y1seq))+1),col = 'red', type = 'l')
    
    lines(x = xseq,y = y0seq,col = 'black', type = 'l',lty = 'dashed')
    points(x = xy0[,1],y = log(exp(xy0[,2])+1))
    lines(x = xseq,y = log(qnbinom(p = prob[1],size = phi0,mu = exp(y0seq))+1),col = 'blue', type = 'l',lty = 'dashed')
    lines(x = xseq,y = log(qnbinom(p = prob[2],size = phi0,mu = exp(y0seq))+1),col = 'red', type = 'l',lty = 'dashed')
    
    ## Plotting counts
    n1 <- 400
    n2 <- 400
    n3 <- 400
    n <- n1+n2+n3
    
    gc <- rbeta(n,alpha,alpha)
    y <- c(rnbinom(n1,size = phi0,mu = exp(a0*gc[1:n1]^2 + b0*gc[1:n1] + c0)),
           rnbinom(n2,size = phi1,mu = exp(a1*gc[(n1+1):(n1+n2)]^2 + b1*gc[(n1+1):(n1+n2)] + c1)),
           rnbinom(n3,size = phi0,mu = exp(a0*gc[(n1+n2+1):n]^2 + b0*gc[(n1+n2+1):n] + c0)))
    
    plot(x = 1:n,y = y,type = 'l',xlab = 'Window',ylab = 'Counts')
    smoothScatter(log(y[-((n1+1):(n1+n2))]+1)~gc[-((n1+1):(n1+n2))],xlim = c(0,1.05),ylim = c(0,5.1),xlab = 'GC',ylab = 'log(Counts+1)',main = 'Background')
    smoothScatter(log(y[(n1+1):(n1+n2)]+1)~gc[(n1+1):(n1+n2)],xlim = c(0,1.05),ylim = c(0,5.1),xlab = 'GC',ylab = 'log(Counts+1)',main = 'Enrichment')
  }
  
  return(list(coeff = cbind('Background' = c(a0,b0,c0),
                            'Enrichment' = c(a1,b1,c1)),
              disp = c(phi0,phi1),
              gc_beta = alpha))
}

pareq <- function(depth,direction,...){
  
  depth <- ifelse(depth=='low',0.9,
                  ifelse(depth=='medium',1,ifelse(depth=='high',1.1,NA)))
  
  x <- c(0.3,0.5,0.7)
  if(direction == 'up'){
    y0 <- c(1.25,1.75,2.25)
    y1 <- c(1.50,2.75,3.25)
  }
  if(direction == 'down'){
    y0 <- c(2.25,1.75,1.25)
    y1 <- c(3.25,2.75,1.50)
  }
  if(direction == 'null'){
    y0 <- c(1.75,1.75,1.75)
    y1 <- c(2.5,2.5,2.5)
  }
  
  xy0 <- cbind(x,y0*depth)
  xy1 <- cbind(x,y1*depth)
  
  eq(xy0 = xy0,xy1 = xy1,...)
}

### Simulation scenarios

scenario = as.data.table(expand.grid(G1_Depth=as.numeric(c('1','2','3')),
                                     G2_Depth=as.numeric(c('1','2','3')),
                                     G1_Direction=as.numeric(c('1','2','3')),
                                     G2_Direction=as.numeric(c('1','2','3')),
                                     Nwindow=25000,Nsim=100,Groups = 2,Replicates =2))

scenario <- scenario[(G2_Depth >= G1_Depth),]
scenario <- scenario[(G2_Direction >= G1_Direction),]

scenario$G1_Depth_Lab <- scenario$G1_Depth %>% plyr::mapvalues(from = 1:3, to = c('low','medium','high')) %>% factor(levels = c('low','medium','high')) %>% as.character()
scenario$G2_Depth_Lab <- scenario$G2_Depth %>% plyr::mapvalues(from = 1:3, to = c('low','medium','high')) %>% factor(levels = c('low','medium','high')) %>% as.character()
scenario$G1_Direction_Lab <- scenario$G1_Direction %>% plyr::mapvalues(from = 1:3, to = c('down','null','up')) %>% factor(levels = c('down','null','up')) %>% as.character()
scenario$G2_Direction_Lab <- scenario$G2_Direction %>% plyr::mapvalues(from = 1:3, to = c('down','null','up')) %>% factor(levels = c('down','null','up')) %>% as.character()


scenario$Label = with(scenario,paste0('G1',G1_Depth,G1_Direction,
                                      'G2',G2_Depth,G2_Direction))
print(scenario)

system('mkdir Data/')

## Initial/transition probability (always begin with background)
N.states = 2^2
d1 = 1-1/(131.92967+1)  
d2 = 1-1/(43.75663+1) 

pi = c(1,0)
gamma = matrix(0,N.states,N.states)
diag(gamma) = c(d1,rep(d2,N.states-1))
gamma[1,2:N.states] = (1-d1)/(N.states-1)
gamma[gamma==0] = (1-d2)/(N.states-1)

### Looping through scenarios

for(i in 1:nrow(scenario)){
  scenario.i = scenario[i,]
  
  if(!dir.exists(file.path('Data',scenario.i$Label))){
    cmd = paste('mkdir',file.path('Data',scenario.i$Label))
    cat('Command: ',cmd)
    system(cmd)   
  }
  
  ### Simulation begins
  # Sample size and number of windows
  ngroups = scenario.i$Groups
  nreps = scenario.i$Replicates
  nwindow = scenario.i$Nwindow
  nstates = N.states
  
  #Data frame with differential paths
  difflist = NULL
  difflist[paste0('Group',1:ngroups)] = list(NULL)
  for(i in names(difflist)){difflist[[i]] = c(NA,i)}
  difflist = as.data.frame(expand.grid(difflist))
  difflist = difflist[order(rowSums(is.na(difflist)),decreasing=T),]
  difflist
  
  for(sim in 1:scenario.i$Nsim)
  {
    if(!file.exists(file.path('Data',scenario.i$Label,paste0(scenario.i$Label,'_',sim,'.RData')))){
      cat("\r",paste('Simulation: '),sim)
      
      # Simulate hidden states
      z = epigraHMM:::generateHMM(gamma,nwindow)
      
      # getting parameters
      Group1 <- pareq(depth = scenario.i$G1_Depth_Lab,direction = scenario.i$G1_Direction_Lab)
      Group2 <- pareq(depth = scenario.i$G2_Depth_Lab,direction = scenario.i$G2_Direction_Lab)
      
      # Now simulating gc
      alpha <- Group1$gc_beta
      gc <- rbeta(nwindow,alpha,alpha)
      
      # Now simulating ChIP
      chiplist = NULL
      chiplist[paste0('Group',1:ngroups)] = list(NULL)
      for(h in names(chiplist)){
        pargroup <- get(h)
        for(k in 1:nreps){
          chiplist[[h]][[paste0('Replicate',k)]] = rep(0,nwindow)
          for(j in 1:nstates){
            if(is.na(difflist[[h]][[j]])){
              chiplist[[h]][[paste0('Replicate',k)]][which(z==j)] = rnbinom(sum(z==j),mu=exp(pargroup$coeff[1,'Background']*gc[which(z==j)]^2+
                                                                                               pargroup$coeff[2,'Background']*gc[which(z==j)]+
                                                                                               pargroup$coeff[3,'Background']),size=pargroup$disp[1])
            } else{
              chiplist[[h]][[paste0('Replicate',k)]][which(z==j)] = rnbinom(sum(z==j),mu=exp(pargroup$coeff[1,'Enrichment']*gc[which(z==j)]^2+
                                                                                               pargroup$coeff[2,'Enrichment']*gc[which(z==j)]+
                                                                                               pargroup$coeff[3,'Enrichment']),size=pargroup$disp[2])
            }
          }
        }
      }
      str(chiplist)
      
      dat = list('sim'=sim,'z'=z,'ChIP'=chiplist,'gc' = gc)
      save(dat,file=file.path('Data',scenario.i$Label,paste0(scenario.i$Label,'_',sim,'.RData')),compress = "xz")
      
      df.chip <- as.data.frame(chiplist)
      dt.chip <- as.data.table(df.chip)[,Window := seq_len(.N)] %>% melt(.,id.vars = 'Window',measure.vars = c('Group1.Replicate1','Group1.Replicate2','Group2.Replicate1','Group2.Replicate2'))
      ref <- rowMeans(log1p(df.chip))
      idx <- 1:1000
      
      pdf(file = file.path('Data',scenario.i$Label,paste0(scenario.i$Label,'_',sim,'.pdf')),width = 7,height = 7)
      pairs(cbind(log1p(df.chip),gc), lower.panel = NULL, upper.panel = my_line)
      pairs(cbind(log1p(df.chip),ref), lower.panel = NULL, upper.panel = my_line)
      ggplot(data = dt.chip[Window %in% idx,],aes(x = Window,y = value))+
        facet_grid(rows = vars(variable))+
        geom_line()+theme_bw()+ylab('Simulated ChIP Counts')
      dev.off()
    }
  }
}

