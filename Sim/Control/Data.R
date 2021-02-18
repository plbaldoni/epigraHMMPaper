library(Rcpp)
library(plyr)

### Simulation scnearios
# Heterogeneity: low (L), medium (M), or high (H)
# Sample: 2, 3, 6, or 9
# Random: intercept (I), or slope (S)
# Nwindow: 25,000
# Nsim: 100

scenario = expand.grid(Heterogeneity=c('L'),Sample=c(2),Random=c('I'),Nwindow=25000,Nsim=200)
scenario$Label = with(scenario,paste0('H',Heterogeneity,'N0',Sample,'R',Random))
print(scenario)

system('mkdir Data/')

## Initial probability (always begin with background) and Transition Matrix
pi = c(1,0)
gamma = matrix(0,2,2,byrow=T)
diag(gamma) = c(0.995,0.99)
gamma[1,2] = 1 - gamma[1,1]
gamma[2,1] = 1 - gamma[2,2]

## Control Parameters: mean (mux), dispersion (phix)
# In real data (Helas3 rep 1, chr1, Control), mux = exp(2.287) =  9.845357
# In real data (Helas3 rep 1, chr1, Control), phix = 2.970831525
mux = 0.7560612#9.8453
phix = 3.230474#2.9708

## ChIP Parameters:
genpar = function(random,heterogeneity){
    if(random=='I'){
        sigma2 = ifelse(heterogeneity=='L',1e-8,ifelse(heterogeneity=='M',1e-8,ifelse(heterogeneity=='H',1e-8,NA)))
        # ChIP parameters
        b1 = c(0.7560612,2.4169519)
        b2 = c(1,0.50)
        phi = c(3.230474,2.891908)
    }
    return(list('sigma2'=sigma2,'b1'=b1,'b2'=b2,'phi'=phi))
}

# Zero-inflation
# Real data suggests p(zero inflation | background) = 0.12 and p(zero inflation | enrichment) = 0.02
# I will use p(zero inflation | background) = 0.10 and p(zero inflation | enrichment) = 0
# The following gives a 10% of zero inflation when colntrol=9 (mean)
# Intercept and Slope for logit model for background only
csi = c(-1000,-1000) 

for(i in 1:nrow(scenario)){
    scenario.i = scenario[i,]
    
    cmd = paste0('mkdir Data/',scenario.i$Label)
    cat('Command: ',cmd)
    system(cmd)
    
    ### Simulation begins
    # Sample size and number of windows
    n = scenario.i$Sample
    m = scenario.i$Nwindow
    
    # Heterogeneity
    random = scenario.i$Random
    
    # Model-specifc parameters
    par = genpar(random,scenario.i$Heterogeneity)
    b1 = par$b1
    b2 = par$b2
    phi = par$phi
    sigma2 = par$sigma2
    
    dat = data.frame()
    for(sim in 1:scenario.i$Nsim)
    {
        cat("\r",paste('Simulation: '),sim)
        # Simulate random effects
        B = rnorm(n,mean=0,sd=1)
        
        # Simulate hidden states
        Z = epigraHMM:::generateHMM(gamma,m)-1
        idx.background = which(Z==0)
        idx.enrichment = which(Z==1)
        
        # Simulate means and responses
        MU = matrix(0,nrow=m,ncol=n)
        Y = matrix(0,nrow=m,ncol=n)
        X = matrix(0,nrow=m,ncol=n)
        for(ind in 1:n){
            # Simulate control
            X[,ind] = rnbinom(m,mu=mux,size=phix)
            
            # Changing the effect of control per condition and per sample
            b2.new <- c(0,0)
            b2.new[1] <- runif(1,b2[1]-0.25,b2[1])
            b2.new[2] <- runif(1,b2[2]-0.25,b2[2])
            
            #Simulate mean
            if(random=='I'){
                MU[idx.background,ind] = exp(b1[1]+sqrt(sigma2)*B[ind]+b2.new[1]*log(X[idx.background,ind]+1))
                MU[idx.enrichment,ind] = exp(b1[2]+sqrt(sigma2)*B[ind]+b2.new[2]*log(X[idx.enrichment,ind]+1))
            }
            
            # Simulate read counts
            Y[idx.background,ind] = rnbinom(length(idx.background),mu=MU[idx.background,ind],size=phi[1])
            Y[idx.enrichment,ind] = rnbinom(length(idx.enrichment),mu=MU[idx.enrichment,ind],size=phi[2])
            
            # Adding Zero-inflation
            pziy= 1/(1+exp(-(csi[1]+csi[2]*log(X[idx.background,ind]+1))))
            idx.ziy = (runif(length(idx.background))<=pziy)
            Y[idx.background[idx.ziy==T],ind] = 0
        }
        dat = data.frame(sim=rep(sim,m),y=Y,z=Z,x=X,b=matrix(B,nrow=m,ncol=n,byrow=T),sigma2=sigma2)
        save(dat,file=paste0('./Data/',scenario.i$Label,'/',scenario.i$Label,'_',sim,'.RData'))
    }
}

