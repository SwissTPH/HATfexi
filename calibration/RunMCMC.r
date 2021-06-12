# ######################################################################
# Copyright (C) 2021 Swiss Tropical and Public Health Institute
# 
# Copyright (C) 2021 University of Warwick
# 
# This HAT model is free software; you can redistribute it and/or
# modify it under the terms of version 2 of the GNU General Public
# License as published by the Free Software Foundation.
# 
# #######################################################################

# Authors: Simon Spencer, Soledad Castano, Swati Patel

library(lhs); library(deSolve); library(rootSolve); library(mvtnorm)
# date 08/10/20

idtask=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
#idtask=1
suffix_save <- sprintf("Mu110420_%g", idtask)

set.seed(idtask)

data.AS<-read.table('dataAS_mushie.txt',header=T)
data.PD<-read.table('dataPD_mushie.txt',header=T)
source('new_params_and_priors.R') # this file gives demographic parameter and upper, lower bounds for fitted parameters
source('functions.R')
source('loglik_related_functions.R')  # this file has the ODEs

par(mar=c(2,2,1.5,.75))

overdispersion<-TRUE 
plots<-FALSE
suffix<-sprintf("Mu110320_%g", idtask)  # this is used to load a previously determined empirical covariance, comment out for initial run

# MCMC control parameters
burnin<-2000 # to find the mode
samples<-10000
thinning<- 2
iterations<-burnin + samples*thinning 

# initialise MCMC from random starting location
X<-"initial"
while (length(X)==1) {
  pars <- drop(fFullSampleLHS(NSAMPLES=1, lowerRange=lowerRange, upperRange=upperRange))  # hypercube sampling for initial parameters
  lpars <- length(pars)
  names(pars) <- names(lowerRange)
  pars.mean<-pars
  # fit model to initial parameters
  X<-call_HAT(pars)
               
  if (overdispersion) {
    pars.od<-c(2,2) # starting values
    names(pars.od)<-c("kappa.PD","kappa.AS") # overdispersion parameter
    lowerRange.od<-c(0,0) 
    upperRange.od<-c(Inf,Inf)
    sigma.od<-c(0.2,0.2) 
    accept.od<-rep(0,2) # monitors
    reject.od<-rep(0,2)
    pars.od.stored<-matrix(NA,2,samples) # stored samples
    target.ar.od<-0.44 # target acceptance rate
    x.od<-function(it) {1+20/(20+it)} # function controls speed of adaptation: must tend to one. 
  }
  if (length(X)==1) {
    cat("initialisation:",X,"\n")
  } else {
    if (overdispersion) {
      lpost<-log.posterior.od(X,pars,pars.od)
    } else {
      lpost<-log.posterior(X,pars) # calculate log likelihood plus log prior
    }
    cat("initial log posterior:",lpost,"\n")
  }
}

# algorithmic parameters for proposals and adaptation
if (file.exists(paste0("Sigma",suffix,".Rdata"))) { # if possible, use sigma from a previous run
  load(paste0("Sigma",suffix,".Rdata"))
  cat("using covariance matrix from", suffix, "\n")
  lambda<-1 # joint proposal variance inflation factor
  min.lambda<-1 # don't shrink by too much!
  target.ar<- 0.234 # target acceptance rate for joint proposals
  n0<-100 # weight of initial covariance matrix is n0/(it+n0); weight of empirical covariance matrix is it/(it+n0).
  pars.Sigma.x<-function(it) {1+10/(10+it)} # speed of adaptation function for joint updates (must tend to 1 asymptotically)
} else {
  pars.Sigma<-diag((upperRange-lowerRange)^2/400) # covariance matrix of joint proposal
  lambda<-1 # joint proposal variance inflation factor
  min.lambda<-0.5 # allow the proposal to shrink to prevent getting stuck.
  target.ar<- 0.1 # target acceptance rate for joint proposals
  n0<-10 # weight of initial covariance matrix is n0/(it+n0); weight of empirical covariance matrix is it/(it+n0).
  pars.Sigma.x<-function(it) {1+500/(500+it)} # speed of adaptation function for joint updates (must tend to 1 asymptotically)
}
pars.sigma<-(upperRange-lowerRange)/4 # single site update proposal sd
pars.sigma.x<-2 # speed of adaptation function (constant here during burnin) for single site updates
target.ar.single<-0.44 # target acceptance rate for single proposals

# set counters
accept<-0 # joint counters
reject<-0
accept.reject<-c() # joint counters
pars.accept<-rep(0,lpars) # single site counters
pars.reject<-rep(0,lpars)

# storage
pars.stored<-matrix(NA,lpars,samples)
lambda.stored<-rep(NA,iterations)
X.stored<-array(NA,c(dim(X),samples))
lpost.stored<-rep(NA,samples)
pars.stored.all<-matrix(NA,lpars,1+burnin) # required for estimating Sigma, only initially
pars.stored.all[,1]<-pars # store initial value  

print(Sys.time())
for (it in 1:iterations) {
  # single site updates only during burnin to get to posterior mode
  if (it<=burnin) {
    for (j in 1:lpars) {
      proposal<-pars
      proposal[j]<-rnorm(1,pars[j],pars.sigma[j])
      if (proposal[j]>=lowerRange[j] & proposal[j]<=upperRange[j]) {
        X.proposed<-call_HAT(proposal)
        if (length(X.proposed)==1 | min(X.proposed)<0) {
          cat("Error with par",j,X.proposed,"\n")
          lpost.proposed<--Inf
        } else if (overdispersion) {          
          lpost.proposed <- log.posterior.od(X.proposed,proposal,pars.od)
        } else {
          lpost.proposed <- log.posterior(X.proposed,proposal)
        }               
        lar<-lpost.proposed-lpost
      } else {
        lar<--Inf # always reject
      }
      u<-runif(1)
      if (log(u)<lar) {#accept
        pars[j]<-proposal[j]
        X<-X.proposed
        lpost<-lpost.proposed
        pars.accept[j]<-pars.accept[j]+1
        pars.sigma[j]<-pars.sigma[j]*pars.sigma.x # things are going well - inflate the proposal sd
      } else { # reject
        pars.reject[j]<-pars.reject[j]+1
        pars.sigma[j]<-pars.sigma[j]*pars.sigma.x^(target.ar.single/(target.ar.single-1)) # things are going badly - shrink the proposal sd
      }
    }
  }
  # joint proposals
  proposal<-drop(rmvnorm(1,pars,2.38^2/lpars*lambda^2*pars.Sigma))
  if (prod(proposal>=lowerRange & proposal<=upperRange)) {
    X.proposed <- call_HAT(proposal)
    if (length(X.proposed)==1 | min(X.proposed)<0) {
       cat("Error",X.proposed,"\n")
       lpost.proposed<--Inf
    } else if (overdispersion) {          
      lpost.proposed <- log.posterior.od(X.proposed,proposal,pars.od)
    } else {
      lpost.proposed <- log.posterior(X.proposed,proposal)
    }        
    lar<-lpost.proposed-lpost
  } else {
    lar<--Inf # always reject
  }
  u<-runif(1)
  if (log(u)<lar) {#accept
    pars<-proposal
    X<-X.proposed
    lpost<-lpost.proposed
    accept<-accept+1
    accept.reject<-c(accept.reject, 1)
    lambda<-lambda*pars.Sigma.x(it) # things are going well - inflate all the proposal sds
  } else { # reject
    reject<-reject+1
    accept.reject<-c(accept.reject, 0)
    lambda<-max(min.lambda,lambda*pars.Sigma.x(it)^(target.ar/(target.ar-1))) # things are going badly - shrink all the proposal sds
  }
  if (overdispersion) { # update overdispersion parameters
    for (j in 1:2) {
      proposal<-pars.od
      proposal[j]<-rnorm(1,pars.od[j],sigma.od[j])
      if (lowerRange.od[j]<=proposal[j] & upperRange.od[j]>=proposal[j]) {
        lpost.proposed<-log.posterior.od(X,pars,proposal)
        lar<-lpost.proposed-lpost
      } else {
        lar<--Inf # reject
      }
      u<-runif(1)
      if (log(u)<lar) { # accept
        pars.od[j]<-proposal[j]
        lpost<-lpost.proposed
        accept.od[j]<-accept.od[j]+1
        sigma.od[j]<-sigma.od[j]*x.od(it)
      } else { # reject
        reject.od[j]<-reject.od[j]+1
        sigma.od[j]<-sigma.od[j]*x.od(it)^(target.ar.od/(target.ar.od-1))
      }
    }
  }
  # housekeeping at the end of the iteration
  if (it<=burnin) {pars.stored.all[,1+it]<-pars}
  # update the estimates of the covariance matrix
  if (it<=2*burnin && it%%2==0) { # new observation replaces the oldest
    pars.mean.new<-pars.mean+(pars-pars.stored.all[,floor(it/2)])/(it-floor(it/2)+1) # need to remove burnin everywhere!
    pars.Sigma<-pars.Sigma+(tcrossprod(pars)-tcrossprod(pars.stored.all[,floor(it/2)])+(it-floor(it/2)+1)*(tcrossprod(pars.mean)-tcrossprod(pars.mean.new)))/(it-floor(it/2)+n0)
    pars.mean<-pars.mean.new        
  } else { # gain an observation
    pars.mean.new<-(pars.mean*(it-floor(it/2)) + pars)/(it-floor(it/2)+1)
    pars.Sigma<-((it-1-floor(it/2)+n0)*pars.Sigma+tcrossprod(pars)+(it-floor(it/2))*tcrossprod(pars.mean)-(it-floor(it/2)+1)*tcrossprod(pars.mean.new))/(it-floor(it/2)+n0)
    pars.mean<-pars.mean.new
  }
  # check the covariance matrix calculation went ok in the burnin
  if (it==burnin) {
    xbar.check<-apply(pars.stored.all[,1+floor(it/2):it],1,mean)   # this check only works if Sigma0 is default.
    Sigma.check<-cov(t(pars.stored.all[,1+floor(it/2):it]))*(it-floor(it/2))/(it-floor(it/2)+n0)+diag((upperRange-lowerRange)^2/400)*n0/(it-floor(it/2)+n0)
    cat(it,"Xbar check for iterative formulae:",xbar.check-pars.mean,"\n")
    cat(it,"Sigma check for iterative formulae:",Sigma.check-pars.Sigma,"\n")
  }
  # take a sample after the burnin
  if (it>burnin && (it-burnin)%%thinning==0) {
    pars.stored[,(it-burnin)/thinning]<-pars
    X.stored[,,(it-burnin)/thinning]<-X
    lpost.stored[(it-burnin)/thinning]<-lpost
    if (overdispersion) {pars.od.stored[,(it-burnin)/thinning]<-pars.od}
  }
  lambda.stored[it]<-lambda
  if (overdispersion) {
    cat(it,lpost,lambda,round(accept/(accept+reject),4),";",pars,pars.od,"\n")
  } else {
    cat(it,lpost,lambda,round(accept/(accept+reject),4),";",pars,"\n")
  }
  #cat(lambda,round(accept/(accept+reject),3),";",round(pars.accept/(pars.accept+pars.reject),3),"\n")
 # if (plots && it%%thinning==0) {plot.fit(X,pars)}
}
print(Sys.time())

save.image(file=paste0("output",suffix_save,".Rdata"))
save(pars.Sigma,file=paste0("Sigma",suffix_save,".Rdata")) 

par(mfrow=c(1,1))
matplot(c(1:(1+burnin)), t(pars.stored.all[c(1:lpars),]), t="l", xlab="val", col=c(1:lpars), lty=c(1:lpars))
legend("topright",names(pars),col=c(1:lpars), lty=c(1:lpars))
matplot(c(1:samples), t(pars.stored[c(1:lpars),]), t="l", col=c(1:lpars), lty=c(1:lpars))

plot.alltraces()
pdf(file=paste0("plot.",suffix_save,".traces.pdf"))
plot.alltraces()
dev.off()

plot.traces()

plot.posterior.fit.od(X.stored, pars.stored, pars.od.stored)
pdf(file=paste0("plot.",suffix_save,".postfit.pdf"))
plot.posterior.fit.od(X.stored, pars.stored, pars.od.stored)
dev.off()

pdf(file=paste0("plot.",suffix_save,".stagefit.pdf"))
plot.stage.fit.od(X.stored, pars.stored, pars.od.stored)
dev.off()

pdf(file=paste0("plot.",suffix_save,".scatter.pdf"))
source('function_scatterplots.R')
dev.off()

pdf(file=paste0("plot.",suffix_save,".acceptance.pdf"))
slide.window.ar(accept.reject, 1000)
dev.off()

pdf(file=paste0("plot.", suffix_save, ".density.pdf"))
source('function_densityplots.R')
dev.off()

pdf(file=paste0("plot.", suffix_save, ".logistic.pdf"))
source('logistic_postsamples.R')
dev.off()

pdf(file=paste0("plot.", suffix_save, ".lpost.pdf"))
hist(lpost.stored[(burnin+1):iterations])
dev.off()