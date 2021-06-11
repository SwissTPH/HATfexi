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

# This script produces a parameter table, saved as a .txt file, that 
# contains a subsample of values from the posterior distributions of 
# parameters fitted by MCMC 

## ######################################################################


## print(Sys.time())  
rm(list=ls())
options(width=400)


## There are 10,000 posterior samples. Select a subsample for running stocastic sims
fsample <- function(lpost.stored, nsamples) { ## ll dataframe with 'np' and 'values' columns
  ll=lpost.stored
  sumll <- sum(ll)
  cumll <- cumsum(ll)
  ## Creates indices that correspond to parameters to be used for resampling.
  sampleIndices <- numeric(nsamples)
  for (ii in 1:nsamples) {
    x  <- runif(1)
    sc <- x*sumll ## Scale random number
    ## Find all indices higher than random number & pick first element
    sampleIndices[ii]  <- which(cumll <= sc)[1]  
  }
  return(sampleIndices)
}

invlogit<-function(y){
  ans<-  exp(y)/(1+exp(y))
  return(ans)
}

# Please change the working directory (WD) to where the repository is saved
WD='~/HATfexi/simulations'
setwd(WD)
load('outputMu110420_1.Rdata') ## load fitted parameter sets
ls();Nt
timepre


## NOTE: lpost.stored  matches pars.stored

## matrix with fitted parameters
parsmat=t(pars.stored)
dim(parsmat)
nrow.parsmat=nrow(parsmat)
colnames(parsmat)=names(pars)
parsmat[1:2,]
##(mypar=parsmat[nrow.parsmat,])
##summary(parsmat)

## vector of loglikelihood,  
if(length(lpost.stored)==dim(parsmat)[1]) print('loglik and params match')


## Sub-sample:
indices<-fsample(lpost.stored,nsamples=1000)
head(indices); tail(indices)
length(unique(indices))


## Compare histogram of full sample and subsample:
par(mfrow=c(2,1))
hist(lpost.stored)
hist(lpost.stored[indices])
##hist(loglik[indices[1:100]])


x=parsmat[indices,]
(a=colnames(x))

xVH1=exp(x[,'logvh1'])
xSPEC=invlogit(x[,'logitspec'])
xPROP=exp(x[,'logprop'])
xPROPVH=exp(x[,'logpropvh'])


y=cbind(x[,1],xVH1, x[,3:6],xSPEC,xPROP,xPROPVH,x[,10])
colnames(y)=c(a[1], 'VH1',a[3:6],'spec','prop','propvh',a[10])
head(y);tail(y)
a


getwd()
setwd(WD)
write.table(y, 'parameters1000_mushie.txt')

