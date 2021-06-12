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

library(parallel)

log.prior <- function(pars, pars.od){
  # sums all prior values
  # Needs to be edited if prior distributions change
    dnorm(pars['logvh1'], mean = 1.1, sd=.05, log=TRUE) + 
      dgamma(pars['x0'], shape = 10, scale=.06, log=TRUE) +
           dgamma(pars['logprop'], shape=1, scale=1, log=TRUE) + 
              dgamma(pars['logpropvh'], shape=1, scale=1, log=TRUE) + 
            dgamma(pars.od[1], shape= 23.5, scale=3, log=TRUE) + dgamma(pars.od[2], shape=23.5, scale=3, log=TRUE)
  # Note: does not include priors for uniform
}


log.negative.binomial <- function(X, d, od_shape){
  #X is an array of proportion of expected detected diseases stages 1 and 2, high and low risk, and false positives and negatives; from the ODE output)
  # X has nrows of years
  # X had two columns, one for each stage
  sum(dnbinom(d[,'hat.cases'], mu=X[,1]+X[,2], size=od_shape,log=TRUE) + dbinom(d[,'P1'], d[,'P1']+d[,'P2'], X[,1]/(X[,1]+X[,2]),log=TRUE))
}

log.posterior.od <-function(X,pars,pars.od) {
  # X is matrix of expected total passive (row 1) and active cases (row 2)
  # parameters of model
  # overdispersion parameters
  log.negative.binomial(X[,1:2], data.PD, pars.od[1]) + log.negative.binomial(X[,3:4], data.AS, pars.od[2]) + log.prior(pars, pars.od)  
}

plot.fit<-function(X,pars) {
  # X is matrix 
  par(mfrow=c(2,1))
  yrs<-length(X[,1])
  plot(1:yrs,X[,2],t="l",main="Active screening",ylab="all cases",
       ylim=c(0,max(X[,2],data.AS[,3])))
  points(1:yrs,data.AS[,3])
  plot(1:yrs,X[,1],t="l",main="Passive detection",ylab="all cases",
       ylim=c(0,max(X[,1],data.PD[,3])))
  points(1:yrs,data.PD[,3])
}

plot.alltraces<-function() {
  lpars<-length(pars.stored[,1])
  samples<-length(pars.stored[1,])
  par(mfrow=c(4,2))
  for (j in 1:lpars) {
    plot(1:samples,pars.stored[j,],t="l",xlab="",ylab=round(pars.accept[j]/(pars.accept[j]+pars.reject[j]),3),
         main=names(pars)[j])
  }
  plot(1:length(lambda.stored),lambda.stored,t="l",xlab="",ylab=round(accept/(accept+reject),3),main="lambda")
}

plot.traces<-function(){
  lpars<-length(pars.stored[,1])
  samples<-length(pars.stored[1,])
  par(mfrow=c(1,1))
  for (j in 1:lpars) {
    pdf(file=paste0("plot.", names(pars)[j], ".pdf"))
    plot(1:samples,pars.stored[j,],t="l",xlab="",ylab=round(pars.accept[j]/(pars.accept[j]+pars.reject[j]),3),
       main=names(pars)[j])
    dev.off()
  }
  
  for (jj in 1:length(pars.od)) {
    pdf(file=paste0("plot.", names(pars.od)[jj], ".pdf"))
    plot(1:samples,pars.od.stored[jj,],t="l",xlab="",ylab=round(accept.od[jj]/(accept.od[jj]+reject.od[jj]),3),
         main=names(pars.od)[jj])
    dev.off()
  }
    
}


plot.posterior.fit.od<-function(X.stored,pars.stored,pars.od.stored,n.cores=1,verbose=FALSE) {
  par(mfrow=c(2,1))
  years<-data.AS[,1]
  yrs<-length(years)
  samples<-length(X.stored[1,1,])
  quantiles<-c(0.025,0.25,0.5,0.75,0.975)
  qtl<-matrix(NA,length(quantiles),yrs)
  #cat(date(),"\n")
  if (verbose) {cat("Starting qnbinom:",date(),"\n")}
  for (i in 1:length(quantiles)) {
    qtl[i,]<-apply(matrix(qnbinom(quantiles[i],mu=X.stored[,3,]+X.stored[,4,],size=rep(pars.od.stored[2,],each=yrs)),yrs,samples),1,mean)
  }
  if (verbose) {cat("Completed qbetabinom:",date(),"\n")}   
  plot(years,t="n",main="Active screening",ylab="all cases",ylim=c(0,max(qtl,data.AS[,'hat.cases'])),xlab="year",xlim=c(min(years),max(years)))
  polygon(c(years,rev(years)),c(qtl[1,],rev(qtl[5,])),border=NA,col="lightgrey")
  polygon(c(years,rev(years)),c(qtl[2,],rev(qtl[4,])),border=NA,col="darkgrey")
  lines(years,qtl[3,],lwd=2)
  points(years,data.AS[,'hat.cases'],col="red")

  #
  qtl<-matrix(NA,length(quantiles),yrs) 
  for (i in 1:length(quantiles)) {
    qtl[i,]<- apply(matrix(qnbinom(quantiles[i],mu=X.stored[,1,]+X.stored[,2,],size=rep(pars.od.stored[1,],each=yrs)),yrs,samples),1,mean)
  }  
  plot(years,t="n",main="Passive detection",ylab="all cases",ylim=c(0,max(qtl,data.PD[,'hat.cases'])),xlab="year",xlim=c(min(years),max(years)))
  polygon(c(years,rev(years)),c(qtl[1,],rev(qtl[5,])),border=NA,col="lightgrey")
  polygon(c(years,rev(years)),c(qtl[2,],rev(qtl[4,])),border=NA,col="darkgrey")
  lines(years,qtl[3,],lwd=2)
  points(years,data.PD[,'hat.cases'],col="red")
}

plot.stage.fit.od<-function(X.stored,pars.stored,pars.od.stored,n.cores=1,verbose=FALSE) {
  par(mfrow=c(2,1))
  years<-data.AS[,1]
  yrs<-length(years)
  samples<-length(X.stored[1,1,])
  quantiles<-c(0.025,0.25,0.5,0.75,0.975)
  qtl<-matrix(NA,length(quantiles),yrs)
  for (i in 1:length(quantiles)) {
    qtl[i,]<-apply(matrix(qbinom(quantiles[i],rep(data.AS[,'P1']+data.AS[,'P2'],samples), X.stored[,3,]/(X.stored[,3,]+X.stored[,4,])),yrs,samples),1,mean)
    #  cat(date(),"\n") 
  }
  plot(years,t="n",main="Stage 1 from AS",ylim=c(0,max(qtl,data.AS[,'P1'])),ylab="stage 1 cases",xlab="year",xlim=c(min(years),max(years)))
  polygon(c(years,rev(years)),c(qtl[1,],rev(qtl[5,])),border=NA,col="lightgrey")
  polygon(c(years,rev(years)),c(qtl[2,],rev(qtl[4,])),border=NA,col="darkgrey")
  lines(years,qtl[3,],lwd=2)
  points(years,data.AS[,'P1'],col="red")
  
  qtl<-matrix(NA,length(quantiles),yrs)
  for (i in 1:length(quantiles)) {
    qtl[i,]<-apply(matrix(qbinom(quantiles[i],rep(data.PD[,'P1']+data.PD[,'P2'],samples),X.stored[,1,]/(X.stored[,1,]+X.stored[,2,])),yrs,samples),1,mean) 
  }
  plot(years,t="n",main="Stage 1 from PD",ylim=c(0,max(qtl,data.PD[,'P1'])),ylab="stage 1 cases",xlab="year",xlim=c(min(years),max(years)))
  polygon(c(years,rev(years)),c(qtl[1,],rev(qtl[5,])),border=NA,col="lightgrey")
  polygon(c(years,rev(years)),c(qtl[2,],rev(qtl[4,])),border=NA,col="darkgrey")
  lines(years,qtl[3,],lwd=2)
  points(years,data.PD[,'P1'],col="red")
}

slide.window.ar<- function(accept.reject, window.size){
 win.ar <- c()
 
 for(k in (window.size+1):length(accept.reject)){
 win.ar <- c(win.ar, sum(accept.reject[(k-window.size):k])/window.size)
 }
 plot((window.size+1):length(accept.reject), win.ar)
}

## Function fFullSampleLHS ## useful when testing the functions below
fFullSampleLHS <-  function(NSAMPLES,lowerRange,upperRange) {
  require(lhs)
  llr <- length(lowerRange); lup <- length(upperRange)
  if(llr!=lup) stop('mismatch between lowerRange and upperRange lengths')
  nParameterSample <- length(lowerRange)
  ## Draw uniform samples in the range [0,1] using LHS
  Lhs <- randomLHS(NSAMPLES,nParameterSample)
  ## Rescale according to lower & upper ranges
  ans <- t(apply(Lhs, 1, function(x) {(upperRange - lowerRange) * x + lowerRange}))
}

