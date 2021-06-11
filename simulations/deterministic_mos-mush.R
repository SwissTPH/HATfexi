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

# This script contains the functions for running the deterministic
# component of the model

## ######################################################################
## Define this function to feed the events argument in the ODEs, to reset the cumulatives to zero every year
get.events2<-function(timepre, startdata, cumvariables=c('Incl', 'Inch', 'Diel', 'Dieh', 'reportedA1','reportedA2','reportedP1','reportedP2')){
    tt<-seq(startdata-timepre,startdata,1)
    ll<-length(tt)
    a1<-rep(cumvariables[1],ll)
    a2<-rep(cumvariables[2],ll)
    a3<-rep(cumvariables[3],ll)
    a4<-rep(cumvariables[4],ll)
    a5<-rep(cumvariables[5],ll)
    a6<-rep(cumvariables[6],ll)
    a7<-rep(cumvariables[7],ll)
    a8<-rep(cumvariables[8],ll)
    df.events <- data.frame(var = c(rbind(a1,a2,a3,a4,a5,a6,a7,a8)), 
                            time = c(matrix(rep(tt,8),nrow=8,byrow=TRUE)),
                            value = rep(0,8*ll),
                            method =rep("rep", 8*ll))
    return(df.events)
}


if(FALSE){
    xx<-0:16
    myx<-seq(0,17,1)
    plot(myx,logisticfun(myx,x0,K,r1const),lwd=1, col='gray')
    (y<-takevalues(logisticfun=logisticfun,x0=x0,K=K,r1const=r1const,l=l))
    points(0:16, takevalues(logisticfun,x0,K,r1const,l=17),col='red')
    lines(xx,y,type='s',lwd=1)
    abline(h=min(y))
    abline(h=max(y))  
    cbind(2000:2016,y)
    range(y)
}



## october 2020
## Readapt fHAT_preinterv from the MCMC results folder to run pre-intervention dynamics (only passive detection ongoing)
fHAT_preinterv_deterministic<- function(times, y, parms) {
    with(as.list(c(y,parms)), {     
        ## Detection rates for each group - only PD for this period
        r1l <- r1const
        r1h <- r1const                       
        r2l <- r2const
        r2h <- r2const     
        ##Inlines
        ## Allow for disease induced deaths to cycle back in (e.g., births make up for deaths so population sizes remain constant)
        betal <- mu*(Sl + El + Il1 + Il2 + Rl) + mu_gamma*Il2 + mu_t*Rl
        betah <- mu*(Sh + Eh + Ih1 + Ih2 + Rh) + mu_gamma*Ih2 + mu_t*Rh
        ##
        ## Human population exposed to bites
        Nl_exposed <- Nl-Rl  ## total minus removed (either die or get treated)
        Nh_exposed <- Nh-Rh 
        ##
        ## terms for biting probabilities
        theta_vLhL <- (sigma*Nl_exposed)/(sigma*(Nl_exposed+((1-xi)*Nh_exposed))+sigma_aL*Nal)
        theta_vLhH <- (sigma*(1-xi)*Nh_exposed) / (sigma*(Nl_exposed+((1-xi)*Nh_exposed))+sigma_aL*Nal)
        theta_vLaL <- (sigma_aL*Nal)/(sigma*(Nl_exposed+((1-xi)*Nh_exposed))+sigma_aL*Nal)
        theta_vHhH <- (sigma*xi*Nh_exposed)/((sigma*xi*Nh_exposed) + (sigma_aH*Nah))
        theta_vHaH <- (sigma_aH*Nah)/((sigma*xi*Nh_exposed) + (sigma_aH*Nah))
        ##
        ## HUMANS
        lambdaL <- (b*f*theta_vLhL)/Nl_exposed 
        lambdaH1 <- (b*f*theta_vLhH)/Nh_exposed
        lambdaH2 <- (b*f*theta_vHhH)/Nh_exposed 
        ##
        ##VECTORS
        c1 <- (f*theta_vLhL*c_h)/Nl_exposed
        c2 <- (f*theta_vLhH*c_h)/Nh_exposed
        c3 <- (f*theta_vLaL*c_aL)/Nal
        Lambda_vL <- c1*(Il1 + Il2) + c2*(Ih1 + Ih2) + c3*Ial
        c4 <- (f*theta_vHhH*c_h)/Nh_exposed
        c5 <- (f*theta_vHaH*c_aH)/Nah
        Lambda_vH <- c4*(Ih1 + Ih2) + c5*Iah
        ##
        ##ANIMALS
        lambda_aL <- (b*f*theta_vLaL*c_aL)/Nal 
        lambda_aH <- (b*f*theta_vHaH*c_aH)/Nah 
        ##    
        ##ODES
        dSldt  <- betal + Delta*Rl - mu*Sl - lambdaL*Ivl*Sl
        dEldt  <- lambdaL*Ivl*Sl - (mu+eta)*El
        dIl1dt <- eta*El - (mu+gamma+r1l)*Il1 
        dIl2dt <- gamma*Il1 - (mu+mu_gamma+r2l)*Il2
        dRldt  <- r1l*Il1 + r2l*Il2 - (mu+mu_t+Delta)*Rl
        ##
        ## 2nd host - humans
        dShdt  <- betah + Delta*Rh - mu*Sh - (lambdaH1*Ivl + lambdaH2*Ivh)*Sh
        dEhdt  <- (lambdaH1*Ivl + lambdaH2*Ivh)*Sh - (mu+eta)*Eh
        dIh1dt <- eta*Eh - (mu+gamma+r1h)*Ih1 
        dIh2dt <- gamma*Ih1 - (mu+mu_gamma+r2h)*Ih2
        dRhdt  <- r1h*Ih1 + r2h*Ih2 - (mu+mu_t+Delta)*Rh
        ##    
        ## Store cummulative number of reported cases by stage
        dIncldt <- eta*El
        dInchdt <- eta*Eh
        dDieldt<- mu_gamma*Il2
        dDiehdt<- mu_gamma*Ih2
        dreportedA1dt<- 0 ## rAS*sensitivity*Il1 + rAS*(1-specificity)*(Sl + El + Rl)  ##truePositives and falsePositives
        dreportedA2dt<- 0 ## rAS*sensitivity*Il2
        dreportedP1dt <- r1h*(Il1 + Ih1)
        dreportedP2dt <- r2h*(Il2 + Ih2)
        ##    
        ##  Non-human vertebrate host - Area 1
        dSaldt <- beta_al + Delta_a*Ral - mu_aL*Sal - lambda_aL*Ivl*Sal
        dEaldt <- lambda_aL*Ivl*Sal - (mu_aL+eta)*Eal
        dIaldt <- eta*Eal - (mu_aL+gamma_aL)*Ial
        dRaldt <- gamma_aL*Ial - (mu_aL+Delta_a)*Ral
        ##  Non-human vertebrate hosst - Area 2
        dSahdt <- beta_ah +Delta_a*Rah - mu_aH*Sah - lambda_aH*Ivh*Sah
        dEahdt <- lambda_aH * Ivh*Sah - (mu_aH+eta)*Eah
        dIahdt <- eta*Eah - (mu_aH+gamma_aH)*Iah
        dRahdt <- gamma_aH*Iah - (mu_aH+Delta_a)*Rah
        ##
        ## Vectors - Area 1
        dSvldt <- beta_vl - Lambda_vL*Svl - mu_v*Svl - alpha*Svl   
        dEvldt <- Lambda_vL*Svl - (mu_v+nu)*Evl
        dIvldt <- nu*Evl - mu_v*Ivl
        dUvldt <- alpha*Svl - mu_v*Uvl
        ##
        ## Vectors - Area 2
        dSvhdt <- beta_vh - Lambda_vH*Svh - mu_v*Svh - alpha*Svh   
        dEvhdt <- Lambda_vH*Svh - (mu_v+nu)*Evh
        dIvhdt <- nu*Evh - mu_v*Ivh
        dUvhdt <- alpha*Svh - mu_v*Uvh
        ##    
        list(c(dSldt, dEldt, dIl1dt, dIl2dt, dRldt, dShdt, dEhdt, dIh1dt, dIh2dt, dRhdt,
               dIncldt, dInchdt, dDieldt, dDiehdt, 
               dreportedA1dt,dreportedA2dt, dreportedP1dt, dreportedP2dt,
               dSaldt, dEaldt, dIaldt, dRaldt, dSahdt, dEahdt, dIahdt, dRahdt,
               dSvldt, dEvldt, dIvldt, dUvldt, dSvhdt, dEvhdt, dIvhdt, dUvhdt))
    })
}

## october 2020
## puedo remover data.PD de los argumentos, no lo necesito
## run pre-interventions period, timepre until 2000, using function fHAT_preinterv_deterministic. Ouput from here feed then the IC of stoch sims
runhat_det<-function(Nt, timepre, kappa, VH1, propvh, r1c, r1diff, r2diff, x0, al, prop, spec, data.AS, data.PD, takevalues2, logisticfun2, get.events2){ ## usar myevents from run_mosango.R si se complica
    ## CONSTANTS:
    daysInYear <- 365  
    alpha      <- 0.2 * daysInYear
    AH1        <- 1.35
    AH2        <- 1.5 
    b          <- 0.433
    c_h        <- 0.065
    c_aL       <- 0 
    c_aH       <- 0 
    Delta      <- 0.006*daysInYear
    Delta_a    <- 0.002*daysInYear
    eta        <- 0.085*daysInYear
    f          <- 0.333 * daysInYear
    gamma      <- (1/526)*daysInYear
    gamma_aL   <- (1/365)*daysInYear
    gamma_aH   <- (1/365)*daysInYear
    mu         <- 1/(60.026*daysInYear) * daysInYear
    mu_aL      <- 0.0016 *daysInYear
    mu_aH      <- 0.0019 *daysInYear
    mu_gamma   <- (1/252)*daysInYear
    mu_t       <- 0*daysInYear 
    mu_v       <- 0.03*daysInYear
    ## Nt         <- 77943 ## Mosango population by 2000
    nu         <- 0.034*daysInYear
    sigma      <- 0.326 
    sigma_aL   <- 0.8 
    sigma_aH   <- 0.396
    xi         <- 0.698
    sensitivity <- 0.91
    startdata=head(data.AS[,'year'],1)
    ## These must be within 'parms' argument in the ode() function:
    VH2 <- propvh*VH1
    Nl <- Nt/(1+kappa) ## non-communter
    Nh <- Nl*kappa     ## communter pop
    Nvl <- Nl*VH1
    Nvh <- Nh*VH2
    Nal <- Nl*AH1
    Nah <- Nh*AH2 
    ##  
    ## Set pre-intervention initial conditions
    Sl  <- Nl
    El  <- 0
    Il1 <- 0
    Il2 <- 0
    Rl  <- 0
    ##
    Sh  <- Nh
    Eh  <- 0
    Ih1 <- 0
    Ih2 <- 0
    Rh  <- 0
    ##  
    ## Cumulatives for treated people and annual incidence. If nothing given, this is initialised to zero
    Incl<-0
    Inch<-0
    Diel <-0
    Dieh <-0
    reportedA1<-0
    reportedA2<-0
    reportedP1<-0
    reportedP2<-0
    ##
    ##animals
    Sal <- Nal
    Eal <- 0 
    Ial <- 0
    Ral <-0
    Sah <- Nah
    Eah <- 0
    Iah <- 0 
    Rah <-0
    ## vectors
    Svl <- Nvl*0.95
    Evl <- Nvl*0.02
    Ivl <- Nvl*0.02
    Uvl <- Nvl*0.01
    Svh <- Nvh*0.95
    Evh <- Nvh*0.02
    Ivh <- Nvh*0.02
    Uvh <- Nvh*0.01
    ##  
    y0 <-  c(Sl=Sl, El=El, Il1=Il1, Il2=Il2, Rl=Rl, Sh=Sh, Eh=Eh, Ih1=Ih1, Ih2=Ih2, Rh=Rh,
             Incl=Incl, Inch=Inch, Diel=Diel, Dieh=Dieh, 
             reportedA1=reportedA1, reportedA2=reportedA2, reportedP1=reportedP1, reportedP2=reportedP2,
             Sal=Sal, Eal=Eal, Ial=Ial, Ral=Ral, Sah=Sah, Eah=Eah, Iah=Iah, Rah=Rah,
             Svl=Svl, Evl=Evl, Ivl=Ivl, Uvl=Uvl, Svh=Svh, Evh=Evh, Ivh=Ivh, Uvh=Uvh)
    ## These betas, that compensate deaths to keep population constant, are constant and go to 'parms' argument of ode() function
    beta_al <- mu_aL*Nal
    beta_ah <- mu_aH*Nah
    beta_vl <- Nvl*mu_v
    beta_vh <- Nvh*mu_v
    ## Constant stage-specific annual passive detection rate
    r1const <- r1c*daysInYear            
    r2const <- prop*r1const        
    ## Time
    t0=startdata-timepre
    Times=seq(t0,startdata,1) ## this spans the pre-intervention period (1600-2000), run X00 dynamics years to ensure endemic equilibrium
    ##  
    Parms<-c(daysInYear=daysInYear, alpha=alpha,AH1=AH1,AH2=AH2,b=b,c_h=c_h,c_aL=c_aL,c_aH=c_aH,Delta=Delta,Delta_a=Delta_a,eta=eta,
             f=f,gamma=gamma,gamma_aL=gamma_aL,gamma_aH=gamma_aH,mu=mu,mu_aL=mu_aL,mu_aH=mu_aH,mu_gamma=mu_gamma,
             mu_t=mu_t,mu_v=mu_v,Nt=Nt,nu=nu,sigma=sigma,sigma_aL=sigma_aL,sigma_aH=sigma_aH,xi=xi,
             VH2=VH2, Nl=Nl,Nh=Nh,Nvl=Nvl,Nvh=Nvh,Nal=Nal,Nah=Nah,
             beta_al=beta_al,beta_ah=beta_ah, beta_vl=beta_vl, beta_vh=beta_vh, r1const=r1const, r2const=r2const, sensitivity=sensitivity, spec=spec, startdata=startdata,prop=prop)##, enddata=enddata)
    ##
    myevents=get.events2(timepre, startdata)
    ## Solve the ODEs for pre-endemic period
    ans <- ode(y=y0, times=Times, func=fHAT_preinterv_deterministic, parms=Parms, events=list(data=myevents))## pop es p/detPassive
    return(ans)
}


if(FALSE){
    startdata
    r1now <- r1const*(1 + (K/(1 + exp(-(times -(startdata+x0)))))) ## when data starts (so when r1(t) starts changing, the exponential argument must be -(x0+1))
    curve(r1const*(1 + (K/(1 + exp(-(x -(startdata+x0+1)))))),from=1995, to=2017,xlim=c(1995,2017))
    abline(h=r1const,col='red')
    abline(v=startdata)
    abline(v=startdata+x0)
}



