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

## Write functions for the MCMC with Simon

posfun <- function(t, y, parms, intervTime, rate_AS){
       with(as.list(y), {
             y[which(y<0)] <- 0  
             return(y)
    })
 }

logisticfun2<-function(x,x0,rdiff,rconst, al){
  ans<-rconst + rdiff*(1/(1 + exp(-al*(x-x0))))
}

takevalues2<-function(logisticfun2=logisticfun2,x0,rdiff,rconst, al, l){ #l: length of data (2000-2016 for this analysis)
  x<-numeric()
  for (i in 1:l)
    x[i]<-logisticfun2(x=i-1,x0,rdiff,rconst, al)
  x
}

## Pre-intervention dynamics (only passive detection ongoing)
fHAT_preinterv <- function(times, y, parms) {
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
    
    list(c(dSldt, dEldt, dIl1dt, dIl2dt, dRldt, dShdt, dEhdt, dIh1dt, dIh2dt, dRhdt,
           dSaldt, dEaldt, dIaldt, dRaldt, dSahdt, dEahdt, dIahdt, dRahdt,
           dSvldt, dEvldt, dIvldt, dUvldt, dSvhdt, dEvhdt, dIvhdt, dUvhdt))
  })
}

## This is valid for 2000-2016, the data period. Active screeening is a continuous rate (rAS), and value changes each year 
fHAT_withInterv <- function(times, y, parms, rateAS, r1vec, r2vec) {
  with(as.list(c(y,parms)), {
    
    integertime <- floor(times)
    specificity<-1
    
    if(times<2017) specificity <- spec ## specificity=1 from 2017 onwards

    a <- 1+integertime-startdata ## 1 corresponds to selecting value for dynamics within year 2000
    rAS<- rateAS[a]
    r1now <- r1vec[a]
    r2now <- r2vec[a]
    
    ## Define detection rates for each group
    r1h <- r1now                       ## high risk setting: only PD
    r2h <- r2now 
    r1l <- r1now + rAS*sensitivity     ## low risk setting: both PD and AS
    r2l <- r2now + rAS*sensitivity    
    
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
    
    ## Store cummulative number of reported cases by stage
    
    ## These are unscaled values, so no comparable to data
    dreportedA1dt<- rAS*sensitivity*Il1 + rAS*(1-specificity)*(Sl + El + Rl)  ##truePositives and falsePositives  
    dreportedA2dt<- rAS*sensitivity*Il2
    dreportedP1dt <- r1now*(Il1 + Ih1)
    dreportedP2dt <- r2now*(Il2 + Ih2)  
        
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
    
    list(c(dSldt, dEldt, dIl1dt, dIl2dt, dRldt, dShdt, dEhdt, dIh1dt, dIh2dt, dRhdt,
           dreportedA1dt,dreportedA2dt, dreportedP1dt, dreportedP2dt,
           dSaldt, dEaldt, dIaldt, dRaldt, dSahdt, dEahdt, dIahdt, dRahdt,
           dSvldt, dEvldt, dIvldt, dUvldt, dSvhdt, dEvhdt, dIvhdt, dUvhdt))
  })
}


## ######################################################################
## Define this function to feed the events argument in the ODEs, to reset the cumulatives to zero every year
get.events<-function(data.AS, cumvariables=c('reportedA1','reportedA2','reportedP1','reportedP2')){
  times.event<-data.AS[,'year']+1
  length.intervention<-length(times.event)
  a1<-rep(cumvariables[1],length.intervention)
  a2<-rep(cumvariables[2],length.intervention)
  a3<-rep(cumvariables[3],length.intervention)
  a4<-rep(cumvariables[4],length.intervention)
  ##
  df.events <- data.frame(var = c(rbind(a1,a2,a3,a4)), 
                          time = c(matrix(rep(times.event,4),nrow=length(cumvariables),byrow=TRUE)),
                          value = rep(0,4*length.intervention),
                          method =rep("rep", 4*length.intervention))
  return(df.events)
}


output_HAT<-function(Nt, timepre, kappa, VH1, propvh, r1c, r1diff, r2diff, x0, al, prop, spec, data.AS, data.PD, myevents, takevalues2, logisticfun2, get.events){
  ## CONSTANTS:
  daysInYear <- 365  
  alpha      <- 0.2 * daysInYear
  AH1        <- 1.35 ## from Stone (high risk setting)
  AH2        <- 1.5 ## from Stone (high risk setting)
  b          <- 0.433 ## from Stone (high risk setting). Existe data on this??
  c_h        <- 0.065 ## Proportion of bites on an infective human that lead to a mature infection in flies.Stone, HR:0.003. W:0.065, Y:0.2
  c_aL       <- 0 ## Proportion of bites on an infective animal of type i that lead to a mature infection in flies
  c_aH       <- 0 ## Proportion of bites on an infective animal of type i that lead to a mature infection in flies
  Delta      <- 0.006*daysInYear ## from Mpanya et al. 2012 
  Delta_a    <- 0.002*daysInYear ##assumed
  eta        <- 0.085*daysInYear
  f          <- 0.333 * daysInYear
  gamma      <- (1/526)*daysInYear ## 526 taken from Checchi 2015. This is s_1 in Stone
  gamma_aL   <- (1/365)*daysInYear ## this is s_ai in Stone; he doesnt show results for HR setting so I assume 1 year infective
  gamma_aH   <- (1/365)*daysInYear ## this is s_ai in Stone; he doesnt show results for HR setting so I assume 1 year infective
  mu         <- 1/(60.026*daysInYear) * daysInYear ## life expectancy for 2017: 60.026 years, from https:##data.worldbank.org/country/congo-dem-rep, accesed on October 2018.
  mu_aL      <- 0.0016 *daysInYear ## from Stone (high risk setting)
  mu_aH      <- 0.0019 *daysInYear ## from Stone (high risk setting)
  mu_gamma   <- (1/252)*daysInYear ## stage 2 duration from Checchi 2015 is 252 days
  mu_t       <- 0*daysInYear ## Death rate due to treatment.
  mu_v       <- 0.03*daysInYear ## from Rogers DJ. A general model for the African trypanosomiases. Parasitology. 1988;97(Pt 1):193-212.
  ## Nt         <- 77943 ## Mosango population by 2000
  nu         <- 0.034*daysInYear ## incubation period:29 days, Ravel S, Grebaut P, Cuisance D, Cuny G. Monitoring the developmental status of Trypanosoma brucei gambiense in the tsetse fly by means of PCR analysis of anal and saliva drops. Acta Trop. 2003;88(2):161-5. doi:10.1016/ s0001-706x(03)00191-8.    ## r1l<-r1
  sigma      <- 0.326 ## from Stone high risk setting
  sigma_aL   <- 0.8 ## from Stone high risk setting
  sigma_aH   <- 0.396 ## from Stone high risk setting
  xi         <- 0.698 ## from Stone high risk setting
  sensitivity <- 0.91
  
  startdata=head(data.AS[,'year'],1)
  enddata=tail(data.AS[,'year'],1)
  
  ## These must be within 'parms' argument in the ode() function:
  VH2 <- propvh*VH1
  Nl <- Nt/(1+kappa) ## non-communter
  Nh <- Nl*kappa     ## communter pop
  Nvl <- Nl*VH1
  Nvh <- Nh*VH2
  Nal <- Nl*AH1
  Nah <- Nh*AH2 
  
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
  
  y0.preEndemic <-  c(Sl=Sl, El=El, Il1=Il1, Il2=Il2, Rl=Rl, Sh=Sh, Eh=Eh, Ih1=Ih1, Ih2=Ih2, Rh=Rh,
                      Sal=Sal, Eal=Eal, Ial=Ial, Ral=Ral, Sah=Sah, Eah=Eah, Iah=Iah, Rah=Rah,
                      Svl=Svl, Evl=Evl, Ivl=Ivl, Uvl=Uvl, Svh=Svh, Evh=Evh, Ivh=Ivh, Uvh=Uvh)
  names(y0.preEndemic) <-  c('Sl', 'El', 'Il1', 'Il2', 'Rl',
                  'Sh', 'Eh', 'Ih1', 'Ih2', 'Rh',
                  'Sal', 'Eal', 'Ial', 'Ral', 'Sah', 'Eah', 'Iah', 'Rah',
                  'Svl', 'Evl', 'Ivl', 'Uvl', 'Svh', 'Evh', 'Ivh', 'Uvh')
  
  
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
  Times.preEndemic=seq(t0,startdata,1) ## this spans the pre-intervention period (1600-2000), run X00 dynamics years to ensure endemic equilibrium
  
  Parms<-c(daysInYear=daysInYear, alpha=alpha,AH1=AH1,AH2=AH2,b=b,c_h=c_h,c_aL=c_aL,c_aH=c_aH,Delta=Delta,Delta_a=Delta_a,eta=eta,
           f=f,gamma=gamma,gamma_aL=gamma_aL,gamma_aH=gamma_aH,mu=mu,mu_aL=mu_aL,mu_aH=mu_aH,mu_gamma=mu_gamma,
           mu_t=mu_t,mu_v=mu_v,Nt=Nt,nu=nu,sigma=sigma,sigma_aL=sigma_aL,sigma_aH=sigma_aH,xi=xi,
           VH2=VH2, Nl=Nl,Nh=Nh,Nvl=Nvl,Nvh=Nvh,Nal=Nal,Nah=Nah,
           beta_al=beta_al,beta_ah=beta_ah, beta_vl=beta_vl, beta_vh=beta_vh, r1const=r1const, r2const=r2const, sensitivity=sensitivity, spec=spec, startdata=startdata,prop=prop,
           enddata=enddata)
  
  names(Parms) <- c('daysInYear', 'alpha','AH1','AH2','b','c_h','c_aL','c_aH','Delta','Delta_a','eta',
                    'f','gamma','gamma_aL','gamma_aH','mu','mu_aL','mu_aH','mu_gamma',
                    'mu_t','mu_v','Nt','nu','sigma','sigma_aL','sigma_aH','xi',
                    'VH2', 'Nl','Nh','Nvl','Nvh','Nal','Nah',
                    'beta_al','beta_ah', 'beta_vl', 'beta_vh', 'r1const', 'r2const', 'sensitivity', 'spec', 'startdata','prop',
                    'enddata')
  ## Solve the ODEs for pre-endemic period
  EE <- ode(y=y0.preEndemic, times=Times.preEndemic, func=fHAT_preinterv, parms=Parms) 
  subEE<-tail(EE,2)
  
  ## Calculate the difference in incidence in the last two years, accept with tolerance 1E-2 per year so in 100 years the error would be <1 individual
  cond<-abs(subEE[2,'Il1'] + subEE[2,'Ih1'] - subEE[1,'Il1'] - subEE[1,'Ih1'])/Nt
  if(cond>=1e-6)
    return("equilibrium condition failed, discard this parameter set")
  
  EE <- tail(EE[,-1],1) ## keep last row (remove time value) to feed new initial conditions
  
  ## #########################################################################################################################################################
  ## If the equilibrium condition is TRUE, start simulating the intervention period (2000-2016) where passive detection is improved and active screening added
  ## #########################################################################################################################################################
  
  ## Define cumulatives to be added to the list of variables of the ODEs
  reportedA1<-0
  reportedA2<-0
  reportedP1<-0
  reportedP2<-0
  
  y0.interv <- c(EE[1:10], reportedA1, reportedA2, reportedP1, reportedP2, EE[11:(length(EE))]) ## added variables are cumulative, and reset to zero every new year via 'events' argument in the ode() function
  ## Add names:
  names(y0.interv) <-  c('Sl', 'El', 'Il1', 'Il2', 'Rl', 'Sh', 'Eh', 'Ih1', 'Ih2', 'Rh',
                         'reportedA1', 'reportedA2', 'reportedP1', 'reportedP2',
                         'Sal', 'Eal', 'Ial', 'Ral', 'Sah', 'Eah', 'Iah', 'Rah',
                         'Svl', 'Evl', 'Ivl', 'Uvl', 'Svh', 'Evh', 'Ivh', 'Uvh')
  
  ## The +1 below is because I consider due year to evaluate reported cases. So what happens from t0=2000 to t1=2001 is recorded in row time=2001 in the ode output. 
  ## I correct this just after solving the ode
  Times.interv <- seq(startdata,enddata+1, 1)  
  
  ## ############################################
  ## ACTIVE SCREENING rate vector for data period
  scaled.pop <- data.AS[,'pop']
  screened<-data.AS[,'peopleScreened']
  popLR<- Nl*scaled.pop/Nt ## scaled population in the low risk setting (it depends on kappa parameter via Nl)
  
  rateAS <-as.numeric(-log(1-(screened/popLR)))
  rateAS<-c(rateAS,tail(rateAS,1)) ## repeat last element to avoid bug in ode solver in times=finaltime
  
  ## #######################################################
  ## ANNUAL r1 PASSIVE DETECTION RATE vector for data period
  r1vec<-takevalues2(logisticfun2=logisticfun2, x0=x0, rdiff=r1diff, rconst=r1const, al=al, l=nrow(data.PD))
  
  secondterm <-takevalues2(logisticfun2=logisticfun2, x0=x0,rdiff=r2diff, rconst=r2const, al=al, l=nrow(data.PD))
  r2vec = r1vec + secondterm - r1const
  
  r1vec<-c(r1vec,tail(r1vec,1)) ## repeat last element to avoid bug in ode solver in times=finaltime
  r2vec<-c(r2vec,tail(r2vec,1)) ## repeat last element to avoid bug in ode solver in times=finaltime
  
  # ## ######################################################################
  # ## Define events, an argument to reset the cumulatives to zero every year
  myevents<-get.events(data.AS)
  
  ## Evaluate the dynamics
  ans <- ode(y=y0.interv, times=Times.interv, func=fHAT_withInterv, parms=Parms, rateAS=rateAS, r1vec=r1vec, r2vec=r2vec, events=list(data=myevents))[-1,] 
  ## Correct time values so the year and reported cases from model match data
  ans[,'time'] <- ans[,'time']-1
  
  unscaled.reportedA1 <- ans[,'reportedA1']
  unscaled.reportedA2 <- ans[,'reportedA2']
  unscaled.reportedP1 <- ans[, 'reportedP1']
  unscaled.reportedP2 <- ans[, 'reportedP2']
  
  ## Scale up reported cases:
  scale.factor <- scaled.pop/Nt
  
  scaled.reportedA1 <- tail(unscaled.reportedA1, length(scaled.pop))*scale.factor
  scaled.reportedA2 <- tail(unscaled.reportedA2, length(scaled.pop))*scale.factor
  scaled.reportedP1 <-tail(unscaled.reportedP1, length(scaled.pop))*scale.factor
  scaled.reportedP2 <-tail(unscaled.reportedP2, length(scaled.pop))*scale.factor
  
  return(cbind(Passive1=scaled.reportedP1, Passive2=scaled.reportedP2, Active1=scaled.reportedA1, Active2=scaled.reportedA2))
}


call_HAT<-function(p) {
  #takes in vector of parameters to be estimated; outputs expected number of active and passive cases in each stage
     kappa<-p[1]
     VH1<-exp(p[2])
     r1c<-p[3]
     r1diff<-p[4]
     r2diff <- p[5]
     x0<-p[6]
     spec<-exp(p[7])/(1+exp(p[7])) 
     prop <- exp(p[8])
     propvh <- exp(p[9])
     al <- p[10]
     return(output_HAT(Nt, timepre, kappa, VH1, propvh, r1c, r1diff, r2diff, x0, al, prop, spec, data.AS, data.PD, myevents, takevalues2, logisticfun2, get.events))
}
