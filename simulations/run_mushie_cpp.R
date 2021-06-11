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

# This is the main script for running the HAT model. This script first 
# runs the deterministic dynamics, followed by the stochastic dynamics 
# using the Gillespie direct method, implemented in C++

# Prior to this, it is important that the scripts 
# 'make_fexi_param_matrix.R' and 'subsamplePosteriors.R' have been run,
# as they produce two files that are imported by this script.

## ######################################################################


rm(list=ls())
a<-Sys.time()
options(width=400)


#Set working directory to correct folder
WD <- '~/HATfexi/simulations' #Change this as appropriate


load('outputMu110420_1.Rdata') ## load fitted parameter sets


## Independently of period used to model calibration, here we need the most complete data set available 

# Case data at the Mushie territory level is not publicly available, but a template with missing values is provided
data.AS<-read.table('dataAS_mushie.txt') 
data.PD<-read.table('dataPD_mushie.txt')

## Load the population projections assuming a 3% annual growth rate
source('population.R') 
pop.projections<- pop.mushie
head(pop.projections) ; tail(pop.projections)


## Read subsamble of fitted parameters
setwd(WD)
pars.stored[,1]
parmat=read.table('parameters1000_mushie.txt')

parmat_rep<- matrix( rep( t( parmat ) , 28 ) , ncol = ncol(parmat) , byrow = TRUE )
colnames(parmat_rep)<-colnames(parmat)
parmat <- parmat_rep


## TIME PARAMETERS
tfdet=2000 ## end of deterministic part, pre-endemic solution, use year unit
tfstoch=51 ## end of stochastic sims, which start with t=0 which is equivalent to 2000, use integer unit

(nbyears<-nrow(data.AS))
startdata=head(data.AS[,'year'],1)  ## year 2000


library(deSolve);library(Rcpp)

sourceCpp("hat_gillespie_fexi.cpp")
SaveData<-TRUE


## ###########################################
## Take arguments from array number, if available ##
CA <- commandArgs(TRUE)
ifelse(length(CA)==0, UseScript <- FALSE, UseScript <- TRUE)
ifelse (UseScript, taskID <- as.integer(CA)[2], taskID <- 407)
ifelse (UseScript, jobID <- as.integer(CA)[1], jobID <- 1)
print(paste(jobID, taskID))


setwd(WD)
source('deterministic_mos-mush.R') # this file contains the functions for the deterministic component of the model



## For a given parameter set, first run deterministic dynamics until tfdet. Take output to then feed stochastic part starting in 2000


## #########################################################################################################################
## Define events, an argument to reset the cumulatives to zero every year to solve ODEs (deterministic, pre-endemic period)
## start.years<-startdata-timepre
Tf=tfdet ##2000
## This is to set cummulatives to zero every year:
if(FALSE){ ## si get.events2 funciona, puedo remover esta definicion de myevents
  tt<-seq(startdata-timepre,Tf,1)
  ll<-length(tt)
  a1<-rep("Incl",ll);    a2<-rep("Inch",ll)
  a3<-rep("Diel",ll);    a4<-rep("Dieh",ll)
  a5<-rep("reportedA1",ll);a6<-rep("reportedA2",ll)  
  a7<-rep("reportedP1",ll); a8<-rep("reportedP2",ll) 
  ##
  myevents <- data.frame(var = c(rbind(a1,a2,a3,a4,a5,a6,a7,a8)), 
                         time = c(matrix(rep(tt,8),nrow=8,byrow=TRUE)),
                         value = rep(0,8*ll),
                         method =rep("rep", 8*ll))
}





## RUN STOCHASTIC DYNAMICS

(iipar=taskID)

pars<-parmat[iipar,]
print(paste('Using row', taskID, 'of parmat (matrix from parameters.txt)'))
kappa <- as.numeric(pars['kappa'])
VH1 <- as.numeric(pars['VH1'])
r1c<- as.numeric(pars['r1c'])
r1diff<- as.numeric(pars['r1diff'])
r2diff<- as.numeric(pars['r2diff'])
x0<- as.numeric(pars['x0'])
spec<- as.numeric(pars['spec'])
prop<- as.numeric(pars['prop'])
propvh<- as.numeric(pars['propvh'])
al<- as.numeric(pars['al'])


## ########################################
## Run pre-endemic deterministic dynamics #
out<-runhat_det(Nt, timepre, kappa, VH1, propvh, r1c, r1diff, r2diff, x0, al, prop, spec, data.AS, data.PD, takevalues2, logisticfun2, get.events2)
##head(out);
## tail(out[,1:11]); tail(out[,12:19]); tail(out[,20:35])
values <-tail(out,1)[-1] ##remove time

## Build vector of initial conditions for stochastic simulations, x0, composed of integer values
initvalues<- round(values)
names(initvalues)<-colnames(out)[-1]
initvalues[1:10];initvalues[11:18];initvalues[19:26];initvalues[27:34]

print(paste("Incidence by 2000 is", initvalues['Incl']+initvalues['Inch'])) ## just a check, can commment it 


## #############################################################################################################
## Build parameters vector for stochastic dynamics - this is passed to propensity (transition) rates function ## 
## #############################################################################################################
if(TRUE){  
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
  ## Nt         <- 77943 ## Mushie population by 2000
  nu         <- 0.034*daysInYear ## incubation period:29 days, Ravel S, Grebaut P, Cuisance D, Cuny G. Monitoring the developmental status of Trypanosoma brucei gambiense in the tsetse fly by means of PCR analysis of anal and saliva drops. Acta Trop. 2003;88(2):161-5. doi:10.1016/ s0001-706x(03)00191-8.    ## r1l<-r1
  sigma      <- 0.326 ## from Stone high risk setting
  sigma_aL   <- 0.8 ## from Stone high risk setting
  sigma_aH   <- 0.396 ## from Stone high risk setting
  xi         <- 0.698 ## from Stone high risk setting
  sensitivity <- 0.91
  startdata=head(data.AS[,'year'],1)
  ## Screened and popAS straight from data.AS correspond to data period
  screened <- data.AS[, 'peopleScreened'] 
  popAS <- data.AS[, 'pop']
}


## WARNING: functions for stochastic dynamics ARE NOT ADAPTED to run this before 2000 

## ##########################################
## List of output variables I want to keep ##
variables.out<-c('time', 'Incl', 'Inch', 'Diel', 'Dieh', 'reportedA1', 'reportedA2', 'reportedP1', 'reportedP2')
lvar<-length(variables.out) ## to define number of columns  in the matrix output (1 matrix per state)

## ##################################################################################
##  Define variables used in stoch simulations that depend on calibrated parameters #
VH2 <- propvh*VH1
(Nl <- Nt/(1+kappa)) ## non-communter
(Nh <- Nl*kappa)     ## communter pop
Nvl <- Nl*VH1
Nvh <- Nh*VH2
Nal <- Nl*AH1
Nah <- Nh*AH2 
##
beta_al <- mu_aL*Nal
beta_ah <- mu_aH*Nah
beta_vl <- Nvl*mu_v
beta_vh <- Nvh*mu_v
##
## Constant stage-specific passive detection (annual)  rate (valid during 'timepre')
r1const <- r1c*daysInYear    
r2const <- prop*r1const      

## ###########################################################################
## Build vectors for 2000-2050 active screening and passive detection rates ##
## ###########################################################################

## Here we assume two things: check in assumption 1 and 2 below
## ##################################################################################################################
## Assumption 1: mean number of people actively screened (2014-2018 period) gives projected values for t>data period.
## Note that with a 3% growth rate, the proportion screened/populationLowRiskSetting keeps decreasing because denominator increases
pop4AS <- tail(pop.projections,tfstoch+1) ## population projections for 2000 onwards
head(pop4AS); tail(pop4AS)

mean.screen <- round(mean(tail(data.AS[, 'peopleScreened'],5)))
(max.screen <- max(data.AS[, 'peopleScreened']))
##
screened<-c(data.AS[, 'peopleScreened'], rep(mean.screen, (tfstoch+1-nbyears)))
screened.mitig2021<-c(data.AS[, 'peopleScreened'], rep(mean.screen,2), rep(max.screen, (tfstoch-nbyears-1)))  ## rep() cubre 2019 y 2020
## cbind(pop4AS,screened) ## check: NON SCALED VALUES  ## cbind(pop.mushie[401:451],screened)

## Scaled population in the low risk setting (where AS applies)
popLR<- Nl*pop4AS[,2]/Nt ## (pop4AS[,2]/Nt) is the scaling factor; it gives 1 for year 2015: check here: cbind(pop4AS[,1],pop4AS[,2]/Nt)

## DEFINE NEW VECTOR OF ANNUAL ACTIVE SCREENING RATES FOR 2000 ONWARDS for strategies with and without mitigation
rateAS <-as.numeric(-log(1-(screened/popLR)))
##rateAS<-c(rateAS,tail(rateAS,1)) ## repeat last element to avoid bug in times=finaltime
rateAS.mitigation <-as.numeric(-log(1-(screened.mitig2021/popLR)))

## ## Also, generate ratesAS data to simulate projections using max screened value from data
## screenedMax<-c(data.AS[, 'peopleScreened'], rep(max.screen, (tfstoch-nbyears)))
## rateASmax <-as.numeric(-log(1-(screenedMax/popLR)))

if(FALSE){
  ##Check: until 2020, both matrices should be the same
  cbind(2000:2051,screened, screened.mitig2021)
  cbind(2000:2051,rateAS, rateAS.mitigation)
}

if(FALSE){
  ##Check:
  popProj=data.frame(pop.projections)
  colnames(popProj)=c('year', 'pop')
  cbind(popProj[popProj$year>=2000,],rateAS) ## pop.projections=pop.mushie
}



## #######################################################################################################################################################
## Assumption 2: projected passive detection rate (r1 and r2) takes maximum from past, which here means project r1[by 2018]; same for r2.
## Note in t>Oct2020, there are two independent logistic improvements for r1 and r2 so for t>2000 in the sims, r2 is no longer just a proportion of r1 !!!!!

## note: l=nrow(data.AS) must be the length of data used in the fit
r1.dataperiod <-takevalues2(logisticfun2=logisticfun2, x0=x0, rdiff=r1diff, rconst=r1const, al=al, l=nrow(data.AS))
r1.value<- r1.dataperiod[length(r1.dataperiod)] ## take last element which SHOULD BE THE maximum value from this vector
## Doing the same for r2:
secondterm <-takevalues2(logisticfun2=logisticfun2, x0=x0,rdiff=r2diff, rconst=r2const, al=al, l=nrow(data.AS))
r2.dataperiod = r1.dataperiod + secondterm - r1const
r2.value <- r2.dataperiod[length(r2.dataperiod)] ## take last element which SHOULD BE THE maximum value from this vector

## Define new vector of stage 1 annual passive detection rate for 2000 onwards
r1vec<-c(r1.dataperiod, rep(r1.value,(tfstoch+2-nbyears))) ## vector must span 2000-2051, so vector length should be 53
r2vec<-c(r2.dataperiod, rep(r2.value,(tfstoch+2-nbyears))) ## vector must span 2000-2051


##add in fexinidazole parameters
mat<-read.table('mat_indices_5col.txt', header=TRUE)

## Read parameters that define the dynamics
val<- as.numeric(mat[taskID,])
print(paste('current parameter set given by index ', val[1]))
print(paste('fexinidazole parameters (psi, phi1, phi2, pdfactor) are', val[2], val[3],val[4], val[5])) 
par<-val[1]
psi<-val[2]
phi1<-val[3]
phi2<-val[4]
pdfactor<-val[5]


iiParams <- list()
iiParams <- within(iiParams, {
  #parameters listed in order of which they are imported into the C++ function
  rateAS=rateAS 
  r_pd1=r1vec 
  r_pd2=r2vec 
  #the following are fexi parameters
  psi=psi;  phi1=phi1;  phi2=phi2; pdfactor=pdfactor; 
  
  alpha=alpha
  b=b
  c_h=c_h; c_aL=c_aL; c_aH=c_aH
  delta=Delta; delta_a=Delta_a
  eta=eta
  f=f
  gamma=gamma; gamma_aL=gamma_aL; gamma_aH=gamma_aH
  mu=mu;  mu_aL=mu_aL;mu_aH=mu_aH; mu_gamma=mu_gamma; mu_t=mu_t; mu_v=mu_v
  nu=nu; 
  sensitivity=sensitivity
  sigma=sigma; sigma_aL=sigma_aL; sigma_aH=sigma_aH
  specificity=spec
  xi=xi
  
  maxyears=years4proj
  
  initHuman <- within(list(), {
    ## HUMAN POPULATION
    S_L=initvalues["Sl"] 
    E_L=initvalues["El"]
    I1_L=initvalues["Il1"] 
    I2_L=initvalues["Il2"]
    T_L=initvalues["Rl"]
    Incidence_1L=initvalues["Incl"] 
    Active1=initvalues["reportedA1"]
    Active2=initvalues["reportedA2"]
    Passive_1=initvalues["reportedP1"]
    Passive_2=initvalues["reportedP2"]    
    #Passive_1L=initvalues["reportedP1"]
    #Passive_2L=initvalues["reportedP2"]
    
    S_H=initvalues["Sh"]
    E_H=initvalues["Eh"]
    I1_H=initvalues["Ih1"]
    I2_H=initvalues["Ih2"]
    T_H=initvalues["Rh"]
    Incidence_1H=initvalues["Inch"]
    #Passive_1H=round(treated_Passive1.H[leng])
    #Passive_2H=round(treated_Passive2.H[leng])
  })
  
  initVector <- within(list(), {
    ## Vector POPULATION
    SvL=initvalues["Svl"]
    EvL=initvalues["Evl"]
    IvL=initvalues["Ivl"]
    UvL=initvalues["Uvl"]
    
    SvH=initvalues["Svh"]
    EvH=initvalues["Evh"]
    IvH=initvalues["Ivh"]
    UvH=initvalues["Uvh"]
  })
  
  initAnimal <- within(list(), {
    ## Animal POPULATION
    SaL=initvalues["Sal"]
    EaL=initvalues["Eal"]
    IaL=initvalues["Ial"]
    RaL=initvalues["Ral"]
    
    SaH=initvalues["Sah"]
    EaH=initvalues["Eah"]
    IaH=initvalues["Iah"]
    RaH=initvalues["Rah"]
  }) 
  
  
})

##########################################################################
## Launch C++ stochastic simulations
##########################################################################
## To add column names 

stochColnames <- c('time',
                   'byYear_Incid_1L', 'byYear_Incid_1H',
                   'byYear_Diagnosis_1L','byYear_Diagnosis_2L','byYear_Diagnosis_1H','byYear_Diagnosis_2H',
                   'byYear_Active1_fexi', 'byYear_Active1_curr', 'byYear_Active2_fexi', 'byYear_Active2_curr',
                   'byYear_Passive_1_fexi', 'byYear_Passive_1_curr','byYear_Passive_2_fexi',
                   'byYear_Passive_2_curr')


## ####################################################################
## Set the dir where subdirectories (1 per scenarios) will be stored ##
## ####################################################################
if(SaveData){
  setwd('stoch_outputs')
  dirStoch <- paste0('taskID_',taskID)
  if(file.exists(dirStoch)==FALSE)
    dir.create(dirStoch)
  setwd(dirStoch)
}
getwd()

for (ss in 1:100){ ## 100 seeds per full parameter set
  ## Run stochastic dynamics
  set.seed(ss)
  ans <- as.matrix(hat_gillespie_fexi(iiParams))
  colnames(ans)<-stochColnames
  saveRDS(ans, file=paste0('sim_', ss,'.rds'))
}


b <- Sys.time()

# This gives how long the code took to run
(a)
(b)

