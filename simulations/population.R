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

# This script contains the functions for producing population projections
# for various DRC health zones assuming a 3% annual growth rate

## ######################################################################


project.pop<-function(timepre,y,x.2014){
  a<-b<-numeric()
  ## to project back
  for (ii in 1:(14+timepre))
    a[ii]<- x.2014/1.03^ii
  ##
  ## to project forward
  for (ii in 1:(y-14))
    b[ii]<- x.2014*1.03^ii
  ##  
  pop<-matrix(c(seq(2000-timepre,2000+y,1),round(c(rev(a),x.2014, b))), ncol=2, byrow=FALSE)
  return(pop)
}


## Generate projections for different DRC health zone populations using OCHA 2014 estimations
## and assuming an annual 3% growth rate

years4proj=51 ## starting from 2000


health.zone=c('Bagata', 'Budjala', 'Bominenge', 'Mbaya', 'Mosango', ' Mushie')
pop=c(161155, 125766, 152259, 64522, 117896, 122030)
(pop.2014<-rbind(health.zone, pop))
## startdata=2000


## health.zone=c('Mosango', 'Bagata',  'Mushie')
## pop=c(117896, 161155, 122030)
## (pop.2014<-rbind(health.zone, pop))

pop.list<-list()
for (ii in 1:length(health.zone)){
    pop.list[[ii]]<-project.pop(timepre=timepre, y=years4proj, x.2014=as.numeric(pop.2014[2,ii]))
    attributes(pop.list[[ii]])$healthzone<-health.zone[ii]
}


## bagata
pop.bagata = pop.list[[1]]

## budjala
pop.budjala = pop.list[[2]]

## bominenge
pop.bominenge = pop.list[[3]]

## mbaya
pop.mbaya = pop.list[[4]]

## mosango
pop.mosango = pop.list[[5]]

##mushie
pop.mushie = pop.list[[6]]


