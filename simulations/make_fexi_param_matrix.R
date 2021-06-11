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

# This script makes mat_indices_5col.txt
# This .txt file contains the values of the fexi parameters of interest

## ######################################################################

WD <-'~/HATfexi/simulations'
setwd(WD)

# Reduced set of fexi parameters of interest
reducedpars <- rbind(c(0.25,0.25,0.75),
                     c(0.25,0.75,0.5),
                     c(0.25,1,0.25),
                     c(0.25,1,0.75),
                     c(0.5,0.75,0.5),
                     c(0.75,0.75,0.5),
                     c(1,0.75,0.5))
pdpars <- c(1.20,1.50,2)
m <- length(pdpars)

index.samples<-1:1000
l<-length(index.samples)
c1<-rep(index.samples,times=dim(reducedpars)[1]*m)
c2<-rep(reducedpars[,1],each=l, times=m) 
c3<-rep(reducedpars[,2],each=l, times=m) 
c4<-rep(reducedpars[,3],each=l, times=m)
c5<-rep(pdpars, each=dim(reducedpars)[1]*l)
mat.cases<-cbind(c1,c2,c3,c4,c5)

colnames(mat.cases)<- c("parindex","psi","phi1","phi2","pdfactor")
write.table(mat.cases, 'mat_indices_5col.txt')