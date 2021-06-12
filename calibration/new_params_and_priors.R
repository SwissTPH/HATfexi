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

timepre=200 ## nb of years the pre-endemic intervention period is simulated
Nt= 125691 ## This value is for Mushie health zone. This will change depending on the health zone; I am using estimations by year 2015 in all cases

intervTime <- nrow(data.AS) # 17 years for current analysis

kappa_lb <- 0
kappa_ub <- 1

logvh1_lb <- 0
logvh1_ub <- log(100)

r1c_lb <- 0
r1c_ub <- .001

r1diff_lb <- 0
r1diff_ub <- 2.5

r2diff_lb <- 0
r2diff_ub <- 2.5

x0_lb <- 0
x0_ub <- intervTime

logitspec_lb <- 6.906
logitspec_ub <- 9.21024

logprop_lb <- log(1)
logprop_ub <- log(50)

logpropvh_lb <- log(1)
logpropvh_ub <- log(50)

al_lb <- 0
al_ub <- 1

## ranges for the priors 
lowerRange <- c(kappa_lb, logvh1_lb, r1c_lb, r1diff_lb, r2diff_lb, x0_lb, logitspec_lb, logprop_lb, logpropvh_lb, al_lb)
upperRange <- c(kappa_ub, logvh1_ub, r1c_ub, r1diff_ub, r2diff_ub, x0_ub, logitspec_ub, logprop_ub, logpropvh_ub, al_ub)
names(lowerRange) <- c('kappa', 'logvh1', 'r1c','r1diff', 'r2diff',  'x0', 'logitspec', 'logprop', 'logpropvh', 'al')
names(upperRange) <- names(lowerRange)
if(prod(upperRange>lowerRange)!=1) stop('error')