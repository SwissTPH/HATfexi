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

par(mfrow=c(2,1))

rpoint <- sample((samples/2):samples, size=10, replace=TRUE)

r1c_vec <- pars.stored[3, rpoint]
rdiff_vec <- pars.stored[4, rpoint]
x0_vec <- pars.stored[6, rpoint]
al_vec <- pars.stored[10, rpoint]

r2diff_vec <- pars.stored[5, rpoint]
logprop_vec <- pars.stored[8, rpoint]
prop_vec <- exp(logprop_vec)

maxr = max(r1c_vec*365*prop_vec + rdiff_vec + r2diff_vec)

par(mar=c(1.75,1.75,1.75,1.75))
par(mfrow=c(1,1))

curve(logisticfun2(x, x0=x0_vec[1], rdiff=rdiff_vec[1], rconst=r1c_vec[1]*365, al=al_vec[1]), from=0, to=x0_vec[1], xlim=c(0,17), ylim=c(0,maxr), main=paste("posterior logistic curves", suffix_save))
curve(logisticfun2(x, x0=x0_vec[1], rdiff=rdiff_vec[1], rconst=r1c_vec[1]*365, al=al_vec[1]), from=x0_vec[1], to=17, add=TRUE, lty=2)

for(i in 2:length(rpoint)){
  curve(logisticfun2(x, x0=x0_vec[i], rdiff=rdiff_vec[i], rconst=r1c_vec[i]*365, al=al_vec[i]), from=0, to=x0_vec[i], add=TRUE, col=i)
  curve(logisticfun2(x, x0=x0_vec[i], rdiff=rdiff_vec[i], rconst=r1c_vec[i]*365, al=al_vec[i]), from=x0_vec[i], to=17, add=TRUE, lty=2, col=i)
}

curve(logisticfun2(x, x0=x0_vec[1], rdiff=rdiff_vec[1] + r2diff_vec[1], rconst=r1c_vec[1]*365*prop_vec[1], al=al_vec[1]), from=0, to=x0_vec[1], xlim=c(0,17), ylim=c(0,maxr), main=paste("stage 2", suffix_save))
curve(logisticfun2(x, x0=x0_vec[1], rdiff=rdiff_vec[1]+ r2diff_vec[1], rconst=r1c_vec[1]*365*prop_vec[1], al=al_vec[1]), from=x0_vec[1], to=17, add=TRUE, lty=2)

for(i in 2:length(rpoint)){
  curve(logisticfun2(x, x0=x0_vec[i], rdiff=rdiff_vec[i]+r2diff_vec[i], rconst=r1c_vec[i]*365*prop_vec[i], al=al_vec[i]), from=0, to=x0_vec[i], add=TRUE, col=i)
  curve(logisticfun2(x, x0=x0_vec[i], rdiff=rdiff_vec[i]+r2diff_vec[i], rconst=r1c_vec[i]*365*prop_vec[i], al=al_vec[i]), from=x0_vec[i], to=17, add=TRUE, lty=2, col=i)
}