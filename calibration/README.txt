Copyright (C) 2021 Swiss Tropical and Public Health Institute
 
Copyright (C) 2021 University of Warwick
 
This HAT model is free software; you can redistribute it and/or
modify it under the terms of version 2 of the GNU General Public
License as published by the Free Software Foundation.

_____

RunMCMC.R is the main script for running the adaptive MCMC algorithm.  It requires three other scripts (see below). It also requires two data files (one for active screening and one for passive detection) with following columns: "year" "peopleScreened" "hat.cases" "P1" "P2" "unstaged" "pop".

 
functions.R has functions for computing the log prior of the sampled parameters and the log likelihood, using the number of active and passive cases in each year in the data as well as the output from the ODE model. It also includes functions used for plotting the quantiles from the numerically solved ODEs with the sampled parameters from the MCMC with the true data overlaid.


loglik_related_functions.R has functions for numerically solving the ODE model given a particular set of parameters. It includes a function for solving the ODE model pre-intervention and another for after interventions commenced.


new_params_and_priors.R defines the bounds for the parameters we are fitting through the MCMC algorithm and also defines the number of years with data as well as the population size.