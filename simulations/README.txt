Copyright (C) 2021 Swiss Tropical and Public Health Institute
 
Copyright (C) 2021 University of Warwick
 
This HAT model is free software; you can redistribute it and/or
modify it under the terms of version 2 of the GNU General Public
License as published by the Free Software Foundation.

_____

The scripts in this folder take the calibrated paramter distributions from the calibration folder and 
use them to run HAT simulations using both deterministic and stochastic dynamics. 

The calibration section must be run first, and the main output used is the file called 'outputMu110420_1.RData'.
This output should be put in the simulation folder before running the simulation scripts.

The files 'make_fexi_param_matrix.R' and 'subsamplePosteriors.R' must be run first before running the main file:
'run_mushie_cpp.R'

The other R scripts, and the C++ script, contain functions imported into 'run_mushie_cpp.R'.

The data on people screened and diagnosed is not publicly available, and so 'dataAS_mushie.txt' and 
'dataPD_mushie.txt' are not populated. Instead, they are included as a template for the data files. 