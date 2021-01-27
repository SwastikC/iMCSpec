# MCMCSPEC
iMCSpec is a tool which combines iSpec(https://www.blancocuaresma.com/s/iSpec) and emcee(https://emcee.readthedocs.io/en/stable/) into a single unit to perform Bayesian analysis of spectroscopic data to estimate stellar parameters. For more details on the individual code please refer to the links above. This code have been tested on Syntehtic dataset as well as GAIA BENCHMARK stars (https://www.blancocuaresma.com/s/benchmarkstars). The example shown here is for the grid genarated MARCS.GES_atom_hfs. If you want to use any other grid, just download it from the https://www.cfa.harvard.edu/~sblancoc/iSpec/grid/ and make the necessary changes in the line_regions.

Dependencies :
iSpec : iSpec v2020.10.01
python = 3+ Branch
numpy>=1.18.5
scipy>=1.5.0
matplotlib>=3.2.2
astropy>=4.0.1.post1
lockfile>=0.12.2
Cython>=0.29.21
pandas>=1.0.5
statsmodels>=0.11.1
dill>=0.3.2
emcee>=3.0.2
corner>=2.1.1

Let us import all the necessary packages that are required for this analysis. 
import os
import sys
import numpy as np
import pandas as pd
import emcee
from multiprocessing import Pool
import matplotlib.pyplot as plt

os.environ["OMP_NUM_THREADS"] = "1"
os.environ['QT_QPA_PLATFORM']='offscreen'
os.environ["NUMEXPR_MAX_THREADS"] = "8"      #CHECK NUMBER OF CORES ON YOUR MACHINE AND CHOOSE APPROPRIATELY 

ispec_dir = '/home/swastik/iSpec'               #MENTION YOUR DIRECTORY WHERE iSPEC is present     
sys.path.insert(0, os.path.abspath(ispec_dir))

import ispec
#np.seterr(all="ignore")                     #FOR MCMC THE WARNING COMES FOR RED BLUE MOVES WHEN ANY PARTICULAR WALKER VALUE DONOT LIE IN THE PARAMETER SPACE
