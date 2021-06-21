# iMCSPEC
[![Build Status](https://img.shields.io/badge/release-1.0.0-orange)](https://github.com/SwastikC/iMCSpec)
[![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-371/)

iMCSpec is a tool which combines iSpec(https://www.blancocuaresma.com/s/iSpec) and emcee(https://emcee.readthedocs.io/en/stable/) into a single unit to perform Bayesian analysis of spectroscopic data to estimate stellar parameters. For more details on the individual code please refer to the links above. This code have been tested on Synthetic dataset as well as GAIA BENCHMARK stars (https://www.blancocuaresma.com/s/benchmarkstars). The example shown here is for the grid genarated MARCS.GES_atom_hfs. If you want to use any other grid, just download it from the https://www.cfa.harvard.edu/~sblancoc/iSpec/grid/ and make the necessary changes in the line_regions.

To run the iMCSpec.py file in a server where you can use MPI : Type the following code -- > time mpirun -n (number of cores,say 6) ipython filename.py 


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

Please check iMCspec.ipynb for complete details

For MPI pool and other advanced features the notebook is not yet updated. The code iMCspec.py is the latest one with MPI pool feature, thanks to Sonith LS for the contribution. You can use that feature in notebook as well by just copying that part. You can also give your suggestions and imporvements for the code...
Happy coding !! 
