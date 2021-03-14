import os
import sys
import numpy as np
import pandas as pd
import emcee
from multiprocessing import Pool
import matplotlib.pyplot as plt
from numpy import inf


os.environ["OMP_NUM_THREADS"] = "1"
os.environ['QT_QPA_PLATFORM']='offscreen'

ispec_dir = '/home/swastik/iSpec'      
sys.path.insert(0, os.path.abspath(ispec_dir))

import ispec
np.seterr(all="ignore")

df = pd.read_csv('/home/swastik/HPksiHya.txt', sep ='\s+')

df = df[df.flux != 0] 
x = df['waveobs'].values
y = df['flux'].values
yerr = df['err'].values
df = np.array(df,dtype=[('waveobs', '<f8'), ('flux', '<f8'), ('err', '<f8')])

def synthesize_spectrum(theta):
    teff ,logg ,MH = theta
#    alpha = ispec.determine_abundance_enchancements(MH)
#    microturbulence_vel = ispec.estimate_vmic(teff, logg, MH) 
#    macroturbulence = ispec.estimate_vmac(teff, logg, MH) 
    microturbulence_vel = 1.23
    macroturbulence = 3.70
    alpha = 0
    limb_darkening_coeff = 0.6
    resolution = 47000
    vsini = 2.4
    code = "grid"
    precomputed_grid_dir = ispec_dir + "/input/grid/SPECTRUM_MARCS.GES_GESv6_atom_hfs_iso.480_680nm/"
    grid = ispec.load_spectral_grid(precomputed_grid_dir)

    atomic_linelist = None
    isotopes = None
    modeled_layers_pack = None
    solar_abundances = None
    fixed_abundances = None
    abundances = None
    atmosphere_layers = None

    
    
    if not ispec.valid_interpolated_spectrum_target(grid, {'teff':teff, 'logg':logg, 'MH':MH, 'alpha':alpha, 'vmic': microturbulence_vel}):
        msg = "The specified effective temperature, gravity (log g) and metallicity [M/H] \
                fall out of the spectral grid limits."
        print(msg)
    regions = None
    # Interpolation
    synth_spectrum = ispec.create_spectrum_structure(x)
#    interpolated_spectrum = ispec.create_spectrum_structure(np.arange(wave_base, wave_top, wave_step))
    synth_spectrum['flux'] = ispec.generate_spectrum(synth_spectrum['waveobs'], \
            atmosphere_layers, teff, logg, MH, alpha, atomic_linelist, isotopes, abundances, \
            fixed_abundances, microturbulence_vel = microturbulence_vel, \
            macroturbulence=macroturbulence, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff, \
            R=resolution, regions=regions, verbose=1,
            code=code, grid=grid)
    return synth_spectrum


walkers = eval(input("Enter Walkers: "))
Iter = eval(input("Enter Iterations: "))

def log_likelihood(theta):
    model = synthesize_spectrum(theta)
    sigma2 = yerr ** 2
    p = (y - (model['flux'])) ** 2/sigma2
#    p[p == +inf] = 0
#    p[p == -inf] = 0
#    s = -0.5 * np.sum(p)
#    print(s)
#    return s
    return -0.5 * np.sum((y - (model['flux'])) ** 2/ sigma2)

def log_prior(theta):
    teff, logg, MH  = theta
    if 4000 < teff < 7500 and 1.1 < logg < 4.8 and -3.49 < MH <= 0.49 :
        return 0.0
    return -np.inf

def log_probability(theta):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta)

initial = np.array([5000,3.0,0.1])
#5873	3.98	-0.07
#5044	2.87	0.14
#4983	2.77	0.13
pos = initial + np.array([100,0.1,0.1])*np.random.randn(walkers, 3)
nwalkers, ndim = pos.shape

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability)
sampler.run_mcmc(pos,Iter, progress=True)

fig, axes = plt.subplots(3, figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
accepted = sampler.backend.accepted.astype(bool)

labels = ["teff","logg","MH","vsini"]
for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)
axes[-1].set_xlabel("step number");

plt.savefig('HPksiHyaa.pdf')

fig, ax = plt.subplots(1, figsize=(10, 7), sharex=True)
samples = sampler.flatchain
theta_max  = samples[np.argmax(sampler.flatlnprobability)]
best_fit_model = synthesize_spectrum(theta_max)
ax.plot(x,best_fit_model['flux'])
ax.plot(x,y)
plt.savefig('HPksiHyab.pdf')
print(('Theta max: ',theta_max))

#discard = input("Enter burns: ")

new_samples = sampler.get_chain(discard=150, thin=1, flat=False)
new_samples = new_samples[:,accepted,:]

fig, axes = plt.subplots(3, figsize=(10, 7), sharex=True)
for i in range(ndim):
    ax = axes[i]
    ax.plot(new_samples[:, :, i], "k", alpha=0.3)
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)
axes[-1].set_xlabel("step number")
plt.savefig('HPksiHyac.pdf')
flat_samples = new_samples.reshape(-1,new_samples.shape[2])
np.savetxt("RNHPksiHya.txt",flat_samples)

