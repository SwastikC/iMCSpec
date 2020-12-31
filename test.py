import os
import sys
import numpy as np
import pandas as pd
import emcee
import matplotlib.pyplot as plt
os.environ['QT_QPA_PLATFORM']='offscreen'
ispec_dir = '/home/swastik/mcspec'
sys.path.insert(0, os.path.abspath(ispec_dir))
import ispec

df = pd.read_csv('/home/swastik/ROXH12.txt', sep ='\s+')
x = df['waveobs'].values
y = df['flux'].values
yerr = df['err'].values

def synthesize_spectrum(theta,code="moog"):
    teff ,logg ,MH,vsini = theta
    resolution = 48000


    microturbulence_vel = ispec.estimate_vmic(teff, logg, MH)
    macroturbulence = ispec.estimate_vmac(teff, logg, MH)

    limb_darkening_coeff = 0.6
    regions = None
    wave_base = 609.0
    wave_top = 620.0

    model = ispec_dir + "/input/atmospheres/ATLAS9.Castelli/"

    atomic_linelist_file = ispec_dir + "/input/linelists/transitions/VALD.300_1100nm/atomic_lines.tsv"
    #atomic_linelist_file = ispec_dir + "/input/linelists/transitions/GESv5_atom_hfs_iso.420_920nm/atomic_lines.tsv"

    isotope_file = ispec_dir + "/input/isotopes/SPECTRUM.lst"
    alpha = ispec.determine_abundance_enchancements(MH)

    atomic_linelist = ispec.read_atomic_linelist(atomic_linelist_file, wave_base=wave_base, wave_top=wave_top)
    atomic_linelist = atomic_linelist[atomic_linelist['theoretical_depth'] >= 0.01]

    isotopes = ispec.read_isotope_data(isotope_file)

    solar_abundances_file = ispec_dir + "/input/abundances/Asplund.2009/stdatom.dat"
    
    solar_abundances = ispec.read_solar_abundances(solar_abundances_file)
    fixed_abundances = None
    modeled_layers_pack = ispec.load_modeled_layers_pack(model)
    atmosphere_layers = ispec.interpolate_atmosphere_layers(modeled_layers_pack, {'teff':teff, 'logg':logg, 'MH':MH, 'alpha':alpha}, code=code)
    synth_spectrum = ispec.create_spectrum_structure(x)
    synth_spectrum['flux'] = ispec.generate_spectrum(synth_spectrum['waveobs'],
            atmosphere_layers, teff, logg, MH, alpha, atomic_linelist, isotopes, solar_abundances,
            fixed_abundances, microturbulence_vel = microturbulence_vel,
            macroturbulence=macroturbulence, vsini=vsini, limb_darkening_coeff=limb_darkening_coeff,
            R=resolution, regions=regions, verbose=0,
            code=code)

    return synth_spectrum


walkers = eval(input("Enter Walkers: "))
Iter = eval(input("Enter Iterations: "))
# lT,uT = input("Enter The Temperature Bounds: ")
# lg,ug = input("Enter log g bounds: ")
# lm,um = input("Enter [M/H] bounds: ")
# lv,uv = input("Enter Vsini : ")

def log_likelihood(theta, x, y, yerr):
    model = synthesize_spectrum(theta,code="moog")
    sigma2 = yerr ** 2
    return -0.5 * np.sum((y - (model['flux'])) ** 2/ sigma2)

def log_prior(theta):
    teff, logg, MH , vsini = theta
    if 3500 < teff < 4900 and 3.2 < logg < 4.8 and -1.49 < MH <= 0.49  and 1 < vsini < 10 :
        return 0.0
    return -np.inf

def log_probability(theta, x, y, yerr):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, yerr)

initial = np.array([4100,4.0,0.20,7])

pos = initial + np.array([100,0.5,0.5,0.3])*np.random.randn(walkers, 4)
nwalkers, ndim = pos.shape

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(x, y, yerr))
sampler.run_mcmc(pos,Iter, progress=True)

fig, axes = plt.subplots(4, figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
accepted = sampler.backend.accepted.astype(bool)

labels = ["teff","logg","MH","vsini"]
for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)

axes[-1].set_xlabel("step number");
plt.savefig('Roxs1.pdf')
fig, ax = plt.subplots(1, figsize=(10, 7), sharex=True)

samples = sampler.flatchain
theta_max  = samples[np.argmax(sampler.flatlnprobability)]
best_fit_model = synthesize_spectrum(theta_max)
ax.plot(x,best_fit_model['flux'])
ax.plot(x,y)
plt.savefig('Roxs2.pdf')
print(('Theta max: ',theta_max))


#discard = input("Enter burns: ")

new_samples = sampler.get_chain(discard=150, thin=1, flat=False)
new_samples = new_samples[:,accepted,:]
fig, axes = plt.subplots(4, figsize=(10, 7), sharex=True)
for i in range(ndim):
    ax = axes[i]
    ax.plot(new_samples[:, :, i], "k", alpha=0.3)
    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)
axes[-1].set_xlabel("step number")
plt.savefig('Roxs3.pdf')
flat_samples = new_samples.reshape(-1,new_samples.shape[2])
np.savetxt("RNRoxsH.txt",flat_samples)

