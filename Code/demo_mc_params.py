
import mbb_emcee
"""
Look at resulting fit from mbb_emcee
Simpliest now, basically just the demo from mbb_emcee README.md
"""
res = mbb_emcee.mbb_results(h5file="BLOOM.h5")
T_val = res.par_cen('T')
print("Temperature/(1+z): {:0.2f}+{:0.2f}-{:0.2f} [K]".format(*T_val))
b_val = res.par_cen('beta')
print("Beta: {:0.2f}+{:0.2f}-{:0.2f}".format(*b_val))
res.compute_peaklambda()
p_val = res.peaklambda_cen()
print("Peak Obs wavelength: {:0.1f}+{:0.1f}-{:0.1f} [um]".format(*p_val))

import matplotlib.pyplot as plt
import numpy as np
wave, flux, flux_unc = res.data
f = plt.figure()
ax1 = f.add_subplot(221)
p_data = ax1.errorbar(wave, flux, yerr=flux_unc, fmt='ro')
p_wave = np.linspace(wave.min() * 0.5, wave.max() * 1.5, 200)
p_fit = ax1.plot(p_wave, res.best_fit_sed(p_wave), color='blue')
ax1.set_xlabel('Wavelength [um]')
ax1.set_ylabel('Flux Density [mJy]')
ax1.set_title('Best-fitted SED')
ax2 = f.add_subplot(222)
h1 = ax2.hist(res.parameter_chain('T/(1+z)'))
ax2.set_xlabel('T / (1 + z)')
ax3 = f.add_subplot(223)
h2 = ax3.hist(res.parameter_chain('beta'))
ax3.set_xlabel(r'$\beta$')
ax4 = f.add_subplot(224)
h3 = ax4.hist(res.lir_chain)
ax4.set_xlabel(r'$L_{IR}$')
f.show()
