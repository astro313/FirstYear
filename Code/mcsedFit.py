#python /Users/admin/mbb_emcee/mbb_emcee/run_mbb_emcee.py -b 100 --initT 11.2 --initBeta 1.5 --initAlpha -2.0 --initLambda0 1.0 --initFnorm 268 --get_dustmass --get_lir --get_peaklambda -v -z 2.221 -p ./ --upT 60.0 --upBeta 2.2 SMGdata.txt BLOOM.h5

#lambda0: The observer frame wavlength where the dust becomes optically thick (lambda0rest * (1+z)) in [um]

--opthin
--noalpha

import mbb_emcee
"""
Look at resulting fit from mbb_emcee
"""
res = mbb_emcee.mbb_results(h5file="BLOOM.h5")
T_val = res.par_cen('T')
print("Temperature/(1+z): {:0.2f}+{:0.2f}-{:0.2f} [K]".format(*T_val))
b_val = res.par_cen('beta')
print("Beta: {:0.2f}+{:0.2f}-{:0.2f}".format(*b_val))
res.compute_peaklambda()
p_val = res.peaklambda_cen()
print("Peak Obs wavelength: {:0.1f}+{:0.1f}-{:0.1f} [um]".format(*p_val))

# Fit results:
# T/(1+z): 18.00 +1.34 -1.59 (low lim: 1.00 upper lim: 60.00) [K]
# beta: 1.28 +0.53 -0.48 (low lim: 0.10 upper lim: 2.20) 
# fnorm: 256.57 +16.86 -16.22 (low lim: 0.03) [mJy]
# lambda0 * (1+z): 573.16 +306.99 -406.96 (low lim: 1.00 upper lim: 3000.00) [um]
# alpha: 2.65 +0.22 -0.22 (low lim: 0.10 upper lim: 20.00)
# Lambda peak: 259.0 +6.4 -6.1 [um]
# L_IR(8.0 to 1000.0um): 89.45 +2.73 -2.71 [10^12 L_sun]
# M_d(kappa=2.64, lam=125.0um): 46.28 +17.15 -15.67 [10^8 M_sun]
# Number of data points: 7
# ChiSquare of best fit point: 4.02
# Saving results to BLOOM

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




# can also read in result and comupute numbers after