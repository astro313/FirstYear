# Radio core fit 


"""
Replaced by SEDfit.py


Extrapolate radio galaxy steep spectrum
and flat spectrum of core using upper limit 0.17 mJy @ 4.86GHz (mullin 2006) and 0.8mJy @ 9 GHz (Haas 2014)

Radio galaxy has lobe synchrotron, cold and warm flit, see Figure 13 in Cleary 2007

CARMA continuum @ 104.2106 GHz = 7.18 for integrated (total contribution)
= 5.56 mJy for peak
= 1.62 mJy = integrated - peak
"""


from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import k, h, c 
import matplotlib

matplotlib.rc('xtick', labelsize=18) 
matplotlib.rc('ytick', labelsize=18) 
font = {'family' : 'normal', 
        'weight' : 'bold', 
        'size'   : 16}
matplotlib.rc('font', **font)


# SMG Data
wavelg_micron=np.array([3.6, 4.5, 70, 100, 160, 250, 350, 500, 1000])
wavelg = wavelg_micron * 1.e-6	# meters
freq_SI=c/wavelg 	# Hz
flux=np.array([0.275, 0.283, 29.5, 102, 289, 440, 403, 268, 51])		# mJy
flux_SI = flux/1000.*10**(-26)

# Radio data
# Added two points from Haas Table 1
# 0.702, 2.124 micron = 4.273e14, 1.412e14
# 0.0253, 0.193 mJy 

# different from the SED fitted in VLA proposal, here: removed 1.49e10 Hz point to fit lobe synchrontron, 
# Cleary 2007: warm dust T = 185 K and cold dust fit = 31 K 

# Added Our 104.2106 GHz point = 7.18 mJy
# Total radio contribution
NED_freq = np.array([4.61e14, 4.273e14, 1.412e14, 1.27e13, 4.20e12, 1.92e12, 1.042106e11, 1.07e10, 1.07e10, 5.00e9, 5.00e9, 2.7e9, 2.7e9, 1.4e9, 1.4e9, 1.4e9, 7.5e8, 7.5e8, 7.5e8, 3.52e8, 3.52e8, 1.78e8, 1.78e8, 1.52e8, 1.52e8, 8.6e7, 7.38e7, 3.8e7, 3.78e7]) # Hz
NED_Jy = np.array([2.60e-5, 2.53e-5, 1.193e-4, 4.8e-3, 2.07e-2, 1.36e-1, 7.18e-3, 2.7e-1, 2.53e-1, 6.4e-1, 6.36e-1, 1.33, 1.34, 2.95, 2.99, 2.80, 5.94, 5.94, 5.6, 1.13e1, 1.16e1, 1.57e1, 1.71e1, 2.26e1, 2.25e1, 5.16e1, 3.75e1, 4.96e1, 6.07e1])	# Jy

# 3C220.3 Core
# 4.86 GHz: no core detected, lobes (steep spectrum dominated), upper limit of core = 0.17 mJy
# 9 GHz: Core and lobes detected: core @ 0.8 mJy.. Based on Haas VLA
# 104.2106 GHz: Core = peak = 5.56 mJy
Core_3c220dot3_Jy = np.array([5.56e-3, 0.8e-3, 0.17e-3])		# Jy
Core_3c220dot3_freq = np.array([1.042106e11, 9.e9, 4.86e9])	# Hz

# compare to similiar FR II galaxies
# 3C6.1: 
# 7.3 mJy @ 8GHz and 10.14mJy at 1.4 GHz
Core_3c6dot1_Jy = np.array([7.3e-3, 10.14e-3]) 		# Jy
Core_3c6dot1_freq = np.array([8e9, 1.4e9])  # Hz
# 3C 263: 
# 161.1 mJy @ 5 GHz and 0.119 Jy @ 8GHz
Core_3c263_Jy = np.array([161.1e-3, 0.119])
Core_3c263_freq = np.array([5e9, 8e9]) 		# Hz

# =======================================================
# #                 _          __       _       _        
#   ___ _ __   __| |   ___  / _|   __| | __ _| |_ __ _ 
#  / _ \ '_ \ / _` |  / _ \| |_   / _` |/ _` | __/ _` |
# |  __/ | | | (_| | | (_) |  _| | (_| | (_| | || (_| |
#  \___|_| |_|\__, _|  \___/|_|    \__, _|\__, _|\__\__, _|
# =======================================================                                                     


# Fit power law to core and get spectral index
# extrapolated to 107 Ghz and 24 micron and 70 micron ??

def BB_freq(freq, T):
	boltz=k*T
	b = 2*h/c**2 * freq**3/(np.exp(h*freq/boltz)-1)
	return b

def graybody_freq(freq, T, beta, freq_0):
	boltz=k*T
	b_nu = 2*h/c**2 * freq**3/(np.exp(h*freq/boltz)-1)
	tau = (freq/freq_0)
	MBB = (1-np.exp(tau**(-beta))) * b_nu
	return MBB


def B_nu(freq, freq_0, C):
	T=36
	beta=1.5
	boltz=k*T
	b_nu = 2*h/c**2 * freq**3/(np.exp(h*freq/boltz)-1)
	tau = (freq/freq_0)**(beta)
	return (1-np.exp(-tau)) * b_nu * C 

def synch(freq, alpha=0.8):
#	alpha=0.75
	del_al=0.5
	freq_b = 1.5e9
	S_ = freq**(-alpha)/ (1+(freq/freq_b)**del_al)
	#S_em = S * (1+z)**(alpha-1)		#  k_correct  to emitted frequency
	return S_

def powerlaw(freq, a, alpha):
	"""
	a is scaling factor
	"""
	powerfit = a*(freq**(alpha))
	return powerfit

def lobeSyn(freq, a):
	"""
	ref: equation 1 in CLEARY 2007
	input: freq in Hz
	a is scaling factor
	"""
	beta = 0.19
	logFluxDensity = (-beta) * (np.log10(freq) - np.log10(3e6))**2 + np.log10(np.exp(-freq/3e12))		# Table 5 equation from CLEARY 2007
	return a*10.**logFluxDensity


#Radio galaxy BB = dust portion part
#radio_coeff1, radio_err1 = curve_fit(, NED_freq[:6], NED_Jy[:6]*10**(-26))			# Dust 

#radio_fit1 = BB_freq(NED_freq[:6], radio_coeff1[0])


FluxDensity = lobeSyn(NED_freq[6:], 65000.)	# Jy
radio_dot75 = powerlaw(NED_freq[6:]/1e9, 5203.6435, -0.75)		# using spectral index 0.75


plt.loglog(NED_freq[6:]/1e9, FluxDensity, 'r--', label='lobe syn fit')
plt.loglog(NED_freq[6:]/1e9, radio_dot75, 'b--', label=r'$\alpha$=0.75')
plt.loglog(NED_freq[6:]/1e9, NED_Jy[6:]*1000, 'k.', label='extended')
plt.loglog(NED_freq[6:7]/1e9, 1.62, 'k.') 		# extended component from CARMA continuum
plt.loglog(NED_freq[6:7]/1e9, NED_Jy[6:7]*1000., '*', label='Our point')
plt.errorbar(NED_freq[6:7]/1e9, NED_Jy[6:7]*1000., yerr=0.5, marker='.', color='k', linestyle='none')
plt.loglog(Core_3c220dot3_freq/1e9, Core_3c220dot3_Jy*1000, 'o-', label='3c220.3 core data')
plt.loglog(Core_3c263_freq/1e9, Core_3c263_Jy*1000, 'o-', label='3c263 core data')
plt.loglog(Core_3c6dot1_freq/1e9, Core_3c6dot1_Jy*1000, 'o-', label='3c6.1 core data points')


plt.xlabel('freq [GHz]', fontsize=20,  fontweight='bold')
plt.ylabel('Flux Density [mJy]', fontsize=20,  fontweight='bold')
plt.grid()
#plt.legend(loc='best')
plt.show()
import sys
sys.exit()




#Power Law synchrotron SMG
freq_ = np.array([27e9, 35.7706e9])
syn_ = synch(freq_)*10**8
#plt.loglog(freq_/1e9, syn_, 'o')
#freq__ = np.linspace(0.1e9, 100.e9, 1000)
#syn__ = synch(freq__)*10**8
#plt.loglog(freq__/1e9, syn__, '-.', label='Synchrotron', linewidth=4)


# MBB
fitParams, fitCovariances  =  curve_fit(B_nu, freq_SI[-4:], flux_SI[-4:], [1.5e12, 1])
#print fitParams
y_SMG = B_nu(freq_SI[-4:], fitParams[0], fitParams[1]) #*1e26*1000.
x_extend = np.linspace(10e9, 300e9, 1000)
y_SMG_extend = B_nu(x_extend, fitParams[0], fitParams[1])
y_SMG_100ghz = B_nu(104.2106e9, *fitParams)
print y_SMG_100ghz*1e26*1000			# for comparison with radio galaxy at 104GHz
plt.loglog(freq_SI[:]/1e9, flux[:], 'o', linewidth=4)
plt.loglog(freq_SI[-4:]/1e9, y_SMG*1e26*1000, 'k-', linewidth=5, label='MBB')
plt.loglog(x_extend/1e9, y_SMG_extend*1e26*1000, 'k--', linewidth=3)


fitpowerSMG, fitpowerSMGCov = curve_fit(powerlaw, freq_SI[2:-3]/1e9, flux_SI[2:-3]*1e26*1000)		# fitting to GHz and mJY
print fitpowerSMG

SMG_powerfit = powerlaw(freq_SI[2:-3]/1e9, *fitpowerSMG)
plt.loglog(freq_SI[2:-3]/1e9, SMG_powerfit, label='power blue')		# Ghz and mJy

plt.xlabel('freq [GHz]', fontsize=20,  fontweight='bold')
#plt.ylim(50**-1, 10**3)
plt.title('SMM J0939+8315', fontsize=20,  fontweight='bold')
#plt.legend(loc='best')
plt.ylabel('Flux Density [mJy]', fontsize=20,  fontweight='bold')
plt.grid()
plt.show()

