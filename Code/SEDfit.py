from scipy.optimize import curve_fit
import scipy.optimize as opt
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import k, h, c 
import matplotlib
import astropy.io.ascii

matplotlib.rc('xtick', labelsize=18) 
matplotlib.rc('ytick', labelsize=18) 
font = {'family' : 'normal', 
        'weight' : 'bold', 
        'size'   : 12}
matplotlib.rc('font', **font)


SMGfilename = 'SMGdata.txt'
RadioFilename = 'Radiodata.txt'

def read_data_fromtxt(filename):
	data = astropy.io.ascii.read(filename, comment='^#')
	if len(filename) == 0:
		errstr = "No data read from %s" % filename
		raise IOError(errstr)
	return data

def data_readin_unit_conversion(micron, flux_mJy, err_mJy):
	"""
	Convert data wavelength in micron, flux in mJy and flux error in mJy to SI units

	Input:
	---------
	int, floats

	Output:
	--------
	SI units
	"""

	wavelg_SI = micron * 1.e-6	# meters
	freq_Hz = c / wavelg_SI 	# Hz
	mJy2SIfac = 1000. * 10**(-26)
	flux_SI = flux_mJy / mJy2SIfac
	error_SI = err_mJy / mJy2SIfac
	return freq_Hz, flux_SI, error_SI

SMG_data = read_data_fromtxt(SMGfilename)
SMGwavelg_micron, SMG_mJy, SMG_error_mJy = ([dat[0] for dat in SMG_data],
		 									[dat[1] for dat in SMG_data],
		 									[dat[2] for dat in SMG_data])

freq_Hz, SMGflux_SI, SMGerr_SI = \
	data_readin_unit_conversion(np.asarray(SMGwavelg_micron),
								np.asarray(SMG_mJy),
								np.asarray(SMG_error_mJy))

Radio_data = read_data_fromtxt(RadioFilename)
Radio_Hz, Radio_Jy, Radio_error =\
	 ([dat[0] for dat in Radio_data],
	  [dat[2] for dat in Radio_data],
	  [dat[3] for dat in Radio_data])	# Jy
Radio_Hz, Radio_Jy, Radio_error = (np.asarray(Radio_Hz),
								   np.asarray(Radio_Jy),
								   np.asarray(Radio_error))	#Jy

# =======================================================
# 3C220.3 Core
# 4.86 GHz: no core detected, lobes (steep spectrum dominated), upper limit of core = 0.17 mJy
# 9 GHz: Core and lobes detected: core @ 0.8 mJy.. Based on Haas VLA
# Our point: 104.2106 GHz: peak = 5.56 mJy
Core_3c220dot3_Jy = np.array([0.8e-3, 0.17e-3])		# Jy
Core_3c220dot3_Hz = np.array([9.e9, 4.86e9])	# Hz
Point_Continuum_Jy = np.array([5.56e-3, 7.18e-3, 1.62e-3])		# Jy: peak, integrated, difference
Point_Cont_Hz = np.array([1.042e11, 1.042e11, 1.02e11])		# Hz
Point_error_Jy = np.array([5.0647e-4, 5.0647e-4, 5.0647e-4])	# Jy
# =======================================================
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

def B_nu(freq, freq_0, C):
	"""
	Simple MBB fit using Haas parameters
	"""
	z = 2.221
	T = 36 / (1+z)
	beta = 1.5
	boltz = k * T
	b_nu = 2.*h/c**2. * freq**3./(np.exp(h*freq/boltz)-1)
	tau = (freq/freq_0)**(beta)
	return (1-np.exp(-tau)) * b_nu * C 


def b_nuAndpower(freq, T, beta, alpha, freq_crit, dummyB, dummyP):
	"""
	MBB + powerlaw

	Input:
	---------
	freq: Hz

	Parameters: 
	------------
	T
	beta
	alpha
	freq_crit
	dummyB
	dummyP


	Returns:
	-----------
	a functional form with power law joint with MBB
	"""

	T = 18
	beta = 1.28
	alpha = -2.67
	freq_crit = 600.e9
	boltz = k * T
	b_nu = 2.*h/c**2 * freq**3. / (np.exp(h*freq/boltz)-1)
	tau = (freq/freq_crit) ** (beta)
	Planck = (1-np.exp(-tau)) * b_nu * dummyB
	power = (freq ** alpha) * dummyP
	f_com = Planck + power
	return f_com

def chi2CompositeFunc(theta, freq, flux, err):
	"""
	Input:
	---------
	theta = [T, beta, alpha, freq_crit, dummyB, dummyP]		# model parameters
	freq: Hz
	flux: SI

	"""

	model = b_nuAndpower(freq, *theta)		
	chisq = np.sum(((flux - model)/err)**2)
	print "chisq: %.2f " %chisq
	return chisq


def powerlaw(freq, dummyC, alpha):
	"""
	Puropse:
	--------
	get spectral index

	Input:
	-------
	freq: Hz
	dummyC: scaling factor
	alpha: spectral index

	Returns:
	--------
	A powerlaw function
	"""
	powerfit = dummyC * (freq**(alpha))
	return powerfit


def lobeSyn(freq, dummyA, beta, freq_t, freq_lobe):
	"""
	Ref: 
	-------
	equation 1 in CLEARY 2007


	Input: 
	--------
	freq: Hz
	dummyA is scaling factor
	beta = a parameter to the fit (beta_Cleary3C220dot3 = 0.19)
	freq_t = freq @ optical depth equals unity
	freq_lobe = cutoff freq for plasma energy (= 3e12 for 3C220.3)

	Returns:
	----------
	Flux density
	"""

	logFluxDensity = (-beta) * (np.log10(freq) - np.log10(freq_t))**2 + np.log10(np.exp(-freq/freq_lobe))		# Table 5 equation from CLEARY 2007
	return dummyA * 10.**logFluxDensity


def chi2func(theta, freq, flux, err):
	"""
	Input:
	---------
	theta = [dummyA, beta, freq_t, freq_lobe]		# model parameters
	freq: Hz
	flux

	Purpose:
	---------
	calculate the chi2 of radio galaxy synchrotron fit

	"""

	model = lobeSyn(freq, *theta)		
	chisq = np.sum(((flux - model)/err)**2)
	print "chisq: %.2f " %chisq
	return chisq


theta_guess = [65000., 0.19, 3.e6, 3.e12]
# result = opt.minimize(chi2func, theta_guess, args=(Radio_Hz[7:], Radio_Jy[7:]*1000., Radio_error[7:]*1000.))		# Hz and mJy of non-thermal data
# if result.success == True:
# 	print "scaling = %2.f, beta = %.2f, freq_t = %.4f GHz, freq_lobe = %.2f GHz"  %(result.x[0], result.x[1], result.x[2]/1.e9, result.x[3]/1.e9)
# else:
# 	print "**** Caution: failed to minimize chi2"

theta_best = opt.fmin(chi2func, theta_guess, args=(Radio_Hz[7:], Radio_Jy[7:]*1000., Radio_error[7:]*1000.))
print "scaling = %2.f, beta = %.2f, freq_t = %.4f GHz, freq_lobe = %.2f GHz"  %(theta_best[0], theta_best[1], theta_best[2]/1.e9, theta_best[3]/1.e9)
print 
print "Compare guess with CLEARY 2007 0.19, 3.e6, 3.e12 since we are using more points"
print 

# plot data and fit and error:
plt.errorbar(Radio_Hz/1.e9, Radio_Jy*1000, Radio_error*1000., fmt='.k', ecolor='lightgray', label='NED')		# mJy
plt.errorbar(Point_Cont_Hz/1.e9,  Point_Continuum_Jy*1000, Point_error_Jy*1000.,fmt='.r', ecolor='darkgrey', label='Our continuum')

yfit = lobeSyn(Radio_Hz[7:], *theta_best)
plt.plot(Radio_Hz[7:]/1.e9, yfit, label='Lobe fit to NED data')
extrapolate_flux = lobeSyn(104.21e9, *theta_best) 	# extrapolate the fit to 104 GHz
plt.plot(104.21e9/1.e9, extrapolate_flux, '*g', label='Extrapolate overall @104.21GHz')

core_fit_3c220dot3, core_fit_3c220dot3_covar = curve_fit(powerlaw, Core_3c220dot3_Hz, Core_3c220dot3_Jy, [7.728e-29, 2.5135477602042937])
# extrapolate to 104.21Ghz
extra104GHz3c220dot3_Jy = powerlaw(104.21e9, 1., -0.415)

core_fit_3c6dot1, core_fit_3c6dot1_covar = curve_fit(powerlaw, Core_3c6dot1_freq, Core_3c6dot1_Jy, [0.53746,-0.18853])
core_fit_3c263, core_fit_3c263_covar = curve_fit(powerlaw, Core_3c263_freq, Core_3c263_Jy, [283943., -0.6444])

print "Spectral index of 3C220.3: %.2f, 3C6.1: %.2f, 3C263: %.2f" %(core_fit_3c220dot3[1], core_fit_3c6dot1[1], core_fit_3c263[1])
print 

plt.plot(Core_3c220dot3_Hz[0]/1.e9, Core_3c220dot3_Jy[0]*1000, 'or', label='3c220.3 core data')
plt.errorbar(Core_3c220dot3_Hz[1]/1.e9, Core_3c220dot3_Jy[1]*1000, Core_3c220dot3_Jy[1]*1000, uplims=True, label='upper limit no core')
plt.plot(104.213e9/1.e9, extra104GHz3c220dot3_Jy*1000, 'k+', label='3c220.3 Extrapolated Core from 9GHz @ 104GHz') 	# using alpha = -0.451 (average of 3c6.1 and 3c263; assuming unity scaling factor)
plt.plot(Core_3c263_freq/1.e9, Core_3c263_Jy*1000, 'sc', label='3c263 core data')
plt.plot(Core_3c6dot1_freq/1.e9, Core_3c6dot1_Jy*1000, 'bx', label='3c6.1 core data points')


###########
# SMG 
###########
ComposTheta_guess = [18.2, 1.28, -2.65, 1200.e9, 4.6264e-12, 1.e-10]
ComposTheta_best = opt.fmin(chi2CompositeFunc, ComposTheta_guess, args=(freq_Hz[2:], SMGflux_SI[2:], SMGerr_SI[2:]))	# MBB and power
print ComposTheta_best
SMG_SED_fit = b_nuAndpower(freq_Hz[2:], *ComposTheta_best)
plt.plot(freq_Hz[2:]/1.e9, SMG_SED_fit*1.e26*1000., '-r', label='MBB and power fit')

fitParams, fitCovariances  =  curve_fit(B_nu, freq_Hz[-4:], SMGflux_SI[-4:], [1.0e12, 1])		# simple MBB fit using known params
y_SMG = B_nu(freq_Hz[-4:], *fitParams)
x_extend = np.linspace(70e9, 300e9, 1000)
y_SMG_extend = B_nu(x_extend, *fitParams)		# extrapolate to longer wavelengths
plt.errorbar(freq_Hz[2:]/1e9, SMG_mJy[2:], SMG_error_mJy[2:], fmt='ko')
plt.errorbar(freq_Hz[:2]/1e9, SMG_mJy[:2], SMG_error_mJy[:2], fmt='r*')
plt.plot(freq_Hz[-4:]/1e9, y_SMG*1e26*1000, 'k-', linewidth=2, label='simple MBB fit using Haas Params')

plt.plot(x_extend/1.e9, y_SMG_extend*1.e26*1000, 'k--', linewidth=2, label='extrapolate simple MBB')
plt.xlabel('freq [GHz]', fontsize=20,  fontweight='bold')
plt.title('SMM J0939+8315', fontsize=20,  fontweight='bold')
plt.ylabel('Flux Density [mJy]', fontsize=20,  fontweight='bold')
plt.yscale('log')
plt.xscale('log')
#plt.legend(loc='best')
plt.grid()
plt.show()


