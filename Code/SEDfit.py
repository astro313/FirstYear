from scipy.optimize import curve_fit
import scipy.optimize as opt
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import k, h, c
import matplotlib
import astropy.io.ascii     # for .txt format

matplotlib.rc('xtick', labelsize=18)
matplotlib.rc('ytick', labelsize=18)
font = {'family': 'sans-serif',
        'weight': 'bold',
        'size': 12}
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

    wavelg_SI = micron * 1.e-6  # meters
    freq_Hz = c / wavelg_SI     # Hz
    mJy2SIfac = 1000. * 10 ** (26)
    flux_SI = flux_mJy / mJy2SIfac
    error_SI = err_mJy / mJy2SIfac
    return freq_Hz, flux_SI, error_SI

SMG_data = read_data_fromtxt(SMGfilename)
SMGwavelg_micron, SMG_mJy, SMG_error_mJy = ([dat[0] for dat in SMG_data],
                                            [dat[1] for dat in SMG_data],
                                            [dat[2] for dat in SMG_data])

# np.asarray converts input to array
freq_Hz, SMGflux_SI, SMGerr_SI = \
    data_readin_unit_conversion(np.asarray(SMGwavelg_micron),
                                np.asarray(SMG_mJy),
                                np.asarray(SMG_error_mJy))

Radio_data = read_data_fromtxt(RadioFilename)
Radio_Hz, Radio_Jy, Radio_error = \
    ([dat[0] for dat in Radio_data],
     [dat[2] for dat in Radio_data],
     [dat[3] for dat in Radio_data])   # Jy
Radio_Hz, Radio_Jy, Radio_error = (np.asarray(Radio_Hz),
                                   np.asarray(Radio_Jy),
                                   np.asarray(Radio_error))  # Jy

# =======================================================
# 3C220.3 Core
# 4.86 GHz: no core detected, lobes (steep spectrum dominated), upper limit of core = 0.17 mJy
# 9 GHz: Core and lobes detected: core @ 0.8 mJy.. Based on Haas VLA
# Our point: 104.2106 GHz: peak = 5.56 mJy
Core_3c220dot3_Jy = np.array([0.8e-3, 0.17e-3])     # Jy
Core_3c220dot3_Hz = np.array([9.e9, 4.86e9])    # Hz

# Our Data: [peak, integrated, difference]
Point_Continuum_Jy = np.array([5.56e-3, 7.18e-3, 1.62e-3])
Point_Cont_Hz = np.array([1.042e11, 1.042e11, 1.02e11])     # Hz
Point_error_Jy = np.array([5.0647e-4, 5.0647e-4, 5.0647e-4])    # Jy
# =======================================================
# compare to similiar FR II galaxies
# 3C6.1:
# 7.3 mJy @ 8GHz and 10.14mJy at 1.4 GHz
Core_3c6dot1_Jy = np.array([7.3e-3, 10.14e-3])      # Jy
Core_3c6dot1_freq = np.array([8e9, 1.4e9])  # Hz
# 3C 263:
# 161.1 mJy @ 5 GHz and 0.119 Jy @ 8GHz
Core_3c263_Jy = np.array([161.1e-3, 0.119])
Core_3c263_freq = np.array([5e9, 8e9])      # Hz
# =======================================================


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
    powerfit = dummyC * (freq ** (alpha))
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
    Flux density: mJy
    """

    logFluxDensity = (-beta) * (np.log10(freq) - np.log10(freq_t)) ** 2 + \
        np.log10(
            np.exp(-freq / freq_lobe))       # Table 5 equation from CLEARY 2007
    return dummyA * 10. ** logFluxDensity


def chi2func(theta, freq, flux, err):
    """
    Input:
    ---------
    theta = [dummyA, beta, freq_t, freq_lobe]       # model parameters
    freq: Hz
    flux: mJy

    Purpose:
    ---------
    calculate the chi2 of radio galaxy lobe synchrotron fit

    """

    model = lobeSyn(freq, *theta)
    chisq = np.sum(((flux - model) / err) ** 2)
    print "chisq: {:.2f} \n".format(chisq)
    return chisq


def B_nu(freq, freq_0, C, beta=1.5):
    """
    Simple MBB fit using Haas parameters
    Fixed z to z_spec


    """
    z = 2.221
    T = 36 / (1 + z)
    boltz = k * T
    b_nu = 2. * h / c ** 2. * freq ** 3. / (np.exp(h * freq / boltz) - 1)
    tau = (freq / freq_0) ** (beta)
    return (1 - np.exp(-tau)) * b_nu * C


def b_nuAndpower(freq, dummyB):
    """
    MBB + powerlaw

    Input:
    ---------
    freq: Hz

    Parameters:
    ------------
    dummyB
    dummyP


    Returns:
    -----------
    a functional form with power law joint with MBB
    flux in SI
    """

    T = 18
    beta = 1.28
    alpha = -2.67
    freq_crit = 600.e9
    boltz = k * T
    b_nu = 2. * h / c ** 2 * freq ** 3. / (np.exp(h * freq / boltz) - 1)
    tau = (freq / freq_crit) ** (beta)
    Planck = (1 - np.exp(-tau)) * b_nu  # * dummyB
    power = (freq ** alpha)  # * dummyP
    f_com = (Planck + power) * dummyB
    return f_com


def chi2CompositeFunc(theta, freq, flux, err):
    """
    Input:
    ---------
    theta = [dummyB]        # model parameters
    freq: Hz
    flux: SI

    """

    model = b_nuAndpower(freq, *theta)
    chisq = np.sum(((flux - model) / err) ** 2)
    print"chisq: {:.2f}".format(chisq)
    return chisq


def Covar2Sig(Covar):
    """
    Purpose:
    ----------
    Take in a covariance matrix (from curve fit in this script) and return sigma, an array containing the diagonal elements in the order of the parameters

    Input:
    -------
    Covariance matrix

    Returns:
    ---------
    sigma: a list, error associated with each parameters
    None if covariance matrix contains inf or nan
    """

    for i in range(len(Covar)):
        diag = Covar[i][i]
        checkInvalid = (
            np.isinf(diag).any() == True or np.isnan(diag).any() == True)
        print i, checkInvalid
        if checkInvalid != True:
            try:
                sigma
            except NameError:
                sigma = []
            sigma.append(np.abs(diag) ** 0.5)
            return sigma
        else:
            return None


def plot_fit_err(x, y, function, **kwargs):
    """
    Purpose:
    ----------
    A function to plot fitted line with the propagated errors associated along the x array..

    Status:
    --------
    Not ready

    Todo:
    -------
    implement fill error with alpha
    """
    # plt.plot(t, fitFunc(t, fitParams[0], fitParams[1], fitParams[2]),\
    #       t, fitFunc(t, fitParams[0] + sigma[0], fitParams[1] - sigma[1], fitParams[2] + sigma[2]),\
    #       t, fitFunc(t, fitParams[0] - sigma[0], fitParams[1] + sigma[1], fitParams[2] - sigma[2])\
    #       )
    return None


Jy2mJy = 1000.
Hz2GHz = 1.e-9
SI2Jy = 1.e26


# dummy, (CLEARY 2007 0.19, 3.e6, 3.e12)
theta_guess = [1., 0.19, 3.e6, 3.e12]
theta_best = opt.fmin(chi2func, theta_guess,
                      args=(Radio_Hz[7:], Radio_Jy[7:] * Jy2mJy, Radio_error[7:] * Jy2mJy))
lobe_fitrestr = "scaling = {:.2f}, beta = {:.2f}, freq_t = {:.4f} Hz, freq_lobe = {:.2f} Hz"
print lobe_fitrestr.format(*theta_best)
import sys
sys.exit()

# plot data and fit and error:
plt.errorbar(Radio_Hz * Hz2GHz,
             Radio_Jy * Jy2mJy,
             Radio_error * Jy2mJy,
             fmt='.k', ecolor='lightgray', label='NED')      # mJy
plt.errorbar(Point_Cont_Hz * Hz2GHz,
             Point_Continuum_Jy * Jy2mJy, Point_error_Jy * Jy2mJy,
             fmt='.r', ecolor='darkgrey', label='Our continuum', ms=10)


# Testing:


# _lobe_fit_params, lobe_cov = \
#     curve_fit(lobeSyn, Radio_Hz[7:], (Radio_Jy[
#               7:] * Jy2mJy), theta_best)   # to get covariance matrix from the best chi2 params
# print Covar2Sig(lobe_cov)
yfit = lobeSyn(Radio_Hz[7:], *theta_best)       # just fitting a curve
plt.plot(Radio_Hz[7:] * Hz2GHz, yfit,
         label='Lobe fit to NED data')


x_syn_extend = np.linspace(10.e9, 104.21e9)
# extrapolate the fit to 104 GHz
extrapolate_mJy = lobeSyn(x_syn_extend, *theta_best)
plt.plot(x_syn_extend * Hz2GHz, extrapolate_mJy, ':',
         label='Extrapolate to 104.21GHz')

core_fit_3c220dot3, core_fit_3c220dot3_covar = \
    curve_fit(powerlaw, Core_3c220dot3_Hz, Core_3c220dot3_Jy,
              [7.728e-29, 2.5135477602042937])
import pdb
pdb.set_trace()

# extrapolate to 104.21Ghz using average spectral index
extra104GHz3c220dot3_Jy = powerlaw(104.21e9, 1., -0.415)

core_fit_3c6dot1, core_fit_3c6dot1_covar = curve_fit(
    powerlaw, Core_3c6dot1_freq, Core_3c6dot1_Jy, [0.53746, -0.18853])

core_fit_3c263, core_fit_3c263_covar = curve_fit(powerlaw,
                                                 Core_3c263_freq,
                                                 Core_3c263_Jy,
                                                 [283943., -0.6444])

print"Spectral index of 3C220.3: {:.2f}, 3C6.1: {:.2f}, 3C263: {:.2f}"\
    .format(core_fit_3c220dot3[1], core_fit_3c6dot1[1], core_fit_3c263[1])
print

plt.plot(Core_3c220dot3_Hz[0] * Hz2GHz,
         Core_3c220dot3_Jy[0] * Jy2mJy,
         'or',
         label='3c220.3 core data')

plt.errorbar(Core_3c220dot3_Hz[1] * Hz2GHz,
             Core_3c220dot3_Jy[1] * Jy2mJy,
             Core_3c220dot3_Jy[1] * Jy2mJy,
             uplims=True,
             label='upper limit no core')

plt.plot(104.213e9 * Hz2GHz,
         extra104GHz3c220dot3_Jy * Jy2mJy,
         'k+',
         markersize=10,
         label='3c220.3 Extrapolated Core from 9GHz @ 104GHz')   # using alpha = -0.451 (average of 3c6.1 and 3c263; assuming unity scaling factor)
plt.plot(Core_3c263_freq * Hz2GHz,
         Core_3c263_Jy * Jy2mJy,
         'sc',
         label='3c263 core data')
plt.plot(Core_3c6dot1_freq * Hz2GHz,
         Core_3c6dot1_Jy * Jy2mJy,
         'bx',
         label='3c6.1 core data points')


###########
# SMG
###########
ComposTheta_guess = [4.6264e-7]  # , 1.e-9]
ComposTheta_best = opt.fmin(chi2CompositeFunc,
                            ComposTheta_guess,
                            args=(freq_Hz[:], SMGflux_SI[:], SMGerr_SI[:]))  # MBB and power
print ComposTheta_best

SMG_SED_fit = b_nuAndpower(freq_Hz[:], *ComposTheta_best)
plt.plot(freq_Hz[:] * Hz2GHz, SMG_SED_fit * SI2Jy * Jy2mJy, '-c',
         label='MBB and power fit')

fitParams, fitCovariances = curve_fit(B_nu, freq_Hz[-4:],
                                      SMGflux_SI[-4:],
                                      [1.0e12, 1])       # simple MBB fit using known params
y_SMG = B_nu(freq_Hz[-4:], *fitParams)
x_extend = np.linspace(70.e9, 300.e9)
# extrapolate to longer wavelengths
y_SMG_extend = B_nu(x_extend, *fitParams)
plt.errorbar(freq_Hz * Hz2GHz, SMG_mJy, SMG_error_mJy, fmt='ko')
plt.plot(freq_Hz[-4:] * Hz2GHz,
         y_SMG * SI2Jy * Jy2mJy,
         'k-',
         linewidth=2,
         label='simple MBB fit using Haas Params')

plt.plot(x_extend * Hz2GHz, y_SMG_extend * SI2Jy * Jy2mJy, 'b--', alpha=0.5,
         linewidth=2, label='extrapolate simple MBB')


plt.xlabel('Freq [GHz]', fontsize=20,  fontweight='bold')
plt.title('SMM J0939+8315', fontsize=20,  fontweight='bold')
plt.ylabel('Flux Density [mJy]', fontsize=20,  fontweight='bold')
plt.loglog()
plt.legend(loc='best')
plt.grid()
plt.show()

# Question is why chi2 not similar to mbb-emcee?
# Still giving me chi2  ~ 60...
# won't integrate to get FIR using this for now, since normalization seems weird
