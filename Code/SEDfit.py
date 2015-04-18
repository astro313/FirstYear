import scipy.optimize as opt
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import k, h, c
import matplotlib


matplotlib.rc('text', usetex=True)
matplotlib.rc('xtick', labelsize=15)
matplotlib.rc('ytick', labelsize=15)
font = {'family': 'sans-serif',
        'weight': 'bold',
        'size': 20}
matplotlib.rc('font', **font)
# monospace sans-serif


def savefigure(f, verbose=True, dirStr='../Figure/'):
    """
    Save figure in designatied directory and make directory if it doesn't alraedy exist

    Inputs:
    -------
    f: str
        filename without extension
    dirStr: str
        default to the ../Figure
    """
    import os
    if verbose is True:
        if not os.path.isdir(dirStr):
            print "... Making directory {:s}...".format(dirStr)
            os.mkdir(dirStr)
        else:
            print "... Directory {:s} exists...".format(dirStr)
    plt.savefig(dirStr + f + '.png')


def read_data_fromtxt(filename):
    """
    Read in rows .txt that are not commented out

    Input:
    ------
    filename: str
        filename with extension '.txt'

    Returns:
    --------
    data: astropy.table.table.Table
        content inside .txt
    """
    import astropy.io.ascii     # for .txt format
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


########################################
# user define area
########################################
SMGfilename = 'SMGdata.txt'
RadioFilename = 'Radiodata.txt'
Jy2mJy = 1000.
Hz2GHz = 1.e-9
SI2Jy = 1.e26
Avg_cont_freq = 104.21e9

# =======================================================
# 3C220.3 Core
# 4.86 GHz: no core detected, lobes (steep spectrum dominated), upper limit of core = 0.17 mJy
# 9 GHz: Core and lobes detected: core @ 0.8 mJy.. Based on Haas VLA
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


########################################
# Read in Data
########################################
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
Radio_Hz, Radio_Jy, Radio_error = ([dat[0] for dat in Radio_data],
                                   [dat[2] for dat in Radio_data],
                                   [dat[3] for dat in Radio_data])   # Jy

Radio_Hz, Radio_Jy, Radio_error = (np.asarray(Radio_Hz),
                                   np.asarray(Radio_Jy),
                                   np.asarray(Radio_error))  # Jy


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
#    powerfit = dummyC * (freq ** (alpha))
    powerfit = (freq ** (alpha))
    return powerfit


def lobeSyn(freq, dummyA, beta, freq_t, freq_lobe):
    """
    Ref:
    -------
    equation 1 in CLEARY 2007


    Input:
    --------
    freq: Hz
    dummyA is scaling factor that is calculating flux in [mJy]
    beta = a parameter to the fit (beta_Cleary3C220dot3 = 0.19)
    freq_t = freq @ optical depth equals unity
    freq_lobe = cutoff freq for plasma energy (= 3e12 for 3C220.3)
    empirical fit

    Returns:
    ----------
    Flux density: mJy
    """

    FluxDensity = dummyA * \
        ((10 ** (-beta)) ** (np.log10(freq / freq_t)) ** 2) * \
        np.exp(-freq / freq_lobe)

    return dummyA * FluxDensity


def chi2func(theta, freq, flux, err, verbose=False):
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
    if verbose == 'True':
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
    b_nu = 2. * h / c ** 2. * freq ** 3. / np.expm1(h * freq / boltz)
    tau = (freq / freq_0) ** (beta)
    return (1 - np.exp(-tau)) * b_nu * C


class MBBpower(object):

    def __init__(self, T, beta, alpha, lambda0=None, opthin=True):
        self.T = T
        self.beta = beta
        self.alpha = alpha
        self.c = 299792458.0
        self.h = 6.62606957e-34
        self.k = 1.3806488e-23
        self._opthin = opthin
        if self._opthin is not True:
            try:
                self.freq_crit = self.c / lambda0      # in Hz, lambda0 in metres
            except Exception, e:
                raise e + ('must supply lambda0 for optically thick')

    def b_nuAndpower(self, freq, dummyB, dummyP):
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

        boltz = self.k * self.T
        b_nu = 2. * self.h / self.c ** 2 * freq ** 3. / np.expm1(self.h * freq / boltz)
        if self._opthin is True:
            Planck = freq ** self.beta * b_nu * dummyB
        else:
            tau = (freq / self.freq_crit) ** (self.beta)
            Planck = (1 - np.exp(-tau)) * b_nu * dummyB
        power = (freq ** self.alpha) * dummyP
        f_com = (Planck + power)
        return f_com


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

    sigma = []
    for i in range(len(Covar)):
        diag = Covar[i][i]
        checkInvalid = (
            np.isinf(diag).any() == True or np.isnan(diag).any() == True)
#        print i, checkInvalid
        if checkInvalid != True:        # non Invalid values
            sigma.append(np.abs(diag) ** 0.5)
        else:
            print "Variance not valid: {.2f}".format(np.abs(diag))
    return sigma


def errorfill(x, y, yerr, color=None, alpha_fill=0.3, ax=None):
    """
    A function to plot fitted line with the propagated errors associated along the x array..
    """
    ax = ax if ax is not None else plt.gca()
    if color is None:
        color = ax._get_lines.color_cycle.next()
    # different errors for each point
    if np.isscalar(yerr) or len(yerr) == len(y):
        ymin = y - yerr
        ymax = y + yerr
    # same errors for all
    elif len(yerr) == 2:
        ymin, ymax = yerr
    ax.plot(x, y, color=color)
    ax.fill_between(x, ymax, ymin, color=color, alpha=alpha_fill)


def GHz2um(GHz):
    """
    Convert GHz x-array to um for plotting

    Inputs:
    GHz: array
        frequency array in GHz
    Returns:
    um: array
        wavelength array in um
    """
    m2um = 1.e6
    um = c / (GHz / Hz2GHz) * m2um
    return um


############################################
# Initialize Twin x axes Plot
############################################
from matplotlib import transforms as mtransforms
import astropy.units as un
import numpy as np
from mpl_toolkits.axes_grid.parasite_axes import SubplotHost


class Freq2WavelengthTransform(mtransforms.Transform):
    input_dims = 1
    output_dims = 1
    is_separable = False
    has_inverse = True

    def __init__(self):
        mtransforms.Transform.__init__(self)

    def transform_non_affine(self, fr):
        return (fr*un.GHz).to(un.mm, equivalencies=un.spectral()).value

    def inverted(self):
        return Wavelength2FreqTransform()


class Wavelength2FreqTransform(Freq2WavelengthTransform):
    input_dims = 1
    output_dims = 1
    is_separable = False
    has_inverse = True

    def __init__(self):
        mtransforms.Transform.__init__(self)

    def transform_non_affine(self, wl):
        """
        z: float
            redshift
        """
        x = (wl*un.um).to(un.GHz, equivalencies=un.spectral()).value
        z = 2.221
        x *= (1. + z)
        return x

    def inverted(self):
        return Freq2WavelengthTransform()

fig = plt.figure(1, figsize=(10, 8))
ax = SubplotHost(fig, 1, 1, 1)
fig.add_subplot(ax)
print "*** !! Redshifting axis for z = 2.221, if undesired, change in function"

aux_trans = mtransforms.BlendedGenericTransform(Wavelength2FreqTransform(), mtransforms.IdentityTransform())
ax_top = ax.twin(aux_trans)
ax_top.set_xlabel(r'$\nu_{Rest} [GHz]$', size=15)
ax_top.set_viewlim_mode("transform")
ax_top.axis["right"].toggle(ticklabels=False)


############################################
# Plot foreground galaxy data and fit
############################################
theta_guess = [1e3, 0.19, 3.e6, 3.e12]      # dummy, (CLEARY 2007: 0.19, 3.e6, 3.e12)
theta_best = opt.fmin(chi2func, theta_guess, args=(Radio_Hz[7:], Radio_Jy[7:] * Jy2mJy, Radio_error[7:] * Jy2mJy))

if False:
    lobe_fitrestr = "scaling = {:.2f}, beta = {:.2f}, freq_t = {:.4f} Hz, freq_lobe = {:.2f} Hz"
    print lobe_fitrestr.format(*theta_best)

ax.errorbar(GHz2um(Radio_Hz[7:] * Hz2GHz),
            Radio_Jy[7:] * Jy2mJy,
            Radio_error[7:] * Jy2mJy,
            fmt='.k', ecolor='darkgray', label='3C220.3 NED Data')      # mJy

ax.errorbar(GHz2um(Point_Cont_Hz * Hz2GHz),
            Point_Continuum_Jy * Jy2mJy, Point_error_Jy * Jy2mJy,
            fmt='.r', ecolor='darkgrey', label='Our continuum', ms=10)
ax.plot(GHz2um(Core_3c263_freq * Hz2GHz), Core_3c263_Jy * Jy2mJy,
        'sc',
        label='3C263 Core')
ax.plot(GHz2um(Core_3c6dot1_freq * Hz2GHz), Core_3c6dot1_Jy * Jy2mJy,
        'bx',
        label='3C6.1 Core')
ax.plot(GHz2um(Core_3c220dot3_Hz[0] * Hz2GHz), Core_3c220dot3_Jy[0] * Jy2mJy,
        'or',
        label='3C220.3 Core')

ax.errorbar(GHz2um(Core_3c220dot3_Hz[1] * Hz2GHz),
            Core_3c220dot3_Jy[1] * Jy2mJy,
            Core_3c220dot3_Jy[1] * Jy2mJy,
            uplims=True,
            fmt='rv')   # label='upper limit on core')

yfit = lobeSyn(Radio_Hz[7:], *theta_best)
ax.plot(GHz2um(Radio_Hz[7:] * Hz2GHz), yfit, 'k-.', linewidth=2, label='Lobe fit')  # 'Lobe fit to NED data')

epsilon = 50.e9     # extra buffer
start_Hz = 10.e9
end_Hz = Avg_cont_freq+epsilon
x_syn_extend = np.linspace(start_Hz, end_Hz)
syn_extend_mJy = lobeSyn(x_syn_extend, *theta_best)
ax.plot(GHz2um(x_syn_extend * Hz2GHz), syn_extend_mJy, 'k-.', linewidth=2)   # label='Extrapolate Lobe fit')

# ------ Core extrapolate -----------
# ------ Don't think worth doing -------
# core_fit_3c220dot3, core_fit_3c220dot3_covar = curve_fit(
#     powerlaw, Core_3c220dot3_Hz, Core_3c220dot3_Jy,
#     [1, 2.57])      # calculated by hand
# # extrapolate to 104.21Ghz using average spectral index
# extra104GHz3c220dot3_Jy = powerlaw(Avg_cont_freq, 1, -0.15212785)
# core_fit_3c6dot1, core_fit_3c6dot1_covar = curve_fit(
#     powerlaw, Core_3c6dot1_freq, Core_3c6dot1_Jy, [1, -0.18853])
# core_fit_3c263, core_fit_3c263_covar = curve_fit(powerlaw,
#                                                  Core_3c263_freq,
#                                                  Core_3c263_Jy,
#                                                  [1., -0.6444])
# print"Spectral index of 3C220.3: {:.2f}, 3C6.1: {:.2f}, 3C263: {:.2f}"\
#     .format(core_fit_3c220dot3[1], core_fit_3c6dot1[1], core_fit_3c263[1])
# plt.plot(Avg_cont_freq * Hz2GHz,
#          extra104GHz3c220dot3_Jy * Jy2mJy,
#          'k+',
#          markersize=10,
#          label='3c220.3 Extrapolated Core from 9GHz @ 104GHz')
#          'k+')

############################################
# SMG Plot MBB results
############################################
import mbb_emcee
filename_thick = 'thick_500_500.h5'
res_thick = mbb_emcee.mbb_results(h5file=filename_thick)
wave, flux, flux_unc = res_thick.data
p_data = ax.errorbar(wave, flux, yerr=flux_unc, fmt='go', label="J0939 IR Data")
p_wave = np.linspace(wave.min() * 0.5, wave.max() * 1.5, 200)
p_fit_thick = ax.plot(p_wave, res_thick.best_fit_sed(p_wave), 'b:', lw=2)
filename_thin = 'thin_500_500.h5'
res_thin = mbb_emcee.mbb_results(h5file=filename_thin)
p_fit_thin = ax.plot(p_wave, res_thin.best_fit_sed(p_wave), '--', color='magenta', lw=2)

ax.set_xlabel(r'$\lambda_{obs} [\mu m]$', size=20)
# ax.set_title('SMM J0939+8315', fontsize=20,  fontweight='bold')
ax.set_ylabel('Flux Density [mJy]', fontsize=20,  fontweight='bold')
plt.legend(loc='best', fontsize='10')
ax.set_xscale('log')
ax.set_yscale('log')
ax_top.set_xscale('log')
ax_top.set_yscale('log')

############################################
# Save Figure
############################################
try:
    response = raw_input('Save fig?: (y/n)')
except EOFError:
    raise ValueError(response)

if response.lower() in ['y', 'yes']:
    filename = '3C220.3FullSED'
    savefigure(filename)
elif response.lower() in ['n', 'no']:
    pass
else:
    print "Response {0} is not valid! ".format(response)
plt.show()
