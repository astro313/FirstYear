import mbb_emcee
import matplotlib.pyplot as plt
import sys
from pylab import savefig
import numpy as np

"""
Purpose:
--------
Look at resulting fits from mbb_emcee
plot probabilities of parameters
"""


def convergence(plotTitle='ln prob', saveFig=False):
    """
    Plot the convergence profile.  I.e., Max(lnprob) - lnprob as a function of
    iteration number.
    """

    lnprob = res.lnprobability.flatten()
    Plotlnprob = max(lnprob) - lnprob
    plt.clf()
    plt.plot(Plotlnprob)  # , ',', alpha=0.5)
    plt.xlabel('iteration')
    plt.ylabel('max(lnprob) - lnprob')
    plt.title(plotTitle)
    plt.semilogy()
    if saveFig == True:
        outfile = 'convergence__'
        savefig('../Figure/' + outfile + filename.replace('.h5', '.png'))
    else:
        plt.show()


def Ramdon():
    """
    Grab a set of values for the parameters randomly from the sampled result

    """
    return res.choice()


def get_Paramvalues(percentile=68.3):
    """
    Input:
    ------
    percentile: +/- return of all parameters


    Print:
    --------
    5 element list of parameters (T, beta, lambda0, alpha, fnorm)
    list of 3, [mean = central confidence, +errorUp, -errorLow]
    default percentile =  68.3
    The returned is not the max. ln prob value, only the mean
    For best value: see BestFit()

    Returns:
    --------
    None at the moment
    """

    paramsDict = {'T': 'K',
                  'beta': ' ',
                  'lambda0': 'um',
                  'alpha': ' ',
                  'fnorm': 'mJy'}

    for k, v in paramsDict.iteritems():
        param = k
        param_val = res.par_cen(param, percentile=percentile)
        unit = v
        print(
            "{0:s}: Mean {2:0.2f}+{3:0.2f}-{4:0.2f} [{1:s}]".format(param, unit, *param_val))
    # Do I want to save this?, also accessible from .log output from
    # mbb_emcee, probably useful if I want a different percentile, this is not the bestFit value anyway...


def get_computation_value(FIR=False, percentile=68.3):
    """
    Purpose:
    --------
    Print peak lambda (p_val [mean, +errorUp, -errorLow])
    Returns, LIR, (LFIR), dustMass

    Input:
    ------
    Keywords:
    FIR for computing FIR luminosity; default: wavemin=42.5, wavemax=122.5
    percentile for computing the value at a different confidence interval


    Returns:
    ---------
    LIR
    LFIR    # if FIR == True
    Dust Mass


    Note:
    ------
    Most quantities are computed before saving results to .h5

    However, if one wish to get an integrate Luminosity of different range of wavelengths, that will require computing in the function, not available from result.
    """

    p_val = res.peaklambda_cen(percentile=percentile)
    print("Peak Obs wavelength: {:0.1f}+{:0.1f}-{:0.1f} [um]".format(*p_val))
    lir = res.lir_cen(percentile=percentile)     # fetch from computed, integrated
    args = (res._lir_min, res._lir_max) + tuple(lir)
    lirstr = "L_IR({:0.1f} to {:0.1f}um): {:0.2f} "\
        "+{:0.2f} -{:0.2f} [10^12 L_sun]\n"
    print(lirstr).format(*args)
    IR_ = {'lir': lir, 'lirChain': res.lir_chain}    # Avoid over-written by FIR

    if FIR == True:
        res.compute_lir(wavemin=42.5, wavemax=122.5)    # compute from chain
        RangeMSG = "Current range of wavelength for computing LIR Luminosity: {0:.2f} - {1:.2f} um"
        print RangeMSG.format(*res.lir_wavelength)
        FIR_val = res.lir_cen(percentile=percentile)
        # _lir_min changed to 42.5
        args = (res._lir_min, res._lir_max) + tuple(FIR_val)
        lirstr = "L_FIR({:0.1f} to {:0.1f}um): {:0.2f} "\
            "+{:0.2f} -{:0.2f} [10^12 L_sun]\n"
        print(lirstr).format(*args)
        IR_['FIRchain'], IR_['LFIR'] = res.lir_chain, FIR_val
    return IR_


def Plot_Standard(saveFig=False):
    """
    Purpose:
    --------
    Make 4 panels plot of
    - Data + errorbar + best-fit SED (unnormalized)
    - Likelihood of T [Temperature rest Frame]
    - Likelihood of beta
    - Likelihood of L_IR
    """
    wave, flux, flux_unc = res.data
    redshiftZ = res.redshift
    f, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))

    ax1 = axs[0, 0]
    p_data = ax1.errorbar(wave, flux, yerr=flux_unc, fmt='ro')
    p_wave = np.linspace(wave.min() * 0.5, wave.max() * 1.5, 200)   # range of x
    # plot  best fit using the params in res
    p_fit = ax1.plot(p_wave, res.best_fit_sed(p_wave), color='blue')
    ax1.set_xlabel('Wavelength [um]')
    ax1.set_ylabel('Flux Density [mJy]')
    ax1.set_title('Best-fitted SED Unnormalized')

    ax2 = axs[0, 1]
    h1 = ax2.hist(res.parameter_chain('T/(1+z)') * (1 + redshiftZ))
    ax2.set_xlabel(r'$T_{d, rest}$')
    ax2.set_ylabel(r'$\cal L$')

    ax3 = axs[1, 0]
    h2 = ax3.hist(res.parameter_chain('beta'))
    ax3.set_xlabel(r'$\beta$')
    ax3.set_ylabel(r'$\cal L$')

    ax4 = axs[1, 1]
    h4 = ax4.hist(res.dustmass_chain)
    ax4.set_xlabel('Dust Mass  ' + r'$[10^8 M_{\odot}]$')
    ax4.set_ylabel(r'$\cal L$')

    if saveFig == True:
        outfile = 'SED_ParamProb__'
        savefig('../Figure/' + outfile + filename.replace('.h5', '.png'))
    else:
        f.show()


def Plot_Chain(IR_chain, saveFig=False, FIR=False):

    f, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))

    ax1 = axs[0, 0]
    h1 = ax1.hist(res.peaklambda_chain)
    ax1.set_xlabel('Peak Obs Wavelength [um]')
    ax1.set_ylabel(r'$\cal L$')

    ax2 = axs[0, 1]
    h2 = ax2.hist(IR_chain['lirChain'])
    ax2.set_xlabel(r'$L_{IR} [10^{12}L_{\odot}]$')
    ax2.set_ylabel(r'$\cal L$')

    # ax3 = axs[1, 0]
    # h3 = ax3.hist(IR_chain['lirChain'])
    # ax3.set_xlabel(r'$L_{IR}\ \ [10^{12}L_{\odot}]$')
    # ax3.set_ylabel(r'$\cal L$')

    if FIR == True:
        ax4 = axs[1, 1]
        h3 = ax4.hist(IR_chain['FIRchain'])
        ax4.set_xlabel(r'$L_{FIR} [10^{12}L_{\odot}]$')
        ax4.set_ylabel(r'$\cal L$')

    if saveFig == True:
        outfile = 'ParamCHAINProb__'
        savefig('../Figure/' + outfile + filename.replace('.h5', '.png'))
    else:
        f.show()


def MCMC_scatterDots(Xparam='T', Yparam='beta', saveFig=False):
    """
    Purpose:
    --------
    Plot scattered chain samples of 2 parameters from MCMC results
    Similar to hist2D

    Input:
    ------
    Xparam: Observed frame dust temeprature
    Yparam: beta
    options: 'fnorm', 'beta', 'alpha', 'lambda_0', 'T'

    """
    X = res.parameter_chain('T')
    Y = res.parameter_chain('beta')

    # Define plot ranges
    Xrange = [min(X) - 0.2, max(X) + 0.2]
    Yrange = [min(Y) - 0.2, max(Y) + 0.2]

    # Define figure size and formatting
    fig = plt.figure(figsize=(7, 7))
    fig.subplots_adjust(left=0.10, bottom=0.09, top=0.98, right=0.98)
    plt.plot(X, Y, 'o', ms=4, alpha=.1, color='b')

    # set plot range, axes ticks, and axes labels
    plt.xlim(Xrange)
    plt.ylim(Yrange)
    plt.xlabel(Xparam)
    plt.ylabel(Yparam)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    # plt.gcf()
    if saveFig == True:
        outfile = 'ScatterPlot' + Xparam + Yparam + '__'
        savefig('../Figure/' + outfile + filename.replace('.h5', '.png'))
    else:
        plt.show()


def PanelPlot(Xparam='beta', Yparam='T/(1+z)', maxLev=False, saveFig=False):
    """
    Take MCMC samples parameters
    and make plots with sides panels = projected distribution of params

    Main panel as Contours of confidence levels and Make 2D histogram using colormesh and hist2d

    For single points , see MCMC_scatterDots()

    Input:
    ------
        keyword:
        - maxLev = using percentage of the maximum value instead of confidence interval

    """
    from matplotlib.gridspec import GridSpec
    X = res.parameter_chain(Xparam)
    Y = res.parameter_chain(Yparam)
    if Yparam.find('z') != -1:
        z = res.redshift
        Y *= (1 + z)
    Xrange = [min(X), max(X)]
    Yrange = [min(Y), max(Y)]

    fig = plt.figure(figsize=(7, 7))           # window 2
    fig.subplots_adjust(
        hspace=0.001, wspace=0.001, left=0.10, bottom=0.095, right=0.98)
    # gridspec allows you to assign different formats to panels in one plot
    # getting plot grid  2 X 2
    gs = GridSpec(2, 2, width_ratios=[1, 4], height_ratios=[4, 1])

    # Main panel top tight contains full 2d histogram
    plt.subplot(gs[1])
    Bins = 24
    hist2D, xedges, yedges = np.histogram2d(
        X, Y, bins=[Bins, Bins], range=[Xrange, Yrange], normed=False)
    # Be Aware: numpy switches axes, so switch back
    hist2D = np.transpose(hist2D)

    # smoothed distribution instead of signal points
    plt.pcolormesh(xedges, yedges, hist2D, cmap=plt.cm.rainbow)

    # Overplot with error contours
    if maxLev == True:
        maximum = np.max(hist2D)
        # maximum
        [L1, L2, L3] = [0.5 * maximum, 0.25 * maximum, 0.125 * maximum]
        fmtdict = {L1: '0.5Max', L2: '0.25Max', L3: '0.125Max'}

    else:   # using sigma as levels
        conf = np.percentile(hist2D, 68.3)
        conf2 = np.percentile(hist2D, 95.4)
        conf3 = np.percentile(hist2D, 99.7)
        [L1, L2, L3] = [conf, conf2, conf3]
        fmtdict = {L1: r'$\sigma$', L2: r'$2\ \sigma$', L3: r'$3\ \sigma$'}

    # Use bin edges to restore extent
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    cs = plt.contour(hist2D, extent=extent, levels=[L1, L2, L3], linestyles=[
        '--', '--', '--'], colors=['orange', 'orange', 'orange'], linewidths=3)

    # colour label for contour levels
    plt.clabel(cs, fmt=fmtdict, inline=True, fontsize=20)

    # Plot the sameple points on top, gets really messy...
    # plt.plot(X, Y ,'bo', markersize=5)

    plt.xlim(Xrange)
    plt.ylim(Yrange)

    # Add side panels showing the projected distributions of X and Y:
    # Bin X, Y separately. In 1D bin, can use more bins now...
    S = 101
    LX = np.histogram(X, bins=S, range=Xrange, normed=True)[0]
    LY = np.histogram(Y, bins=S, range=Yrange, normed=True)[0]

    # restore positions lost by binning, similiar to making arrays in IDL
    X = Xrange[0] + (Xrange[1] - Xrange[0]) * \
        np.array(range(0, len(LX))) / float(len(LX) - 1)
    Y = Yrange[0] + (Yrange[1] - Yrange[0]) * \
        np.array(range(0, len(LY))) / float(len(LY) - 1)

    # bottom right panel: projected density of x.
    plt.subplot(gs[3])
    plt.plot(X, LX, '-', lw=3, color='black')

    # Adjust tick mark size and take out yticks
    plt.xticks(fontsize=16)
    # take out tick marks, replace with simply curly L symbol for likelihood
    plt.yticks([])
    plt.xlabel(Xparam, fontsize=24)
    plt.ylabel(r'$\cal L$', fontsize=24)          # likelihood
    plt.xlim(Xrange)
    plt.ylim(0.0, 1.1 * np.max(LX))

    # Top left panel: Projected density of Y
    plt.subplot(gs[0])
    # Likelihood on horizontal axis
    plt.plot(LY, Y, '-', lw=3, color='k')

    plt.yticks(fontsize=16)
    plt.xticks([])
    plt.xlabel(r'$\cal L$', fontsize=24)
    plt.ylabel(Yparam, fontsize=24)
    plt.xlim(0.0, 1.1 * np.max(LY))
    plt.ylim(Yrange)
    if saveFig == True:
        outfile = 'mcmcSamplesPanels' + Xparam + 'T' + '__'
        savefig('../Figure/' + outfile + filename.replace('.h5', '.png'))
    else:
        plt.show()

def BestFit(verbose=True):
    """
    Grab best fit parameters, e.g. chi_square, best_set of param values
    """
    BestFitIndex = res._best_fit[2]
    best_fit_chisq = -2.0 * res._best_fit[1]
    if verbose == True:
        print "BestFitParams {:s}".format(res._best_fit[0])
        print "BestFitLogLike %.2f" %(res._best_fit[1])
        print "Best Chi Square: %.2f" %best_fit_chisq
    return best_fit_chisq


if __name__ == '__main__':
    """
    run script.py filename True
    sys.argv[2] determines whether or not to get values and plots for FIR
    """
    if len(sys.argv) < 3:
        errmsg = "Invalid number of arguments: {0:d}\n  run script.py filename True"
        raise IndexError(errmsg.format(len(sys.argv)))

    try:
        filename = sys.argv[1]
        if not filename.endswith('.h5'):
            filename += '.h5'
    except IOError:
        print"'{0:s}' is not a valid file".format(filename)

    print "... Retriving data from {0:s} ...".format(filename)
    global filename
    global res
    res = mbb_emcee.mbb_results(h5file=filename)
    print "... Done reading file ..."
    fir_op = True if sys.argv[2].lower() == 'true' else False
    # convergence()
    get_Paramvalues()
    IR = get_computation_value(FIR=fir_op)
    Plot_Standard(saveFig=True)
    Plot_Chain(IR, FIR=fir_op, saveFig=True)
#    MCMC_scatterDots()
    PanelPlot(saveFig=True)
    chi2 = BestFit()






def PlotSED_twoXAxes():
    """
    Plot SED with top showing wavelength in observed frame, bottom x axis showing observed freqeucny.

    Need to fix a bunch of stuff in here
    """
    # look into res.parameter_chain('fnorm') for normalization
    #Normalization factor
    if self._xnorm > self._xmerge:
        self._normfac = self._fnorm * self._xnorm**self._alpha / \
            self._kappa
    else:
        expmfac = math.expm1(-(self._xnorm / self._x0)**self._beta)
        self._normfac = -self._fnorm * math.expm1(self._xnorm) / \
            (self._xnorm**3 * expmfac)





    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.xlabel('Frequency')      # bottom

    # Make a second axis to plot the Top axis = wavelength
    ax1 = plt.twiny(ax)

    # Convert Freq to wavelength
    frequecny_arr = np.linspace(100.e9, 1.e14, 1000)
    c = 3.e8        # replace with accurate speed of light
    # is a function of frequency_arr; use to intepolate the desire
    # wavelength_arr to plot
    _wavelength = [c / f for f in frequecny_arr]
    wavelength_arr = np.arange(0, 4000, 500)    # range of wavelength to show

    # Then interpolate to the (freq) values at which we want ticks.
    WavelgTicks = np.interp(wavelength_arr, _wavelength, frequecny_arr)
    ax1.set_xticks(WavelgTicks)
    ax1.set_xticklabels([str(v) for v in WavelgTicks])

    # Make both axes have the same start and end point.
    x0, x1 = ax.get_xlim()
    ax1.set_xlim(x0, x1)
    ax1.set_xlabel('Observed wavelength')

    plt.show()
