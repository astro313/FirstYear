"""
Return scientific values using derived values + extra Parmas (e.g. Line )

Under development
"""


import numpy as np
import scipy.constants as scConst

"""
optionparser
r31 = CO (3-2) (1-0)= 0.7 or 1?
L_CO = ?
L_FIR =
L_IR =
D_L =
M_dust =
alpha_conversion = 0.8 = "mass-to-light ratio" or "CO-to-H2 conversion factor"
beta =
"""


class prettyGalaxies():
    def __init__(self):
        self._lir = blah
        self._D_L = blah

    def __str__(self):
        print '\n', '=' * 40
        print "All the results..."

    def I2T(self):
        """
        Convert velocity integrated flux density to brightness temperature in RJ

        Input:
        I_nu: flaot
            surface brightness of line in Jy
        freq_obs: float
            observed frequency in Hz

        Returns:
        --------
        T_b: float
            brightness temperature in observed frequency in K
        """
        I_nu = I_nu * 1.e-26
        T_b = I_nu * c**2 / (2 * scConst.k * freq_obs ** 2)

        return T_b

    def LineRatioCO(self):
        """
        Calcaulate CO line ratio in the RJ limit
        e.g. R_31 = T_32 / T_10
        If R = 1, corresponds to LTE
        If R < 1, gas is sub-thermally excited to the high J

        Inputs:
        -------
        freq_high: rest-frame CO of the higher J transition
        freq_low: rest-frame CO of the lower J transition


        Returns:
        --------
        R_ul_T: using brightness temperatures
        R_ul_I: using veloicity integrated line Intensities

        """
        R_ul_T = T_u / T_l
        R_ul_I = I_u / I_l * (freq_low / freq_high) ** 2
        return R_ul_I

    def L_transition(self):
        """
        Ref: Solomon & Vanden Bout 2005

        Input:
        ------
        D_L = luminosity distance in Mpc
        freq_res_line = line frequency in GHz in rest-frame
        I = velocity-integrated line flux in [Jy km/s] of CO(1-0)
        R = CO line ratio = I_{1-0} / I_{3-2}
        mu = magnification


        Returns:
        --------
        L: in [L_sun]
        """
        L = 1.04e-03 * (I * R) * D_L ** 2 * freq_res_line / (1 + z) / mu
        return L

    def Lprime_transition(self):
        """
        Ref: e.g., Solomon et al. 1992, 1997
        Calculate luminosity of line transition using flux of higher-J CO line

        Inputs:
        -------
        D_L: float
            luminosity distance in Mpc
        freq_res_line: float
            line frequency in GHz in rest-frame
        I: float
            velocity-integrated line flux in [Jy km/s] of CO(higher-J) (not 1-0)
        R: float
            CO line ratio = I_{1-0} / I_{3-2} = R_lu
        mu: flaot
            magnification

        Returns:
        L_line': float
            in [K km/s pc^2]
        """
        L_prime = 3.25e7 / freq_res_line ** 2 * (D_L ** 2) * (I * R) / (1 + z) / mu
        return L_prime

    def L_linePrime2line(self):
        """
        Convert from L' to L

        Input:
        ------
        freq_rest: float
            rest frequency of line in GHz
        Lprime: float
            L' of line [K km/s pc^2]
        Output:
        -------
        L: float
            L of line [L_sun]
        """
        L = 3.e-11 * freq_rest ** 3 * Lprime

        return L

    def X2alpha(self):
        """
        Converting from X_CO to alpha_CO for computing M_gas
        X Does not include contribution of Helium
        X_HE includes contriubtion of He in the molecular gas mass

        Inputs:
        -------
        X: float     [cm^-2 /K km/s]
        X_HE: float     [cm^-2 /K km/s]


        Returns:
        --------
        Default:
            alpha: float            [ M_sun pc^-2 (K km/s)^-1 ]
        """

        alpha = 6.3e-19 * X
        alpha_HE = 4.65e-19 * X_HE
        return alpha

    def M_gas(self):
        """
        ref: Downes & Solomon 1998

        Calculate cold H2 molecular gas masses from L_CO (lensing corrected)

        alpha_disk = 3.6    # z ~ 1.5 disk galaxies (Daddi et al. 2010)
        alpha_BzK_low = 2.5
        alpha_BzK_high = 4.0    # BzK with high uncertainties of +/- 1.4
        alpha = 0.8     # local ULIRGs, typically adopted for SMGs
        alpha = 4.6     # MW, optically thick ISM

        Returns:
        --------
        M_gas ()
        """

        M_gas = alpha * L_primeOne2zero
        return M_gas

    def SFR(self, IMF='Chabrier'):
        """
        Calculate the SFR using L_IR (8-1000) or (42.5-122.5) avoid AGN
        ref: Kennicutt 1998 conversion using 1-100 M_sun IMF of some sort

        Inputs:
        -------
        L_IR: [L_sun]

        Returns:
        --------
        SFR_IR: [M_sun / yr]

        """

        if IMF == 'Chabrier':
            SFR = 1.0e-10 * L_IR
        if IMF == 'Salpeter':
            SFR = 1.71e-10 * L_IR
        return SFR


    def M_dust(self, kappa=2.64):
        """
        Compute dust mass

        ref: Greve 2012 for M_dust
        ref: mbb_emcee uses different formular (Riechers 2013), calculated using chain in mcmc, Kappa default = 2.64 m^2/kg^-1 at 125 micron

        Inputs:
        -------
        v_rest: float
            GHz

        Returns:
        --------
        M_dust
        """
        def planck(T):
            x = scConst.h * freq / scConst.k / T
            B = 2 * scConst.h * freq ** 3 / c**2 / (np.exp(x) - 1)
            return B

        kappa_greve = 0.045 / ((v_rest/250.)**beta)
        B_T = planck(T_d)
        T_CMB_z = 2.73 * (1+z)
        B_CMB = planck(T_CMB_z)
        M_dust = D_L ** 2 * I * (B_T - B_CMB) ** (-1) / (1+z) / kappa_greve / mu
        return M_dust

    def SFE(self):
        """
        Forms:
        Some author uses L_FIR / M_gas, L_IR / L'_CO , SFR / M_H2
        We use the more direct more without conversion factor alpha and IMF assumption: L_IR / L'_CO

        SFE values for ULIRGs and distant SMGs typically exceed 100 L⊙ (K km s−1 pc^2)−1
        (e.g., Yao et al. 2003; Neri et al. 2003; Greve et al. 2005; Bouch´e et al. 2007).


        Inputs:
        -------
        L_IR: float
            Luminosity in FIR or IR?
        L_prime: float
            CO luminosity
        Returns:
        --------
        SFE: float
        """

        SFE = L_IR / L_prime
        return SFE


    def depleTime_Gas(self):
        """
        Compute depletion time scale assuming some value for alpha in computing M_gas = alpha * L_CO

        Returns:
        --------
        t: float        [years]

        """
        t = M_gas / SFR
        return t

    def depleTime_ISM(self):
        """
        Computing the depletion time using the mass of the ISM equation from Scoville 2014, where M_ISM = M_HI + M_H2

        Inputs:
        ------
        mu: float
            magnification factor
        D_L: float
            Luminosity distance Mpc
        beta: float
        freq_obs: float
            GHz
        I_obs: float
            flux density in observed frequency

        Returns:
        --------
        tau: float

        """

        M_ISM = 1.2e4 * D_L ** 2 * (350./ freq_obs) ** beta * (1 + z)**(-(1+beta)) * I_obs / mu
        tau = M_ISM / SFR
        return tau

    def f_mol_total():
        """
        compute the ratio of molecular gas to total gas
        """
        f = M_gas / M_ISM
        return f


    def M_SF(self):
        """
        ref: Scoville 2004
        Compute lower limit for M_sf using the Eddington-limited star foramtion efficiency

        SFE_max = L_IR / M_SF < 500 [L_sun/M_sun] in units of L_sun [K km/s pc^2]^{-1}
        L_IR = lensing-corrected IR luminosity
        """
        SFE_Max = 500       # [L_sun / M_sun]
        M_SF = L_IR / SFE_max

        return M_SF


    def qFactor(self, FIR, S1dot4GHz, unit='flux'):
        """
        Compute the FIR-radio correlation q to distinguish star forming from AGN dominated regions

        ref: Bell (2003)

        Input:
        -------
        flux_FIR = luminosity in the fa infared
        Flux: radio 1.4GHz

        either in luminosity [L_sun] or flux units SI
        """

        a = np.log10()


        Lradio = 4*pi*D_L**2 * flux
        q = blah if unit =='flux' else (np.log10(FIR/9.8e15) - np.log10(S1dot4GHz))
        return q


def mat2LaTex(arr):
    """
    Convert numpy array into latex table format
    """
    tab = " \\\\\n".join([" & ".join(map(str, line)) for line in arr])
    return tab

def me(theta_Es, e_theta_Es, zlenses, zsources):
    """
    Created 2013 March 13
    Author: Shane Bussmann
    Purpose: Returns mass and velocity dispersion inside Einstein radius given an Einstein radius (in arcsec), source redshift, and lens redshift
    Inputs: must be in numpy array format!
    """
    from math import pi
    from astropy.cosmology import WMAP9 as cosmo
    import numpy

    ntarg = theta_Es.size
    M_E = numpy.zeros(ntarg)
    e_M_E = numpy.zeros(ntarg)
    vdisp = numpy.zeros(ntarg)
    e_vdisp = numpy.zeros(ntarg)

    for targ in numpy.arange(ntarg):

        zsource = zsources[targ]
        zlens = zlenses[targ]
        theta_E = theta_Es[targ]
        e_theta_E = e_theta_Es[targ]

        # luminosity distances
        d_LS = cosmo.luminosity_distance(zsource).value * 3.08e24
        d_LL = cosmo.luminosity_distance(zlens).value * 3.08e24

        # comoving distances
        d_MS = d_LS / (1 + zsource)
        d_ML = d_LL / (1 + zlens)

        # angular diameter distances
        d_ALS = 1 / (1 + zsource) * ( d_MS - d_ML )
        d_AL = d_LL / (1 + zlens)**2
        d_AS = d_LS / (1 + zsource)**2

        # einstein radius in cm (7.1 kpc/" at z=0.7)
        theta_E_cm = theta_E / 206265. * d_AL
        e_theta_E_cm = e_theta_E / 206265. * d_AL

        # get a distribution of Einstein radii
        niters = 1e3
        theta_E_iters = numpy.random.normal(loc=theta_E_cm, scale=e_theta_E_cm, size=niters)

        # compute the mass enclosed within the Einstein radius
        c = 3e10
        G = 6.67e-8
        sigma_crit = c**2 / 4 / pi / G * d_AS / d_AL / d_ALS
        M_E_iters = pi * sigma_crit * theta_E_iters**2 / 2e33
        M_E[targ] = numpy.mean(M_E_iters)
        e_M_E[targ] = numpy.std(M_E_iters)

        vdisp2 = theta_E_iters / d_AL / 4 / pi * c**2 * d_AS / d_ALS
        vdisp[targ] = numpy.mean(numpy.sqrt(vdisp2) / 1e5)
        e_vdisp[targ] = numpy.std(numpy.sqrt(vdisp2) / 1e5)

    return M_E, e_M_E, vdisp, e_vdisp
