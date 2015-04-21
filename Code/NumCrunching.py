"""
Return scientific values using derived values + extra Parmas (e.g. Line )

"""
import numpy as np
import scipy.constants as scConst

class prettyGalaxies:
    def __init__(self):
        self.LIR = None
        self.LFIR = None
        self.beta = None
        self.R_ul = 1.         # Riechers 2013 QSO host
        self.z = 0.
        self.mu = 0.
        self.alphaCO = 0.8
        self.alpha_radio = -0.80
        self.Sradio_obs = None
        self.FWHM_Line = None
        self.I_line = None
        self.M_dust = None

    def __str__(self):
        print '\n', '=' * 40
        print "All the results..."
        print 'LIR: {:.3f} * 10^12 [L sun]'.format(self.LIR/1e12)
        print 'LFIR: {:.3f} * 10^12 [L sun]'.format(self.LFIR/1e12)
        print 'beta: {:.2f}'.format(self.beta)
        print 'line ratio upper-lower: {:.2f}'.format(self.R_ul)
        print 'z: ', self.z
        print 'mu: ', self.mu
        print 'alpha_co conversion: ', self.alphaCO
        print 'alpha_radio spectral for q', self.alpha_radio
        print 'flux density observed frame 1.4Ghz: [mJy]', self.Sradio_obs*1.e3
        print 'line velocity: {:.3f} [km/s]'.format(self.FWHM_Line)
        print 'Line integrated Intensity: {:.2f} [Jy km/s /beam]'.format(self.I_line)
        print 'gas mass: {:.3f} * 10^10 [M_sun]'.format(self.M_gas/1.e10)
        print "Luminosity Distance: {:.3f} Mpc Using WMAP9 Cosmo".format(self.lum_dist)
        print "L'_co (1-0) {:.2f} * 10^10 [K km/s pc^2]".format(self.L_prime/1.e10)
        return '=' * 40

    def D_L(self):
        """
        calculate luminosity distance based on redshift
        """
        from astropy.cosmology import WMAP9 as cosmo
        self.lum_dist = cosmo.luminosity_distance(self.z).value

    def I2T(self, freq_obs):
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
        I_line_SI = self.I_line * 1.e-26
        T_b = I_line_SI * scConst.c ** 2 / (2 * scConst.k * freq_obs ** 2)
        return T_b

    def LineRatioCO(self, upper, lower, upper_freq, lower_freq, UseTemp=None):
        """
        Calcaulate CO transition lines ratio in the RJ limit
        e.g. R_31 = T_32 / T_10
        If R = 1, corresponds to LTE
        If R < 1, gas is sub-thermally excited to the high J

        Inputs:
        -------
        upper: float
            upper transition velocity integrated line Intensity or upper transition brightness temperature
        lower: float
            lower transition velocity integrated line Intensity or lower transition brightness temperature
        freq_high: float [Hz]
            rest-frame CO of the higher J transition
        freq_low: float  [Hz]
            rest-frame CO of the lower J transition
        UseTemp: bool
            compute R using brightness temperature ratios instead


        Returns:
        --------
        R_ul: using brightness temperatures
        R_ul: using veloicity integrated line Intensities

        """
        if UseTemp:
            R_ul = T_u / T_l
        else:
            R_ul = I_u / I_l * (freq_low / freq_high) ** 2
        return R_ul

    def L_transition(self, freq_res_line):
        """
        Ref: Solomon & Vanden Bout 2005

        Input:
        ------
        freq_res_line: float
            line frequency in Hz in rest-frame
        I: float
            velocity-integrated line flux in [Jy km/s] of CO(3-2)
        R_ul: float
            CO line ratio = I_{3-2} / I_{1-0}
        mu: float

        Returns:
        --------
        L: in [L_sun]
        """

        freq_res_line /= 1e9
        self.L = 1.04e-03 * (self.I_line / self.R_ul) * self.lum_dist ** 2 * freq_res_line / (1 + self.z) / self.mu

    def Lprime_transition(self, freq_res_line):
        """
        Ref: e.g., Solomon et al. 1992, 1997
        Calculate luminosity of line transition using flux of higher-J CO line

        Inputs:
        -------
        D_L: float
            luminosity distance in Mpc
        freq_res_line: float
            line frequency in Hz in rest-frame
        I: float
            velocity-integrated line flux in [Jy km/s] of CO(higher-J) (not 1-0)
        R: float
            CO line ratio = I_{3-2} / I_{1-0} = R_ul
        mu: flaot
            magnification

        Returns:
        L_line': float
            in [K km/s pc^2]
        """

        freq_res_line /= 1e9
        self.L_prime = 3.25e7 / freq_res_line ** 2 * (self.lum_dist ** 2) * (self.I_line / self.R_ul) / (1 + self.z) / self.mu

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

    def X2alpha(self, X):
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
        # alpha_HE = 4.65e-19 * X_HE
        return alpha

    def M_gas(self, alpha=None):
        """
        ref: Downes & Solomon 1998

        Calculate cold H2 molecular gas masses from L_CO (lensing corrected)
        alpha [M_sun K km/s pc^2]-1
        alpha_disk = 3.6    # z ~ 1.5 disk galaxies (Daddi et al. 2010)
        alpha_BzK_low = 2.5
        alpha_BzK_high = 4.0    # BzK with high uncertainties of +/- 1.4
        alpha = 0.8     # local ULIRGs, typically adopted for SMGs
        alpha = 4.6     # MW, optically thick ISM

        Returns:
        --------
        M_gas [M_sun]
        """
        if alpha is None:
            alpha = 0.8
        self.M_gas = alpha * self.L_prime

    def SFRFunc(self, L_IR, IMF='Chabrier'):
        """
        Calculate the SFR using L_IR (8-1000) or (42.5-122.5) avoid AGN dust heating
        ref: Kennicutt 1998 conversion using 1-100 M_sun IMF of some sort

        Inputs:
        -------
        L_IR: [L_sun]

        Returns:
        --------
        SFR_IR: [M_sun / yr]

        """

        if IMF == 'Chabrier':
            SFR = 1.0e-10 * L_IR / self.mu
        if IMF == 'Salpeter':
            SFR = 1.71e-10 * L_IR / self.mu
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

    def SFEFunc(self, L_IR):
        """
        Forms:
        Some author uses L_FIR / M_gas, L_IR / L'_CO , SFR / M_H2
        We use the more direct more without conversion factor alpha and IMF assumption: L_IR / L'_CO

        SFE values for ULIRGs and distant SMGs typically exceed 100 L_sun [K km/s pc^2]-1
        (e.g., Yao et al. 2003; Neri et al. 2003; Greve et al. 2005; Bouche et al. 2007).


        Inputs:
        -------
        L_IR: float
            Luminosity in FIR or IR [L_sun]
        L_prime: float
            CO luminosity [K km/s pc^2]-1
        Returns:
        --------
        SFE: float
        """

        SFE = (L_IR / self.mu) / self.L_prime
        return SFE

    def depleTime_Gas(self, SFR):
        """
        Compute depletion time scale assuming some value for alpha in computing M_gas = alpha * L_CO

        Inputs:
        -------
        SFR: float   [M_sun/yr]

        Returns:
        --------
        t: float        [years]

        """
        t = self.M_gas / SFR
        return t

    def depleTime_ISM(self, freq_obs, SFR):
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

        self._M_ISM = 1.2e4 * self.lum_dist ** 2 * (350. / freq_obs) ** self.beta * (1 + self.z)**(-(1 + self.beta)) * I_obs / self.mu
        tau = self._M_ISM / SFR
        return tau

    def Mdyn(self, r):
        """
        Compute dynamical mass using "isotropic virial estimator"

        Inputs:
        --------
        FWHM_line: float
            [km/s] of line
        r: float
            half light radius [kpc], major axis of line (resolved)
        """
        self.M_dyn = 2.8e5 * (self.FWHM_Line) ** 2 * r

    def f_molGas_dyn(self):
        """
        Compute gas mass fraction = M_gas / M_dyn, M_dyn includes stellar mass (see notes)

        Inputs:
        -------
        M_gas: float
            molecular gas mass
        M_dyn: float
            Dynamical gas mass

        Returns:
        --------
        molecular gas mass fraction
        """
        return self.M_gas / self.M_dyn

    def f_mol_total(self):
        """
        compute the ratio of molecular gas mass to total gas mass, different from ratio of gas mass to dynamical gas mass

        Inputs:
        -------
        M_ISM: float
            ISM that includes mass of HI and H2
        M_gas: float
            H2 molecular gas mass
        Returns:
        --------
        f: float
        """
        M_ISM = self._M_ISM
        f = self.M_gas / M_ISM
        return f

    def f_gas_dust(self):
        """
        compute the ratio of molecular gas mass to dust mass

        Inputs:
        -------
        M_gas: float
            H2 molecular gas mass [M_sun]
        M_dust: float
            dust mass [M_sun]
        Returns:
        --------
        f: float
        """
        f = self.M_gas / self.M_dust / self.mu
        return f

    def M_SF(self):
        """
        ref: Scoville 2004
        Compute lower limit for M_sf using the Eddington-limited star foramtion efficiency

        Inputs:
        -------
        Default -- SFE_max = L_IR / M_SF < 500 [L_sun/M_sun] in units of L_sun [K km/s pc^2]^{-1}
        L_IR: float
            lensing-corrected IR luminosity

        Returns:
        --------
        M_SF: float
            lower limit of star forming mass
        """
        SFE_Max = 500.       # [L_sun / M_sun]
        M_SF = L_IR / self.mu / SFE_max
        return M_SF

    def qFactor_Helou(self, FIR, S_radio):
        """
        Compute the FIR-radio correlation q to distinguish star forming from AGN dominated regions

        q = log (S_FIR / 3.75 e 12 [W m^-2]) - log (S_{1.4 GHz} / 1 [W m^-2 Hz^-1] )
        where S_{1.4GHz} is the monochromatic rest-frame 1.4 GHz luminosity (Helou et al. 1985; Condon 1992) and S_FIR is luminosity which was originally computed from the rest-frame 60- and 100- micron IRAS fluxes (Helou et al. 1985) under the assumption of a typical dust temperature of ~ 30 K.

            However, these definitions are not practical for ULIRGs that have higher dust temperatures, and therefore use the integrated S_FIR of the warm dust component.

        Input:
        -------
        S_radio: float
            observed frame

        """
        q = np.log10(S_FIR / 3.75e12) - np.log10(S_radio)
        return q

    def qFactor_normL(self, L_FIR):
        """
        Compute the q factor using integrated FIR luminosity, rest-frame 1.4 GHz

        ref: Riecher et al. 2013

        q = log10(L_FIR / 9.8e-15 [L_sun]) - log10(L_1.4 /[W Hz^-1])
        Input:
        -------
        L_FIR: float [L_sun]
            luminosity in the far infared (integrated L), some use 8-1000 micron, most use 40-1000micron for cold dust, star formation
        L_radio: float
            luminosity at 1.4GHz, k-corrected (rest frame) [W/Hz]

        self.alpha_radio: negative float


        """
        Mpc2m = 3.08567758e22
        Jy2SI = 1.e-26
        L_radio = 4 * np.pi * (self.lum_dist * Mpc2m) ** 2 * (1+self.z) ** (-(1 + self.alpha_radio)) * (self.Sradio_obs * Jy2SI)
        q = np.log10(L_FIR / self.mu / 9.8e-15) - np.log10(L_radio)
        return q


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