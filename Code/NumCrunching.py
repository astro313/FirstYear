import numpy as np

"""
L_CO = ?
L_CO' = 
L_FIR = 
L_IR = 
M_dust =
alpha_conv = 
beta = 
"""


class prettygalaxies(object):
	def __init__(self):


	def M_dyn(FWHMvel, halfLightradius):
		"""
		halfLightradius ~ major axis FWHM of line (very rough)

		"""
		m = 2.8e5 * (FWHMvel)**2 * halfLightradius
		return m
	
	def qFactor(FIR, S1dot4GHz, unit='flux'):
		"""
		Compute the FIR-radio correlation q
	
		Input:
		-------
		flux_FIR = luminosity in the fa infared
		Flux: radio 1.4GHz 
	
		either in luminosity [L_sun] or flux units SI
		"""

		Lradio = 4*pi*D_L**2 * flux
		q = blah if unit =='flux' else (np.log10(FIR/9.8e15) - np.log10(S1dot4GHz))
		return q
	
	def SFRfromIR(IMF=Chabrier):
		if IMF == 'Chabrier':
			blah
		if IMF == 'Salpeter':
			blah
		alpha = blah
		return SFR
	
	def L_IR(optical=thick):
		"""
		L of 42.5 micron to 122.5 micron
		integrate SED
		"""
		if optical=='thick':
			general SED
		if optical=='thin':
			blah
		return
	
	def fraction(self):
		"""
		ratio of M_dust / M_gas
		"""

		return
	def f_gas():
		"""
		ratio of m_gas to m_dyn
		"""
		fgas = blah
		return fgas

	
	def SFR():
		return


	
	def M_dust(kappa=2.64):
		"""
		Compute dust mass
	
		Kappa default = 2.64 m^2/kg^-1 at 125 micron
	
		"""
		M_dust = 
		return