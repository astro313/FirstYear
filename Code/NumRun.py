from NumCrunching import prettyGalaxies, me
import numpy as np

SMG = prettyGalaxies()
SMG.LIR = 88.52e12      # L_sun, before lensing correction
SMG.LFIR = 53.33e12      # L_sun, before lesning correction
SMG.z = 2.2214
SMG.mu = 9.7368538
SMG.beta = 1.90267276
SMG.Sradio_obs = 0.45655848e-3   # Jy, before lensing correction
SMG.FWHM_Line = 541.65      # km/s
SMG.I_line = 14.6       # CO 3-2 [mJy km/s /beam]
SMG.M_dust = 50.47e8
SMG.D_L()
freq_CO32_rest = 345.651e9
SMG.Lprime_transition(freq_CO32_rest)
SMG.L_transition(freq_CO32_rest)
SMG.M_gas()
SFR_IR = SMG.SFRFunc(SMG.LIR)
SFR_FIR = SMG.SFRFunc(SMG.LFIR)
SFE_IR = SMG.SFEFunc(SMG.LIR)
SFE_FIR = SMG.SFEFunc(SMG.LFIR)
t_IR = SMG.depleTime_Gas(SFR_IR)
t_FIR = SMG.depleTime_Gas(SFR_FIR)
q_IR = SMG.qFactor_normL(SMG.LIR)
q_FIR = SMG.qFactor_normL(SMG.LFIR)
gas_dust = SMG.f_gas_dust()
R_halflight = 0.023484507515718567    # arscec, from lens model of dust continuum
M_dyn = SMG.Mdyn(R_halflight)
f_gas_dyn = SMG.f_molGas_dyn()
print SMG
print "SFR using LIR: {:.2f} [M_sun/yr]".format(SFR_IR)
print "SFR using LFIR: {:.2f} [M_sun/yr]".format(SFR_FIR)
print "SFE using LIR: {:.2f}".format(SFE_IR)
print "SFE using LFIR: {:.2f}".format(SFE_FIR)
print "depletion Time using IR {:.2f} Myr".format(t_IR/1e6)
print "Depletion Time using FIR {:.2f} Myr".format(t_FIR/1e6)
print "q correlation using IR: {:.2f}".format(q_IR)
print "q correlation using FIR: {:.2f}".format(q_FIR)
print "gas_dust fraction: {:.2f}".format(gas_dust)
print "Dynamical mass: {:.2f} M_sun ".format(M_dyn/1e10)
print "gas to dyn mass fraction: {:.2f}".format(f_gas_dyn)

r_E0 = 1.2230522136336528
r_E1 = 0.73374533852753454

r1 = np.array([r_E0])
r2 = np.array([r_E1])
e_r = np.array([0.01012153668])
e_r2 = np.array([0.0150122127])
z_arr = np.array([0.685])
z_source = np.array([2.2214])
x = me(r1, e_r, z_arr, z_source)
x2 = me(r2, e_r2, z_arr, z_source)

print x, x2
