from NumCrunching import prettyGalaxies
import numpy as np

SMG = prettyGalaxies('numinput.yaml')
SMG.Sradio_obs = 0.45655848e-3   # Jy, before lensing correction
SMG.D_L()
SMG.Lprime_transition()
SMG.L_transition()
SMG.M_gas()
SFR_IR = SMG.SFRFunc(SMG.LIR)
SFR_FIR = SMG.SFRFunc(SMG.LFIR)
SFE_IR = SMG.SFEFunc(SMG.LIR)
SFE_FIR = SMG.SFEFunc(SMG.LFIR)
t_IR = SMG.depleTime_Gas(SFR_IR)
t_FIR = SMG.depleTime_Gas(SFR_FIR)
# q_IR, Lradio = SMG.qFactor_normL(SMG.LIR)
# q_FIR, Lradio = SMG.qFactor_normL(SMG.LFIR)
gas_dust = SMG.f_gas_dust()
M_dyn = SMG.Mdyn()
f_gas_dyn = SMG.f_molGas_dyn()
print SMG
# print "SFR using LIR: {:.2f} [M_sun/yr]".format(SFR_IR)
print "SFR using LFIR: {:.2f} [M_sun/yr]".format(SFR_FIR)
# print "SFE using LIR: {:.2f}".format(SFE_IR)
print "SFE using LFIR: {:.2f}".format(SFE_FIR)
# print "depletion Time using IR {:.2f} Myr".format(t_IR/1e6)
print "Depletion Time using FIR {:.2f} Myr".format(t_FIR/1e6)
# print "rest frame 1.4 GHz Luminosity: {:.2f} [W/Hz]".format(Lradio)
# print "q correlation using IR: {:.2f}".format(q_IR)
# print "q correlation using FIR: {:.2f}".format(q_FIR)
print "gas_dust fraction: {:.2f}".format(gas_dust)
print "Dynamical mass: {:.2f} x 10^10 [M_sun] ".format(M_dyn/1e10)
print "gas to dyn mass fraction: {:.2f}".format(f_gas_dyn)

r_E0 = 1.218  # 1.2230522136336528 (Best)
r_E1 = 0.745  # 0.73374533852753454 (Best)

r1 = np.array([r_E0])
r2 = np.array([r_E1])
e_r = np.array([0.01012153668])
e_r2 = np.array([0.0150122127])
z_arr = np.array([0.685])
z_source = np.array([2.221])
x = SMG.me(r1, e_r, z_arr)
x2 = SMG.me(r2, e_r2, z_arr)

print x, x2
