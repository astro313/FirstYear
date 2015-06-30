freq_CO32 = 345.651     # GHz
dl = 19052          # Mpc at 2.221
z = 2.221
I_line = 10.7
del_I = 2.12    # sigma of CO(3-2) line intensity Jy km/s
mu = 10.1298
del_mu = 1.38
alpha = 0.8
rkpc = 0.105907102426 * 8.406
del_kpc = 0.03309 * 8.406
vel = 541.65
del_vel = 31.9
LIR = 88.52e12    # before lensing correction, from mbb
del_IR_noMu = 2.62e12
M_dust = 50.47e+8  # before lensing correction, from mbb
sig_Mdust_noMu = 20.42e8
M_dyn = 7.35e10
Mgas = 2.326e10    # with lensing correction

#####################
#####################
sig_LprimeCO = (3.25e7 / freq_CO32**2 * dl**2 / (1+z))**(.5) * ((del_I/mu)**2 + (del_mu/mu**2*I_line)**2)**(.5)
del_LprimeCO = 3.25e7 / freq_CO32**2 * dl**2 / (1+z) * ((del_I/mu) + (del_mu/mu**2*I_line))

sig_LprimeCO_noMu = (3.25e7 / freq_CO32**2 * dl**2 / (1+z))**(.5) * (del_I/mu)
del_LprimeCO_noMu = 3.25e7 / freq_CO32**2 * dl**2 / (1+z) * (del_I/mu)

LprimCO_noMu = 3.25e7/freq_CO32**2 * dl**2 / (1+z) * I_line

sig_Mgas = alpha * sig_LprimeCO
del_Mgas = alpha * del_LprimeCO

Mgas_noMu = alpha * LprimCO_noMu
sig_Mgas_noMu = alpha * sig_LprimeCO_noMu
del_Mgas_noMu = alpha * del_LprimeCO_noMu

# With lensing correction M_dust
sig_Mdust = ((M_dust/mu**2 * del_mu)**2 + (sig_Mdust_noMu/mu)**2)**(.5)
del_Mdust = M_dust/mu**2 * del_mu + sig_Mdust_noMu/mu

sig_Mdyn = (2.8e5)**(.5) * ((2. * vel * del_vel * rkpc)**2 + (vel**2 * del_kpc)**2)**(.5)
del_Mdyn = 2.8e5 * (2. * vel * del_vel * rkpc + vel**2 * del_kpc)

sig_SFR = (1.e-10)**(.5) * ((LIR/mu**2*del_mu)**2 + (del_IR_noMu/mu)**2)**(.5)
del_SFR = 1.e-10 * (LIR/mu**2*del_mu + del_IR_noMu/mu)

SFR_noMu = 1.e-10 * LIR
sig_SFR_noMu = (1.e-10)**(.5) * (del_IR_noMu/mu)
del_SFR_noMu = 1.e-10 * (del_IR_noMu/mu)

sig_tauDepl = ((sig_Mgas_noMu/SFR_noMu)**2 + (Mgas_noMu/SFR_noMu**2*sig_SFR_noMu)**2)**(.5)
del_tauDepl = (del_Mgas_noMu/SFR_noMu) + (Mgas_noMu/SFR_noMu**2*del_SFR_noMu)

sig_SFE = ((del_IR_noMu/LprimCO_noMu)**2 + (sig_LprimeCO_noMu/LprimCO_noMu**2*LIR)**2)**(.5)
del_SFE = del_IR_noMu/LprimCO_noMu + del_LprimeCO_noMu/LprimCO_noMu**2 * LIR

sig_fgasdust = ((Mgas_noMu/M_dust**2*sig_Mdust_noMu)**2 + (sig_Mgas_noMu/M_dust)**2)**(.5)
del_fgasdust = Mgas_noMu/M_dust**2*sig_Mdust_noMu + del_Mgas_noMu/M_dust

sig_fgasdyn = ((sig_Mgas/M_dyn)**2 + (Mgas/M_dyn**2*sig_Mdyn)**2)**(.5)
del_fgasdyn = del_Mgas/M_dyn + Mgas/M_dyn**2*del_Mdyn

# print "delta Lprime CO: ", sig_LprimeCO
# print "delta Mgas: ", sig_Mgas
# print "delta_Mdyn: ", sig_Mdyn
# print "delta SFR: ", sig_SFR
# print "delta_TAu: ", sig_tauDepl
# print "delta_SFE: ", sig_SFE
# print "delta f gas dust: ", sig_fgasdust
# print "delta f gas dyn: ", sig_fgasdyn


print "delta Lprime CO: ", del_LprimeCO/1.e10
print "delta Mgas: ", del_Mgas/1.e10
print "delta_Mdyn: ", del_Mdyn/1.e10
print "delta SFR: ", del_SFR
print "delta_TAu: ", del_tauDepl/1.e6
print "delta_SFE: ", del_SFE
print "delta f gas dust: ", del_fgasdust
print "delta f gas dyn: ", del_fgasdyn