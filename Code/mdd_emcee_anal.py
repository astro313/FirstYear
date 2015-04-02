import mbb_emcee
"""
Look at resulting fit from mbb_emcee
"""


class getMBB_Params(object):

    def __str__(self):
        print("blah blah ")
    return none



if parames not fixed:
    msg = "alpha: {:0.2f} +{:0.2f} -{:0.2f}"
    print(msg.format(par_central_values[3][0],      # unpacking from par_cent()
                    par_central_values[3][1],
                    par_central_values[3][2],
        ))

same as the following
T_val = res.par_cen('T')
print("Temperature/(1+z): {:0.2f}+{:0.2f}-{:0.2f} [K]".format(*T_val))
b_val = res.par_cen('beta')
print("Beta: {:0.2f}+{:0.2f}-{:0.2f}".format(*b_val))
res.parameter_chain('T')        #get flattened chain for parameter
res.par_cen('T')        # get the central confidence and +/- 1 sigma
if parames not fixed:
    msg = "alpha: {:0.2f} +{:0.2f} -{:0.2f}"
    print(msg.format(par_central_values[3][0],      # unpacking from par_cent()
                    par_central_values[3][1],
                    par_central_values[3][2],
        ))

# something that requires computation as mcmc runs
res.compute_peaklambda()
p_val = res.peaklambda_cen()        # similiar to par_cen()
print("Peak Obs wavelength: {:0.1f}+{:0.1f}-{:0.1f} [um]".format(*p_val))


# LIR compute, this part is sort of working.....
res.compute_lir()
lir = res.lir_cen()
args = (res._lir_min, res._lir_max) + tuple(lir)     # adding wavelength
lirstr = "L_IR({:0.1f} to {:0.1f}um): {:0.2f} "\
                     "+{:0.2f} -{:0.2f} [10^12 L_sun]\n"
print(lirstr).format(*args)
# L_IR(8.0 to 1000.0um): 89.52 +2.56 -2.59 [10^12 L_sun]


res.compute_lir(wavemin=42.5, wavemax=122.5)    # walkers again using stuff inside
lir = res.lir_cen()
args = (res._lir_min, res._lir_max) + tuple(lir)     # adding wavelength
lirstr = "L_FIR({:0.1f} to {:0.1f}um): {:0.2f} "\
                     "+{:0.2f} -{:0.2f} [10^12 L_sun]\n"
print(lirstr).format(*args)
# L_FIR(42.5 to 122.5um): 52.70 +1.13 -1.12 [10^12 L_sun]


if __name__ == '__main__':
    filename = "BLOOM.h5"
    res =  mbb_emcee.mbb_results(h5file= filename)       # improve using option parser and argvs