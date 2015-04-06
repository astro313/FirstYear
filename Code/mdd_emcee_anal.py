"""
For any random things I wanna test, dump below
"""

import mbb_emcee
import matplotlib.pyplot as plt
import sys
from pylab import savefig
import numpy

if len(sys.argv) < 2:
    errmsg = "Invalid number of arguments: {0:d}"
#    import pdb; pdb.set_trace()
    raise IndexError(errmsg.format(len(sys.argv)))

try:
     filename = sys.argv[1]
     if not filename.endswith('.h5'):
        filename += '.h5'
except IOError:
    print"'{0:s}'' is not a valid file".format(filename)

print "... Retriving data from {0:s} ...".format(filename)
res = mbb_emcee.mbb_results(h5file=filename)
print "... Done reading file ..."





