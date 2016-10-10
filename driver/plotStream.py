#!/usr/bin/env python
import numpy as np
import sys,pylab
import bmxsim as bs
import cPickle as cP
from optparse import OptionParser
import matplotlib.pyplot as plt

parser = OptionParser(usage="Usage: %prog [options] picke_filename")
parser.add_option("-o", "--output", dest="outfile", default=None,
                  help="Output filename, if none plot to screen", type="string")
# default
for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
(o, args) = parser.parse_args()
if (len(args)!=1):
    parser.print_help()
    sys.exit(0)

st=cP.load(open(args[0]))
plt.figure(figsize=(10,10))
plt.ylabel("freq [MHz]")
plt.xlabel("t [h]")
plt.imshow(st.streams[0],interpolation='nearest',extent=(0.,st.tMax_h(),st.nuMin(), st.nuMax()),aspect='auto')
plt.colorbar()
if (o.outfile is None):
    plt.show()
else:
    plt.savefig(o.outfile)
