#!/usr/bin/env python
import numpy as np
import sys,pylab
sys.path.append('py/')
import DataStream
from BMXDemo import BMXDemoSingleDish
import cPickle as cP

from optparse import OptionParser


parser = OptionParser()
parser.add_option("--dt", dest="dt", default=10,
                  help="dt for timesteps (in seconds!)", type="float")
parser.add_option("--tmax", dest="tmax_h", default=1.,
                  help="Length of observation (in hours!)", type="float")
parser.add_option("--field", dest="whichfield",default="cosmo",
                  help="Which field to simulate; to wit, cosmo, ptso, gfree, or egfree", type="string") 
parser.add_option("-o","--output", dest="outfile", default="testStream.pickle",
                  help="Output filename", type="string")
# default
for option in parser.option_list:
    if option.default != ("NO", "DEFAULT"):
        option.help += (" " if option.help else "") + "[default: %default]"
(o, args) = parser.parse_args()

telescope=BMXDemoSingleDish()
tlist=DataStream.getTimeList(dt=o.dt, Ns=int((o.tmax_h*3600.)/o.dt)+1)
stream=DataStream.DataStream(telescope,tlist=tlist)
stream.fillStream(reader=None,whichfield=o.whichfield,Npix=201,Nfwhm=3)
cP.dump(stream,open(o.outfile,'w'))
