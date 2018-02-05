#!/usr/bin/env python
from datetime import datetime, timedelta
import numpy as np
import sys,pylab
import bmxsim as bs
import cPickle as cP
import argparse
from astropy.time import Time
import astropy.units as u

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--dt", dest="dt", default=10,
                  help="dt for timesteps (in seconds!)", type=float)
parser.add_argument("--tmax", dest="tmax_h", default=1.,
                  help="Length of observation (in hours!)", type=float)
parser.add_argument("--field", dest="whichfield",default="cosmo",
                  help="Which field to simulate;  cosmo, ptso, gfree, egfree or none,"+
                  "combine the with +;  ", type=str)
parser.add_argument("--NVSS", dest="nvss", default=False, action="store_true",
                  help="Add NVSS 1Jy+ sources.")
parser.add_argument("--sun", dest="sun", default=False, action="store_true",
                  help="Add Sun")
parser.add_argument("--moon", dest="moon", default=False, action="store_true",
                  help="Add Moon")
parser.add_argument("-p", "--parallel", dest="parallel", default=False, action="store_true", help="Run in parallel")
parser.add_argument("--tag", dest="tag", type=str, default=None, help="BMX reduce file tag")
parser.add_argument("--start", dest="start", type=str, default="2016-08-01 00:00:00", help="observation start time")
parser.add_argument("--serial", dest="serial", type=str, required=True, help="Five digit serial number in the file name")

o = parser.parse_args()

telescope=bs.BMXDemoSingleDish()

if o.tag:
    stream = bs.DataStream(telescope, tag=o.tag, serial=o.serial)
else:
    start_time = Time(o.start) + 4*u.hour
    tlist=bs.getTimeList(tstart=start_time, dt=o.dt, Ns=int((o.tmax_h*3600.)/o.dt)+1)
    stream=bs.DataStream(telescope, tlist=tlist, serial=o.serial)

#point sources
psources=None
if (o.nvss or o.moon or o.sun):
    psources=bs.PointSourceCatalog(o.nvss,o.sun, o.moon)
#fields
field=o.whichfield
if field=="none":
    field=None
stream.fillStream(reader=None,whichfield=field,Npix=201,Nfwhm=3,psources=psources, parallel=o.parallel)
stream.save()
