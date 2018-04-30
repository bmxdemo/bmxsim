#!/usr/bin/env python
from BMXDemo import BMXDemoSingleDish
from DataStream import *
from datetime import datetime, timedelta
import numpy as np
import argparse
from astropy.time import Time
import astropy.units as u

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("--tstart", dest="tstart", default="2016-08-01 00:00:00", 
                    help="observation start time in UTC", type=str)
parser.add_argument("--dt", dest="dt", default=10.,
                    help="dt for timesteps (in seconds!)", type=float)
parser.add_argument("--tmax", dest="tmax", default=1.,
                    help="Length of observation (in hours!)", type=float)
parser.add_argument("--numin", dest="numin", default=800.,
                    help="Min freq in MHz", type=float)
parser.add_argument("--numax", dest="numax", default=1400.,
                    help="Max freq in MHz", type=float)
parser.add_argument("--dnu", dest="dnu", default=1.,
                    help="dnu in MHz", type=float)
parser.add_argument("--nu_f", dest="nu_f", default=800.,
                    help="freq at which aperture is fully illuminated in MHz. At "+
                    "higher frequencies, aperture is considered under-illuminated "+
                    "and beam will not scale with frequency. Set to np.inf or a "+
                    "very high number to turn on beam chromaticity for all frequencies.", type=float)
parser.add_argument("--field", dest="field", default="colore",
                    help="Which field to simulate;  cosmo, ptso, gfree, egfree, colore, hi4pi; "+
                    "combine with +, in which case fields are summed", type=str)
parser.add_argument("--NVSS", dest="nvss", default=False,
                    help="Add NVSS 1Jy+ sources.", action="store_true")
parser.add_argument("--sun", dest="sun", default=False, 
                    help="Add Sun", action="store_true")
parser.add_argument("--moon", dest="moon", default=False, 
                    help="Add Moon", action="store_true")
parser.add_argument("--tag", dest="tag", default=None,
                    help="BMX file tag to simulate. If set, overrides tstart, dt, "+
                    "tmax, numin, numax, and dnu", type=str)
parser.add_argument("--sn", dest="sn", default='xxxxx', 
                    help="Five digit serial number in the file name", type=str)
parser.add_argument("--decdither", dest="decdither", default = 0.0,
                    help="Beam declination dither in degrees", type=float)
parser.add_argument("--Npix", dest="Npix", default = 201,
                    help="Beam postage stamp image will be Npix x Npix", type=int)
parser.add_argument("--Nfwhm", dest="Nfwhm", default = 3.0,
                    help="Beam postage stamp will extend out to Nfwhm", type=float)

o = parser.parse_args()

telescope = BMXDemoSingleDish(nu_f=o.nu_f)

#point sources
psources=None
if (o.nvss or o.moon or o.sun):
    psources = PointSourceCatalog(o.nvss, o.sun, o.moon)

# Frequency parameters
nuparams = (o.numin,o.numax,o.dnu)

# Time list (overridden internally in DataStream if tag is specified)
tlist = getTimeList(tstart=o.tstart, dt=o.dt, tmax=o.tmax)

stream = DataStream(telescope, tag=o.tag, sn=o.sn, tlist=tlist, nuparams=nuparams)
stream.fillStream(whichfield=o.field, Npix=o.Npix, Nfwhm=o.Nfwhm, psources=psources)
stream.save(o)
