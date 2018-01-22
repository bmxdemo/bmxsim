#
# Basic object for having datastream
#
import cPickle
import os

import numpy as np
from scipy.interpolate import interp1d
import astropy.units as u
from CrimeReader import *
from ObserveSky import *

def getTimeList(tstart=Time('2016-08-01 00:00:00')+4*u.hour, dt=1, Ns=3600):
    """ helper routine that gets a list of times
    tstart -- start, default if 1st Aug 2016 midnight, astropy Time object
    dt -- delta t in seconds
    Ns -- number of samples, integer
    """
    ## in august we are 4 hours of UTC, otherwise 5
    return [tstart+i*dt*u.s for i in range(Ns)]

class DataStream(object):
    def __init__ (self, telescope, tlist=getTimeList(), tag=None, serial=None):
        """ Constructor:
            telescope : TelescopeBase object
            tlist : list of times used here in the format
                    returned b getTimeList
        """
        self.telescope=telescope
        self.tlist=tlist
        Nbeams=len(self.telescope.beams)
        dnu=1.0
        self.nulist=np.arange(telescope.numin, telescope.numax, dnu)
        self.streams=[[None]*Nbeams]
        self.serial = serial
        
        self.has_tag = False
        if tag is not None:
            self.has_tag = True
            from bmxreduce import datamanager
            self.dm = datamanager()
            reduced_fname = self.dm.getreducedfname(tag)
            reduced_file = np.load(reduced_fname)
            self.reduced_data = dict(reduced_file)
            reduced_file.close()
            self.tlist = [Time(i, format="mjd") for i in self.reduced_data["mjd"]]

        self.skyc_list = [] # two dimensional array
        for beam in telescope.beams:
            skycs = []
            for t in self.tlist:
                aaz = beam.AltAz(t, telescope.location)
                skyc = aaz.transform_to(apc.ICRS)
                skycs.append(skyc)
            self.skyc_list.append(skycs)


    def nuMin(self):
        return self.nulist[0]
    def nuMax(self):
        return self.nulist[-1]
    def tMax_h(self):
        return (self.tlist[-1]-self.tlist[0]).to(u.h).value



    def setNuList(self, nulist):
        """ sets nulist to nulist and updated parameters in the telescope obj
            nulist : list of frequencies [MHz]
        """
        self.nulist=nulist
        dnu=(nulist[-1]-nulist[0])/(len(nulist)-1)
        self.telescope.setNuParams(nulist[0], nulist[-1])

    def fillStream(self,reader=None,whichfield='cosmo+gfree+gsync',Npix=201, Nfwhm=3,
                   psources=None, parallel=False):
        """
        Fills stream using ObserveSky and crimereader objects.
        reader is a CrimerReader instance, if None, it will take its own
        field specifies which field in CrimeReader to use

        For Npix, Nfwhm, see ObserveSky.py
        if psoruces is not None, add point sources

        """
        if reader is None:
            reader=CrimeReader()
        self.setNuList(reader.freq)
        print "Set frequencies from crime reader:", self.telescope.numin, self.telescope.numax, len(self.nulist)

        ## set the data fields
        if parallel:
            import celery
            from .celery_tasks import get_stream
            task_list = celery.group([get_stream.s(self.telescope, self.skyc_list, nu, i, whichfield, psources) for i, nu in enumerate(self.nulist)])
            task_promise = task_list()
            task_results = task_promise.get()
        else:
            from .celery_tasks import get_stream
            task_results = [get_stream(self.telescope, self.skyc_list, nu, i, whichfield, psources) for i, nu in enumerate(self.nulist)]

        self.streams=[np.zeros((len(self.nulist),len(self.tlist)))]
        for i, perfreqstreams in enumerate(task_results):
            for b, stream in enumerate(perfreqstreams):
                self.streams[b][i, :] = stream

        if self.has_tag:
            new_f_idx = np.bitwise_and(self.reduced_data["f"] <= self.nuMax(), self.reduced_data["f"] >= self.nuMin())
            new_f = self.reduced_data["f"][new_f_idx]
            interpmat = np.zeros_like(self.reduced_data["data"])
            interpmat[:, :, :] = np.nan
            for i in range(len(self.tlist)):
                f = interp1d(self.nulist, self.streams[0][:, i], kind="cubic")
                interpmat[0, i, new_f_idx] = f(new_f)
            interpmat[1] = interpmat[0].copy()
            self.reduced_data["data"] = interpmat

        self.fields = whichfield.split("+")

    def save(self):
        """save the data with name format defined in datamanager.
        This function will create directories if necessary."""
        filename = self.dm.getreducedsimfname(str(self.reduced_data["tag"]), self.serial, self.fields, pickle_file=not self.has_tag)
        dirname = os.path.dirname(filename)
        if not os.path.isdir(dirname): 
            if not os.path.isdir(os.path.dirname(dirname)): 
                os.mkdir(os.path.dirname(dirname)) # Create Serial directory
                os.chmod(os.path.dirname, 0o775)
            os.mkdir(dirname) # Create Serial/YYMM directory
            os.chmod(dirname, 0o775)

        if self.has_tag:
            np.savez(filename, **self.reduced_data)
        else:
            with open(filename, "wb") as f:
                cPickle.dump(self, f)

        os.chmod(filename, 0o664)
