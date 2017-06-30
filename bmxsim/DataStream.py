#
# Basic object for having datastream
#
import numpy as np
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
    def __init__ (self, telescope, tlist=getTimeList()):
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
                   psources=None):
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
        self.streams=[np.zeros((len(self.nulist),len(self.tlist)))]
        for i,nu in enumerate(self.nulist):
            print "Doing ",i,nu
            if whichfield is not None:
                # Read input map
                field=reader.named_slice(whichfield,i)
                # Generate time stream
                perfreqstreams=getIntegratedSignal(self.telescope, self.tlist, field, nu, Npix=201, Nfwhm=3)
            else:
                # Return Zero
                perfreqstreams=[np.zeros(len(self.tlist)) for i in range(len(self.telescope.beams))]
                
            # Add point sources if requested
            if (psources is not None):
                perfreqs=getPointSourceSignal(self.telescope, self.tlist, psources, nu)
                for i,s in enumerate(perfreqs):
                    perfreqstreams[i]+=s

            # Populate stream object
            for b,stream in enumerate(perfreqstreams):
                self.streams[b][i,:]=stream
    
