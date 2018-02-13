#
# Basic object for having datastream
#
import cPickle
import os
import numpy as np
from scipy.interpolate import interp1d
import astropy.units as u
from MapReader import *
from ObserveSky import *
from astropy.coordinates import Angle, SkyCoord

def getTimeList(tstart=Time('2016-08-01 00:00:00')+4*u.hour, dt=1, Ns=3600):
    """ helper routine that gets a list of times
    tstart -- start, default if 1st Aug 2016 midnight, astropy Time object
    dt -- delta t in seconds
    Ns -- number of samples, integer
    """
    ## in august we are 4 hours of UTC, otherwise 5
    return [tstart+i*dt*u.s for i in range(Ns)]


def get_stream(reader, telescope, skyc_list, i, fld, psources, Npix, Nfwhm):

    # Read input map
    hmap = reader.gethmap(fld, i)

    # Generate time stream
    perfreqstreams=getIntegratedSignal(telescope, skyc_list, hmap, reader.freq[fld][i], Npix=Npix, Nfwhm=Nfwhm)

    # Add point sources if requested. This is probably broken for now.
    if (psources is not None):
        perfreqs=getPointSourceSignal(telescope, skyc_list, psources, nu)
        for i,s in enumerate(perfreqs):
            perfreqstreams[i]+=s

    return perfreqstreams


class DataStream(object):
    def __init__ (self, telescope, tlist=getTimeList(), tag=None, sn=None,
                  nuparams=(800,1400,1), decdither=None):
        """ Constructor:
            telescope : TelescopeBase object
            tlist : list of times used here in the format
                    returned b getTimeList
            tag : tag string of reduced data to simulate, overrides tlist if not None
            nuparams: three elment tuple specifying (numin, numax, dnu) in MHz. 
            decdither: dec dither in degrees
            sn: 5 digit serial number used to save simulation results

        """
        self.telescope=telescope
        self.tlist=tlist
        Nbeams=len(self.telescope.beams)
        self.nulist=np.arange(nuparams[0], nuparams[1]+nuparams[2], nuparams[2])*1.0
        self.streams=[[None]*Nbeams]
        self.sn = sn
        
        self.has_tag = False
        if tag is not None:
            # Simulate an actual BMX reduced data file on disk
            self.has_tag = True
            from bmxreduce import datamanager
            self.dm = datamanager()
            reduced_fname = self.dm.getreducedfname(tag)
            reduced_file = np.load(reduced_fname)
            self.reduced_data = dict(reduced_file)
            reduced_file.close()
            # Set gain to 1 since we will be simulating in mK
            self.reduced_data['g'][:] = 1.0
            self.tlist = [Time(i, format="mjd") for i in self.reduced_data["mjd"]]
            # Set nulist, we will interpolate to it later
            self.nulist = self.reduced_data['f']

        # Pre-calculate sky coordinates of pointings to reduce overhead in the
        # per-frequency for loop in fillStream
        self.skyc_list = [] # two dimensional array
        self.skyc_gal_list = [] # two dimensional array
        for beam in telescope.beams:
            skycs = []
            skycs_gal = []
            for t in self.tlist:
                aaz = beam.AltAz(t, telescope.location)
                skyc = aaz.transform_to(apc.ICRS)
                # Dither in Dec
                if decdither is not None:
                    skyc = SkyCoord(ra=skyc.ra, dec=skyc.dec + decdither*u.deg,
                                    frame = 'icrs')
                skyc_gal = skyc.transform_to(apc.Galactic)
                skycs.append(skyc)
                skycs_gal.append(skyc_gal)
            self.skyc_list.append(skycs)
            self.skyc_gal_list.append(skycs_gal)


    def fillStream(self, whichfield='gfree+gsync/colore/hi4pi', mapsn='000', Npix=201,
                   Nfwhm=3, psources=None, parallel=False):
        """
        Fills stream using ObserveSky.

        For Npix, Nfwhm, see ObserveSky.py
        if psoruces is not None, add point sources
        parallel = True/False use celery or not
        mapsn = 3 digit serial number string, or list of strings, designating
                serial number of input maps to read

        """

        # Get input map properties
        reader = MapReader(whichfield, self.nulist, mapsn)

        ## set the data fields
        v = {}
        for fld in reader.fields:
            if parallel:
                import celery
                from bmxsim.celery_tasks import get_stream_celery
                task_list = celery.group([get_stream_celery.s(reader, self.telescope, self.skyc_gal_list, i, fld, psources, Npix, Nfwhm) for i,val in enumerate(reader.freq[fld])])
                task_promise = task_list()
                x = task_promise.get()
            else:
                x = [get_stream(reader, self.telescope, self.skyc_gal_list, i, fld, psources, Npix, Nfwhm) for i,val in enumerate(reader.freq[fld])]
                
            v[fld] = np.squeeze(np.array(x)).T
    
        # Get mean of fine frequency bins, which have now been averaged
        reader.collapsefreq()

        # Interpolate to frequency grid
        v_interp = {}
        for fld in reader.fields:
            vi = np.zeros((len(self.tlist),len(self.nulist)))*np.nan
            for i in range(len(self.tlist)):
                f = interp1d(reader.freq[fld], v[fld][i], kind="cubic", fill_value=np.nan, bounds_error=False)
                vi[i] = f(self.nulist)
            v_interp[fld] = vi

        if 0:
            # Make diagnostic plots
            import matplotlib.pyplot as plt
            nfld = len(reader.fields)
            plt.figure(figsize=(20,12))
            for k,fld in enumerate(reader.fields):
                plt.subplot(nfld,1,k+1)
                plt.plot(reader.freq[fld], v[fld][0])
                plt.plot(self.nulist, v_interp[fld][0],'.')
                plt.xlim(self.nulist.min(),self.nulist.max())
                plt.title(fld)
            plt.savefig('interp_diagnostic.png')

        # Now sum fields. Multiple beams not (re)implemented yet.
        for k,fld in enumerate(reader.fields):
            if k==0:
                self.streams[0] = v_interp[fld]
            else:
                self.streams[0] = np.nansum(np.stack((self.streams[0], v_interp[fld])), 0)

        # Now convert mK -> K and replace data if simulating real BMX tag. Save
        # ra/dec too.
        if self.has_tag:
            for k in range(self.reduced_data['data'].shape[0]):
                self.reduced_data['data'][k] = self.streams[0]*1e-3
            self.reduced_data['ra'] = np.array([c.ra.deg for c in self.skyc_list[0]])
            self.reduced_data['dec'] = np.array([c.dec.deg for c in self.skyc_list[0]])
            self.reduced_data['l'] = np.array([c.l.deg for c in self.skyc_gal_list[0]])
            self.reduced_data['b'] = np.array([c.b.deg for c in self.skyc_gal_list[0]])

        self.fields = whichfield.replace('+',';').replace('/',';').split(';')

    def save(self):
        """save the data with name format defined in datamanager.
        This function will create directories if necessary."""
        filename = self.dm.getreducedsimfname(str(self.reduced_data["tag"]), self.sn, self.fields, pickle_file=not self.has_tag)
        dirname = os.path.dirname(filename)
        if not os.path.isdir(dirname): 
            os.makedirs(dirname) # Create Serial/YYMM directory
            os.chmod(dirname, 0o775)

        if self.has_tag:
            np.savez(filename, **self.reduced_data)
        else:
            with open(filename, "wb") as f:
                cPickle.dump(self, f)

        os.chmod(filename, 0o664)
