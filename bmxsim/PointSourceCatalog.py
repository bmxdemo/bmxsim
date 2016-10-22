#
# A class to deal with point sources, essentially NVSS sources + moon + sun
#
import pyfits, os
import numpy as np
import astropy as ap
import astropy.coordinates as apc
import astropy.units as u
import scipy.constants as const

class PointSourceCatalog(object):
    

    # constructor
    def __init__ (self, have_NVSS=True, have_sun=True, have_moon=True):
        self.have_NVSS=have_NVSS
        self.have_sun=have_sun
        self.have_moon=have_moon
        if have_NVSS:
            path=""
            if os.environ.has_key("BMXSIM_ROOT"):
                path=os.environ["BMXSIM_ROOT"]+"/"
            da=pyfits.open("data/nvss-1Jy+.fits")
            dat=da[1].data
            self.flux20=dat["FLUX_20_CM"]*1e-3## convert mJy to Jy
            self.skyc=apc.SkyCoord(ra=dat["RA"]*u.deg, dec=dat["DEC"]*u.deg)
            self.spec_index=-0.5*np.ones(len(self.skyc))
            self.nvss_nu0=1.4e3
        else:
            self.skyc=None
    def skyCoordList_fixed(self):
        return self.skyc

    def fluxes_mK_fixed(self,nu, beam):
        flux_Jy=self.flux20*(nu/self.nvss_nu0)**self.spec_index
        ## convert Jy to mK
        flux_mK=beam.Jy2mK(flux_Jy,nu)
        return flux_mK
        
    def skyCoordList_movable(self, time, telescope):
        if (not self.have_moon) and (not self.have_sun):
            return None
        ra=[]
        dec=[]
        if (self.have_sun):
            sun=apc.get_sun(time)
            ra.append(sun.ra)
            dec.append(sun.dec)
        if (self.have_moon):
            moon=apc.get_moon(time,telescope.location)
            ra.append(moon.ra)
            dec.append(moon.dec)
        if len(ra)==1:
            ra=ra[0]
            dec=dec[0]
        return apc.SkyCoord(ra=ra, dec=dec)

    def fluxes_mK_movable(self,nu, beam):
        flux_Jy=[]
        if self.have_sun:
            flux_Jy.append(1e6*(nu/1000.)**2)
        if self.have_moon:
            flux_Jy.append(500.*(nu/1000.)**2)
        flux_mK=beam.Jy2mK(np.array(flux_Jy),nu)
        return flux_mK
