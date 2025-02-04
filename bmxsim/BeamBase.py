#
# Base class for a radio beam
#
import numpy as np
from astropy.coordinates import AltAz
import astropy.units as u
import scipy.constants as const

class BeamBase(object):
    """
     This is a base abstract class for a radio telescope's beam, a single beam.
     Beam defines position on the grid (for intefereometers) and pointing and
     actual radio beam.
    """

    def __init__ (self, x,y, alt,az):
        """ Constructor
        x,y - position in meters on a gird
        alt,az - fixed pointing in degrees
        """
        self.x=float(x)
        self.y=float(y)
        self.alt=alt*np.pi/180.
        self.az=az*np.pi/180.

    def beam(self, delta_phi, delta_theta, nu):
        """ returns beam delta_phi and delta_theta from the
            boresight -- these are really coordinates in gnomonomic projection, 
            and not delta_ra, delta_dec
            delta_phi, delta_theta in *radians*
            nu -- observing frequency in MHz
        """
        raise NotImplemented

    
    def hpbw(self,nu):
        """ returns Half Power Beam Width in radians """
        raise NotImplemented

    def fwhm(self,nu):
        """ returns full width half maximum in radians
        """
        raise NotImplemented
    
    def Omega(self,nu):
        """ returns beam area in solid angle (radians^2) 
            By default, we implent https://science.nrao.edu/facilities/vla/proposing/TBconv
            pi/(4*ln(2))=1.13309
        """
        return 1.13309*self.hpbw(nu)**2
    
    def Jy2mK(self,flux,nu):
        """ convert Jy to mK into receiver temperature """
        lamb=const.c/(1e6*nu)
        Kelvin=lamb**2/(2*const.k*self.Omega(nu))*flux*1e-26
        return Kelvin*1e3
    
    def beamImage(self, N, reso, nu):
        """ returns beam image.
        N - number of pixels in x,y directions
        reso - pixel size in radians
        nu - frequency in MHz

        returns NxN np array containing beam image.
        """
        dx = np.arange(-(N-1)/2.,(N-1)/2.+1)*reso
        xx,yy = np.meshgrid(dx,dx)
        b = self.beam(xx,yy,nu)
        if np.any(np.isnan(b)):
            print "Beam is NaN. Not good!"
            stop()
        return b
            
    def AltAz(self, time, location):
        """Returns the astropy AltAz object for
        time : astropy Time object
        location : astropy Location object
        returns: astropy AltAz object
        """
        return AltAz(obstime=time,location=location, alt=self.alt*u.radian, az=self.az*u.radian)

# add more functions to deal with x,t,alt,az manipulations as required

        
