#
# Base class for a radio beam
#
import numpy as np
from astropy.coordinates import AltAz
import astropy.units as u

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

    def fwhm(self,nu):
        """ returns beam FWHM in radians """
        raise NotImplemented

    def beamImage (self, N, reso, nu):
        """ returns beam image.
        N - number of pixels in x,y directions
        reso - pixel size in radians
        nu - frequency in MHz

        returns NxN np array containing beam image.
        """
        toret=np.zeros((N,N))
        for i in range (N):
            for j in range(N):
                b=self.beam(reso*(i-(N-1)/2.), reso*(j-(N-1)/2.),nu)
                if (np.isnan(b)):
                    print "Beam is NaN. Not good!"
                    print i,j, reso*(i-(N-1)/2.), reso*(j-(N-1)/2.),nu
                    stop()
                toret[i,j]=b

        return toret
            
    def AltAz(self, time, location):
        """Returns the astropy AltAz object for
        time : astropy Time object
        location : astropy Location object
        returns: astropy AltAz object
        """
        return AltAz(obstime=time,location=location, alt=self.alt*u.radian, az=self.az*u.radian)

# add more functions to deal with x,t,alt,az manipulations as required

        
