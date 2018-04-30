#
# BMX Demo Telescopes
#
from astropy.coordinates import AltAz, EarthLocation
import astropy.units as u
from TelescopeBase import TelescopeBase
from BeamAiry import BeamAiry
import numpy as np

class BMXDemoSingleDish(TelescopeBase):
    """
    A single dish in the basin
    """
    def __init__(self, airy=True, nu_f=800.0):
        """ Constructor
            if airy==True, use airy disk approximation to the system
        """
        # 3.5 meter dish pointing at zenith 
        beams=[BeamAiry(0, 0, 90., 0, 3.5, nu_f)]
        TelescopeBase.__init__(self, beams, 
                               EarthLocation(lat=40.87792*u.deg, lon=-72.85852*u.deg, height=0*u.m),
                               "BMXDemoSingleDish")

    def BeamAltAz(self, time):
        """ return AltAz object given time """
        return [AltAz(obstime=time,location=self.location(),alt=b.alt, az=b.az)
                for b in self.beams]

