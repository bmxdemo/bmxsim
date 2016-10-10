#
# BMX Demo Telescopes
#
from TelescopeBase import TelescopeBase
from BeamAiry import BeamAiry
from astropy.coordinates import EarthLocation
import astropy.units as u
class BMXDemoSingleDish (TelescopeBase):
    """
    A single dish in the basin
    """
    def __init__(self, airy=True):
        """ Constructor
            if airy==True, use airy disk approximation to the system
        """
        # 4 meter dish pointing at zenith, fully illuminated at 800MHz
        beams=[BeamAiry(0,0,90.,0,4.,800.)]
        TelescopeBase.__init__(self, beams, 800., 1400.,
                               EarthLocation(lat=40.87792*u.deg, lon=-72.85852*u.deg, height=0*u.m),
                               "BMXDemoSingleDish")
                               
    def BeamAltAz(self, time):
        """ return AltAz object given time """
        return [AltAz(obstime=time,location=self.location(),alt=b.alt, az=b.az)
                for b in self.beams]

