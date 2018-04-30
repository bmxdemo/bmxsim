#
# HFSS beam
#
import numpy as np
from BeamBase import BeamBase
from astropy import units as u
from scipy.interpolate import RegularGridInterpolator

replace = np.core.defchararray.replace
find = np.core.defchararray.find

class BeamHFSS(BeamBase):
    """
    Implemenetation of Beam for which the beam is a saved HFSS csv file,
    B(theta,phi,nu)
    """
    def __init__(self, x,y, alt, az, fname):
        """ Constructor
        x,y - position in meters on a grid of the beam coordinate system
              origin. 
        alt,az - fixed pointing in degrees of the north pole (theta=0) of the
                 HFSS coordinate system in which the beam is saved (typically nominal zenith)
        fname - filename of the HFSS csv file
        """
        BeamBase.__init__(self,x,y,alt,az)
        self.fname = fname
        self.readhfss()

    def readhfss(self):
        """Load hfss"""
        x = np.loadtxt(self.fname, delimiter=',', dtype=np.str)
        self.th = x[1:,0].astype(np.float) # theta

        b = replace(replace(x[1:,1:],'i','j'),' ','').astype(np.complex)
        h = x[0,1:] # header information

        # Find stuff
        f = []
        phi = []
        xy = []
        for k,s in enumerate(h):
            xy0, f0, phi0 = self.parseheader(s)
            f.append(f0)
            phi.append(phi0)
            xy.append(xy0)

        f = np.array(f)
        phi = np.array(phi)
        xy = np.array(xy)

        phiu = np.unique(phi)
        fu = np.unique(f)
        xyu = np.unique(xy)

        nfields = len(xyu)
        nth = len(self.th)
        nf = len(fu)
        nphi = len(phiu)

        self.b = np.reshape(b, (nth,nphi,nf,nfields), order='F').T
        self.f = fu
        self.phi = phiu
        self.xy = xyu

        # Setup for interpolation
        self.setupinterp()

    def parseheader(self, s):
        """Parse a header, return X/Y, Freq in MHz, and phi in deg"""

        i = s.find('rEL3X')
        if i>0:
            xy = 'x'
        i = s.find('rEL3Y')
        if i>0:
            xy = 'y'

        i1 = s.find('Freq=')
        i2 = s.find('Phi=')

        if i2>i1:
            f = s[(i1+5):i2]
            phi = s[(i2+4):]
        else:
            f = s[(i1+5):]
            phi = s[(i2+4):i1]
        
        f = f.replace("'","").replace(' ','').replace('"','')
        phi = phi.replace("'","").replace(' ','').replace('"','')

        f = 1.0*u.Unit(f).to(u.MHz)
        phi = 1.0*u.Unit(phi).to(u.deg)

        return xy, f, phi


    def setupinterp(self):
        """Set up the interpolation function"""

        # Initialize interpolation function
        self.bi = []

        # Loop over channels
        for k,val in enumerate(self.b):
            self.bi.append(RegularGridInterpolator((self.f,self.phi,self.th),
                                                   val, bounds_error=False, fill_value=0))
            

    def beam(self, delta_phi, delta_theta, nu):
        """ returns beam delta_phi and delta_theta from the
            boresight -- these are really coordinates in gnomonomic projection, 
            and not delta_ra, delta_dec
            delta_phi, delta_theta in *radians*, can be scalars or equal sized arrays
            nu -- observing frequency in MHz, must be scalar
        """
        
        # Convert to arrays
        delta_phi = np.atleast_1d(delta_phi)
        delta_theta = np.atleast_1d(delta_theta)

        # Convert x/y to polar coordinate. Will want to make this match the
        # gnomoni projection used in the convolution later
        phi = np.arctan2(delta_theta, delta_phi)*180/np.pi
        phi[phi<0] = phi[phi<0] + 360
        phi[phi>360] = phi[phi>360] - 360

        th = np.sqrt(delta_phi**2 + delta_theta**2)*180/np.pi

        # Interpolate
        b_complex = self.bi[0]((nu, phi, th))
        b_real = np.real(b_complex * np.conj(b_complex))
        
        return b_real


    def fwhm(self,nu):
        """ returns full width half maximum in radians
        """
        return 4.0 * np.pi/180


