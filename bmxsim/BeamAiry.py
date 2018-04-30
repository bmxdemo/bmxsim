#
# Idealised airy beam
#
import numpy as np
import scipy.constants as const
from scipy.special import jn
from BeamBase import BeamBase

class BeamAiry(BeamBase):
    """
    Implemenetation of Beam for which the beam is the basic airy disc 
    """
    def __init__ (self, x,y, alt, az, D, nu_f):
        """ Constructor
        x,y - position in meters on a gird
        alt,az - fixed pointing in degrees
        D - mirror diameter (in meters)
        nu_f -- full illuminator frequency (in MHz)
        """
        BeamBase.__init__(self,x,y,alt,az)
        self.D=float(D)
        self.nu_f=float(nu_f)
        self.lambdamax=const.c/(1e6*self.nu_f)

    def DoverLam(self,nu):
        """ returns effective D/L"""
        if (nu>self.nu_f):
            lamb=self.lambdamax
        else:
            lamb=const.c/(1e6*nu)
        return self.D/lamb
    
    def beam(self, delta_phi, delta_theta, nu):
        """ returns beam delta_phi and delta_theta from the
            boresight -- these are really coordinates in gnomonomic projection, 
            and not delta_ra, delta_dec
            delta_phi, delta_theta in *radians*, can be scalar or equal sized arrays
            nu -- observing frequency in MHz, must be scalar
        """

        # Convert to arrays
        delta_phi = np.atleast_1d(delta_phi)
        delta_theta = np.atleast_1d(delta_theta)

        ## note flat sky approx!! FIX!!!
        theta=np.sqrt(delta_phi**2+delta_theta**2) 
        x=np.pi*self.DoverLam(nu)*theta
        beam=(2*jn(1,x))**2/x**2
        beam[np.where(x==0)]=1.0

        return beam

    def hpbw(self,nu):
        """ returns Half Power Beam Width in radians 
            see e.g https://arxiv.org/pdf/1609.09376.pdf"""
        return 1.028/self.DoverLam(nu)

    def fwhm(self,nu):
        """ returns full width half maximum in radians
        """
        return 1.41/self.DoverLam(nu)


