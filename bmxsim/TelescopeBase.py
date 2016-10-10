#
# Base class for telescopes
#

class TelescopeBase(object):
    """
     This is a base abstract class for a radio telescope.
     For the time being, we assume the telescope to contain a fixed number of dishes with fixed alt, ez,
     but 

    """
    def __init__(self, beams, numin, numax, location=None, name="telescope"):
        """ Constructor:
        beams -- list of BeamBase derived objects representing beams of this telescope
        numin,numax -- min and max frequency of the telescope in MHz
        location -- telescope location, astropy Earth Location 
        """
        self.beams=beams
        self.numin=float(numin)
        self.numax=float(numax)
        self.location=location
        self.name=name

    def setNuParams(self, numin, numax):
        self.numin=numin
        self.numax=numax
        
    
