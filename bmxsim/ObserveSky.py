#
# Given a telescope, simulate observations
#

import numpy as np
from astropy.time import Time, TimeMJD
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, ICRS, FK5
from astropy.coordinates import get_sun
import healpy as hp


def getIntegratedSignal (telescope, tlist, sigslice, nu, Npix=201, Nfwhm=3):
    """ Integrate the smooth component of the signal:
        tlist : list of Time object where you want signal integrated
        sigslice : CrimeReader slice
        nu : frequnecy
        Npix : number of pixels in the beam image to integrate over
        Nfwhm : total size of beam image to use
        beams : list of beams to simulate.
        returns: A list of np arrays of length of tslice corresponding to 
                power measured at times tlist. If single beam (beams = integer):
                return that particular np array.
    """

    beams=telescope.beams
    toret=[]
    for beam in beams:
        reso=Nfwhm*beam.fwhm(nu)/Npix
        beam_img=beam.beamImage(Npix, reso, nu)
        beam_img=beam_img**2 ## amplitude to power
        beam_img/=beam_img.sum()
        Nside=int(np.sqrt(len(sigslice)/12))
        ## this defines a lambda vec2pix function to feed to projmap later
        vec2pix=lambda x,y,z:hp.vec2pix(Nside,x,y,z)
        intlist=[]
        for i,t in enumerate(tlist):
            aaz=beam.AltAz(t, telescope.location)
            skyc=aaz.transform_to(FK5)
            rot=(skyc.ra.deg, skyc.dec.deg, 0.)
            proj=hp.projector.GnomonicProj(xsize = Npix, ysize = Npix, rot = rot, reso = reso)
            mp=proj.projmap(sigslice,vec2pix)
            csig=(mp*beam_img).sum()
            print i,csig,'\r',
            intlist.append(np.array(csig))
        toret.append(intlist)
    if (type(beams)==int):
        return toret[beams]
    else:
        return toret


