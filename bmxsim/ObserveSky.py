#
# Given a telescope, simulate observations
#

import numpy as np
from astropy.time import Time, TimeMJD
from astropy import units as u
import astropy.coordinates as apc
import healpy as hp
import astropy.coordinates.sky_coordinate as X
import sys
from matplotlib.pyplot import *
import bmxsim as bs

tlist_cache = None
skyc_cache = None

def getIntegratedSignal(telescope, tlist, sigslice, nu, Npix=201, Nfwhm=3):
    """ Integrate the smooth component of the signal:
        telescope : telescope object to simulate
        tlist : list of Time object where you want signal integrated
        sigslice : CrimeReader slice
        nu : frequency
        Npix : number of pixels in the beam image to integrate over
        Nfwhm : total size of beam image to use
        returns: A list of np arrays of length of tslice corresponding to
                power measured at times tlist. If single beam (beams = integer):
                return that particular np array.
    """
    global tlist_cache, skyc_cache
    ion()
    beams=telescope.beams
    toret=[]
    for beam in beams:
        reso=Nfwhm*beam.fwhm(nu)/Npix
        beam_img=beam.beamImage(Npix, reso, nu)
        beam_img/=beam_img.sum()
        Nside=int(np.sqrt(len(sigslice)/12))
        ## this defines a lambda vec2pix function to feed to projmap later
        vec2pix=lambda x,y,z:hp.vec2pix(Nside,x,y,z)
        intlist=[]
        if tlist_cache is None or len(tlist_cache) != len(tlist) or tlist_cache[0] != tlist[0]:
            tlist_cache = tlist
            skyc_cache = []
            for t in tlist:
                aaz = beam.AltAz(t, telescope.location)
                skyc = aaz.transform_to(apc.ICRS)
                skyc_cache.append(skyc)

        for i,t in enumerate(tlist):
            rot=(skyc_cache[i].ra.deg, skyc_cache[i].dec.deg, 0.)
            proj=hp.projector.GnomonicProj(xsize = Npix, ysize = Npix, rot = rot, reso = reso*180*60/np.pi)
            mp=proj.projmap(sigslice,vec2pix)
            csig=(mp*beam_img).sum()

            # Test plots
            if 0:
                clf()
                subplot(1,2,1)
                imshow(mp); #axis('square')
                subplot(1,2,2)
                imshow(beam_img); #axis('square')
                show()
                pause(0.01)

            print i,csig,'\r',
            sys.stdout.flush()
            intlist.append(np.array(csig))
        toret.append(intlist)
    return toret


def getPointSourceSignal(telescope, tlist, sources, nu):
    """ Integrate the point-source component of the signal:
        telescope : telescope object to simulate
        tlist : list of Time object where you want signal integrated
        sourcelist : an object of type PointSourceList
        nu : frequency
        returns: A list of np arrays of length of tslice corresponding to
                power measured at times tlist. If single beam (beams = integer):
                return that particular np array.
    """
    beams=telescope.beams
    toret=[]
    srcf=sources.skyCoordList_fixed()
    for beam in beams:
        li=np.zeros(len(tlist))
        for i,t in enumerate(tlist):
            aaz=beam.AltAz(t, telescope.location)
            skyc=apc.SkyCoord(aaz).transform_to(apc.ICRS)
            ## there is a stupid bug in astropy -- this seems to be the only way
            ## to make it "forget" its locationa, otherwise spherical_offsets_to
            ## fails claiming frame is wrong.
            skyc=apc.SkyCoord(ra=skyc.ra, dec=skyc.dec)
            sig=0.0
            #
            # For efficiency, let's do fixed and movable separately
            #
            # first fixed sources
            sig=0
            if srcf is not None:
                ofs=skyc.spherical_offsets_to(srcf)
                beamsup=beam.beam(ofs[0].rad,ofs[1].rad,nu)
                # total flux is sum over sources of flux * beam supression
                sig+=(beamsup*sources.fluxes_mK_fixed(nu, beam)).sum()

            srcm=sources.skyCoordList_movable(t,telescope)
            if srcm is not None:
                # next movable sources (moon, sun)
                ofs=skyc.spherical_offsets_to(srcm)
                beamsup=beam.beam(ofs[0].rad,ofs[1].rad,nu)
                # total flux is sum over sources of flux * beam supression
                sig+=(beamsup*sources.fluxes_mK_movable(nu,beam)).sum()
            li[i]=sig
            print i,sig,'\r',
            sys.stdout.flush()
        toret.append(li)
    return toret

