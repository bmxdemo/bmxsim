import healpy as hp
import numpy as np
from astropy.io import fits
from glob import glob
import gc
from datamanager import datamanager

mapdir = 'maps/HI4PI/'

fn = np.sort(glob(mapdir+'raw/HPX_*.fits'))
nmaps = len(fn)
vrng = [-600,600] # km/s
nside = 1024

# There is not enough memory to open up an Nfreq x Npix healpix map on one
# node. Break it up into several frequencies
Nfblk = 2

class stitch_HI4PI_maps(object):

    def __init__(self):
        """Do nothing"""
        return

    def stitch(self):
        """Stitch maps from raw data"""

        for j in range(Nfblk):

            for k,fn_raw in enumerate(fn):

                # Load map
                print('loading map {0} of {1}'.format(k+1,nmaps))
                x = fits.getdata(fn_raw,1)
                map = x['data']
                pix = x['HPXINDEX']

                if k==0:
                    # Get frequency axis
                    nchan = map.shape[1]
                    v = np.linspace(vrng[0], vrng[1], nchan)
                    f0 = 1420.406 # MHz
                    c = 2.99792458e5 # km/s
                    z = np.sqrt( (1+v/c)/(1-v/c) ) - 1
                    f = f0/(1+z)

                    # Save for saving nutable
                    if j==0:
                        f_master = f*1.0

                    # Now cut down on the channels
                    nfreq = nchan / Nfblk
                    ind0 = nfreq*j
                    if j < (Nfblk-1):
                        ind1 = nfreq*(j+1)
                    else:
                        ind1 = nchan

                    f = f[ind0:ind1]
                    nchan = f.size

                    # Create a bunch of healpix maps
                    npix = hp.nside2npix(nside)
                    hmap = np.zeros((nchan,npix))

                # Put map data in right place
                hmap[:,pix] = (map.T)[ind0:ind1,:]

            # Convert to mK and save maps
            for k,val in enumerate(hmap):
                val = val*1e3 # K -> mK
                fn_out = mapdir + 'stitched/HI_mK_f{:03d}.fits'.format(ind0+k)
                hp.write_map(fn_out,val)

            del hmap
            gc.collect()

        # Save nutable
        np.savetxt(mapdir + 'stitched/nuTable.txt', f_master, delimiter='\n')

        return
        

