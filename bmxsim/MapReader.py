#
# Class for reading input maps of all sorts
#
import healpy as hp
import os,sys
import numpy as np

class MapReader:

    def __init__ (self, whichfields, freq, sn='000'):

        self.fields = whichfields.split('/')

        if type(sn) is not list:
            sn = [sn]*len(self.fields)
        self.sn = {}
        for k,val in enumerate(self.fields):
            self.sn[val] = sn[k]
            
        # Get map locations
        self.getmapinfo(freq)
        
    def getmapinfo(self, freq):
        """Get map directories, prefixes, frequencies, and map index
        numbers (zero indexed, so have to add 1 if one indexed in the
        filenames). This will return only maps in the min/max frequency range
        provided in freq, with a buffer of 1 index on either side to help later
        interpolation. HI4PI frequencies, which are very finely spaced, will be
        a list of lists of frequencies. Sub-lists will be coadded prior to beam
        convolution. """ 
        self.mapdir = {}
        self.prefix = {}
        self.freq = {}
        self.mapind = {}
        self.maptype = {}
        for k,val in enumerate(self.fields):
            # Grab first field. This requires that fields concatenated with "+"
            # all be the same type, e.g. CRIME, Colore, HI4PI, and have same
            # serial number.
            x = val.split('+')
            if x[0] in ['cosmo','egfree','gfree','gsync','psources']:
                # CRIME
                self.mapdir[val] = 'maps/CRIME/{:s}/'.format(self.sn[val])
                self.prefix[val] = [k+'_' for k in x]
                self.maptype[val] = 'crime'
                self.freq[val], self.mapind[val] = self.getfreq(val, freq)
            elif x[0] in ['colore']:
                # COLORE
                self.mapdir[val] = 'maps/colore/{:s}/'.format(self.sn[val])
                self.prefix[val] = ['colore_imap_s1_nu']
                self.maptype[val] = 'colore'
                self.freq[val], self.mapind[val] = self.getfreq(val, freq)
            elif x[0] in ['hi4pi']:
                # HI4PI
                self.mapdir[val] = 'maps/HI4PI/stitched/'
                self.prefix[val] = ['HI_mK_f']
                self.maptype[val] = 'hi4pi'
                self.freq[val], self.mapind[val] = self.getfreq(val, freq)
            else:
                raise ValueError('{:s} not a recognized field'.format(val))

    def getfreq(self, fld, freq):
        """Get frequencies from CRIME nuTable.txt"""
        
        # Get frequency list
        if self.maptype[fld] == 'crime':
            x = np.loadtxt(self.mapdir[fld] + 'nuTable.txt')
            f = 0.5 * (x[:,1] + x[:,2])
        elif self.maptype[fld] == 'colore':
            x = np.loadtxt(self.mapdir[fld] + 'colore_nuTable.txt')
            f = 0.5 * (x[:,0] + x[:,1])
        elif self.maptype[fld] == 'hi4pi':
            f = np.loadtxt(self.mapdir[fld] + 'nuTable.txt')            

        # Cut down to only requested frequencies
        df = freq[1] - freq[0]
        ind = np.where( (f>=freq.min()-df/2) & (f<=freq.max()+df/2) )[0]

        # Error if there are no maps in the frequency range. Could just return
        # zeros, but better to make the user eliminate the problematic map
        # type. 
        if len(ind) == 0:
            raise ValueError('{:s} has no maps within the requested frequency range'.format(fld))

        # Put a buffer map on either end to aid in interpolation later
        if ind.min() > 0:
            ind = np.hstack((ind.min()-1,ind))
        if ind.max() < len(f)-1:
            ind = np.hstack((ind, ind.max()+1))

        # Only return needed frequencies
        f = f[ind]

        # Make a list of lists over which to coadd input maps prior to beam
        # convolution 
        if self.maptype[fld] == 'hi4pi':
            # Construct bin edges
            be_lo = freq - df/2
            be_hi = freq + df/2

            ff = []
            indd = []
            for k in range(len(freq)):
                ind2 = np.where( (f>=be_lo[k]) & (f<be_hi[k]) )[0]
                if len(ind2) > 0:
                    ff.append(f[ind2])
                    indd.append(ind[ind2])
            f = ff
            ind = indd

        return f, ind


    def gethmap(self, fld, i):
        """Read field by signal type and slice index. 
           fld = 'cosmo', 'egfree', 'gfree', 'gsync', 'ptso', 'colore', 'hi4pi'
           i = slice index
           
           Returns healpix map 

           Can combine maps of the same type with '+' between
           e.g. field = self.named_slice('cosmo+gsync', 0)
        """

        mt = self.maptype[fld]

        if mt in ['crime','colore']:
            # 1 indexed
            mi = self.mapind[fld][i] + 1
        elif mt in ['hi4pi']:
            # Zero indexed
            mi = self.mapind[fld][i]

        # If only a single map just load it, otherwise loop over maps and
        # add. This is needed for maps with high frequency resolution like
        # HI4PI.
        if not hasattr(mi,'__iter__'):
            mi = [mi]
        else:
            # This is duplicated in collapsefreq. This call is needed when
            # running in parallel so the proper freq can be passed to
            # getIntegratedSignal. But since it is in parallel, the modificaiton
            # of this attribute does not get passed back to fillStream. 
            self.freq[fld][i] = self.freq[fld][i].mean()

        toret=None
        for j in range(len(fld.split('+'))):
            # Loop over foreground and signal components
            hmap = None
            for mi0 in mi:
                # Loop over frequencies
                fname = self.mapdir[fld] + self.prefix[fld][j] + '{:03d}.fits'.format(mi0)
                hmap0 = hp.read_map(fname)
                if hmap is None:
                    hmap = hmap0
                    nhits = np.isfinite(hmap0).astype('float')
                else:
                    hmap = np.nansum(np.vstack((hmap,hmap0)),0)
                    nhits += np.isfinite(hmap).astype('float')

            # Average healpix maps if multiple frequencies of same
            # component. This works if units of map are temperature units.
            hmap = hmap / nhits

            # Simply add different components
            if toret is None:
                toret = hmap
            else:
                toret += hmap

        return toret


    def collapsefreq(self):
        """Get mean frequency in bins for fine frequency maps"""
        # Bin center of mean of maps
        for k,fld in enumerate(self.fields):
            for i in range(len(self.freq[fld])):
                if hasattr(self.freq[fld][i], '__iter__'):
                    self.freq[fld][i] = self.freq[fld][i].mean()

