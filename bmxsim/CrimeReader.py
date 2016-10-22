#
# Class for reading CRIME maps
#
import healpy as hp
import os,sys
import numpy as np

class CrimeReader:
    def __init__ (self):
        try:
            self.dirname=os.environ['CRIME_DATA']
        except:
            print "Set environment CRIME_DATA to point to CRIME root"
            sys.exit(1)

        self.z_tab = np.loadtxt(self.dirname + "/nuTable.txt")
        self.n_slices=len(self.z_tab)
        self.z_mean= 0.5 * (self.z_tab[:,3] + self.z_tab[:,4])
        self.freq= 0.5 * (self.z_tab[:,1] + self.z_tab[:,2])
        print "Have ",self.n_slices,"slices."
        print "z min, max:",self.z_tab[-1,3],self.z_tab[0,4]

        self.cosmo_slices=[None]*self.n_slices
        self.egfree_slices=[None]*self.n_slices
        self.gfree_slices=[None]*self.n_slices
        self.gsync_slices=[None]*self.n_slices
        self.ptso_slices=[None]*self.n_slices

        self.cosmo_prefix='cosmo_'
        self.gfree_prefix='gfree_'
        self.egfree_prefix='egfree_'
        self.gsync_prefix='gsync_'
        self.ptso_prefix='psources_'

        
    def setOxfordPrefixes():
    # sets prefixes to oxford convention
        self.cosmo_prefix='cosmo/sim_2048_'
        self.gfree_prefix='gfree/gfree_'
        self.egfree_prefix='egfree/egfree_'
        self.gsync_prefix='gsync/gsync_'
        self.ptso_prefix='psources/psources_'
        
    def named_slice(self,names,i):
        toret=None
        for name in names.split("+"):
            if name=='cosmo':
                sli=self.cosmo_slice(i)
            elif name=='egfree':
                sli=self.egfree_slice(i)
            elif name=='gfree':
                sli=self.gfree_slice(i)
            elif name=='gsync':
                sli=self.gsync_slice(i)
            elif name=='ptso':
                sli=self.ptso_slice(i)
            else:
                print "bad name slice in named_slice"
                raise NotImplemented
            if toret is None:
                toret=sli
            else:
                toret+=sli
        return toret

    def cosmo_slice(self,i):
        if (self.cosmo_slices[i]==None):
            fname=self.dirname+"/"+self.cosmo_prefix+"%03d.fits" % (i+1)
            print "Reading ",fname
            self.cosmo_slices[i],di= hp.read_map(fname, h=True)
            self.cosmo_d=dict(di)
        return self.cosmo_slices[i]

    def egfree_slice(self,i):
        if (self.egfree_slices[i]==None):
            fname=self.dirname+"/"+self.egfree_prefix+"%03d.fits" % (i+1)
            print "Reading ",fname
            self.egfree_slices[i],di= hp.read_map(fname, h=True)
            self.egfree_d=dict(di)
        return self.egfree_slices[i]

    def gfree_slice(self,i):
        if (self.gfree_slices[i]==None):
            fname=self.dirname+"/"+self.gfree_prefix+"%03d.fits" % (i+1)
            print "Reading ",fname
            self.gfree_slices[i],di= hp.read_map(fname, h=True)
            self.gfree_d=dict(di)
        return self.gfree_slices[i]

    def gsync_slice(self,i):
        if (self.gsync_slices[i]==None):
            fname=self.dirname+"/"+self.gsync_prefix+"%03d.fits" % (i+1)
            print "Reading ",fname
            self.gsync_slices[i],di= hp.read_map(fname, h=True)
            self.gsync_d=dict(di)
        return self.gsync_slices[i]

    def ptso_slice(self,i):
        if (self.ptso_slices[i]==None):
            fname=self.dirname+"/"+self.ptso_prefix+"%03d.fits" % (i+1)

            print "Reading ",fname
            self.ptso_slices[i],di= hp.read_map(fname, h=True)
            self.ptso_d=dict(di)
        return self.ptso_slices[i]

    def plot_slice(self, slice_s, N, reso, rot=(10,0,0)):
        pt=slice_s
        hp.visufunc.gnomview(map = pt, xsize = N, ysize = N,
                             coord = 'G',rot = (10,0,0), reso = reso)
    
