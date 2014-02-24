import numpy as np
from readcol import *

def bp_binning( cls, bin_file, verbose=False):
    code = ' > bp_binning: '
    if verbose:
        print code+'Assuming l=0, len(cl)-1; lmax=len(cl)-1'
        print code+'Reading bin_file %s' %bin_file
    nl = len(cls)
    fl,ll = readcol(bin_file, verbose=verbose)
    nb = len(ll)
    if verbose:
        print nb, ll[nb-1], ll[nb]
    lmax = np.min( [nl-1,ll[nb]] )
    if verbose:
        print 'lmax = %s' %lmax
    #bok = (ll - lmax) > 0
    #bmax = nb-1
    #if bok[0] > 0:
    #    bmax = bok[0]-1
    bin_cls = np.zeros( nb )
# for ib=0,bmax do bin_cls[ib] = total( cls[fl[ib]:ll[ib]] ) / (ll[ib]-fl[ib]+1)
    for ib in np.arange( nb ):
        #bin_cls[ib] = np.sum( cls[fl[ib]:ll[ib]] ) / (ll[ib]-fl[ib]+1)
        bnl = np.sum( np.nonzero( cls[fl[ib]:ll[ib]] ) )
        bin_cls[ib] = np.sum( cls[fl[ib]:ll[ib]] ) / bnl #(ll[ib]-fl[ib]+1)
        if verbose:
            print ib, fl[ib], ll[ib], bin_cls[ib]
    tnl = np.sum( np.nonzero( bin_cls ) )
    return bin_cls[0:tnl]
