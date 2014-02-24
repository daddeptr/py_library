def make_sky_cut(gal_cut, ns):
    import numpy as np
    import healpy as hp
    npix = 12 * ns^2
    cut = np.ones( npix )
    ipix = np.arange( npix, dtype=np.long )
    theta,phi = hp.pix2ang( ns, ipix )
    radeg = 180./np.pi
    theta_deg = ( np.pi/2.-theta ) * radeg
    phi_deg = phi * radeg
    cut[ np.logical_and( theta>gal_cut, theta<-gal_cut )] = 0.
    return cut
