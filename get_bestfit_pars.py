import sys
import os
import numpy
import matplotlib.pyplot as plt

def read_pars(root,filename='results/xf_chains_TP.likestats',dir='/global/scratch2/sd/dpietrob/Software/XFaster/outputs/'):
    print dir+root+'/'+filename
    f = open(dir+root+'/'+filename, 'r')
    lines = f.readlines()[3:]
    pars = {}
    for line in lines:
        line = line.strip('\n')
        entries = line.split()
        #print entries
        key = ''.join( entries[6:] )
        #print key
        val = entries[1]
        #print val
        if key not in pars:
            pars[key] = val
    #print pars
    return pars
    
def write_pars(pars,root='',dir='/global/scratch2/sd/dpietrob/Software/XFaster/outputs/'):
    tags = { \
        '\Omega_bh^2': 'ombh2', \
        '\Omega_ch^2': 'omch2', \
        'H_0': 'hubble',        \
        '10^9A_s': 'scalar_amp(1)', \
        'n_s': 'scalar_spectral_index(1)', \
        '\\tau': 're_optical_depth', \
        'A^{tSZ}': 'AtSZ', \
        'A^{kSZ}': 'AkSZ', \
        'A^{CIB}': 'Acib', \
        'A^{PS}_{TT}': 'ApsTT', \
        'A^{PS}_{EE}': 'ApsEE', \
        'A^{PS}_{TE}': 'ApsTE' \
        }

    #print pars['10^9 A_s']
    pars['10^9A_s'] = str( float( pars['10^9A_s'] ) * 1.e-9 )
    #print pars['10^9 A_s']

    file = dir+'/'+root+'/'+'xf_camb.ini'
    print file
    f = open(file, 'w')
    f.write( 'DEFAULT(/global/scratch2/sd/dpietrob/Software/XFaster/src_v1.9/parfiles/camb_par.base)\n' )
    f.write( 'output_root = '+dir+root+'/xf_chains\n' )
    for key in tags:
        print tags[key], pars[key]
        f.write( '%s = %s \n' % (tags[key], pars[key]) )
    f.close()

def run_camb(root, cambdir = '/global/scratch2/sd/dpietrob/Software/cosmomc_cosmoslik/camb/', dir='/global/scratch2/sd/dpietrob/Software/XFaster/outputs/', cambfile='xf_camb.ini'):
    parfile = dir + root + '/' + cambfile
    cmd = cambdir+'camb '+parfile
    os.system(cmd)

def compute_bestfit(root, datadir = '/global/scratch2/sd/dpietrob/Software/cosmomc_cosmoslik/data/fg_temp', dir='/global/scratch2/sd/dpietrob/Software/XFaster/outputs/', cambparfile='xf_camb.ini', \
                    clsfile='xf_chains_lensedCls.dat'):
    ksz_file = datadir + '/' + 'dls_ksz_148_trac.txt'
    tsz_file = datadir + '/' + 'dls_tsz_143_eps0.50.txt'
    ps_file  = datadir + '/' + 'dls_poisson.txt'

    print ksz_file
    ksz = numpy.loadtxt( ksz_file, usecols=([1]) )
    tsz = numpy.loadtxt( tsz_file, usecols=([1]) )
    ps  = numpy.loadtxt( ps_file, usecols=([1])  )

    file = dir + '/' + root + '/' + clsfile  
    print file
    l, tt, ee, bb, te = numpy.loadtxt( file, usecols=(0,1,2,3,4), unpack=True )

    print tt
    
    f = open( dir+root+'/'+cambparfile, 'r' ).readlines()
    fl = f.pop(0)
    pars = {}
    for line in f:
        line = line.strip('\n')
        key, value = line.split('=')
        key = key.strip()
        value = value.strip()
        #print value
        pars[key] = value

    print pars

    cib = ( numpy.arange(4999) +2 ) / 3000.

    lmax = min( [len(tt), len(cib), len(ksz), len(tsz), len(ps)] )
    print lmax

    tot_tt = tt[0:lmax] + \
             ps[0:lmax]  * float( pars['ApsTT']) + \
             cib[0:lmax] * float( pars['Acib'] ) + \
             ksz[0:lmax] * float( pars['AkSZ'] ) + \
             tsz[0:lmax] * float( pars['AtSZ'] )

    tot_ee = ee[0:lmax] + \
             ps[0:lmax]  * float( pars['ApsEE'] )

    tot_te = te[0:lmax] + \
             ps[0:lmax]  * float( pars['ApsTE'] )

### --- Foreground best fit only
    fg_tt =  ps[0:lmax]  * float( pars['ApsTT']) + \
             cib[0:lmax] * float( pars['Acib'] ) + \
             ksz[0:lmax] * float( pars['AkSZ'] ) + \
             tsz[0:lmax] * float( pars['AtSZ'] )

    fg_ee =  ps[0:lmax]  * float( pars['ApsEE'] )

    fg_te =  ps[0:lmax]  * float( pars['ApsTE'] )

    win = 1
    width = 1

    #print tot_tt, tot_ee, tot_te
    x = numpy.transpose( [ l, tot_tt, tot_ee, tot_te] ) 
    file = dir+root+'/'+ 'bestfit.txt'
    numpy.savetxt( file, x , fmt='%12i %15.8f %15.8f %15.8f', newline='\n')
    
    x = numpy.transpose( [ l, fg_tt, fg_ee, fg_te] ) 
    file = dir+root+'/'+ 'fg_bestfit.txt'
    numpy.savetxt( file, x , fmt='%12i %15.8f %15.8f %15.8f', newline='\n')
    
