import re
import sys
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from readcol import *
from bp_binning import *

def read_xfaster_newdata_output( file, pol=False, win=1, old=False, otit="XF T-P spectra", res=False, xpol=False, lmax=2000, chkbb=False, verbose=False ):
    code = ' > read_xfaster_newdata_output: '
    xfdir = '/global/scratch2/sd/dpietrob/Software/XFaster/'
    if verbose:
        print code+'XFaster folder = '+xfdir
    tcl = hp.fitsfunc.read_cl(xfdir+'data/planck_lcdm_cl_uK_xf1.e-3.fits')
    tcl[4] = 0.
    tcl[5] = 0.
    nl = len( tcl[0] )
    l = np.arange( len( tcl[0] ) )
    ll = l * (l+1)/2./np.pi
### NB consistent with my definition of .newdat
    if old:
        nbtt,nbee,nbbb,nbtb,nbte,nbeb = readcol(file, format=np.arange(6)+1, skipline=1, numline=1, verbose=verbose)
    else:
        nbtt,nbee,nbbb,nbte,nbtb,nbeb = readcol(file, format=np.arange(6)+1, skipline=1, numline=1, verbose=verbose)
    
    print ' TT # bin = ', nbtt
    print ' EE # bin = ', nbee
    print ' TE # bin = ', nbte
    print ' BB # bin = ', nbbb
    print ' TB # bin = ', nbtb
    print ' EB # bin = ', nbeb

    #if verbose:
    #print code+'binning theoretical power spectrum...'
    btcltt = bp_binning( tcl[0]*ll, xfdir+'data/bins/ctp/CTP_bin_TT', verbose=verbose )
#    btclee = bp_binning( tcl[1]*ll, xfdir+'data/bins/ctp/CTP_bin_EE' )
#    btclte = bp_binning( tcl[3]*ll, xfdir+'data/bins/ctp/CTP_bin_TE' )
#    btclbb = bp_binning( tcl[2]*ll, xfdir+'data/bins/ctp/CTP_bin_EE' )
#    btcltb = bp_binning( tcl[4]*ll, xfdir+'data/bins/ctp/CTP_bin_TE' )
#    btcleb = bp_binning( tcl[5]*ll, xfdir+'data/bins/ctp/CTP_bin_EE' )

"""pro read_xfaster_newdata_output, file, init=init, pol=pol, win=win, old=old,otit=otit, res=res, xpol=xpol, lmax=lmax, chkbb=chkbb

    xfdir = '/global/scratch2/sd/dpietrob/Software/XFaster/'
;##    fits2cl, tcl, xfdir+'data/lcdm_cl_uK.fits'
    fits2cl, tcl, xfdir+'data/planck_lcdm_cl_uK_xf1.e-3.fits'
    tcl[*,4] = 0.
    tcl[*,5] = 0.
    l=findgen(n_elements(tcl[*,0]))
    ll=l*(l+1)/2./!pi

;## NB consistent with my definition of .newdat
    if keyword_set(old) then readcol,file,nbtt,nbee,nbbb,nbtb,nbte,nbeb, format='f,f,f,f,f,f',skipline=1,numline=1 else $
      readcol,file,nbtt,nbee,nbbb,nbte,nbtb,nbeb, format='f,f,f,f,f,f',skipline=1,numline=1


;    print, 6+2*nbtt+2*nbee+2*nbbb

; ------- Binning
;    if not keyword_set(bincl) then begin
;    endif else begin
;        btcltt = bp_binning(tcl[*,0], xfdir+'data/bins/ctp/CTP_bin_TT') * bp_binning(ll, xfdir+'data/bins/ctp/CTP_bin_TT')
;        btclee = bp_binning(tcl[*,1], xfdir+'data/bins/ctp/CTP_bin_EE') * bp_binning(ll, xfdir+'data/bins/ctp/CTP_bin_EE')
;        btclte = bp_binning(tcl[*,3], xfdir+'data/bins/ctp/CTP_bin_TE') * bp_binning(ll, xfdir+'data/bins/ctp/CTP_bin_TE')        
;    endelse

    readcol,file,cltt,cltter,lmntt,lmxtt,format='x,f,f,x,x,f,f', skipline=6, numline=nbtt[0]

    if (not keyword_set(pol)) then begin
        il = findgen(max(lmxtt))+2
        !p.multi = 0
        window, win, xsize=720, ysize=450,tit=otit
        plot, l[il], tcl[il,0]*ll[il], chars=1.5, xtit='!6l', ytit='!6D!dl!uTT!n [!7l!6K!u2!n]', xr=[1,lmax]
        oplot, (lmntt+lmxtt)/2,cltt, psym=4
        errplot, (lmntt+lmxtt)/2,cltt-cltter, cltt+cltter
        oplot, l[il], tcl[il,0]*ll[il], col=245
        xyouts, max(il), 5500, 'TT Spectrum', chars=1.5

    endif else begin
        il = findgen(max(lmxtt))+2
        !p.multi = [0,3,1]
        if (keyword_set(res)) then !p.multi=[0,3,2]
        window, win, xsize=720*1.8, ysize=450*1.25,tit=otit
        print, lmntt[0], lmxtt[0]
        plot, l[il], tcl[il,0]*ll[il], chars=3, xtit='!6l', ytit='!6D!dl!uTT!n [!7l!6K!u2!n]', xr=[1,lmax]
        oplot, (lmntt+lmxtt)/2,cltt, psym=4
        errplot, (lmntt+lmxtt)/2,cltt-cltter, cltt+cltter
        oplot, l[il], tcl[il,0]*ll[il], col=245

        readcol,file,clee,cleeer,lmnee,lmxee,format='x,f,f,x,x,f,f', skipline=2*nbtt[0]+6+1, numline=nbee[0]
        print, lmnee[0], lmxee[0]
        il = findgen(max(lmxee))+2
        ll=l*(l+1)/2./!pi
        plot, l[il], tcl[il,1]*ll[il], chars=3, xtit='!6l', ytit='!6D!dl!uEE!n [!7l!6K!u2!n]', xr=[1,lmax]
        oplot, (lmnee+lmxee)/2,clee, psym=4
        errplot, (lmnee+lmxee)/2,clee-cleeer, clee+cleeer
        oplot, l[il], tcl[il,1]*ll[il], col=245

        readcol,file,clte,clteer,lmnte,lmxte,format='x,f,f,x,x,f,f', skipline=6+2*nbtt[0]+2*nbee[0]+2*nbbb[0]+1+1+1, numline=nbte[0]
        print, lmnte[0], lmxte[0]
        il = findgen(max(lmxte))+2
        ll=l*(l+1)/2./!pi
        plot, l[il], tcl[il,3]*ll[il], chars=3, xtit='!6l', ytit='!6D!dl!uTE!n [!7l!6K!u2!n]', xr=[1,lmax]
        oplot, (lmnte+lmxte)/2,clte, psym=4
        errplot, (lmnte+lmxte)/2,clte-clteer, clte+clteer
        oplot, l[il], tcl[il,3]*ll[il], col=245
; ------ Residuals
        if (keyword_set(res)) then begin
            plot, (lmntt+lmxtt)/2,cltt-btcltt, psym=4, xtit='!6l', ytit='!6Residuals [!7l!6K!u2!n]', chars=3, yr=[-250,250], xr=[1,lmax]
            errplot, (lmntt+lmxtt)/2, cltt-btcltt-cltter, cltt-btcltt+cltter
            oplot, (lmntt+lmxtt)/2,cltt*0., thick=0.5, col=245

            plot, (lmnee+lmxee)/2,clee-btclee, psym=4, xtit='!6l', ytit='!6Residuals [!7l!6K!u2!n]', chars=3, xr=[1,lmax], yr=[-5,5] 
            errplot, (lmnee+lmxee)/2, clee-btclee-cleeer, clee-btclee+cleeer
            oplot, (lmnee+lmxee)/2,clee*0., thick=0.5, col=245

            plot, (lmnte+lmxte)/2,clte-btclte, psym=4, xtit='!6l', ytit='!6Residuals [!7l!6K!u2!n]', chars=3, xr=[1,lmax], yr=[-15,15] 
            errplot, (lmnte+lmxte)/2, clte-btclte-clteer, clte-btclte+clteer
            oplot, (lmnte+lmxte)/2,clte*0., thick=0.5, col=245
        endif

    endelse

    if (keyword_set(xpol)) then begin
        il = findgen(max(lmxtt))+2
        !p.multi = [0,3,2]
        window, win, xsize=720*1.8, ysize=450*1.25,tit=otit
; --- TT
        print, lmntt[0],lmxtt[0]
        plot, l[il], tcl[il,0]*ll[il], chars=3, xtit='!6l', ytit='!6D!dl!uTT!n [!7l!6K!u2!n]', xr=[1,lmax]
        oplot, (lmntt+lmxtt)/2,cltt, psym=4
        errplot, (lmntt+lmxtt)/2,cltt-cltter, cltt+cltter
        oplot, l[il], tcl[il,0]*ll[il], col=245
; --- EE
        readcol,file,clee,cleeer,lmnee,lmxee,format='x,f,f,x,x,f,f', skipline=2*nbtt[0]+6+1, numline=nbee[0]
        print, lmnee[0],lmxee[0]
        il = findgen(max(lmxee))+2
        ll=l*(l+1)/2./!pi
        plot, l[il], tcl[il,1]*ll[il], chars=3, xtit='!6l', ytit='!6D!dl!uEE!n [!7l!6K!u2!n]', xr=[1,lmax]
        oplot, (lmnee+lmxee)/2,clee, psym=4
        errplot, (lmnee+lmxee)/2,clee-cleeer, clee+cleeer
        oplot, l[il], tcl[il,1]*ll[il], col=245

; --- BB
        readcol,file,clbb,clbber,lmnbb,lmxbb,format='x,f,f,x,x,f,f', skipline=2*nbtt[0]+2*nbee[0]+6+1+1, numline=nbbb[0]
        print, lmnbb[0],lmxbb[0]
        if (keyword_set(chkbb)) then begin
            oplot, (lmnee+lmxee)/2,clee-clbb, col=90, psym=6
            legend, ['C!dl!uEE!n-C!dl!uBB!n'], psym=6, col=90, /top, /left, chars=1.8
        endif
        il = findgen(max(lmxbb))+2
        ll=l*(l+1)/2./!pi
        plot, l[il], tcl[il,2]*ll[il], chars=3, xtit='!6l', ytit='!6D!dl!uBB!n [!7l!6K!u2!n]', xr=[1,lmax], yr=[0,max(tcl[where(il lt lmax),1]*ll[where(il lt lmax)])]
        oplot, (lmnbb+lmxbb)/2,clbb, psym=4
        errplot, (lmnbb+lmxbb)/2,clbb-clbber, clbb+clbber
        oplot, l[il], tcl[il,2]*ll[il], col=245
        if (keyword_set(chkbb)) then begin
            oplot, (lmnee+lmxee)/2,clee-btclee, col=90, psym=6
            legend, ['C!dl!uEE!n-C!dl!uEE,Th!n'], psym=6, col=90, /top, /right, chars=1.8
        endif
; --- TE
        readcol,file,clte,clteer,lmnte,lmxte,format='x,f,f,x,x,f,f', skipline=6+2*nbtt[0]+2*nbee[0]+2*nbbb[0]+1+1+1, numline=nbte[0]
        print, lmnte[0],lmxte[0]
        il = findgen(max(lmxte))+2
        ll=l*(l+1)/2./!pi
        plot, l[il], tcl[il,3]*ll[il], chars=3, xtit='!6l', ytit='!6D!dl!uTE!n [!7l!6K!u2!n]', xr=[1,lmax]
        oplot, (lmnte+lmxte)/2,clte, psym=4
        errplot, (lmnte+lmxte)/2,clte-clteer, clte+clteer
        oplot, l[il], tcl[il,3]*ll[il], col=245
; --- TB
        readcol,file,cltb,cltber,lmntb,lmxtb,format='x,f,f,x,x,f,f', skipline=6+2*nbtt[0]+2*nbee[0]+2*nbbb[0]+2*nbte[0]+1+1+1+1, numline=nbtb[0]
        print, lmntb[0],lmxtb[0]
        il = findgen(max(lmxtb))+2
        ll=l*(l+1)/2./!pi
        plot, l[il], tcl[il,4]*ll[il], chars=3, xtit='!6l', ytit='!6D!dl!uTB!n [!7l!6K!u2!n]', xr=[1,lmax],yr=[-50,50]
        oplot, (lmntb+lmxtb)/2,cltb, psym=4
        errplot, (lmntb+lmxtb)/2,cltb-cltber, cltb+cltber
        oplot, l[il], tcl[il,4]*ll[il], col=245
; --- EB
        readcol,file,cleb,cleber,lmneb,lmxeb,format='x,f,f,x,x,f,f', skipline=6+2*nbtt[0]+2*nbee[0]+2*nbbb[0]+2*nbte[0]+2*nbtb[0]+1+1+1+1+1, numline=nbeb[0]
        print, lmneb[0],lmxeb[0]
        il = findgen(max(lmxeb))+2
        ll=l*(l+1)/2./!pi
        plot, l[il], tcl[il,5]*ll[il], chars=3, xtit='!6l', ytit='!6D!dl!uEB!n [!7l!6K!u2!n]', xr=[1,lmax], yr=[-1,1]
        oplot, (lmneb+lmxeb)/2,cleb, psym=4
        errplot, (lmneb+lmxeb)/2,cleb-cleber, cleb+cleber
        oplot, l[il], tcl[il,5]*ll[il], col=245
    endif

!p.multi = 0
stop
end"""
