import re
import sys
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp

# ----------------------------------------------------------------------
#def ctp_binning(cls):
#    intype = type(cls)
#    if intype.upper() == "STRING":
    
# ----------------------------------------------------------------------
def bp_binning( cls, bin_file, verbose=False):
    code = ' > bp_binning: '
    if verbose:
        print code+'Assuming l=0, len(cl)-1; lmax=len(cl)-1'
        print code+'Reading bin_file %s' %bin_file
    nl = len(cls)
    fl,ll = readcol(bin_file)
    #nb = len(ll)
    lmax = np.min( [nl-1, ll[-1]] )
    #print nl-1, ll[-1]
    if verbose:
        print 'lmax = %s' %lmax
    bok = np.sum( ( fl < lmax ) )
    nb = bok
    #print nb
    if verbose:
        print code+'Maximum bin = ', (fl[nb-1],ll[nb-1])
    bin_cls = np.zeros( nb )
# for ib=0,bmax do bin_cls[ib] = total( cls[fl[ib]:ll[ib]] ) / (ll[ib]-fl[ib]+1)
    for ib in np.arange( nb ):
        #bin_cls[ib] = np.sum( cls[fl[ib]:ll[ib]] ) / (ll[ib]-fl[ib]+1)
        if ll[ib] <= lmax:
            bnl = ll[ib]-fl[ib]+1
        else:
            bnl = lmax-fl[ib]+1
        #print fl[ib], ll[ib], cls[fl[ib]:ll[ib]+1], bnl
        bin_cls[ib] = np.sum( cls[fl[ib]:ll[ib]+1] ) / bnl #(ll[ib]-fl[ib]+1)
        if verbose:
            print ib, fl[ib], ll[ib]+1, bin_cls[ib], bnl
    #tnl = np.sum( ( bin_cls > 0.) )
    return bin_cls

# ----------------------------------------------------------------------
def plot_xfaster_newdata_output( file, pol=False, win=1, old=False, otit="XF T-P spectra", res=False, xpol=False, lmax=2000, chkbb=False, verbose=False, eelim=4, ttlim=250, telim=15 ):
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
        nbtt,nbee,nbbb,nbtb,nbte,nbeb = readcol(file, format=np.arange(6)+1, skipline=2, numline=1, verbose=verbose)
    else:
        nbtt,nbee,nbbb,nbte,nbtb,nbeb = readcol(file, format=np.arange(6)+1, skipline=2, numline=1, verbose=verbose)
    
    print ' TT # bin = ', nbtt
    print ' EE # bin = ', nbee
    print ' TE # bin = ', nbte
    print ' BB # bin = ', nbbb
    print ' TB # bin = ', nbtb
    print ' EB # bin = ', nbeb

    #if verbose:
    #print code+'binning theoretical power spectrum...'
    btcltt = bp_binning( tcl[0]*ll, xfdir+'data/bins/ctp/CTP_bin_TT' )
    btclee = bp_binning( tcl[1]*ll, xfdir+'data/bins/ctp/CTP_bin_EE' )
    btclte = bp_binning( tcl[3]*ll, xfdir+'data/bins/ctp/CTP_bin_TE' )
    btclbb = bp_binning( tcl[2]*ll, xfdir+'data/bins/ctp/CTP_bin_EE' )
    btcltb = bp_binning( tcl[4]*ll, xfdir+'data/bins/ctp/CTP_bin_TE' )
    btcleb = bp_binning( tcl[5]*ll, xfdir+'data/bins/ctp/CTP_bin_EE' )

    cltt,cltter,lmntt,lmxtt = readcol( file, format=[2,3,6,7], skipline=7, numline=nbtt, verbose=verbose )
#    print cltt[0:10]

    width = 10
    if (not pol) and (not xpol):
        il = np.arange(np.max(lmxtt),dtype=np.int16) + 2
        #print code+'Here %s' %lmax
        #fig = plt.figure(num=win, figsize=(width,450./720*width) )
        #fig.figsize=(12,12*450/720)
        #print code+'Here'
        if not res:
            fig = plt.figure(num=win, figsize=(width,450./720*width) )
            fig.set_size_inches(width,450./720*width)
            ax = fig.add_subplot(111)
        else:
            fig = plt.figure(num=win, figsize=(width,450./720*width*1.25) )
            fig.set_size_inches(width,450./720*width*1.25)
            ax = fig.add_subplot(211)
        yl = tcl[0]*ll
        ax.plot( l[il], yl[il], "r", label="Planck best fit")
        ax.plot( (lmntt+lmxtt)/2, cltt, 'k.' )
        ax.errorbar( (lmntt+lmxtt)/2, cltt, yerr=cltter, fmt='k.', label="TT spectrum")
        tit=file.split('.newdat')
        tit=tit[0].split('/')
        ax.set_title( tit[-1] )
        ax.set_xlabel('$\ell$')
        ax.set_ylabel('$\ell(\ell+1) C_\ell /2\pi \; [\mu K^2]$')
        ax.set_xlim([1,lmax*1.05])
        ax.set_ylim([-500,6500])
        # legend
        leg = plt.legend(frameon=True)
        # remove box around legend
        leg.get_frame().set_edgecolor("white")
        leg.get_frame().set_alpha(.8)

        if res:
            ax = fig.add_subplot(212)
            ax.plot( l[il], yl[il]*0., "r", label="Planck best fit")
            ax.plot( (lmntt+lmxtt)/2, cltt-btcltt[0:len(cltt)], "k." )
            ax.errorbar( (lmntt+lmxtt)/2, cltt-btcltt[0:len(cltt)], yerr=cltter, fmt='k.', label="TT spectrum")
            ax.set_xlabel('$\ell$')
            ax.set_ylabel('$\ell(\ell+1) (C_\ell-C_\ell^th) /2\pi \; [\mu K^2]$')
            ax.set_xlim([1,lmax*1.05])
            ax.set_ylim([-ttlim,ttlim])

    if pol:
        clee,cleeer,lmnee,lmxee = readcol( file, format=[2,3,6,7], skipline=7+2*nbtt+1, numline=nbee, verbose=verbose )
        clbb,clbber,lmnbb,lmxbb = readcol( file, format=[2,3,6,7], skipline=7+2*nbtt+2*nbee+1+1, numline=nbbb, verbose=verbose )
        clte,clteer,lmnte,lmxte = readcol( file, format=[2,3,6,7], skipline=7+2*nbtt+2*nbee+2*nbbb+1+1+1, numline=nbte, verbose=verbose )
        #print clee[0:10], btclee[0:10]
        #print clte[0:10], btclte[0:10]
        #print clbb[0:10], btclbb[0:10]

        il = np.arange(np.max(lmxtt),dtype=np.int16) + 2
        #fig = plt.figure(win)
        if not res:
            fig = plt.figure(num=win, figsize=(width*1.8,450./720*width) )
            fig.set_size_inches(width*1.8,450./720*width)
            ax = fig.add_subplot(131)
        else:
            fig = plt.figure(num=win, figsize=(width*1.8,450./720*width*1.25) )
            fig.set_size_inches(width*1.8,450./720*width*1.25)
            ax = fig.add_subplot(231)
        tit=file.split('.newdat')
        tit=tit[0].split('/')
        #plt.title( otit )
        yl = tcl[0]*ll
        ax.plot( l[il], yl[il], "r", label="Planck best fit")
        ax.plot( (lmntt+lmxtt)/2, cltt, "k." )
        ax.errorbar( (lmntt+lmxtt)/2, cltt, yerr=cltter, fmt='k.', label="TT spectrum")
        ax.set_xlabel('$\ell$')
        ax.set_ylabel('$\ell(\ell+1) C_\ell^{TT} /2\pi \; [\mu K^2]$')
        ax.set_xlim([1,lmax*1.05])
        ax.set_ylim([-500,6500])
        leg = plt.legend(frameon=True)
        leg.get_frame().set_edgecolor("white")
        leg.get_frame().set_alpha(.8)
        if res:
            ax = fig.add_subplot(234)
            ax.plot( l[il], yl[il]*0., "r", label="Planck best fit")
            ax.plot( (lmntt+lmxtt)/2, cltt-btcltt[0:len(cltt)], "k." )
            ax.errorbar( (lmntt+lmxtt)/2, cltt-btcltt[0:len(cltt)], yerr=cltter, fmt='k.', label="TT spectrum residuals")
            ax.set_xlabel('$\ell$')
            ax.set_ylabel('$\ell(\ell+1) (C_\ell^{TT}-C_\ell^{th}) /2\pi \; [\mu K^2]$')
            ax.set_xlim([1,lmax*1.05])
            ax.set_ylim([-ttlim,ttlim])

# EE
        if not res:
            axee = fig.add_subplot(132)
        else:
            axee = fig.add_subplot(232)
        yl = tcl[1]*ll
        axee.plot( l[il], yl[il], "r", label="Best fit")
        axee.plot( (lmnee+lmxee)/2, clee, "k." )
        axee.errorbar( (lmnee+lmxee)/2, clee, yerr=cleeer, fmt='k.', label="EE spectrum")
        axee.set_xlabel('$\ell$')
        axee.set_ylabel('$\ell(\ell+1) C_\ell^{EE} /2\pi \; [\mu K^2]$')
        axee.set_xlim([1,lmax*1.05])
        eerange = np.max(yl[0:lmax])
        axee.set_ylim( [-0.1*eerange,eerange] )
        leg = plt.legend(frameon=True)
        leg.get_frame().set_edgecolor("white")
        leg.get_frame().set_alpha(.8)
        if res:
            ax = fig.add_subplot(235)
            ax.plot( l[il], yl[il]*0., "r", label="Planck best fit")
            ax.plot( (lmnee+lmxee)/2, clee-btclee[0:len(clee)], "k." )
            ax.errorbar( (lmnee+lmxee)/2, clee-btclee[0:len(clee)], yerr=cleeer, fmt='k.', label="EE spectrum residuals")
            ax.set_xlabel('$\ell$')
            ax.set_ylabel('$\ell(\ell+1) (C_\ell^{EE}-C_\ell^{th}) /2\pi \; [\mu K^2]$')
            ax.set_xlim([1,lmax*1.05])
            ax.set_ylim([-eelim,eelim])

# TE
        if not res:
            axte = fig.add_subplot(133)
        else:
            axte = fig.add_subplot(233)
        yl = tcl[3]*ll
        axte.plot( l[il], yl[il], "r", label="Planck best fit")
        axte.plot( (lmnte+lmxte)/2, clte, "k." )
        axte.errorbar( (lmnte+lmxte)/2, clte, yerr=clteer, fmt='k.', label="TE spectrum")
        axte.set_xlabel('$\ell$')
        axte.set_ylabel('$\ell(\ell+1) C_\ell^{TE} /2\pi \; [\mu K^2]$')
        axte.set_xlim([1,lmax*1.05])
        axte.set_ylim([-150,150])
        leg = plt.legend(frameon=True)
        leg.get_frame().set_edgecolor("white")
        leg.get_frame().set_alpha(.8)
        if res:
            axte = fig.add_subplot(236)
            axte.plot( l[il], yl[il]*0., "r", label="Planck best fit")
            axte.plot( (lmnte+lmxte)/2, clte-btclte[0:len(clte)], "k." )
            axte.errorbar( (lmnte+lmxte)/2, clte-btclte[0:len(clte)], yerr=clteer, fmt='k.', label="TE spectrum residuals")
            axte.set_xlabel('$\ell$')
            axte.set_ylabel('$\ell(\ell+1) (C_\ell^{TE}-C_\ell^{th}) /2\pi \; [\mu K^2]$')
            axte.set_xlim([1,lmax*1.05])
            axte.set_ylim([-telim,telim])



"""pro read_xfaster_newdata_output, file, init=init, pol=pol, win=win, old=old,otit=otit, res=res, xpol=xpol, lmax=lmax, chkbb=chkbb

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

        readcol,file,clte,clteer,lmnte,lmxte,format='x,f,f,x,x,f,f', skipline=7+2*nbtt[0]+2*nbee[0]+2*nbbb[0]+1+1+1, numline=nbte[0]
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
        readcol,file,clte,clteer,lmnte,lmxte,format='x,f,f,x,x,f,f', skipline=7+2*nbtt[0]+2*nbee[0]+2*nbbb[0]+1+1+1, numline=nbte[0]
        print, lmnte[0],lmxte[0]
        il = findgen(max(lmxte))+2
        ll=l*(l+1)/2./!pi
        plot, l[il], tcl[il,3]*ll[il], chars=3, xtit='!6l', ytit='!6D!dl!uTE!n [!7l!6K!u2!n]', xr=[1,lmax]
        oplot, (lmnte+lmxte)/2,clte, psym=4
        errplot, (lmnte+lmxte)/2,clte-clteer, clte+clteer
        oplot, l[il], tcl[il,3]*ll[il], col=245
; --- TB
        readcol,file,cltb,cltber,lmntb,lmxtb,format='x,f,f,x,x,f,f', skipline=7+2*nbtt[0]+2*nbee[0]+2*nbbb[0]+2*nbte[0]+1+1+1+1, numline=nbtb[0]
        print, lmntb[0],lmxtb[0]
        il = findgen(max(lmxtb))+2
        ll=l*(l+1)/2./!pi
        plot, l[il], tcl[il,4]*ll[il], chars=3, xtit='!6l', ytit='!6D!dl!uTB!n [!7l!6K!u2!n]', xr=[1,lmax],yr=[-50,50]
        oplot, (lmntb+lmxtb)/2,cltb, psym=4
        errplot, (lmntb+lmxtb)/2,cltb-cltber, cltb+cltber
        oplot, l[il], tcl[il,4]*ll[il], col=245
; --- EB
        readcol,file,cleb,cleber,lmneb,lmxeb,format='x,f,f,x,x,f,f', skipline=7+2*nbtt[0]+2*nbee[0]+2*nbbb[0]+2*nbte[0]+2*nbtb[0]+1+1+1+1+1, numline=nbeb[0]
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

# ----------------------------------------------------------------------
def readcol(file,skipline=0l,numline=-1l,format=[1,2],verbose=False, twod=False):
    """Routine to emulate IDL readcol.pro"""
    code = ' > readcol: '
    if verbose:
        print ' > Readcol: being talkative ', verbose
    if verbose:
        print ' > Readcol routine. D.Pietrobon Oct 2013'
    if verbose:
        print ' > Readcol: reading file "%s"' %file
    nlines   = int( get_lineNum(file,verbose=verbose) )
    skipline = int( skipline )
    numline  = int( numline )
    #print nlines, skipline, numline
    if verbose:
        print ' > Readcol: total number of lines: %s' %nlines
    try:
        file_content = open(file, 'r')
        
        if skipline > 0:
            if verbose:
                print ' > Readcol: skipping %s lines...' %skipline
        if numline < 0:
            numline = nlines-skipline
        else:
            if numline > (nlines-skipline):
                print ' > Readcol: WARNING, not enough lines. Reading %s lines:' %nlines-skipline
                numline = nlines-skipline
        # Getting format
        ifields = np.asarray( format, np.int16 )
        nfields = len( format )
        if verbose:
            print ' > readcol: number of fields to be read = %s ' %nfields
        if verbose:
            print ' > readcol: fields read = %s' %ifields
        ifields[:] -= 1
        values = np.zeros( (numline, nfields) )
#        for iline in np.arange(np.max( [nlines,numline] ))+skipline:
        cnt = 0
        for iline,line in enumerate( file_content ):
            if (iline >= skipline) and (iline < skipline+numline):
                if verbose:
                    print ' > line %s' %iline
                    print ' >      %s' %line
                line = line.strip('\n')
                entries = np.asarray( line.split() )
                #print ifields
                #print entries[ ifields[0] ]
                #print entries[ ifields ]
                #if verbose:
                #    print entries[ [ifields] ]
                #print entries[ifields]
                entries = np.asarray( entries[ifields] )
                #if verbose:
                #    print entries
                #values.append(entries)
                values[cnt,:] = entries
                if verbose:
                    print values[cnt,:]
                cnt += 1


        file_content.close()
        #print [ values[:,icol] for icol in np.arange(nfields) ]
        if twod:
            if verbose:
                print code+'Exit-'
            return values
        else:
            if verbose:
                print code+'Exit-'
            return [ values[:,icol] for icol in np.arange(nfields) ]

    except IOError:
        print "File not found: %s" %file
        
# ----------------------------------------------------------------------
def get_lineNum(file,verbose=False):
    if verbose:
        print ' > get_lineNum routine. D.Pietrobon Oct 2013'
        nlines = -1
    try:
        f = open(file, 'r')
        nlines=0
        for iline in f:
            nlines+=1
        if verbose:
            print " > get_lineNum: line = %s" %nlines
        return nlines
    except IOError:
        print "File not found: %s" %file

# ----------------------------------------------------------------------
def extract_xfaster_newdata_output(file, ncl='EE', xfdir='/global/scratch2/sd/dpietrob/Software/XFaster/', tclfile='data/planck_lcdm_cl_uK_xf1.e-3.fits', old=False):
    xfdata = {}
    code = 'extract_xfaster_newdata_output: '
    print code + 'WARNING - new xfaster output assumed.'

    tcl = hp.fitsfunc.read_cl( tclfile )
    tcl[4] = 0.
    tcl[5] = 0.
    l = np.arange( len( tcl[0]))
    ll = l*(l+1)/2./np.pi

### NB consistent with my definition of .newdat
    if (ncl == '1') or (ncl.upper() == 'TT'):
        nbtt=readcol(file,format=[1],skipline=2,numline=1)
        print ' TT # bin = ', nbtt
    else:
        if old:
            nbtt,nbee,nbbb,nbtb,nbte,nbeb = readcol(file, format=np.arange(6)+1,skipline=2,numline=1)
        else:
            nbtt,nbee,nbbb,nbte,nbtb,nbeb = readcol(file, format=np.arange(6)+1,skipline=2,numline=1)
        
        print ' EE # bin = ', nbee
        print ' TE # bin = ', nbte
        print ' BB # bin = ', nbbb
        print ' TB # bin = ', nbtb
        print ' EB # bin = ', nbeb

# ------- Binning
    if (ncl == '1') or (ncl.upper() == 'TT'):
        cl,cler,lmn,lmx = readcol(file,format=[2,3,6,7], skipline=7, numline=nbtt[0])
        btcl = bp_binning(tcl[0]*ll, xfdir+'data/bins/ctp/CTP_bin_TT')

    if (ncl == '2') or (ncl.upper() == 'EE'):
        cl,cler,lmn,lmx = readcol(file,format=[2,3,6,7], skipline=2*nbtt[0]+7+1, numline=nbee[0])
        btcl = bp_binning(tcl[1]*ll, xfdir+'data/bins/ctp/CTP_bin_EE')
        
    if (ncl == '3') or (ncl.upper() == 'BB'):
        cl,cler,lmn,lmx = readcol(file,format=[2,3,6,7], skipline=2*nbtt[0]+2*nbee[0]+7+1+1, numline=nbbb[0])
        btcl = bp_binning(tcl[2]*ll, xfdir+'data/bins/ctp/CTP_bin_TE')

    if (ncl == '4') or (ncl.upper() == 'TE'):
        cl,cler,lmn,lmx = readcol(file,format=[2,3,6,7], skipline=7+2*nbtt[0]+2*nbee[0]+2*nbbb[0]+1+1+1, numline=nbte[0])
        btcl = bp_binning(tcl[3]*ll, xfdir+'data/bins/ctp/CTP_bin_EE')

    if (ncl == '5') or (ncl.upper() == 'TB'):
        cl,cler,lmn,lmx = readcol(file,format=[2,3,6,7], skipline=7+2*nbtt[0]+2*nbee[0]+2*nbbb[0]+2*nbte[0]+1+1+1+1, numline=nbtb[0])
        btcl = bp_binning(tcl[4]*ll, xfdir+'data/bins/ctp/CTP_bin_TE')

    if (ncl == '6') or (ncl.upper() == 'EB'):
        cl,cler,lmn,lmx = readcol(file,format=[2,3,6,7], skipline=7+2*nbtt[0]+2*nbee[0]+2*nbbb[0]+2*nbte[0]+2*nbtb[0]+1+1+1+1+1, numline=nbeb[0])
        btcl = bp_binning(tcl[5]*ll, xfdir+'data/bins/ctp/CTP_bin_EE')

    lcen = (lmn+lmx)/2.
    res = cl-btcl

    xfdata.append(lcen)
    xfdata.append(cl)
    xfdata.append(res)
    xfdata.append(btcl)

    return xfdata

# ----------------------------------------------------------------------
def make_sky_cut(gal_cut, ns, verbose=False):
    import numpy as np
    import healpy as hp
    npix = 12 * ns**2
    cut = np.ones( npix )
    ipix = np.arange( npix, dtype=np.long )
    theta,phi = hp.pix2ang( ns, ipix )
    radeg = 180./np.pi
    theta_deg = ( np.pi/2.-theta ) * radeg
    phi_deg = phi * radeg
    cut[ np.logical_and( theta_deg<gal_cut, theta_deg>-gal_cut )] = 0.
    if verbose:
        hp.mollview(cut)
    return cut
#  Simple ILC code -- Translate from IDL
"""def ilc(mapfiles=maps, maskfile=maskfile, gal_cut=0, silent=False, check=False, show_mask=False, double=False, md_rem=False, verbose=False, show_map=False):

    import healpy as hp
    import numpy as np

    if (len(maskfile) != 0) and ( gal_cut == 0) and verbose ):
        print ' --> ILC performed on the whole sky.'

    if (not keyword_set(maps)) then begin
    nfreq = len( mapfiles )
    input_type = type( mapfiles[0] )
    if input_type == 'str':
        if verbose:
            print ' - determining map size...:'    
        map = hp.read_map( mapfiles[0] )
        npix = len( map )
        ns = np.sqrt( npix/12 )
        if verbose:
            print ' - determining map size...: %s' %nside    

        maps = np.zeros( npix, nfreq, dtype=np.float )
        print 'reading maps...'

        psmask = np.ones( npix )
        for ifreq in range(nfreq):
            print ' - map ', ifreq+1, ' of', nfreq, ' :'+mapfiles[ifreq]
            temp = hp.read_map( mapfiles[ifreq] )
            psmask[ int (np.logical_not( (temp == -1.6375e30) or (np.isfinite(temp)) ) ) ] = 0.
            maps[*,ifreq] = temp
            if check:
                hp.mollview(temp, min=-300, max=300., fig=ifreq+1)
            

    else:
        npix = n_elements( maps[*,0] )
        ns = np.sqrt( npix/12 )
        psmask = np.ones( npix )
        for ifreq in range(nfreq):
            print ' - map ', ifreq+1, ' of', nfreq, ' :'+mapfiles[ifreq]
            temp = ( mapfiles[ifreq] )
            psmask[ int (np.logical_not( (temp == -1.6375e30) or (np.isfinite(temp)) ) ) ] = 0.
            if check:
                hp.mollview(temp, min=-300, max=300., fig=ifreq+1)

    mask = np.ones( npix )

    if len( maskfile ) gt 0:
        print, ' - reading '+maskfile+'...'
        gmask = hp.read_map( maskfile )
        mask_ns = len(gmask)
        if mask_ns != ns:
            gmask = hp.ud_grade( gmask, order_in='RING', nside_out=ns)
            thres = int( np.logical_not( gmask >= 0.5) )
            gmask[ thres ] = 0
            thres = int( np.logical_not( gmask < 0.5) )
            gmask[ thres ] = 1
        mask = mask * gmask

    if gal_cut>0.:
        dpfunc.make_sky_cut( gal_cut, ns ) 
        mask = mask * sky_cut
    if (keyword_set(show_mask) or keyword_set(check)) then mollview, mask, px=600, win=0, chars=1.5, tit='!6Mask used'
    
;for ifreq = 0, nfreq-1 do maps[*,ifreq] = maps[*,ifreq] * cmb2sz[ifreq]

    if (keyword_set(md_rem)) then begin
        print, 'removing mono/dipole...'
        for ifreq=0,nfreq-1 do begin
            temp = maps[*,ifreq]
            remove_dipole, temp, mask, ordering='ring', nside=ns
        endfor
    endif

gpix = where( (mask[*,0] gt 0.) and (pmask gt 0.) )
ngpix = n_elements(gpix)
bpix = where(mask[*,0] eq 0.)
nbpix = n_elements(bpix)


gR = dblarr(nfreq, nfreq)
bR = dblarr(nfreq, nfreq)
print, ' - :DP - ILC: computing correlation matrix...'
;------ Pedestrian
if (False) then begin
    ave = dblarr(nfreq,2)

    for ifreq = 0, nfreq-1 do begin
        gavei = mean(maps[gpix, ifreq])
        ave[ifreq,0] = gavei
        if (bpix[0] ne -1) then begin
            bavei = mean(maps[bpix, ifreq])
            ave[ifreq,1] = bavei
        endif
        for jfreq = ifreq, nfreq-1 do begin
            avej = mean(maps[gpix,jfreq])
            gR[ifreq, jfreq] = total( (maps[gpix,ifreq]-gavei) * (maps[gpix,jfreq]-avej) ) / ngpix
;        gR[ifreq, jfreq] = total( (maps[gpix,ifreq]) * (maps[gpix,jfreq]) ) / ngpix
            gR[jfreq, ifreq] = gR[ifreq, jfreq]
            
            if (bpix[0] ne -1) then begin        
                avej = mean(maps[bpix,jfreq])
                bR[ifreq, jfreq] = total( (maps[bpix,ifreq]-bavei) * (maps[bpix,jfreq]-avej) ) / nbpix
;        bR[ifreq, jfreq] = total( (maps[bpix,ifreq]) * (maps[bpix,jfreq]) ) / nbpix
                bR[jfreq, ifreq] = bR[ifreq, jfreq]
            endif
        endfor
    endfor
endif
;------
for ifreq = 0, nfreq-1 do begin
    if not ( keyword_set(silent) ) then print, ifreq
    for jfreq = ifreq, nfreq-1 do begin
        gR[ifreq, jfreq] = correlate( maps[gpix,ifreq], maps[gpix,jfreq], /covariance )
        gR[jfreq, ifreq] = gR[ifreq, jfreq]

        if (bpix[0] ne -1) then begin        
            bR[ifreq, jfreq] = correlate( maps[bpix,ifreq], maps[bpix,jfreq], /covariance )
            bR[jfreq, ifreq] = bR[ifreq, jfreq]
        endif
    endfor
endfor
; ------


gRm1 = invert(gR, /double, status)
print, status
if (bpix[0] ne -1) then begin
    bRm1 = invert(bR, /double, status)
    print, status
endif
a = findgen(nfreq) * 0. + 1.

gw = dblarr(nfreq)
gw = total(gRm1,2) / total(gRm1)
if (bpix[0] ne -1) then bw = total(bRm1,2) / total(bRm1)

print, ' - gw: ', gw
if (bpix[0] ne -1) then print, ' - bw: ', bw

ilc = fltarr(npix,3)
if (keyword_set(double)) then ilc = dblarr(Npix, 3)

for ifreq = 0, nfreq-1 do begin
;##	ilc[gpix] = ilc[gpix] + maps[gpix,ifreq] * gw[ifreq] 
;##        ilc[bpix] = ilc[bpix] + maps[bpix,ifreq] * bw[ifreq]
	ilc[*,0] = ilc[*,0] + maps[*,ifreq] * gw[ifreq] 
        if (bpix[0] ne -1) then ilc[*,1] = ilc[*,1] + maps[*,ifreq] * bw[ifreq]
endfor

ilc[gpix,2] = ilc[gpix,0]
ilc[bpix,2] = ilc[bpix,1]

print, 'STDDEV ILC gp (fs,out,in)       = ', stddev(ilc[*,0]), stddev(ilc[gpix,0]), stddev(ilc[bpix,0])
print, 'STDDEV ILC bp (fs,out,in)       = ', stddev(ilc[*,1]), stddev(ilc[gpix,1]), stddev(ilc[bpix,1])
print, 'STDDEV ILC combined (fs,out,in) = ', stddev(ilc[*,2]), stddev(ilc[gpix,2]), stddev(ilc[bpix,2])

;write_fits_map, 'ffp4_01a_ILCout.fits', ilc[*,0], /ring, units='!7l!6K'
;write_fits_map, 'ffp4_01a_ILCin.fits', ilc[*,1], /ring, units='!7l!6K'


if (not keyword_set(silent)) then mollview, ilc[*,0], chars=1.5, tit='!6ILC: weights outside mask. N!dside!n='+string(ns,format='(i4.4)'), grat=[10,10], px=650;, no_monopole=true, gal_cut=40, min=-300, max=300
if (not keyword_set(silent) and (bpix[0] ne -1) ) then mollview, ilc[*,1], chars=1.5, tit='!6ILC: weights inside mask. N!dside!n='+string(ns,format='(i4.4)'), grat=[10,10], px=650;, no_monopole=true, gal_cut=40, min=-300, max=300
if (not keyword_set(silent) and (bpix[0] ne -1) ) then mollview, ilc[*,2], chars=1.5, tit='!6ILC: combined', grat=[10,10], px=650;, no_monopole=true, gal_cut=40, min=-300, max=300
if (not keyword_set(silent) and (bpix[0] ne -1) ) then mollview, ilc[*,0]-ilc[*,1], chars=1.5, tit='!6ILC Difference', grat=[10,10], px=650;, no_monopole=true, gal_cut=40, min=-15, max=15

if (do_png) then begin
    restore, 'chains/pix_01a_v2.res.sav'
    mollview, cmb, min=-300, max=300, chars=1.5, win=4, tit='!6Commander', no_monopole=true, gal_cut=40
    mollview, cmb-ilc[*,0], min=-30, max=30, chars=1.5, win=5, tit='!6Commander-ILC!dout!n', no_monopole=true, gal_cut=40
    mollview, cmb-ilc[*,1], min=-30, max=30, chars=1.5, win=6, tit='!6Commander-ILC!din!n', no_monopole=true, gal_cut=40
    
    read_fits_map, 'ffp4_scalar_cmb_ns128_60arcmin_uK.fits', inp
    mollview, cmb-inp, min=-30, max=30, chars=1.5, win=7, tit='!6Commander-Input', no_monopole=true, gal_cut=40
    mollview, ilc[*,0]-inp, min=-30, max=30, chars=1.5, win=8, tit='!6ILC!dout!n-Input', no_monopole=true, gal_cut=40
    mollview, ilc[*,1]-inp, min=-30, max=30, chars=1.5, win=9, tit='!6ILC!din!n-Input', no_monopole=true, gal_cut=40

    mollview, ilc[*,0], min=-300, max=300, chars=1.5, win=-1, tit='!6ILC: weights outside mask', no_monopole=true, gal_cut=40, png='ffp4_01a_ILCout.png'
    mollview, ilc[*,1], min=-300, max=300, chars=1.5, win=-2, tit='!6ILC: weights inside mask', no_monopole=true, gal_cut=40, png='ffp4_01a_ILCin.png'
    mollview, ilc[*,0]-ilc[*,1], min=-30, max=30, chars=1.5, win=-3, tit='!6ILC Difference', no_monopole=true, gal_cut=40, png='ffp4_01a_ILCout-in.png'

    mollview, cmb, min=-300, max=300, chars=1.5, win=-4, tit='!6Commander', no_monopole=true, gal_cut=40, png='ffp4_01a_CMD.png'
    mollview, cmb-ilc[*,0], min=-30, max=30, chars=1.5, win=-5, tit='!6Commander-ILC!dout!n', no_monopole=true, gal_cut=40, png='ffp4_01a_CMD-ILCout.png'
    mollview, cmb-ilc[*,1], min=-30, max=30, chars=1.5, win=-6, tit='!6Commander-ILC!din!n', no_monopole=true, gal_cut=40, png='ffp4_01a_CMD-ILCin.png'

    mollview, cmb-inp, min=-30, max=30, chars=1.5, win=-7, tit='!6Commander-Input', no_monopole=true, gal_cut=40, png='ffp4_01a_CMD-INP.png'
    mollview, ilc[*,0]-inp, min=-30, max=30, chars=1.5, win=-8, tit='!6ILC!dout!n-Input', no_monopole=true, gal_cut=40, png='ffp4_01a_ILCout-INP.png'
    mollview, ilc[*,1]-inp, min=-30, max=30, chars=1.5, win=-9, tit='!6ILC!din!n-Input', no_monopole=true, gal_cut=40, png='ffp4_01a_ILCin-INP.png'
endif

print, ' --- End of Program ---'

return, ilc
;stop

; matrix multiplication failures
gw = grm1##a 
n = reform(a##(grm1##a))
gw[*] = gw[*] / n[0]

print, gw

bw = brm1##a
n = reform(a##(brm1##a))
bw[*] = bw[*] / n[0]


stop


end"""
