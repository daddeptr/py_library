pro merge_xfaster_newdat_files, file1, file2, l=l, bin=bin, xfdir=xfdir, old=old

    if (not keyword_set(xfdir)) then xfdir = '/global/scratch2/sd/dpietrob/Software/XFaster/'
    if (not keyword_set(bin)) and (not keyword_set(l)) then stop, 'ERROR: Neither l nor bin specified.'

;## NB consistent with my definition of .newdat
    if (ncl eq 1) or (strupcase(ncl) eq 'TT') then begin
        readcol,file,nbtt, format='f',skipline=2,numline=1
        print, ' TT # bin = ', nbtt
    endif else begin
        if keyword_set(old) then readcol,file,nbtt,nbee,nbbb,nbtb,nbte,nbeb, format='f,f,f,f,f,f',skipline=2,numline=1 else $
          readcol,file,nbtt,nbee,nbbb,nbte,nbtb,nbeb, format='f,f,f,f,f,f',skipline=2,numline=1
        
        print, ' EE # bin = ', nbee
        print, ' TE # bin = ', nbte
        print, ' BB # bin = ', nbbb
        print, ' TB # bin = ', nbtb
        print, ' EB # bin = ', nbeb
    endelse

; ------- Binning
    if (ncl eq 1) or (strupcase(ncl) eq 'TT') then begin
        readcol,file,cl,cler,lmn,lmx,format='x,f,f,x,x,f,f', skipline=7, numline=nbtt[0]
        btcl = bp_binning(tcl[*,0]*ll, xfdir+'data/bins/ctp/CTP_bin_TT')
    endif
    if (ncl eq 2) or (strupcase(ncl) eq 'EE') then begin
        readcol,file,cl,cler,lmn,lmx,format='x,f,f,x,x,f,f', skipline=2*nbtt[0]+7+1, numline=nbee[0]
        btcl = bp_binning(tcl[*,1]*ll, xfdir+'data/bins/ctp/CTP_bin_EE')
    endif
    if (ncl eq 3) or (strupcase(ncl) eq 'BB') then begin
        readcol,file,cl,cler,lmn,lmx,format='x,f,f,x,x,f,f', skipline=2*nbtt[0]+2*nbee[0]+7+1+1, numline=nbbb[0]
        btcl = bp_binning(tcl[*,2]*ll, xfdir+'data/bins/ctp/CTP_bin_TE')
    endif
    if (ncl eq 4) or (strupcase(ncl) eq 'TE') then begin
        readcol,file,cl,cler,lmn,lmx,format='x,f,f,x,x,f,f', skipline=7+2*nbtt[0]+2*nbee[0]+2*nbbb[0]+1+1+1, numline=nbte[0]
        btcl = bp_binning(tcl[*,3]*ll, xfdir+'data/bins/ctp/CTP_bin_EE')
    endif
    if (ncl eq 5) or (strupcase(ncl) eq 'TB') then begin
        readcol,file,cl,cler,lmn,lmx,format='x,f,f,x,x,f,f', skipline=7+2*nbtt[0]+2*nbee[0]+2*nbbb[0]+2*nbte[0]+1+1+1+1, numline=nbtb[0]
        btcl = bp_binning(tcl[*,4]*ll, xfdir+'data/bins/ctp/CTP_bin_TE')
    endif
    if (ncl eq 6) or (strupcase(ncl) eq 'EB') then begin
        readcol,file,cl,cler,lmn,lmx,format='x,f,f,x,x,f,f', skipline=7+2*nbtt[0]+2*nbee[0]+2*nbbb[0]+2*nbte[0]+2*nbtb[0]+1+1+1+1+1, numline=nbeb[0]
        btcl = bp_binning(tcl[*,5]*ll, xfdir+'data/bins/ctp/CTP_bin_EE')
    endif

    lcen = (lmn+lmx)/2.
    res = cl-btcl
    return, cl
end
