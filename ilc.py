#  Simple ILC code -- Translate from IDL
def ilc(mapfiles=maps, maskfile=maskfile, gal_cut=0, silent=False, check=False, show_mask=False, double=False, md_rem=False, verbose=False, show_map=False)

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


end
