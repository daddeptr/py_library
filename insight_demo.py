import dpfunc
import healpy as hp
import numpy as np

map100 = hp.read_map('../../XFaster/data/maps/ffp7/ffp7_IQUmap_100_ns2048_uK_hrhs.fits', field=(0,1,2))
hp.mollview(map100, min=-500, max=500, unit='$\mu K$', title='100GHz Temperature Map')

map143 = hp.read_map('../../XFaster/data/maps/ffp7/ffp7_IQUmap_143_ns2048_uK_hrhs.fits', field=(0,1,2))
hp.mollview(map143, min=-500, max=500, unit='$\mu K$', title='143GHz Temperature Map')

map217 = hp.read_map('../../XFaster/data/maps/ffp7/ffp7_IQUmap_217_ns2048_uK_hrhs.fits', field=(0,1,2))
hp.mollview(map217, min=-500, max=500, unit='$\mu K$', title='217GHz Temperature Map')

map353 = hp.read_map('../../XFaster/data/maps/ffp7/ffp7_IQUmap_353_ns2048_uK_hrhs.fits', field=(0,1,2))
hp.mollview(map353[0], min=-500, max=10000, unit='$\mu K$', title='353GHz Temperature Map')
hp.mollview(np.sqrt(map353[1]**2+map353[2]**2), min=0, max=1000, unit='$\mu K$', title='353GHz Polarization Map')

map143undust = hp.read_map('../../XFaster/data/maps/ffp7/ffp7_IQUmap_143_undusted_insideM_ns2048_uK_hrhs.fits', field=(0,1,2))
hp.mollview(map143undust[0], min=-500, max=500, unit='$\mu K$', title='CMB Temperature Map')
hp.mollview(np.sqrt(map143undust[1]**2+map143undust[2]**2), min=0, max=100, unit='$\mu K$', title='CMB Polarization Map')

map143undust_noise = hp.read_map('../../XFaster/data/maps/ffp7/ffp7_IQUmap_143_undusted_insideM_ns2048_uK_hrhd.fits', field=(0,1,2))
hp.mollview(map143undust_noise, min=-500, max=500, unit='$\mu K$', title='CMB Temperature Noise Map')

xffile_cmb='../../XFaster/outputs/ffp7_IQUcmb_map_143_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
xffile_raw='../../XFaster/outputs/ffp7_IQUmap_143_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'
xffile_clean='../../XFaster/outputs/ffp7_IQUmap_143_undusted_insideM_ns2048_uK_hrhs_xfcl_l2000_dx9_common_xGal06_mask_ns2048_X_QUmask_DX10_2048_T5.0_60_X_ps100-353_l3000.newdat'

dpfunc.read_xfaster_newdata_output(xffile_clean, pol=True, res=True)
