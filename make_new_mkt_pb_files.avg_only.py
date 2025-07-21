import numpy as np
from astropy.io import fits

wget https://archive-gw-1.kat.ac.za/public/repository/10.48479/wdb0-h061/data/MeerKAT_L_band_primary_beam.npz

mdv = np.load("MeerKAT_L_band_primary_beam.npz")
bm = mdv['beam']
ants = mdv['antnames']
freqs = mdv['freq_MHz']*1e+6
degs = mdv['margin_deg']
i0 = (len(degs) // 2) - 1

xx = bm[0,-1,:,1:,1:]
xy = bm[1,-1,:,1:,1:]
yx = bm[2,-1,:,1:,1:]
yy = bm[3,-1,:,1:,1:]
# normalize by beam gain at centre, since the bandpass will do that for us
xx_norm = xx / xx[:,i0,i0][:,np.newaxis,np.newaxis]
xy_norm = xy / yy[:,i0,i0][:,np.newaxis,np.newaxis]
yx_norm = yx / xx[:,i0,i0][:,np.newaxis,np.newaxis]
yy_norm = yy / yy[:,i0,i0][:,np.newaxis,np.newaxis]
xx_re = np.real(xx_norm)
xx_im = np.imag(xx_norm)
xy_re = np.real(xy_norm)
xy_im = np.imag(xy_norm)
yx_re = np.real(yx_norm)
yx_im = np.imag(yx_norm)
yy_re = np.real(yy_norm)
yy_im = np.imag(yy_norm)

xx_mag = np.sqrt(xx_re**2 + xx_im**2)
xy_mag = np.sqrt(xy_re**2 + xy_im**2)
yx_mag = np.sqrt(yx_re**2 + yx_im**2)
yy_mag = np.sqrt(yy_re**2 + yy_im**2)
components = {"xx_re": xx_re, "xx_im":xx_im, "xy_re":xy_re, "xy_im":xy_im, "yx_re":yx_re, "yx_im":yx_im, "yy_re":yy_re, "yy_im":yy_im}
for key, value in components.items():
    hdu = fits.PrimaryHDU(value[:,:,:])
    hdr = hdu.header
    hdr['CRPIX1'] = i0+1
    hdr['CRPIX2'] = i0+1
    hdr['CRPIX3'] = 1
    hdr['CRVAL1'] = 0
    hdr['CRVAL2'] = 0
    hdr['CRVAL3'] = freqs[0]
    hdr['CDELT1'] = degs[1] - degs[0]
    hdr['CDELT2'] = degs[1] - degs[0]
    hdr['CDELT3'] = freqs[1] - freqs[0]
    hdr['CTYPE1'] = 'X'
    hdr['CTYPE2'] = 'Y'
    hdr['CTYPE3'] = 'FREQ'
    hdr['CUNIT1'] = 'deg'
    hdr['CUNIT2'] = 'deg'
    hdr['CUNIT3'] = 'Hz'
    hdu.writeto('new_mkt_pb_avg_reverse_%s.fits' % key, overwrite=True)


