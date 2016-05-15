from astropy.io import fits
from astropy import wcs
import numpy as np
#import matplotlib.pyplot as plt
from scipy.ndimage.interpolation import zoom

'''
Given a fits image and fiber positions and diameters,
this function will return the integrated flux through each fiber.

'''
def get_flux(sdss_fits_fname,fiber_ra,fiber_dec,fiber_dia,zoom_factor=1,errors=False) :

    with fits.open(sdss_fits_fname) as h:
        img_data = h[0].data
        w = wcs.WCS(h[0].header)
    pixcrd = w.wcs_world2pix(np.array(zip(fiber_ra,fiber_dec),dtype='float64'),1)

    #plt.scatter(pixcrd[:,0],pixcrd[:,1])
    #plt.show()
    #img_data[np.round(pixcrd).astype(np.int)] = 0x000000

    #plt.imshow(img_data,interpolation='nearest',cmap='Greys_r',origin='lower')
    #plt.show()

    rp = w.wcs_world2pix([[fiber_ra[0]+fiber_dia,fiber_dec[0]]],1)
    rp -= w.wcs_world2pix([[fiber_ra[0],fiber_dec[0]]],1)
    rp = 1 # TODO fix this!!!

    # we'll take a sum over each fiber's pixels.
    fibsum = np.zeros(pixcrd.shape[0])
    if errors:
        fiberr = np.zeros(pixcrd.shape[0])

    zoomed_img = zoom(input=img_data,zoom=zoom_factor,order=1)
    effr = zoom_factor*rp-0.5
    nx,ny = zoomed_img.shape

    for i in range(pixcrd.shape[0]):
        a,b = zoom_factor*pixcrd[i]
        # take a circular mask around the fiber's center
        y,x = np.ogrid[-a:nx-a,-b:ny-b]
        mask = x*x+y*y <= effr*effr
        fibsum[i] = np.sum(zoomed_img[mask])/(zoom_factor**2)
        if errors:
            fiberr[i] = np.abs(fibsum[i] - np.sum(zoomed_img[x*x+y*y <= rp*rp*zoom_factor**2])/zoom_factor**2)
            # this error is roughly the possible missed flux from pixels on the edge
    if errors:
        return fibsum,fiberr
    else: return fibsum

if __name__ == '__main__':
    with open('./offsets/ifuPos.txt','r') as f:
        f.readline()
        fiber_dia = float(f.readline())
        fiber_ra , fiber_dec = [],[]
        for line in f:
            a,b = line.split()
            fiber_ra.append(float(a))
            fiber_dec.append(float(b))
    print("Flux at the first 30 fibers: ")
    print(get_flux("imaging/frame-g-002326-3-0078.fits",fiber_ra,fiber_dec,fiber_dia,zoom_factor=1,errors=False))
