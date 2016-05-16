from astropy.io import fits
from astropy import wcs
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.interpolation import map_coordinates,zoom

'''
Given a fits image and fiber positions and diameters,
this function will return the integrated flux through each fiber.

'''
def get_flux(sdss_fits_fname,fiber_ra,fiber_dec,fiber_dia,zoom_factor,errors=False) :

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
    rp = 5 # TODO fix this!!!

    # we'll take a sum over each fiber's pixels.
    fibsum = np.zeros(pixcrd.shape[0],dtype='float32')
    if errors:
        fiberr = np.zeros(pixcrd.shape[0],dtype='float32')

    # we're going to find the region of the image we want to sum over for each fiber
    # and 'zoom' into it (i.e., make a 3x3 region a 3*zoom_factor x 3*zoom_factor region,
    #   interpolating to fill in the gaps)
    # 'zoom'-ing like this reduces the error from summing square pixels in a circular region
    effr = zoom_factor*rp-0.5

    rad_to_zoom = rp+3 # keep a few extra pixels for when we interpolate
    dzoom = rad_to_zoom*zoom_factor
    y,x = np.ogrid[-dzoom:dzoom,-dzoom:dzoom]
    mask = x*x+y*y <= effr*effr
    errmask = x*x+y*y <= rp*rp*zoom_factor**2

    for i in range(pixcrd.shape[0]):
        a,b = pixcrd[i]
        # zoom in to the cropped region around our fiber
        #plt.imshow(img_data[int(a-rad_to_zoom):int(a+rad_to_zoom),int(b-rad_to_zoom):int(b+rad_to_zoom)],interpolation='nearest',cmap='Greys_r',origin='lower')
        #plt.show()
        zoomedregion=zoom(input=img_data[int(a-rad_to_zoom):int(a+rad_to_zoom),int(b-rad_to_zoom):int(b+rad_to_zoom)]
                ,zoom=zoom_factor,order=3)
        #plt.imshow(zoomedregion,interpolation='nearest',cmap='Greys_r',origin='lower')
        #plt.show()
        fibsum[i] = np.sum(zoomedregion[mask])/(zoom_factor**2)
        if errors:
            fiberr[i] = np.abs(fibsum[i] - np.sum(zoomedregion[errmask])/zoom_factor**2)
            # this error is roughly the possible missed flux from approximating the circle
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
    print(get_flux("imaging/frame-g-002326-3-0078.fits",fiber_ra[:30],fiber_dec[:30],fiber_dia,zoom_factor=2.0,errors=False))
