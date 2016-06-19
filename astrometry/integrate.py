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
    #img_data[np.round(pixcrd).astype(np.int)] = 0x000000

    #plt.imshow(img_data,interpolation='nearest',cmap='Greys_r',origin='lower')
    #plt.show()
    rp = w.wcs_world2pix([[fiber_ra[0],fiber_dec[0]+float(fiber_dia)/3600.]],1)
    rp -= w.wcs_world2pix([[fiber_ra[0],fiber_dec[0]]],1)
    rp = rp[0,0]
    #print rp
    pixcrd = w.wcs_world2pix(np.array(zip(fiber_ra,fiber_dec),dtype='float64'),1)
    # we'll take a sum over each fiber's pixels.
    fibsum = np.zeros(pixcrd.shape[0],dtype='float32')
    if errors:
        fiberr = np.zeros(pixcrd.shape[0],dtype='float32')

    # we're going to find the region of the image we want to sum over for each fiber
    # and 'zoom' into it (i.e., make a 3x3 region a 3*zoom_factor x 3*zoom_factor region,
    #   interpolating to fill in the gaps)
    # 'zoom'-ing like this reduces the error from summing square pixels in a circular region
    effr = zoom_factor*rp

    rad_to_zoom = np.ceil(rp+3) # keep a few extra pixels for when we interpolate
    dzoom = rad_to_zoom*zoom_factor
    y,x = np.ogrid[-dzoom:dzoom,-dzoom:dzoom]
    mask = x*x+y*y <= effr*effr
    errmask = x*x+y*y <= rp*rp*zoom_factor**2
    for i in range(pixcrd.shape[0]):
        a,b = pixcrd[i]
        # zoom in to the cropped region around our fiber
        #plt.imshow(img_data[int(a-rad_to_zoom):int(a+rad_to_zoom),int(b-rad_to_zoom):int(b+rad_to_zoom)],interpolation='nearest',cmap='Greys_r',origin='lower')
        #plt.show()
        zoomedregion=zoom(input=img_data[int(b-rad_to_zoom):int(b+rad_to_zoom),int(a-rad_to_zoom):int(a+rad_to_zoom)]
                ,zoom=zoom_factor,order=3)
        #plt.imshow(zoomedregion,interpolation='nearest',cmap='Greys_r',origin='lower')
        #plt.show()
        fibsum[i] = np.sum(zoomedregion[mask])/(zoom_factor**2)
        if errors:
            fiberr[i] = np.abs(fibsum[i] - np.sum(zoomedregion[errmask])/zoom_factor**2)
            # this error is roughly the possible missed flux from approximating the circle
    '''
    Visualization
    '''
    ######
    #fibsum=np.arange(len(fibsum))
    xmin = int(np.floor(np.min(pixcrd[:,0])))
    xmax = int(np.ceil(np.max(pixcrd[:,0])))
    ymin = int(np.floor(np.min(pixcrd[:,1])))
    ymax = int(np.ceil(np.max(pixcrd[:,1])))

    test = np.mean(fibsum)*np.ones((xmax-xmin,ymax-ymin))
    i=0
    for x,y in pixcrd:
        test[int(np.floor(y-ymin-2)):int(np.ceil(y-ymin+2)),int(np.floor(x-xmin-2)):int(np.ceil(x-xmin+2))] = fibsum[i]
        i+=1
    
    hdu = fits.PrimaryHDU(test) #create new hdu
    hdulist = fits.HDUList([hdu]) #create new hdulist
    hdulist.writeto('CsdssDR12g.fits',clobber =True)

    virusflux = np.loadtxt('virusflux.txt')
    '''
    Visualization
    '''


    ######
    xmin = int(np.floor(np.min(pixcrd[:,0])))
    xmax = int(np.ceil(np.max(pixcrd[:,0])))
    ymin = int(np.floor(np.min(pixcrd[:,1])))
    ymax = int(np.ceil(np.max(pixcrd[:,1])))
    #plt.imshow(img_data[xmin:xmax,ymin:ymax],interpolation='nearest',cmap='greys_r',origin='lower')
    #plt.scatter(pixcrd[:,0]-xmin,pixcrd[:,1]-ymin)
    #plt.show()

    test = np.mean(virusflux)*np.ones((xmax-xmin,ymax-ymin))
    i=0
    for x,y in pixcrd:
        test[int(np.floor(y-ymin-2)):int(np.ceil(y-ymin+2)),int(np.floor(x-xmin-2)):int(np.ceil(x-xmin+2))] = virusflux[i]
        i+=1
    hdu = fits.PrimaryHDU(test) #create new hdu
    hdulist = fits.HDUList([hdu]) #create new hdulist
    hdulist.writeto('VsdssDR12g.fits',clobber =True)


    h=fits.open(sdss_fits_fname)
    x=xmin
    y=ymin
    '''
    for i in range(160):
        h[0].data[y+160,x+i]=0
        h[0].data[y,x+i]=0
    for j in range(160):
        h[0].data[y+j,x+160]=0
        h[0].data[y+j,x]=0
    '''
    h[0].data=h[0].data[y:160+y,x:160+x]
    h.writeto('VCsdssDR12g.fits',clobber =True)
    ###########

    if errors:
        return fibsum,fiberr
    else: return fibsum

if __name__ == '__main__':
    with open('./offsets/ifuPos094.txt','r') as f:
        fiber_ra , fiber_dec = [],[]
        for line in f:
            a,b,c,d = line.split()
            fiber_ra.append(float(a))
            fiber_dec.append(float(b))
        fiber_dia = c
    #print fiber_ra
    fiber_ra,fiber_dec=np.array(fiber_ra)+0.0053840509754081722,np.array(fiber_dec)+0.00027915552097468321
    get_flux("imaging/sdssDR12g.fits",fiber_ra,fiber_dec,fiber_dia,zoom_factor=4.0,errors=False)
