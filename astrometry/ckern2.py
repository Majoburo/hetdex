import sys
import os
import glob
import numpy as np
import argparse
import os.path as op
import pyfits
from scipy.signal import convolve2d, correlate2d
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage import zoom,interpolation,filters
from scipy import interpolate
from numpy.linalg import inv
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.io import fits 
gfilter = 'filters/g_sdss.dat'
sdss_image = 'imaging/sdssDR12g.fits'
#sdss_image = 'imaging/J134225.00+282357.0-g.fits'
fwhm_vir  = 1.6 #DUMMY VALUE!!!
fwhm_sdss = 1.53 # ARCSEC FWHM
parangle = 255.976452 #HARCODED FROM FE FILE FOR NOW...
sigma = pow(fwhm_vir*fwhm_vir-fwhm_sdss*fwhm_sdss,0.5)*0.4247/3600
#ifu_ra,ifu_dec= 205.63873144, 28.41886915

def cirkernel(radious):
    y,x = np.ogrid[-radious:radious+1,-radious:radious+1]
    mask = x*x+y*y<=radious*radious
    return mask*1

def degtopix(hdu,degrees):

    #transforming degrees to pixels
    mat = np.round([[hdu[0].header['CD1_1'],hdu[0].header['CD1_2']],[hdu[0].header['CD2_1'],hdu[0].header['CD2_2']]],10)
    CDELT1 = pow(abs(np.linalg.det(mat)),0.5) #Coordinate increment per pixel in DEGREES/PIXEL
    #print np.arcsin(np.trunc(1000*mat/CDELT1)/1000)*180./np.pi
    #vec = np.dot(inv(mat),[[degrees],[0]])
    #sdssangle = np.arcsin(np.round(hdu[0].header['CD1_1'],12)/CDELT1)*180./np.pi
    #print sdssangle
    #print 'determinant'
    return degrees/CDELT1
    #print vec
    #pix = pow(vec[0]*vec[0]+vec[1]*vec[1],0.5)[0]
    #return pix

def filterweights():
    """
    FOR NOW USING HARDCODED VALUES SINCE CUBE HAS NO HEADERS!!!

    OUTPUT:
        int
        x = [filtered weights]
    """

    xmin = 3500
    x_range = 1054
    step = 1.9
    xmax = step*(x_range-1) + xmin
    wl = np.linspace(xmin,xmax,x_range)
    
    ifilter = np.loadtxt(gfilter)
    weightf=interpolate.interp1d(ifilter[:,0],ifilter[:,1],fill_value=0,bounds_error=False)
    weight=[]
    for wavelenght in wl:
        weight.append(weightf(wavelenght))
    
    return weight

def rotate(hdu):
    
    mat = [[hdu[0].header['CD1_1'],hdu[0].header['CD1_2']],[hdu[0].header['CD2_1'],hdu[0].header['CD2_2']]]
    CDELT1 = pow(abs(np.linalg.det(mat)),0.5) #Coordinate increment per pixel in DEGREES/PIXEL
    sdssangle = np.arcsin(hdu[0].header['CD1_1']/CDELT1)*180./np.pi
    #pixwidth=1./3600 #cube has one arsec pixel width; this is PIXEL/DEGREE
    #resizefactor = degtopix(hdu,pixwidth)
    #CDELT1= CDELT1*resizefactor
    
    #print degtopix(hdu,0.75/3600)
    hdu[0].data = interpolation.rotate(hdu[0].data,(parangle-sdssangle))
    hdu[0].header['CD1_1']= np.sin(parangle*np.pi/180.)*CDELT1
    hdu[0].header['CD1_2']= np.cos(parangle*np.pi/180.)*CDELT1
    hdu[0].header['CD2_1']= np.cos(parangle*np.pi/180.)*CDELT1
    hdu[0].header['CD2_2']= -np.sin(parangle*np.pi/180.)*CDELT1
    hdu[0].header['CRPIX1']= hdu[0].data.shape[1]/2+0.5
    hdu[0].header['CRPIX2']= hdu[0].data.shape[0]/2+0.5
    hdu.writeto('imaging/LsdssDR12g.fits',clobber =True)
    #print degtopix(hdu,0.75/3600)
    return hdu

def crop(hdu,ifu_ra,ifu_dec):
    from astropy.nddata import Cutout2D
    from astropy.coordinates import SkyCoord
    c=SkyCoord(ifu_ra,ifu_dec, unit="deg")
    size = (124*2, 124*2)     # cutting a box twice the size of the ifu
    w=WCS('imaging/LsdssDR12g.fits')
    cutout = Cutout2D(hdu[0].data, c, size, wcs=w)
    plt.imshow(hdu[0].data, origin='lower')
    cutout.plot_on_original(color='white')
    plt.show()

    hdu = fits.PrimaryHDU(cutout.data,header=cutout.wcs.to_header()) #create new hdu
    hdulist = fits.HDUList([hdu]) #create new hdulist
    hdulist.writeto('imaging/CCCsdssDR12g.fits',clobber =True)
    #CRVAL1  = 2.05574173860000E+02 # RA at Reference Pixel                          
    #CRVAL2  = 2.84984933500000E+01 # DEC at Reference Pixel
    
    return hdu

def flatten_wspectracube(datacube,weight):
    '''
    Weighting each spectra slice of the cube
    '''
    i=0
    for wv in datacube:
        wv = wv*weight[i]
        i+=1
    datacube = datacube.sum(axis=0)

    return datacube

#def resize(hdu):
def resize(hdu,cube):
    pixwidth=1./3600 #cube has one arsec pixel width; this is pixels per degree
    resizefactor = degtopix(hdu,pixwidth)
    # CDELT1=1./resizefactor
    cube = zoom(cube,resizefactor)
    #hdu[0].data = zoom(hdu[0].data,CDELT1)
    #cube = interpolation.rotate(cube,-(parangle-10.1137)) #10.1137 is the rotation angle of the sdss image

    return cube

def parseargs(argv=None):
    '''
    Options parser
    '''
    #parser = OptionParser()

    parser = argparse.ArgumentParser(description = "Calculate astrometry based on Imaging")


    parser.add_argument("basename",nargs='?', type=str,action="store",
                              help="basename of fits FILEs to extract spectra from")
    
    args = parser.parse_args(args=argv)

    if args.basename is None:
        msg = 'The base name was not provided'
        parser.error(msg)
    else:

        args.fitsfiles = []
        args.position = []
        #Searching for files with that name...
        searchname = args.basename + '*.fits'
        filenames = glob.glob(searchname)
        if not filenames:
            msg = 'No files found searching for: {:s}'.format(searchname)
            parser.error(msg)
        else:
            for i in xrange(len(filenames)):
                args.fitsfiles.append(pyfits.open(filenames[i])[0])
                args.position.append(filenames[i][-12:-9])

    return args



def main():
    args = parseargs()
    fitsfiles = args.fitsfiles
    print fitsfiles,args.position

    hdu = pyfits.open(sdss_image)
    #resizing sdss image to cube's size
    #integrating over wavlengths
    weight = filterweights()
    i=0
    for fitsfile in fitsfiles:#for now just doing one ifu
        virus_data = flatten_wspectracube(fitsfile.data,weight)
        fitsfile.data = virus_data
        #plt.imshow(fitsfile.data)
        #plt.show()
        virus_data = resize(hdu,virus_data)
        args.position[i]
        fitsfile.writeto('imaging/cube'+args.position[i]+'.fits',clobber =True)
        i+=1

    #Convolving the image with a gaussian instead of the kernel, since the kernel is to small.
    sigmapix = degtopix(hdu,sigma)
    hdu[0].data = gaussian_filter(hdu[0].data,sigmapix)
    ckern = cirkernel(1.89553)
    hdu[0].data = convolve2d(hdu[0].data,ckern,mode='same')
    ifu_cen=np.loadtxt('../greg/ifuPos.txt')
    print ifu_cen
    hdu = rotate(hdu)
    sdsshdu=hdu
    for ifu in ifu_cen:
        print ifu
        ifu_id,ifu_ra, ifu_dec = ifu[0],ifu[1],ifu[2]
        hdu = crop(hdu,ifu_ra,ifu_dec)
        a=hdu.data
        #print np.std(a)
        #print np.max(a)
        #print np.min(a)
        sdssimg1 = hdu.data
        hdu.data=(a<np.std(a)*10)*a+(a>np.std(a)*10)*a[a<np.std(a)*10].mean()
        hdu.writeto('imaging/Rsdssimg.fits',clobber =True)
        sdssimg = hdu.data
        #virus_data = sdssimg[10:40,40:60]
        #print virus_data.shape
        virus_data = virus_data- virus_data.mean()
        sdssimg =sdssimg - sdssimg.mean()
        
       # hdu[0].data = sdssimg
       # hdu.writeto('imaging/sdssimg.fits',clobber =True)
       # hdu[0].data = virus_data
       # hdu.writeto('imaging/imgcut.fits',clobber =True)
        
        hdu.data = correlate2d(sdssimg,virus_data,boundary='symm',mode='same')
        y, x = np.unravel_index(np.argmax(hdu.data), hdu.data.shape)
        print y,x
        from astropy import wcs
        w = wcs.WCS(hdu.header)
        ra,dec= w.wcs_pix2world([[x,y]],1)[0]
        print 'IFU %d IS CENTERED AT:'%ifu_id
        print ra,dec
        print 'The offset is:'
        print(ra-ifu_ra,dec-ifu_dec)
        '''
        #Plotting
        
        fig, (ax_orig, ax_template) = plt.subplots(1, 2)
        ax_orig.imshow(sdssimg)
        ax_orig.set_title('SDSS-image')
        ax_orig.set_axis_off()
        ax_template.imshow(virus_data)
        ax_template.set_title('IFU 094')
        ax_template.set_axis_off()
        ax_orig.plot(x, y, 'ro')
        plt.show()
        '''
        hdu.data=sdssimg1
        hdu.data[y,x]=0
        for i in range(124):
            hdu.data[y+124/2,x+i-124/2]=0
            hdu.data[y-124/2,x+i-124/2]=0
        for j in range(124):
            hdu.data[y+j-124/2,x+124/2]=0
            hdu.data[y+j-124/2,x-124/2]=0
        hdu.writeto(('imaging/GBsdssDR12g%d.fits')%ifu_id,clobber =True)
        hdu=sdsshdu
if __name__ == "__main__":
    main()
