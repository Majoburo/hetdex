from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import sys
import os
import glob
import argparse
import os.path as op

import numpy as np
from scipy.signal import convolve2d, correlate2d
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage import zoom, interpolation, filters
from scipy import interpolate
from numpy.linalg import inv
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import wcs
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord

gfilter = 'filters/g_sdss.dat'
sdss_image = 'imaging/J134225.00+282357.0-g.fits'
fwhm_vir  = 1.6  # DUMMY VALUE!!!
fwhm_sdss = 1.53  # ARCSEC FWHM
#parangle = 255.976452  # HARCODED FROM FE FILE FOR NOW...
parangle = 256.054951
sigma = pow(fwhm_vir*fwhm_vir-fwhm_sdss*fwhm_sdss,0.5)*0.4247/3600


def cirkernel(radius):
    """Makes a circular mask of a given radius*

    Parameters
    ----------
    radius : float
        radius, in pixels

    Returns
    -------
    2darray
        mask with ones inside the circle and zeros outside
    """
    y,x = np.ogrid[-radius:radius+1,-radius:radius+1]
    mask = x*x+y*y<=radius*radius
    return mask*1

def degtopix(hdu,degrees):
    """Converts degrees to pixels

    Parameters
    ----------
    degrees : float
    hdu : fits file

    Returns
    -------
    float
        pixels corresponding to the given degrees
    """
    mat = np.round([[hdu[0].header['CD1_1'],hdu[0].header['CD1_2']],[hdu[0].header['CD2_1'],hdu[0].header['CD2_2']]],10)
    CDELT1 = pow(abs(np.linalg.det(mat)),0.5) #Coordinate increment per pixel in DEGREES/PIXEL

    return degrees/CDELT1


def filterweights(filter,xmin,x_range,step):
    """Interpolates a given set of spectral weights for given spectral range.
    *For now keeping this as a function. Could really just be a list of values.*

    Parameters
    ----------
    filter : 1d array
        Weights for given spectra (in Angstroms)
    xmin : float
        Value of the minimun spectral weight to calculate (in Angstroms)
    xrange : int
        Number of spectral values to calculate.
    xstep : float
        Step of the gives values to calculate (in Angstroms)

    Returns
    -------
    1darray
        Array of spectral weights at the desired wavelenghts.

    """

    #xmin = 3500
    #x_range = 1054
    #step = 1.9
    xmax = step*(x_range-1) + xmin
    wl = np.linspace(xmin,xmax,x_range)
    ifilter = np.loadtxt(filter)
    weightf=interpolate.interp1d(ifilter[:,0],ifilter[:,1],fill_value=0,bounds_error=False)
    weight=[]
    for wavelenght in wl:
        weight.append(weightf(wavelenght))

    return weight

def rotate(parangle,hdu):
    """Rotates fits file to match the rotation of the data cube.

    Parameters
    ----------
    parangle : float
        Paralactic angle of the cube.
    hdu : fits file

    Returns
    -------
    hdu : fits file
        Rotated fits file.
    """
    mat = [[hdu[0].header['CD1_1'],hdu[0].header['CD1_2']],[hdu[0].header['CD2_1'],hdu[0].header['CD2_2']]]
    CDELT1 = pow(abs(np.linalg.det(mat)),0.5) #Coordinate increment per pixel in DEGREES/PIXEL
    sdssangle = np.arcsin(hdu[0].header['CD1_1']/CDELT1)*180./np.pi

    hdu[0].data = interpolation.rotate(hdu[0].data,(parangle-sdssangle))
    hdu[0].header['CD1_1']= np.sin(parangle*np.pi/180.)*CDELT1
    hdu[0].header['CD1_2']= np.cos(parangle*np.pi/180.)*CDELT1
    hdu[0].header['CD2_1']= np.cos(parangle*np.pi/180.)*CDELT1
    hdu[0].header['CD2_2']= -np.sin(parangle*np.pi/180.)*CDELT1
    hdu[0].header['CRPIX1']= hdu[0].data.shape[1]/2+0.5
    hdu[0].header['CRPIX2']= hdu[0].data.shape[0]/2+0.5
    hdu.writeto('imaging/LsdssDR12g.fits',clobber =True)

    return hdu

def leastsq(sdsscrop,virusdata):
    leastsq = np.zeros((124,124))
    for i, row in enumerate(sdsscrop[62:-62,62:-62]):
        for j, element in enumerate(row):
            for k,rowk in enumerate(virusdata):
                for l, elementk in enumerate(rowk):
                    leastsq[i,j]= leastsq[i,j]+(sdsscrop[i-62+k,j-62+l] - virusdata[k,l])**2
    return leastsq

def crop(hdu,ifu_ra,ifu_dec,size):
    """Crops fits file by twice the size of the data cube around a given ra and dec.

    Parameters
    ----------
    hdu : fits file
        Fits file to crop.
    ifu_ra : float
        RA center of the square to crop
    ifu_dec : float
        DEC center of the square to crop
    size : int
        Size in pixels

    Returns
    -------
    hdu : fits file
        Cropped fits file.
    """

    c=SkyCoord(ifu_ra,ifu_dec, unit="deg")
    size = (size,size)     # cutting a box twice the size of the ifu
    w=wcs.WCS('imaging/LsdssDR12g.fits')
    cutout = Cutout2D(hdu[0].data, c, size, wcs=w)
    #plt.imshow(hdu[0].data, origin='lower')
    #cutout.plot_on_original(color='white')
    #plt.show()
    hdu = fits.PrimaryHDU(cutout.data,header=cutout.wcs.to_header()) #create new hdu
    hdulist = fits.HDUList([hdu])
    return hdulist

def flatten_wspectracube(datacube,weight):
    """Weights each spectra slice of the cube by a given array of weigths.

    Parameters
    ----------
    datacube : 3darray
        Cube with the first dimension being the spectra.
    weight : 1darray
        Array of weights, it's dimension must be the same as that of the first dimension of the datacube.

    Returns
    -------
    datacube : 3darray
        Weighted datacube.
    """
    i=0
    for wv in datacube:
        wv = wv*weight[i]
        i+=1
    datacube = datacube.sum(axis=0)

    return datacube


def resize(hdu,cube):
    """Resizes cube to fit for pixel size of hdu.

    Parameters
    ----------
    hdu : fits file
        Fits file with a given pixel with to resize the cube.
    cube : 2darray
        Cube. The standard pixels with is assumed here. 1''

    Returns
    -------
    cube : 2darray
        Zoomed cube.
    """
    pixwidth=1./3600 #cube has one arsec pixel width; this is pixels per degree
    resizefactor = degtopix(hdu,pixwidth)
    cube = zoom(cube,resizefactor)

    return cube



def parseargs(argv=None):
    '''
    Options parser
    '''

    parser = argparse.ArgumentParser(description = "Calculate astrometry based on Imaging")


    parser.add_argument("basename",
                        help="basename of fits FILEs to extract spectra from")
    parser.add_argument("--sig",dest ="sig",default=8,type=float, help="sigma value to cut the stars")
    args = parser.parse_args(args=argv)

    args.fitsfiles = []
    args.position = []
    #Searching for files with that name...
    searchname = args.basename + '*.fits'
    filenames = glob.glob(searchname)
    if not filenames:
        msg = 'No files found searching for: {:s}'.format(searchname)
        parser.error(msg)
    else:
        for fn in filenames:
            args.fitsfiles.append(fits.open(fn)[0])
            args.position.append(fn[-12:-9])

    return args



def main():
    args = parseargs()
    fitsfiles = args.fitsfiles
    #sdss_image = args.sdss_image
    hdu = fits.open(sdss_image)
    sigmapix = degtopix(hdu,sigma)
    
    output=open('positions.dat','w')
    output.write('{0:14} {1:14} {2:14} \n'.format("#IFU","RA","DEC"))
    ifu_cen=np.loadtxt('../greg/ifuPos.txt')

    #integrating over wavelengths
    filter = 'filters/g_sdss.dat'
    xmin = 3500
    x_range = 1054
    step = 1.9
    weight = filterweights(filter,xmin,x_range,step)

    sig= args.sig
    signoise = 0.3
    #resizing sdss image to cube's size
    virus_data = []
    for i, fitsfile in enumerate(fitsfiles):
        virus_data.append(flatten_wspectracube(fitsfile.data,weight))
        fitsfile.data = virus_data[i]
        virus_data[i] = resize(hdu,virus_data[i])
        #virus_data[i] = gaussian_filter(virus_data[i],sigmapix*2)
        fitsfile.writeto('imaging/cube'+args.position[i]+'.fits',clobber =True)

    #Convolving the image with a gaussian instead of convolving the kernel with a gaussian, since the kernel is to small.
    hdu[0].data = gaussian_filter(hdu[0].data,sigmapix)

    ckern = cirkernel(1.89553)
    hdu[0].data = convolve2d(hdu[0].data,ckern,mode='same')
    hdu = rotate(parangle, hdu)
    sdsshdu = [i.copy() for i in hdu]



    k=0
    for ifu in ifu_cen:
        virusdata = virus_data[k]
        ifu_id,ifu_ra, ifu_dec = ifu[0],ifu[1],ifu[2]
        
        hducropped = crop(sdsshdu, ifu_ra, ifu_dec, 124*2)
        #hdu = crop(hdu, ifu_ra, ifu_dec, 124*2)
        #hdu=hdu[0]
        #virusdata = hdu.data[124/2+20:-124/2+20,124/2+20:-124/2+20]
        sdssimgcrop = hducropped[0].data#+np.random.randn(*hdu.data.shape) * 1
        

        
        
        #substracting by the mean
        virusdata = virusdata - virusdata.mean()
        sdssimgcrop = sdssimgcrop - sdssimgcrop.mean()
        #virusdata = virusdata + np.amin(virusdata)
        #sdssimgcrop = sdssimgcrop + np.amin(sdssimgcrop)
        sdssimgcrop = sdssimgcrop/np.sum(np.absolute(sdssimgcrop))
        virusdata = virusdata/np.sum(np.absolute(virusdata))
        #sdssimg = np.copy(sdssimgcrop)
        #sdssimg_mean = sdssimg[below_threshold].mean()
        above_threshold = sdssimgcrop > np.std(sdssimgcrop)*sig
        sdssimgcrop[above_threshold] = np.std(sdssimgcrop)*sig
        below_threshold = virusdata < np.std(virusdata)*signoise
        virusdata[below_threshold] = 0
        #correlating
        sdssshift = correlate2d(sdssimgcrop**2, np.ones(virusdata.shape),
                                boundary='fill', mode='same')
        virusshift = correlate2d(sdssimgcrop, virusdata, boundary='fill', mode='same')
        fit = virusshift-sdssshift
        fit = fit[62:-62,62:-62]
        y, x = np.unravel_index(np.argmax(fit),fit.shape)
        y, x = y + 62, x + 62
        print(y,x)
        w = wcs.WCS(hducropped[0].header)
        new_ra, new_dec= w.wcs_pix2world([[x,y]],1)[0]
        #hdu[0].data=sdssimgcrop
        hducropped2 = crop(sdsshdu, new_ra, new_dec,int(124*2))
        sdssimgcrop = hducropped2[0].data
        '''
        rotations=[]
        
        #for angle in np.linspace(-4,4,num=9):
        for angle in [0]:
        #plt.imshow(sdssshift)
            #plt.show()
            #sdssimgcrop_rot = interpolation.rotate(sdssimgcrop,angle)
            #sdssshift_rot = interpolation.rotate(sdssshift,angle)
            cs = len(sdssshift_rot)/2
            cv = len(sdssimgcrop_rot)/2
            #cut=(len(sdssshift_rot)-len(sdssimgcrop_rot))/2
            sdssshift_rot = sdssshift_rot[cs-cv:cs+cv,cs-cv:cs+cv]
            #plt.imshow(sdssshift_rot)
            #plt.show()
            virusshift_rot = correlate2d(sdssimgcrop_rot, virusdata, boundary='fill',fillvalue=1, mode='same')
            #virusshift = Cutout2D(virusshift,(cv[0]/2,cv[0]/2),int(124*1.5), wcs=w).data
            #sdssshift_rot = Cutout2D(sdssshift_rot,(cs[0]/2,cs[0]/2),int(124*1.5), wcs=w).data
            #plt.imshow(sdssshift_rot)
            #plt.show()
            
            #diffvi = len(virusshift)-int(124*1.5)
            #diffsd = len(sdssshift_rot)-int(124*1.5)
            #fit = virusshift[cv-diffvi/2:cv+diffvi/2,cv-diffvi/2:cv-diffvi/2]-sdssshift_rot[cs+diffsd/2:cs-diffsd/2,cs+diffsd/2:cs-diffsd/2]
            fit = virusshift_rot - sdssshift_rot
            #fit = fit[62:-62,62:-62]
            print(fit.shape)
            center = len(fit)
            y, x = np.unravel_index(np.argmax(fit[center-62:center+62,center-62:center+62]),fit[center-62:center+62,center-62:center+62].shape)
            y, x = y + 62, x + 62
            
            rotations.append([np.amax(fit),angle,y,x])
        rotations = np.array(rotations)
        maxrot = np.argmax(rotations[:,0])
        rotation = rotations[maxrot, 1]
        y, x = rotations[maxrot, 2], rotations[maxrot, 3]
        print(y, x, rotation)
        ra, dec= w.wcs_pix2world([[x,y]],1)[0]
        #hdu[0].data=sdssimgcrop
        hducropped = crop(hdu, ra, dec,124*2)
        
        i#plt.imshow(hducropped[0].data)
        #plt.show()
        #least = leastsq(sdssimgcrop,virusdata)
        #y, x = np.unravel_index(np.argmax(least),(124,124))
        #print(y,x)
        #y,x= y+62,x+62
        '''

        print('IFU %d IS CENTERED AT:' % ifu_id)
        print(new_ra,new_dec)
        print('The offset is:')
        print(new_ra-ifu_ra,new_dec-ifu_dec)
        hducropped[0].data=sdssimgcrop
        hducropped[0].data[y,x]=np.amax(hducropped[0].data)
        for i in range(124):
            hducropped[0].data[y+62,x+i-62]=np.amax(hducropped[0].data)
            hducropped[0].data[y-62,x+i-62]=np.amax(hducropped[0].data)
        for j in range(124):
            hducropped[0].data[y+j-62,x+62]=np.amax(hducropped[0].data)
            hducropped[0].data[y+j-62,x-62]=np.amax(hducropped[0].data)
        #plt.imshow(hducropped[0].data)
        #plt.show()
        hducropped.writeto(('imaging/GBsdssDR12g%d.fits')%ifu_id,clobber =True)
        k += 1
        output.write('{0:14} {1:14} {2:14} \n'.format(ifu_id,(new_ra-205.658155675)*3600,(new_dec-28.3641961267)*3600))
        #output.write('{0:14} {1:14} {2:14} \n'.format(ifu_id,x,y))
        sdsshdu = [i.copy() for i in hdu]

    #hdu.writeto('imaging/test.fits',clobber =True)
    #print(corr)
    output.close()


if __name__ == "__main__":
    main()
