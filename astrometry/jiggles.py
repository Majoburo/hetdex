#!/usr/bin/env python

from __future__ import division
import sys
import os
import glob
import subprocess
from telarchive import fetchsdss
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os.path as op
import tools
from scipy import interpolate
from astropy import wcs
from astropy.io import fits
ppath,f = os.path.split( os.path.realpath(__file__) )
sys.path.append(ppath)
DEBUG = True
sdss_fits_fname = 'imaging/sdssDR12g.fits'
pixCrd = 'offsets/pixCrd.txt'
if(os.getenv('YODASRC')==None):
    os.environ['YODASRC'] = os.environ['WORK']+"/yoda/src"
if(os.getenv('FETCHSDSS')==None):
    os.environ['FETCHSDSS'] = os.environ['WORK']+"/fetchsdss"

def findchi2(vifu,dssifu,evifu,edssifu):

    '''
    Evaluate chi-squared between two distributions.
    '''
    chi2 = 0.0
    for n in range(len(vifu)):
        residual = (vifu[n] - dssifu[n])
        chi2 = chi2 + residual*residual/pow((evifu[n]*edssifu[n]),0.5)
    return chi2

def wcs2pix(fiber_ra,fiber_dec,fitsfile):

    '''
    Converting ra and dec into pixels (mainly to fit yoda's input)

    '''
    with fits.open(fitsfile) as h:
        img_data = h[0].data
        w = wcs.WCS(h[0].header)
    pixcrd = w.wcs_world2pix(np.array(zip(fiber_ra,fiber_dec),dtype='float64'),1)
    pixcrd = [[i+1,pixcrd[i].tolist()[0],pixcrd[i].tolist()[1]] for i in range(len(pixcrd[:,0]))]
    np.savetxt(pixCrd,pixcrd)

    return

def zoom(positions,virus_flux,zooms=10):

    '''
    Routine to zoom in the jiggles.

    '''
    for z in range(zooms):

        z=z+1
        ddec=round(0.001000/z,6)
        dra=round(0.001000/z,6)
        chi2min = jiggle(positions,virus_flux,ddec,dra)
        positions = np.array([chi2min[1],chi2min[2]])

    return chi2min

def jiggle(positions,virus_flux,ddec=0.001000,dra=0.001000,steps = 5):
    '''
    Routine that jiggle the ifu and calculates the minimun displacement

    INPUT:
        positions   (numpy.array[2]) Position of fibers     (units:sexagesimal)
        virus_flux  (numpy.array[1]) Flux of each fiber in hetdex IFU.
        steps       (int)            Steps to jiggle in.
        ddec        (float)          Lenght of steps in dec (units:sexagesimal)
        dra         (float)          Lenght of steps in ra  (units:sexagesimal)

    OUTPUT:
        chi2min     ([chi2,ratemp,dectemp,dss_flux,s1,s2])
            
            int             chi2           Minimun chi 2 computed
            numpy.array[1]  ratemp         Ras for the fibers with min chi2
            numpy.array[1]  dectemp        Decs for the fibers with min chi2
            rumpy.array[1]  dss_flux       Flux for the fibers with min chi2
            int             s1             Step in the grid of decs with min chi2
            int             s2             Step in the grid of ras with min chi2


    Shots will be jiggled over a range of steps
    with 5 steps : (where x marks the shot center)
    Initial:     Movements:
    o o o o o -> 1 2 3 4 5
    o o o o o -> 6 7 8 9 10
    o o x o o -> and so on...
    o o o o o ->
    o o o o o ->
    '''
    
    chi2min=[1000000000000000000]  ##Gotta find a better way to enforce this...
    ra0=positions[0]
    dec0=positions[1]
    chi2pos=[]
    dss_flux=[]

    '''
    For debugging purposes and due to normalization problems, i'll compare photometry of in sdss itself
    
### This starts debbuging
    #print ra0[0],dec0[0]
    wcs2pix(ra0,dec0,sdss_fits_fname)
    dssifu = photometry()
    virus_flux = dssifu[0]
    evirus_flux = dssifu[1]
### This ends debbuging
    '''
    evirus_flux = np.ones(len(virus_flux)) #just to test, without errors
    #f, axarr = plt.subplots(5, 5)

    dectemp=dec0-ddec*steps/2 #min value of dec to scan
    for s1 in range(steps):
        ratemp=ra0-dra*steps/2 #min value of ra to scan
        for s2 in range(steps):
        #MAIN ROUTINE
            #print "For %3.6f %2.6f:"%(ratemp[0]+dra,dectemp[0]+ddec)
            wcs2pix(ratemp,dectemp,sdss_fits_fname) #converting ra and dec to pixels in sdss image for yoda input
            
            #print "     Calculating photometry."
            dssifu = tools.photometry(pixCrd,sdss_fits_fname)
            dss_flux = dssifu[0]
            edss_flux = dssifu[1]

            #print ratemp[0],dectemp[0]        
            #print dss_flux.sum()
            #print virus_flux.sum()
            chi2 = findchi2(virus_flux,dss_flux,evifu=evirus_flux,edssifu=edss_flux)

            #np.savetxt('debug/jiggled_data_%s_%s.cat'%(s1,s2),map(list,zip(*[ratemp,dectemp])),fmt=['%3.6f','%2.6f'] )
            
            #print('chi2 = %f, virus_flux = %f,dss_flux= %f, s1=%d,s2=%d'%(chi2,virus_flux.sum(),dss_flux.sum(),s1,s2))
            #print "     Searching for min chi2"
            if (chi2min[0] > chi2):
                chi2min = [chi2,ratemp,dectemp,dss_flux,s1,s2]
            '''    #print('chi2min = %f'%(chi2min[0]))
            axarr[s1, s2].scatter(virus_flux,dss_flux)
            axarr[s1, s2].set_title('offset %d %d'%(s1,s2))
            axarr[s1, s2].axis(xmin=-0.1,xmax=0.1,ymax=0.5,ymin=-0.5)
            '''
            ratemp=ratemp+dra #step in ra
        dectemp=dectemp+ddec #step in theta
    #plt.show()
    #import scipy.optimize as optimization
    #def linfun(x,a,b):
    #    return a+b*x
    #optimization.curve_fit(linfun, vifu, ifu,x0, sigma)
    return chi2min

def wavelenghtrange(hdulist):
    """
    INPUT:
        Fits file

    OUTPUT:
        x = [Wavelength range]
    """
    xmin = hdulist.header['CRVAL1']
    xdelta = hdulist.header['CDELT1']
    x_range = hdulist.header['NAXIS1']
    xmax = xmin + x_range*xdelta
    wl = np.linspace(xmin,xmax,x_range)
    return wl

def weigthspectra(data,x):
    '''
    Make a weight based on some filter file and return filter spectra
    '''
    ifilter = np.loadtxt(ppath+'/filters/g_sdss.dat')
    weightf=interpolate.interp1d(ifilter[:,0],ifilter[:,1],fill_value=0,bounds_error=False)
    weight=[]
    for wavelenght in x:
        weight.append(weightf(wavelenght))
    wdata=np.zeros((len(data),len(x)))

    j=0
    for fiber in data:
        i=0
        for wavelenght in x:
            wdata[j,i]=fiber[i]*weight[i]
            i=i+1
        j=j+1
    return wdata


def parseargs(argv=None):
    '''
    Options parser
    '''
    #parser = OptionParser()

    parser = argparse.ArgumentParser(description = "Calculate astrometry based on Imaging")


    parser.add_argument("basename",nargs='?', type=str,action="store",
                              help="basename of fits FILEs to extract spectra from")
    parser.add_argument("ifuPos",nargs='?', type=str,
                        help='File with initial guess of ifu positions in the sky')
    parser.add_argument("--xmin",dest ="xmin",default=4800.13,type=float, help="xmin value (in A)")
    parser.add_argument("--xmax",dest ="xmax",default=5469.87,type=float, help="xmax value (in A)")
    args = parser.parse_args(args=argv)

    if args.basename is None:
        msg = 'The base name was not provided'
        parser.error(msg)
    else:

        args.fitsfiles = []
        searchname = args.basename + '*.fits'
        filenames = glob.glob(searchname)
        if not filenames:
            msg = 'No files found searching for: {:s}'.format(searchname)
            parser.error(msg)
        else:
            args.fitsfiles = [pyfits.open(filenames[i])[0] for i in xrange(len(filenames))] 
    
    return args



def main():
    args = parseargs()

    print """Jiggles BETA version
             --Algorithm used in VENGA astrometry--"""

    fitsfiles = args.fitsfiles
    print "Calculating wavelength range from virus files."
    wl = wavelenghtrange(fitsfiles[0])
    
    print 'Weighting spectral data from IFU using sdss g filter and integrating over wavelenght range for each fiber...'
    data=[]
    data = [spectra for fitsfile in fitsfiles for spectra in fitsfile.data]
    wdata = weigthspectra(data,wl)
    virus_flux = [x.sum() for x in data]

    print 'Retrieving IFU position that were created by Greg\'s visualization_tool.py.'

    #Getting ifu positions for all the fibers in one ifu. The input file is created by greg's program.
    positions = tools.get_txt_data(args.ifuPos,[0,1])

    print 'Using fetchsdss to get image data from SDSS-III DR12.'
    tools.getsdssimg(positions[0,0],positions[1,0])

    jiggled_data_min = zoom(positions,virus_flux)
    
    print "RA is off by %3.6f"%(positions[0,0]-jiggled_data_min[1][0])
    print "DEC is off by %2.6f"%(positions[1,0]-jiggled_data_min[2][0])
    #print('Best fit is plot (%d,%d)'%(jiggled_data_min[4],jiggled_data_min[5]))
    np.savetxt('debug/jiggled_data_min_%d_%d.cat'%(jiggled_data_min[4],jiggled_data_min[5]),map(list,zip(*[jiggled_data_min[1],jiggled_data_min[2]])))
    #dss_flux = jiggled_data_min[3]

if __name__ == "__main__":
    main()
