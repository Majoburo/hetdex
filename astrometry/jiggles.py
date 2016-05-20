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
from scipy import interpolate
from astropy import wcs
from astropy.io import fits
ppath,f = os.path.split( os.path.realpath(__file__) )
sys.path.append(ppath)
DEBUG = True
sdss_fits_fname = 'imaging/frame-g-002326-3-0078.fits'
#sdss_fits_fname = 'imaging/J134225.00+282357.0-g.fits'
pixCrd = 'offsets/pixCrd.txt'
if(os.getenv('YODASRC')==None):
    os.environ['YODASRC'] = os.environ['WORK']+"/yoda/src"
if(os.getenv('FETCHSDSS')==None):
    os.environ['FETCHSDSS'] = os.environ['WORK']+"/fetchsdss"

def getsdssimg(ra,dec):
    '''
    Getting imaging from SDSS DR12.
    
    INPUT:

        float   ra      RA to get image from
        float   dec     DEC to get image from

    '''
    subprocess.call([os.environ['FETCHSDSS']+'/do_fetchsdss.py',"g", "--coords=%s %s"%(ra,dec),"--nodr7","--getdr12","--nosuffix","--output=imaging/sdssDR12"])
    cmd = "bunzip2 imaging/sdssDR12g.fits.bz2 "
    os.system(cmd)
    sdss_fits_fname = 'imaging/sdssD12g.fits'


    return


def photometry():

    '''
    Calling YODA to do photometry on the fibers.
    
    Parameters:
        - A: pixel diameter of the fibers. Hetdex has fibers of 2.5 arsec, this converts in pixel space roughly to 5 pixels.

    '''
    cmd = "$YODASRC/yoda -P --no-kron-ap -p imaging/image.phot -M %s -A 5  %s &> /dev/null" % (pixCrd,sdss_fits_fname)
    #cmd = "$YODASRC/yoda -P --no-kron-ap -p imaging/image.phot -M %s -A 5  %s" % (pixCrd,sdss_fits_fname)
    #print "> " + cmd
    os.system(cmd)

    dssifu= get_txt_data('imaging/image.phot',[12,13])
    return dssifu

def get_txt_data(txtfile,columns):

    '''
    Get data from txt file.
    INPUT:
        txtfile: location of txt file (type:string)
        columns: columns to get  (type:ints)
    OUTPUT:
        np.array[2] data    [[column],[column]...]

    '''
    data=[]
    f = np.loadtxt(txtfile)
    for i in columns:
        data.append(f[:,i])

    data = np.array(data)
    return  data


def findchi2(vifu,dssifu,evifu,edssifu):

    '''
    Evaluate chi-squared between two distributions.
    '''
    chi2 = 0.0
    for n in range(len(vifu)):
        #print('%f,%f,%f,%f'%(vifu[n],dssifu[n],evifu[n],edssifu[n]))
        residual = (vifu[n] - dssifu[n])/(evifu[n]*edssifu[n]) #is this right???
        #residual = vifu[n] - dssifu[n] #without errors
        chi2 = chi2 + residual*residual
    return chi2

def wcs2pix(fiber_ra,fiber_dec,fitsfile=sdss_fits_fname):

    '''
    Converting ra and dec into pixels (mainly to fit yoda's input)

    '''
    with fits.open(fitsfile) as h:
        img_data = h[0].data
        w = wcs.WCS(h[0].header)
    pixcrd = w.wcs_world2pix(np.array(zip(fiber_ra,fiber_dec),dtype='float64'),1)
    #pixcrd = map(list,enumerate(pixcrd))
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
    '''
### This starts debbuging
    #print ra0[0],dec0[0]
    wcs2pix(ra0,dec0)
    dssifu = photometry()
    virus_flux = dssifu[0]
    evirus_flux = dssifu[1]
### This ends debbuging

    #f, axarr = plt.subplots(5, 5)

    dectemp=dec0-ddec*steps/2 #min value of dec to scan
    for s1 in range(steps):
        ratemp=ra0-dra*steps/2 #min value of ra to scan
        for s2 in range(steps):
        #MAIN ROUTINE
            wcs2pix(ratemp,dectemp)
            dssifu = photometry()
            dss_flux = dssifu[0]
            edss_flux = dssifu[1]

            #print ratemp[0],dectemp[0]
            
            #print dss_flux.sum()
            #print virus_flux.sum()
            chi2 = findchi2(virus_flux,dss_flux,evifu=evirus_flux,edssifu=edss_flux)

            np.savetxt('debug/jiggled_data_%s_%s.cat'%(s1,s2),map(list,zip(*[ratemp,dectemp])),fmt=['%3.6f','%2.6f'] )
            
            #print('chi2 = %f, virus_flux = %f,dss_flux= %f, s1=%d,s2=%d'%(chi2,virus_flux.sum(),dss_flux.sum(),s1,s2))
            
            if (chi2min[0] > chi2):
                chi2min = [chi2,ratemp,dectemp,dss_flux,s1,s2]
            '''    #print('chi2min = %f'%(chi2min[0]))
            axarr[s1, s2].scatter(virus_flux,dss_flux)
            axarr[s1, s2].set_title('offset %d %d'%(s1,s2))
            axarr[s1, s2].axis(xmin=-0.1,xmax=0.1,ymax=0.5,ymin=-0.5)
            '''
            ratemp=ratemp+dra #step in ra
        dectemp=dectemp+ddec #step in theta
    #  plt.show()
    #dssifu = get_flux(sdss_fits_fname,phi,theta,1,1) #calling function that given an ra and dec will give u flux in sloan image centered there
    #import scipy.optimize as optimization
    #def linfun(x,a,b):
    #    return a+b*x
    #optimization.curve_fit(linfun, vifu, ifu,x0, sigma)
    return chi2min

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

def getsdssimageold():
    '''
    Get SDSS-II image for the field using ds9 tools. Will have to be replaced by a SQL query for SDSS-III DR12
    '''

    f = open(ifuPos,'r')
    f.readline()
    ra,dec = f.readline().split()
    f.close()
    import pyds9
    ds9 = pyds9.DS9()
    ds9.set('dsseso size %f %f degrees' % (40./60., 40./60.))
    ds9.set('frame delete all')
    ds9.set('frame new')
    ds9.set('dsseso coord %f %f degrees' % (ra, dec))
    ds9.set('dsseso close')
    ds9.save('sdss2.fits')

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


def main():
    args = parseargs()

    print """Jiggles BETA version"""

    fitsfiles = args.fitsfiles
  # if args.verbose:
  #     print "Calculating wavelenght range \n"
    wl = wavelenghtrange(fitsfiles[0])
    
    'Handeling spectral data from IFU'
    data=[]
    data = [spectra for fitsfile in fitsfiles for spectra in fitsfile.data]
    wdata = weigthspectra(data,wl)
    virus_flux = [x.sum() for x in data]
    #print np.sum(virus_flux)


    'Handeling imaging data'

    #Getting ifu positions for all the fibers in one ifu. The input file is created by greg's program.
    positions = get_txt_data(args.ifuPos,[0,1])
    getsdssimg(positions[0,0],positions[1,0])
    jiggled_data_min = zoom(positions,virus_flux)
    print('Best fit is plot (%d,%d)'%(jiggled_data_min[4],jiggled_data_min[5]))
    np.savetxt('debug/jiggled_data_min_%d_%d.cat'%(jiggled_data_min[4],jiggled_data_min[5]),map(list,zip(*[jiggled_data_min[1],jiggled_data_min[2]])))
    dss_flux = jiggled_data_min[3]

if __name__ == "__main__":
    main()
