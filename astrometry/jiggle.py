#!/usr/bin/env python

import sys
import os
import glob

import pyfits
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser
import os.path as op
from scipy import interpolate
from astropy import wcs
from astropy.io import fits

ppath,f = os.path.split( os.path.realpath(__file__) )
sys.path.append(ppath)
DEBUG = True
sdss_fits_fname = 'imaging/frame-g-002326-3-0078.fits'
ifuPos = 'offsets/ifuPos075.txt'
pixCrd = 'offsets/pixCrd.txt'
os.environ['YODASRC'] = "../yoda/src"

def photometry():
    cmd = "$YODASRC/yoda -P --no-kron-ap -p imaging/image.phot -M %s -A 1.8   %s " % (pixCrd,sdss_fits_fname)
    #print "> " + cmd
    os.system(cmd)

    dssifu=[]
    phot=np.loadtxt('imaging/image.phot')
    dssifu=[[atribute[12],atribute[13]] for atribute in phot]
    return dssifu

def getifuPos(ifuPos):

    f = np.loadtxt(ifuPos)
    positions=[f[:,0].tolist(),f[:,1].tolist()]
    return  positions


def findchi2(vifu,dssifu,evifu=1,edssifu=1):
# Evaluate chi-squared.
    chi2 = 0.0
    for n in range(len(vifu)):
        #residual = (vifu[n] - dssifu[n])/(evifu[n]*edssifu[n]) #is this right???
        residual = vifu[n] - dssifu[n] #without errors
        chi2 = chi2 + residual*residual
    
    return chi2

def wcs2pix(fiber_ra,fiber_dec,fitsfile=sdss_fits_fname):

    with fits.open(fitsfile) as h:
        img_data = h[0].data
        w = wcs.WCS(h[0].header)
    pixcrd = w.wcs_world2pix(np.array(zip(fiber_ra,fiber_dec),dtype='float64'),1)
    pixcrd = [[i+1,pixcrd[i].tolist()[0],pixcrd[i].tolist()[1]] for i in range(len(pixcrd[:,0]))]
    np.savetxt(pixCrd,pixcrd)

    return



def jiggle(positions,virus_flux,steps = 5,ddec=0.001,dra=0.005):

   # Shots will be jiggled over a range of steps
   # with 5 steps : (where x marks the shot center)
   # Initial:     Movements:
   # o o o o o -> 1 2 3 4 5
   # o o o o o -> 6 7 8 9 10
   # o o x o o -> and so on...
   # o o o o o ->
   # o o o o o ->
    chi2min=0
    ra0=np.array(positions[0])
    dec0=np.array(positions[1])
    chi2pos=[]
    dectemp=dec0-ddec*steps/2 #min value of dec to scan
    for s1 in range(steps):
        ratemp=ra0-dra*steps/2 #min value of ra to scan
        for s2 in range(steps):
        #MAIN ROUTINE
            wcs2pix(ratemp,dectemp)
            dssifu = photometry()
            dss_flux = np.array(dssifu)[:,0]
            '''
            For debugging purposes and due to normalization problems, i'll compare photometry of in sdss itself
            '''
            print dss_flux
            virus_flux = dss_flux
            chi2 = findchi2(virus_flux,dss_flux)

            np.savetxt('debug/jiggled_data_%s_%s.cat'%(s1,s2),map(list,zip(*[ratemp,dectemp])),fmt=['%3.6f','%2.6f'] )
           
            if (chi2min == 0):
                chi2min = [chi2,ratemp,dectemp,dss_flux,s1,s2]
            elif (chi2min[0] > chi2):
                chi2min = [chi2,ratemp,dectemp,dss_flux,s1,s2]


            ratemp=ratemp+dra #step in ra
        dectemp=dectemp+ddec #step in theta
          

    #dssifu = get_flux(sdss_fits_fname,phi,theta,1,1) #calling function that given an ra and dec will give u flux in sloan image centered there
    #chi2list = [chi2(vifu,dssifu),phi,theta]
    #chi2,phi,theta=np.amin(chi2list,0)
    #optimizelin()
    #with  phi theta as new centers repeat
    #optimizelin()
    #when finished fit a line (or maybe every time)
    #import scipy.optimize as optimization
    #def linfun(x,a,b):
    #    return a+b*x
    #optimization.curve_fit(linfun, vifu, ifu,x0, sigma)
    return chi2min

def weigthspectra(data,x):
    #make a weight based on some filter file and return filter spectra
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

def getsdssimage():

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

def parseargs():
    '''
    Options parser
    '''
    parser = OptionParser()

    description = "Calculate astrometry based on Imaging"



    parser.add_option("--base", dest="basename",action="store",
                              help="basename of fits FILEs to extract spectra from")
    parser.add_option("--xmin",dest ="xmin",default=4800.13,type=float, help="xmin value (in A)")
    parser.add_option("--xmax",dest ="xmax",default=5469.87,type=float, help="xmax value (in A)")

    (options, args) = parser.parse_args()

    if options.basename is None:
        msg = 'The base name was not provided'
        parser.error(msg)
    else:

        fitsfiles = []
        searchname = options.basename + '*.fits'
        filenames = glob.glob(searchname)
        if not filenames:
            msg = 'No files found searching for: {:s}'.format(searchname)
            parser.error(msg)
        else:
            fitsfiles = [pyfits.open(filenames[i])[0] for i in xrange(len(filenames))] 
    
    return fitsfiles

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

    fitsfiles = parseargs()
    
    wl = wavelenghtrange(fitsfiles[0])
    
    'Handeling spectral data from IFU'
    data=[]
    data = [spectra for fitsfile in fitsfiles for spectra in fitsfile.data]
    wdata = weigthspectra(data,wl)
    virus_flux = [x.sum() for x in data]
    #print np.sum(virus_flux)


    'Handeling imaging data'
  #  getsdssimage()
    positions = getifuPos(ifuPos)
   #   for i in range(4)
    jiggled_data_min = jiggle(positions,virus_flux)
    
    np.savetxt('debug/jiggled_data_min_%d_%d.cat'%(jiggled_data_min[4],jiggled_data_min[5]),map(list,zip(*[jiggled_data_min[1],jiggled_data_min[2]])))
    dss_flux = jiggled_data_min[3]
#    plt.scatter(virus_flux,dss_flux)
#    plt.show()
if __name__ == "__main__":
    main()
