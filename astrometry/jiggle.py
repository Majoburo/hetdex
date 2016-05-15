#!/usr/bin/env python

import pyfits
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser
import h5py
from integrate import get_flux
import sys
import os
import os.path as op
from scipy import interpolate
ppath,f = os.path.split( os.path.realpath(__file__) )
sys.path.append(ppath)

sdss_fits_fname = 'imaging/frame-g-002326-3-0078.fits'
ifuPos = 'offsets/ifuPos.txt'
f = open(ifuPos,'r')
f.readline()
f.readline()
ra0 , dec0 = [],[]
for line in f:
    a,b = line.split()
    ra0.append(float(a))
    dec0.append(float(b))




def chi2(vifu,dssifu,evifu=1,edssifu=1):
# Evaluate chi-squared.
    chi2 = 0.0
    for n in range(len(vifu)):
        residual = (vifu[n] - dssifu[n])/(evifu[n]*edssifu[n])
        chi2 = chi2 + residual*residual
    
    return chi2


def optimizelin(vifu,dssifu):

   # Shots will be jiggled over a range of steps
   # with 5 steps : (where x marks the shot center)
   # Initial:     Movements:
   # o o o o o -> 1 2 3 4 5
   # o o o o o -> 6 7 8 9 10
   # o o x o o -> and so on...
   # o o o o o ->
   # o o o o o ->
#   (sdss_fits_fname,fiber_ra,fiber_dec,fiber_dia,zoom_factor,errors=False)
   steps = 5
   ddec=0.01
   dra=0.01
   theta=dec0-ddec*steps/2 #min value of theta to scan
   for s1 in range(steps):
       theta=theta+ddec #step in theta
       phi=ra0-dra*steps/2 #min value of phi to scan
       for s2 in range(steps):
           phi=phi+dra #step in phi
           print theta,phi
           dssifu,edssifu = get_flux(sdss_fits_fname,phi,theta,1,1) #calling function that given an ra and dec will give u flux in sloan image centered there
           chi2list = [chi2(vifu,dssifu),phi,theta]
       chi2,phi,theta=np.amin(chi2list,0)
    #optimizelin()
    #with  phi theta as new centers repeat
    #optimizelin()
    #when finished fit a line (or maybe every time)
    #import scipy.optimize as optimization
    #def linfun(x,a,b):
    #    return a+b*x
    #optimization.curve_fit(linfun, vifu, ifu,x0, sigma)



    


def weigthspectra(data,x):
    #make a weight based on some filter file and return filter spectra
    ifilter = np.loadtxt(ppath+'/filters/g_sdss.dat')
    weightf=interpolate.interp1d(ifilter[:,0],ifilter[:,1],fill_value=0,bounds_error=False)
    weight=[]
    for wavelenght in x:
        weight.append(weightf(wavelenght))
    wdata=np.zeros((len(data[:,0]),len(x)))

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
    #options parser
    parser = OptionParser()
    parser.add_option("-base", dest="basename",action="store",
                              help="basename of fits FILEs to extract spectra from")
    parser.add_option("--xmin",dest ="xmin",default=4800.13,type=float, help="xmin value (in A)")
    parser.add_option("--xmax",dest ="xmax",default=5469.87,type=float, help="xmax value (in A)")

    (options, args) = parser.parse_args()

    if options.basename is None:
        msg = 'The base name was not provided'
        parser.error(msg)
    else:
        hdulist = []
        searchname = args.basename + '*.fits'
        filenames = glob.glob(searchname)
        if not filenames:
            msg = 'No files found searching for: {:s}'.format(searchname)
            parser.error(msg)
        else:
            for i in xrange(len(filenames)):
                hdulist[i]= pyfits.open(filenames[i])
#for now running only for one ifu
    return hdulist[0]

def getfitsdata(hdulist):
    #gets data fits file and returns wavelenght range plus the data (it assumes it's in log scale and changes it back)

  #  for i in hdulist[:,0]:
    data = hdulist[0].data
    xmin= np.exp(hdulist[0].header['CRVAL1'])
    xdelta= hdulist[0].header['CDELT1']
    x_range= hdulist[0].header['NAXIS1']
    xmax=xmin*np.exp(x_range*xdelta)
    x = np.linspace(xmin,xmax,x_range)
    hdulist.close()
    return data,x


def main():
    hdulist = parseargs()
    data,x = getfitsdata(hdulist)
    #getsdssimage()

    wdata = weigthspectra(data,x)
    vifu = [x.sum() for x in data]
    print data.shape,len(vifu)
if __name__ == "__main__":
    main()
