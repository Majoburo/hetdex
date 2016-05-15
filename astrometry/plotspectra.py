#!/usr/bin/env python

import pyfits
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser
import h5py

import sys
import os
from scipy import interpolate
ppath,f = os.path.split( os.path.realpath(__file__) )
sys.path.append(ppath)

def weigthspectra(data,x):
    #make a weight based on some filter file and return filter spectra
    ifilter = np.loadtxt(ppath+'/filters/g_sdss.dat')
    weightf=interpolate.interp1d(ifilter[:,0],ifilter[:,1],fill_value=0,bounds_error=False)
    weight=[]
    for wavelenght in x:
        weight.append(weightf(wavelenght))
    fdata=np.zeros((len(data[:,0]),len(x)))

    j=0
    for fiber in data:
        i=0
        for wavelenght in x:
            fdata[j,i]=fiber[i]*weight[i]
            i=i+1
        j=j+1
    return fdata

def getsdssimage(ra0,dec0):

    f = open('ifuPos.txt','w')
    f.readline()
    ra0,dec0 = f.readline().split()
    f.close()
    import pyds9
    ds9 = pyds9.DS9()
    ds9.set('dsseso size %f %f degrees' % (40./60., 40./60.))
    ds9.set('frame delete all')
    ds9.set('frame new')
    ds9.set('dsseso coord %f %f degrees' % (ra0, dec0))
    ds9.set('dsseso close')
    ds9.save('sdss2.fits')

def parseargs():
    #options parser
    parser = OptionParser()
    parser.add_option("-f", dest="filename",action="store",
                              help="fits FILE to extract spectra from")
    parser.add_option("--xmin",dest ="xmin",default=4800.13,type=float, help="xmin value (in A)")
    parser.add_option("--xmax",dest ="xmax",default=5469.87,type=float, help="xmax value (in A)")

    (options, args) = parser.parse_args()
    
    return options

def getdata(options):
    #gets data fits file and returns wavelenght range plus the data (it assumes it's in log scale and changes it back)

    hdulist = pyfits.open(options.filename)
    data = hdulist[0].data
    xmin= np.exp(hdulist[0].header['CRVAL1'])
    xdelta= hdulist[0].header['CDELT1']
    x_range= hdulist[0].header['NAXIS1']
    xmax=xmin*np.exp(x_range*xdelta)
    x = np.linspace(xmin,xmax,x_range)
    hdulist.close()
    
    return data,x

def main():
    options = parseargs()
    data,x = getdata(options)
    getsdssimage()

    print weigthspectra(data,x)
if __name__ == "__main__":
    main()
