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



#class shot(object):
#
#    def __init__(self,ifupositions,fitsfile):
#        self.ra0 = []
#        self.dec0 =[]
#    
#    center
#    wavelength range
#     xmin
#     xmax
#    fiber
#     flux
#     position
#    @classmethod
#    def center():
#    
#    @classmethod
#    def fiber(fibernumber):
#
#
#    @classmethod
#    def xrange():

#global parameters

j_steps=5 #jiggle steps

def jiggle(shot,dra,ddec,step):

   # Shots will be jiggled over a range of steps
   # with 5 steps : (where x marks the shot center)
   # Initial:     Movements:
   # o o o o o -> 1 2 3 4 5
   # o o o o o -> 6 7 8 9 10
   # o o x o o -> and so on...
   # o o o o o ->
   # o o o o o ->


    #initial jiggle parameters:


    tshot=shot

    for fiber in tshot:
        fiber[0]=fiber[0]+dra*step
        fiber[1]=fiber[1]+ddec*step
    
    return tshot

def chi2(vwifu,dssifu,sigma):
# Evaluate chi-squared.
    chi2 = 0.0
    for n in range(len(vwifu)):
        residual = (vwifu[n] - dssifu[n])/(vwsigma[n]dsssigma[n])
        chi2 = chi2 + residual*residual
    
    return chi2


def optimizelin():

    theta=fiber[1]-ddec*step/2 #min value of theta to scan
    for s1 in range(j_steps):
        theta=theta+ddec*s2 #step in theta
        phi=fiber[0]-dra*step/2 #min value of phi to scan
        for s2 in range(j_steps):
            phi=phi+dra*s1 #step in phi
            dssifu=getflux(phi,theta) #calling function that given an ra and dec will give u flux in sloan image centered there
            chi2list=chi2(vwifu,dssifu,sigma),phi,theta
    chi2min(chi2list)
    #with  phi theta as new centers repeat
    #optimizelin()
    #when finished fit a line (or maybe every time)
    import scipy.optimize as optimization
    def linfun(x,a,b):
        return a+b*x
    optimization.curve_fit(linfun, vwifu, ddsifu,x0, sigma)



    


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
    #getsdssimage()

    print weigthspectra(data,x)
if __name__ == "__main__":
    main()
