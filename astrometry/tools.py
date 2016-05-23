
#!/usr/bin/env python

import sys
import os
import glob
import subprocess

import numpy as np
import pyfits

from telarchive import fetchsdss
if(os.getenv('YODASRC')==None):
    os.environ['YODASRC'] = os.environ['WORK']+"/yoda/src"
if(os.getenv('FETCHSDSS')==None):
    os.environ['FETCHSDSS'] = os.environ['WORK']+"/fetchsdss"


def photometry(pixCrd,sdss_fits_fname):

    '''
    Calling YODA to do photometry on the fibers.
    
    Parameters:
        - A: pixel diameter of the fibers. Hetdex has fibers of 2.5 arsec, this converts in pixel space roughly to 5 pixels.

    '''
    cmd = "$YODASRC/yoda -P --no-kron-ap -p imaging/image.phot -M %s -A 5  %s &> /dev/null" % (pixCrd,sdss_fits_fname)
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


def getsdssimg(ra,dec):
    '''
    Getting imaging from SDSS DR12.
    
    INPUT:

        float   ra      RA to get image from
        float   dec     DEC to get image from

    '''
    print "The RA,DEC positions to search SDSS-III for imaging are: %3.6f, %2.6f"%(ra,dec)
    print "Calling fetchsdss to retrieve the field, please check that your IFU lies within it."
    print('Hit ENTER...')
    sys.stdin.readline()
    subprocess.call([os.environ['FETCHSDSS']+'/do_fetchsdss.py',"g","--nosuffix","--output=imaging/sdssDR12","--coords=%s %s"%(ra,dec),"--notables"])
    cmd = "gunzip -f imaging/sdssDR12g.fits.gz"
    os.system(cmd)
    import pyds9
    ds9 = pyds9.DS9()
    ds9.set('frame delete all')
    ds9.set('frame new')
    ds9.set('fits imaging/sdssDR12g.fits')
    ds9.set('regions load all ../greg/regions/temp.reg')
    print "Please verify the regions file falls entirely within the SDSS field."
    ref = sys.stdin.readline()
    ds9.set('exit')

    return
