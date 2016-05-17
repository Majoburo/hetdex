# -*- coding: utf-8 -*-
"""
@author: gregz
"""
from __future__ import print_function, absolute_import
import argparse as ap
import os.path as op
import os
import numpy
import sys
import pyfits
import glob
import pyhetdex.het.fplane as fplane_parser
from pyhetdex.het import ifu_centers as ifc
import astropy.units as units
from astropy.coordinates import SkyCoord
from astropy.coordinates import FK5

class regioninfo(object):
    @classmethod
    def writeHeader(cls, f):
        """Write the header to file ``f``

        Parameters
        ----------
        f : file-like object
            where to write to; must have a ``write`` method
        """
        s = []
        s.append("# global color=green dashlist=8 3 width=1 font=\"helvetica 8 normal roman\" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1")
        s.append("fk5")
        f.write('\n'.join(s) + "\n")
    
    @classmethod
    def writeifufiberpos(cls,ifuPos,ra,dec,size,status=1):
        """Write the fiber positions for each IFU to file ``ifuPos``

        Parameters
        ----------
        ifuPos : file-like object
            where to write to; must have a ``write`` method
        ra     : ra of fiber
        dec    : dec of fiber
        size   : diameter of fiber
        status : mask for each fiber
        """
        ifuPos.write('{:3.6f} {:2.6f} {} {} \n'.format(float(ra),float(dec),size,status))
        ifuPos.flush()

    @classmethod
    def writeregion(cls, f, ra, dec, size, name):
        """Write something to file ``f``

        Parameters
        ----------
        f : file-like object
            where to write to; must have a ``write`` method
        """

        s = ("circle(%3.6f,%2.6f,%1.2f\") # text={%s}" %
             (ra, dec, size, name))
        f.write(s)
        f.write("\n")
        f.flush()

def load_config(configfile):
    """Read configuration file

    Parameters
    ----------
    configfile : string
        name of the configuration file

    Returns
    -------
    config : :class:`~ConfigParser.ConfigParser` instance
    """
    import ConfigParser
    config = ConfigParser.ConfigParser()
    if not config.read(configfile):
        msg = ("ERROR: file '{}' does not exist. Please provide a valid"
               " configuration file")
        raise ap.ArgumentTypeError(msg.format(configfile))
    return config


def parse_args(argv=None):
    """Parse the command line arguments

    Parameters
    ----------
    argv : list of string
        arguments to parse; if ``None``, ``sys.argv`` is used

    Returns
    -------
    Namespace
        parsed arguments
    """
    description = "Shuffle the HETDEX shots"
    parser = ap.ArgumentParser(description=description,
                            formatter_class=ap.ArgumentDefaultsHelpFormatter)

    parser.add_argument("basename", nargs='?', type=str, 
                        help='''Filename to be mapped to the sky.
                        Ex. "Folder_Location/sci/Fepses20160409T090248.0"''')

    parser.add_argument("output", nargs='?', type=str, 
                        help='''Name of the region file (automatically appends .reg)
                        Ex. "standard_star".
                        Default: "temp"''', default='temp')

    parser.add_argument("-c", "--config", help="""Name of the configuration
                        file. When parsing the command line, the file is loaded
                        into a configuration object""",
                        default="./shuffle.cfg", type=load_config)

    args = parser.parse_args(args=argv)

    # get the positional arguments and check that they are either all ``None``
    # or not ``None``
    if args.basename is None:
        msg = 'The base name was not provided'
        parser.error(msg)
    else:
        searchname = args.basename + '*.fits'
        filenames = glob.glob(searchname)
        if not filenames:
            msg = 'No files found searching for: {:s}'.format(searchname)
            parser.error(msg)
        else:
            args.ifuslot = [op.basename(f).split('_')[1] for f in filenames]
            args.side = [op.basename(f).split('_')[-1].split('.fits')[0] 
                         for f in filenames]
            args.pa = []
            args.ra = []
            args.dec = []
            for i in xrange(len(filenames)):
                p = pyfits.open(filenames[i])
                rastr = p[0].header['TELRA']
                decstr = p[0].header['TELDEC']
                equinox = 2016 + (p[0].header['MJD'] - 57388)/365.24
                args.pa.append(p[0].header['PARANGLE'])
                ratemp = (float(rastr.split(':')[0]) * 15.
                           + float(rastr.split(':')[1]) * 15./60.
                           + float(rastr.split(':')[2]) * 15./3600.)
                dectemp = (float(decstr.split(':')[0])
                           + float(decstr.split(':')[1]) * 1./60.
                           + float(decstr.split(':')[2]) * 1./3600.)
                gc = SkyCoord(ratemp*units.degree, dectemp*units.degree,
                              frame='fk5', equinox='J{:0.4f}'.format(equinox))

                s = gc.transform_to(FK5(equinox='J2000.0'))
                args.ra.append(s.ra.deg)
                args.dec.append(s.dec.deg)
    return args

    
    
def get_ifuslot_list(args):
    """
    Get the ifu centers and ids of the VIRUS+LRS IFUs for the adjusted fplane file.
    """
    fplane_adj_file = args.config.get("General", "fplane_file")
    fplane = fplane_parser.FPlane(fplane_adj_file)
    ifu_centers = numpy.array([[ifu.x, ifu.y] for ifu in fplane.ifus])
    ifu_centers = ifu_centers / 3600.
    ifu_id = numpy.array([[ifu.ihmpid] for ifu in fplane.ifus])
    return ifu_centers , ifu_id    

def main():

    args = parse_args()

    if args.ra is None:
        # read the observing schedule from the file specified in configuration
        # file
        print ( '''No input arguments were given.
                    Example:
                    python visualize_observation.py FOLDER/sci/Fe''')
        return None
    else:
        #ra0 = args.ra[0]
        #dec0 = args.dec[0]
        #import pyds9
        #ds9 = pyds9.DS9()
        #ds9.set('height 1000')
        #ds9.set('width 1200')
        #ds9.set('dsseso size %f %f degrees' % (40./60., 40./60.))
        #ds9.set('frame delete all')
        #ds9.set('frame new')
        #ds9.set('dsseso coord %f %f degrees' % (ra0, dec0))
        #ds9.set('dsseso close')
        fiberdiam = 1.5 / 2.
        if not op.exists('regions'):
            os.mkdir('regions')
        if args.output[-4:]=='.reg':
            regionfile_fn = op.join('regions',args.output)
        else:
            regionfile_fn = op.join('regions',args.output+'.reg')

        print("Creating {:s}".format(regionfile_fn))
        regionfile = open(regionfile_fn, 'w')
        regioninfo.writeHeader(regionfile)  
       
        print("Creating ifuPos*.txt file")
        ifuPos = open("ifuPos%s.txt" %(args.ifuslot[0]),'w')

        ifu_centers, ifu_id = get_ifuslot_list(args)
        fiber_center = ifc.IFUCenter('IFUcen_HETDEX.txt')
        for i in xrange(len(args.ra)):
            ifu_ind    = numpy.array ( [ind for ind,ifuname in enumerate(ifu_id) if ifuname == args.ifuslot[i]], dtype = numpy.int)
            rpa = numpy.pi / 180. * (args.pa[i]-180.) 
            for j in xrange(len(fiber_center.xifu[args.side[i]])):
                dx = -1. * (ifu_centers[ifu_ind,0] + fiber_center.xifu[args.side[i]][j] / 3600.)
                dy = -1. * (ifu_centers[ifu_ind,1] + fiber_center.yifu[args.side[i]][j] / 3600.)
                dxr =       numpy.cos(rpa) * dx + numpy.sin(rpa) * dy
                dyr = -1. * numpy.sin(rpa) * dx + numpy.cos(rpa) * dy
                dra = dxr / numpy.cos(args.dec[i]*numpy.pi/180.)
                ddec = dyr
                ra_center  = args.ra[i]  + dra  # degrees
                dec_center = args.dec[i] + ddec # degrees
                regioninfo.writeregion(regionfile, ra_center,dec_center,fiberdiam,
                                       "{:03d}".format(int(fiber_center.fib_number[args.side[i]][i])))
                
                regioninfo.writeifufiberpos(ifuPos,ra_center,dec_center,fiberdiam)

        #regionfile.close()
        #ds9.set('scale asinh')
        #ds9.set('scale mode zmax')
        #ds9.set('regions system wcs')
        #ds9.set('regions sky fk5')
        #print('Hit ENTER...')
        #sys.stdin.readline()
        
if __name__ == '__main__':
    main()
