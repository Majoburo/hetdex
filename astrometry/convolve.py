import numpy as np
import pyfits
from scipy.ndimage.filters import gaussian_filter
from scipy import misc
import matplotlib.pyplot as plt
#from scipy.interpolate import interp2di
from scipy.signal import convolve2d
from scipy.interpolate import RectBivariateSpline
import tools

ifuPos = 'offsets/pixCrd.txt'
sdss_image= 'imaging/sdssDR12g.fits'

def cmask(position,radious,image):
    '''Pixelated circular mask'''
    a,b = position
    lenx,leny = image.shape
    y,x = np.ogrid[-a:lenx-a,-b:leny-b]
    mask = x*x+y*y <= radious*radious
    
    return mask

def drawtotalmask(img,positions):
    radious=5
    xindex,yindex=np.rint(positions)
    mask0 = np.zeros(np.shape(img))!=0
    for xind,yind in zip(xindex,yindex):
        mask = cmask([int(xind),int(yind)],radious,img)
        mask0 = mask0+mask

    return mask0 
                    
def cirkernel(radious):
    y,x = np.ogrid[-radious:radious+1,-radious:radious+1]
    mask = x*x+y*y<=radious*radious
    return mask*1


#face = misc.face(gray=True)
#face[drawtotalmask(face,pos)]=0
#plt.imshow(face)
#plt.show()
'''
print RectBivariateSpline(range(face.shape[0]),range(face.shape[1]),face).ev([50,51],[50,51])
sigma=3
blurred_face = gaussian_filter(face, sigma=3)
pos=[150.5,130.6]
radious=5
xindex,yindex=np.rint(pos)
x,y=pos
'''
#Square mask
#smask = np.ogrid[int(xindex)-radious:int(xindex)+radious,int(yindex)-radious:int(yindex)+radious]
#blurred_face[smask]=0

'''
    Iterate over positions to get mask of a bunch of circular regions
'''






pos=tools.get_txt_data(ifuPos,[0,1])
with pyfits.open(sdss_image) as hdu:
    #plt.imshow(hdu[0].data)
    #plt.show()
    #plt.imshow(drawtotalmask(hdu[0].data,pos)*1)
    #plt.show()
    ckern=cirkernel(2.5)
    hdu[0].data = convolve2d(hdu[0].data,ckern)
    
    #hdu[0].data[drawtotalmask(hdu[0].data,pos)]=1000000000
    #hdu[0].data = gaussian_filter(hdu[0].data,sigma=10)
    plt.imshow(hdu[0].data)
    plt.show()

    hdu.writeto('BsdssDR12g.fits',clobber =True)

