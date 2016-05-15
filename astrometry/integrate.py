from astropy.io import fits
from astropy import wcs
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage.interpolation import map_coordinates

with open('../greg/ifuPos.txt','r') as f:
    f.readline()
    fiber_dia = float(f.readline())
    fiber_ra , fiber_dec = [],[]
    for line in f:
        a,b = line.split()
        fiber_ra.append(float(a))
        fiber_dec.append(float(b))
with fits.open("sdss2.fits") as h:
    img_data = h[0].data
    w = wcs.WCS(h[0].header)

pixcrd = w.wcs_world2pix(np.array(zip(fiber_ra,fiber_dec),dtype='float64'),1)

#plt.scatter(pixcrd[:,0],pixcrd[:,1])
#plt.show()
#img_data[np.round(pixcrd).astype(np.int)] = 0x000000

#plt.imshow(img_data,interpolation='nearest',cmap='Greys_r',origin='lower')
#plt.show()

fiber_r_pix = w.wcs_world2pix([[fiber_ra[0]+fiber_dia,fiber_dec[0]]],1)
fiber_r_pix -= w.wcs_world2pix([[fiber_ra[0],fiber_dec[0]]],1)
fiber_r_pix = 1 # TODO fix this!!!

# take sum over each fiber's pixels.
fibsum = np.zeros(pixcrd.shape[0])
# should also take fiber_r into account somehow. At the moment fiber_r is ~1 so this kind of works
map_coordinates(img_data,pixcrd.T,fibsum,order=3, mode='nearest')

print(fibsum)
