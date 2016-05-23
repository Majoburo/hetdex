To run jiggles you will need:

FECHSDSS (to fetch SDSS images corresponding to your shot)
YODA (for photometry on SDSS images)

To install FECHTSDSS:


To install yoda:
svn co svn://luna.mpe.mpg.de/yoda/trunk yoda
cd yoda
./bootstrap
./configure
make

Make sure to make and environment variable for their paths:

export YODASRC="your_path_to_yoda/yoda/src"

Then:

$ . runjiggles.sh

Setting files:
- ifuPos*.txt : Text file with the positions of each ifu.

Disclaimer:
Algorithm based on Guille's astrometry code for VENGA (found under the folder scripts).
This python version of that is in very beta version. Pretty much a squeleton of what will be.

    It runs photometry using yoda, an old (yet very fast) program Niv stopped updating a decade ago.

TO DO (big ones I can think of):

- Include errors on flux.
- Normalize both flux distributions.
- Mask bad fibers.


