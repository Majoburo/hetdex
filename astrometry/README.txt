To run jiggles:

First install yoda:
svn co svn://luna.mpe.mpg.de/yoda/trunk yoda
cd yoda
./bootstrap
./configure
make

Then make an environment variable for it:

export YODASRC="your_path_to_yoda/yoda/src"

Then:

$ . runjiggles.sh

Setting files:
- ifuPos*.txt : Text file with the positions of each ifu.

Disclaimer:

Jiggles is in very very beta version. Pretty much a squeleton of what will be.

    It runs photometry using yoda, an old (yet very fast) program Niv stopped updating a decade ago.
    Currently it's not using the flux calculated by virus because the normalization is off by a lot. It works comparing self-dithers on the image from sloan though.

TO DO (big ones I can think of):

- Zoom in (Currently it only dithers).
- Include errors on flux.
- Normalize both flux distributions.
- Mask bad fibers.


