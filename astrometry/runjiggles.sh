#!/bin/bash
# set you path to greg's dir. Now it's set to settings in my work dir
cd ../greg
python visualize_observation.py ../astrometry/spectra/Fepses20160409T100026.3_094_sci
mv ifuPos* ../astrometry/offsets
cd ../astrometry
python jiggles.py spectra/Fepses20160409T100026.3_094_sci_ offsets/ifuPos094.txt 
#python jiggles.py spectra/FeSpses20160512T045012.9_075_sci_ offsets/ifuPos075.txt
