#!/bin/bash
# set you path to greg's dir. Now it's set to settings in my work dir
python ../greg/visualize_observation.py Fepses20160409T100026.3_085_sci
mv ifuPos* ../astrometry/offsets
python jiggles.py /work/03946/hetdex/maverick/20160409/virus/virus0000301/exp01/virus/20160409T100026.3_075
#python jiggles.py spectra/FeSpses20160512T045012.9_075_sci_ offsets/ifuPos075.txt
