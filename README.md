ZZ Ceti pipeline written by JT Fuchs and J Meza.

spectools.py - A module containing commonly done processes on 1D spectra, including reading in arrays, reading headers, rebinning, summing.

spec_sens.py - This program encompasses IRAF's standard, sensfunc, and calibrate tasks into one program. Uses tools written by Ian Crossfield for rebinning the spectra.

ReduceSpec.py - Written primarily by J Meza. This performs the inital calibration on the 2D images, including combining biases and flat, and applying those to the spectra.

spec_extract.py - This program calls superextract.py. All the setting up for the spectral extraction and saving the file are done here.

superextract.py - Written mostly by Ian Crossfield and available here: http://www.lpl.arizona.edu/~ianc/python/ This program uses Horne (1986) and Marsh (1989) to perform optimal extraction of 2D data. Updated by JT Fuchs.

superextract_tools.py - Tools used by superextract.py, mostly for tracing the spectrum.

The reduction flow is a follows: Bias-subtract, flat-field, trim (ReduceSpec.py). Extract spectrum (spec_extract.py). Flux calibration (spec_sens.py).

Dependencies:
- pyfits
- numpy
- scipy
- pylb
- mpfit (can be found at http://code.google.com/p/astrolibpy/source/browse/trunk/)



