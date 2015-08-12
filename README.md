ZZ Ceti pipeline written by JT Fuchs and J Meza.

spectools.py - A module containing commonly done processes on 1D spectra, including reading in arrays, reading headers, rebinning, summing.

spec_sens.py - This program encompasses IRAF's standard, sensfunc, and calibrate tasks into one program. 

ReduceSpec.py - Written primarily by J Meza. This performs the inital calibration on the 2D images, including combining biases and flat, and applying those to the spectra.

superextract.py - Written mostly by Ian Crossfield and available here: http://www.lpl.arizona.edu/~ianc/python/ This program uses Horne (1986) and Marsh (1989) to perform optimal extraction of 2D data. Updated by JT Fuchs.

superextract_tools.py - Tools used by superextract.py, mostly for tracing the spectrum.



Dependencies:
- pyfits
- numpy
- scipy
- pylb
- pysynphot 0.9.6 is available from https://pypi.python.org/pypi/pysynphot - pysynphot does the rebinning of the spectra. When installing pysynphot, you might have to update spectrum.py and observations.py to 'import exceptions as exceptions'. 
