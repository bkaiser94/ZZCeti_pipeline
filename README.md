ZZ Ceti pipeline written by JT Fuchs and J Meza.

spectools.py - A module containing commonly done processes on 1D spectra, including reading in arrays, reading headers, rebinning, summing.

spec_sens.py - This program encompasses IRAF's standard, sensfunc, and calibrate tasks into one program. 

Dependencies:
- pyfits
- numpy
- pysynphot 0.9.6 is available from https://pypi.python.org/pypi/pysynphot - pysynphot does the rebinning of the spectra. When installing pysynphot, you might have to update spectrum.py and observations.py to 'import exceptions as exceptions'. 
