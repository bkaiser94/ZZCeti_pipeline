ZZ Ceti pipeline written by JT Fuchs and J Meza.

spectools.py - A module containing commonly done processes on 1D spectra, including reading in arrays, reading headers, rebinning, summing.

spec_sens.py - This program encompasses IRAF's standard, sensfunc, and calibrate tasks into one program. Uses tools written by Ian Crossfield for rebinning the spectra.

ReduceSpec.py - Written primarily by J Meza. This performs the inital calibration on the 2D images, including combining biases and flats, trimming, and applying those to the spectra.

spec_extract.py - This program calls superextract.py. All the setting up for the spectral extraction and saving the file are done here.

superextract.py - Written mostly by Ian Crossfield and available here: http://www.lpl.arizona.edu/~ianc/python/ This program uses Horne (1986) and Marsh (1989) to perform optimal extraction of 2D data. Updated by JT Fuchs.

superextract_tools.py - Tools used by superextract.py, mostly for tracing the spectrum.

WaveCal.py - Written by J Meza. Fits the grating equation to a spectral lamp and applies it to a 1D spectrum

plotspec.py - Uses the save grating equation to plot the spectrum. You choose which extension.


The reduction flow is a follows: Bias-subtract, flat-field, trim (ReduceSpec.py). Extract spectrum (spec_extract.py). Flux calibration (spec_sens.py).

python ReduceSpec.py listZero, listFlat, listSpec, listFe
       - The final output files are tfb*fits.

python spec_extrac.py tfb.CD-32_9927_930_blue_series.fits
       -Will be asked if you want to extract a lamp too. If so, you must enter the lamp name. Output file will be a 1D spectrum with the same base name and .ms.fits at the end.

python WaveCal.py lamp_spec.fits
       - You will first shift the guess and can see if it works, but usually will have to refit the lines, selecting ones you want to use. After you are happy, you can save the parameters to the fit to the header and are given the option to apply it to another spectrum. If you do this, the program will fit a skyline to determine the new zero point. Output files have w added to the beginning of the filename. 

python spec_sens.py liststandard listflux liststar
       - You have to click the regions you want to exclude from the polynomial fit. 

Dependencies:
- pyfits
- numpy
- scipy
- pylb
- mpfit (can be found at http://code.google.com/p/astrolibpy/source/browse/trunk/)



