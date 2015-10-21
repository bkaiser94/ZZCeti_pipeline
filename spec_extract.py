'''
Spectral Extraction from a 2D image.
Uses superextract written by Ian Crossfield and modified by JTF.
superextract is based on optimal spectral extraction as detailed by Marsh (1989) and Horne (1986)
Crossfield's version can be found at www.lpl.arizona.edu/!ianc/python/index.html
Dependencies: superextract.py and superextrac_tools.py

For best results, first bias-subract and flat-field the 2D image before running this. It is also best to set the extract_radius to be large as this helps the profiles fit better.

To Do:
- Make it possible to read in a file with different parameters

'''

import numpy as np
import pyfits as fits
import matplotlib.pyplot as plt

import spectools as st
import superextract
from pylab import *

#SOAR parameters (Assumes 200 kHz, ATTN 0 - standard for ZZ Ceti mode)
gain = 1.4 #electrons/ADU
rdnoise = 4.74 #electrons

specfile = 'tnb.0526.WD1422p095_930_blue.fits'
datalist = fits.open(specfile)
data = datalist[0].data
data = data[0,:,:]
data = np.transpose(data)

#Calculate the variance of each pixel in ADU
varmodel = (rdnoise**2. + np.absolute(data)*gain)/gain


output_spec = superextract.superExtract(data,varmodel,gain,rdnoise,pord=2,tord=2,bord=2,bkg_radii=[40,60],bsigma=3,extract_radius=5,dispaxis=1,verbose=True,csigma=5.,polyspacing=1)
#pord = order of profile polynomial. Default = 2. This seems appropriate, no change for higher or lower order.
#tord = degree of spectral-trace polynomial, 1 = line
#bord = degree of polynomial background fit
#bkg_radii = inner and outer radii to use in computing background. Goes on both sides of aperture.  
#bsigma = sigma-clipping thresholf for computing background
#extract_radius: radius for spectral extraction. Want this value to be around the FWHM as this will optimize the S/N of the 1D spectrum.
#csigma = sigma-clipping threshold for cleaning & cosmic-ray rejection. Default = 5.
#qmode: how to compute Marsh's Q-matrix. 'fast-linear' default and preferred.
#nreject = number of outlier-pixels to reject at each iteration. Default = 100
#polyspacing = Marsh's S: the spacing between the polynomials. This should be <= 1. Default = 1. Best to leave at 1. S/N decreases dramatically if greater than 1. If less than one, slower but final spectrum is the same. Crossfield note: A few cursory tests suggests that the extraction precision (in the high S/N case) scales as S^-2 -- but the code slows down as S^2.

###########
# In superextract, to plot a 2D frame at any point, use the following
#   plt.clf()
#   plt.imshow(np.transpose(frame),aspect='auto',interpolation='nearest')
#   plt.show()
##########

print 'Done extracting!'

sigSpectrum = np.sqrt(output_spec.varSpectrum)
plt.clf()
#plt.imshow(data)
plt.plot(output_spec.spectrum,'b')
#plt.plot(output_spec.raw,'g')
#plt.plot(output_spec.varSpectrum,'r')
#plt.plot(sigSpectrum,'r')
#plt.plot(output_spec.trace,'m')
#plt.plot(output_spec.background,'k')
plt.show()

#Get the image header and add keywords
header = st.readheader(specfile)
header.set('BANDID1','Optimally Extracted Spectrum')
header.set('BANDID2','Raw Extracted Spectrum')
header.set('BANDID3','Background at trace')
header.set('BANDID4','Sigma Spectrum')
header.set('DISPCOR',0) #Dispersion axis of image

#Save the extracted image
Ni = 4. #Number of extensions
Nx = 1. #All 1D spectra
Ny = len(output_spec.spectrum[:,0])

spectrum = np.empty(shape = (Ni,Nx,Ny))
spectrum[0,:,:] = output_spec.spectrum[:,0]
spectrum[1,:,:] = output_spec.raw[:,0]
spectrum[2,:,:] = output_spec.background
spectrum[3,:,:] = sigSpectrum[:,0]

#newim = fits.PrimaryHDU(data=spectrum,header=header)
#newim.writeto('test_0344_rad7.ms.fits')


