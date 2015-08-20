'''
Testing spectral extraction

'''

import numpy as np
import pyfits as fits
import matplotlib.pyplot as plt

import superextract
from pylab import *

#SOAR parameters
gain = 1.4
rdnoise = 4.74

specfile = '0526.WD1422p095_930_blue.fits'
datalist = fits.open(specfile)
data = datalist[0].data
data = data[0,:,:]
data = np.transpose(data)


#Calculate the variance of each pixel
#From Horne (1986) eq. 12
varmodel = rdnoise**2. + np.absolute(data)/gain


output_spec = superextract.superExtract(data,varmodel,gain,rdnoise,pord=2,tord=2,bord=2,bkg_radii=[50,60],bsigma=3,extract_radius=7,dispaxis=1,verbose=True)
#pord = order of profile polynomial, Default = 2
#tord = degree of spectral-trace polynomial, 1 = line
#bord = degree of polynomial background fit
#bkg_radii = inner and outer radii to use in computing background. Goes on both sides of aperture.
#bsigma = sigma-clipping thresholf for computing background
#extract_radius: radius for spectral extraction
#csigma = sigma-clipping threshold for cleaning & cosmic-ray rejection. Default = 5.
#qmode: how to compute Marsh's Q-matrix. 'fast-linear' default and preferred.
#nreject = number of outlier-pixels to reject at each iteration. Default = 100

print 'Done extracting!'

plt.clf()
#plt.imshow(data)
plt.plot(output_spec.spectrum-1000.,'b')
plt.plot(output_spec.raw,'g')
plt.plot(output_spec.varSpectrum,'r')
plt.plot(output_spec.trace,'m')
plt.plot(output_spec.background,'k')
plt.show()
