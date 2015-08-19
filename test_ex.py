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


output_spec = superextract.superExtract(data,varmodel,gain,rdnoise,pord=1,tord=2,bkg_radii=[50,60],extract_radius=3,dispaxis=1,verbose=True)
#tord = degree of spectral-trace polynomial, 1 = line

print 'Done extracting!'

plt.clf()
#plt.imshow(data)
plt.plot(output_spec.spectrum,'b')
plt.plot(output_spec.raw,'g')
plt.plot(output_spec.varSpectrum,'r')
plt.plot(output_spec.trace,'m')
plt.plot(output_spec.background,'k')
plt.show()
