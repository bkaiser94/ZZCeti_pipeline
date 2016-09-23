'''
Plot 1D spectrum. If wavelength solution exists from grating equation, use that. If not, plot by pixel number. 

Written by Josh Fuchs, UNC. June 2016.

:OPTIONAL INPUTS:
       specname: string, filename of 1D spectrum to plot. If not supplied, user will be prompted. 
'''

import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt
from sys import argv

def DispCalc(Pixels, alpha, theta, fr, fd, fl, zPnt):
    # This is the Grating Equation used to calculate the wavelenght of a pixel
    # based on the fitted parameters and angle set up.
    # Inputs: 
    # Pixels= Vector of Pixel Numbers
    # alpha=  of Grating Angle
    # aheta=  Camera Angle
    # fr= fringe density of grating
    # fd= Camera Angle Correction Factor
    # zPnt= Zero point pixel 
    Wavelengths= [] # Vector to store calculated wavelengths 
    for pix in Pixels:    
        beta = np.arctan( (pix-zPnt)*15./fl ) + (fd*theta*np.pi/180.) - (alpha*np.pi/180.) 
        wave = (10**6.)*( np.sin(beta) + np.sin(alpha*np.pi/180.) )/fr
        Wavelengths.append(wave)
    return Wavelengths

if len(argv) == 2:
    script, specname = argv[0], argv[1]
elif len(argv) == 1:
    specname = raw_input('Name of spectrum to plot: ')
else:
    print 'Too many inputs. Please try again.'
    exit()

#specname = 'wtfb.wd1425-811_930_red_flux.ms.fits'
spec_data= pf.getdata(specname)
spec_header= pf.getheader(specname)

#See if wavelength solution exists. If so, use it. Otherwise, use pixel numbers
try:
    alpha = float(spec_header['GRT_TARG'])
    theta = float(spec_header['CAM_TARG'])
    fr = float(spec_header['LINDEN'])
    fd = float(spec_header['CAMFUD'])
    fl = float(spec_header['FOCLEN'])
    zPnt = float(spec_header['ZPOINT'])
    
    trim_sec= spec_header["CCDSEC"]
    trim_offset= float( trim_sec[1:len(trim_sec)-1].split(':')[0] )-1
    bining= float( spec_header["PARAM18"] ) 
    nx= np.size(spec_data[0])
    Pixels= bining*(np.arange(0,nx,1)+trim_offset)
    
    WDwave = DispCalc(Pixels, alpha, theta, fr, fd, fl, zPnt)
except:
    WDwave = np.arange(len(spec_data[0,0,:]))

#np.savetxt('OWJ1818-2434_930_blue.txt',np.transpose([WDwave2,spec_data[0,0,:],spec_data[3,0,:],spec_data[2,0,:]]),header='wavelengths, optimally extracted spectrum, sigma spectrum, sky spectrum')


if '_fe_' in specname.lower():
    plt.clf()
    plt.plot(WDwave,spec_data[0,:])
    plt.show()
else:
    plotagain = 'y'
    print 'Extension options:'
    print '0: Optimally extracted spectrum'
    print '1: Raw extracted spectrum'
    print '2: Sky spectrum'
    print '3: Sigma spectrum'
    while plotagain == 'y':
        extension = raw_input('Which extenstion? >> ')
        plt.clf()
        plt.plot(WDwave,spec_data[extension,0,:])
        plt.show()
        plotagain = raw_input('Would you like to plot a different extension? (y/n)')