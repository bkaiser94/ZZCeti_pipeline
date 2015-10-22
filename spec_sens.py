"""
Written by JT Fuchs in July 2015
Based off pySALT redution routine specsens.py by S. Crawford
And reading the darned IRAF documentation

spec_sens.p calculates the calibration curve given an
observation and a standard star.

Inputs: 1D standard star spectrum, 1D spectrum of observed star, filename of standard details, order of polynomial fit for sensitivity function, output filename for flux calibrated spectrum

To do:
- Make the sensitivity funtion fitting interactive. Need to be able to test different polynomials.

"""

import os
import sys
import time
import numpy as np
import pyfits as pf
import spectools as st

import matplotlib.pyplot as plt

#These functions are to help with excluding regions from the sensitivity function
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def onclick(event):
    global ix,iy
    ix, iy = event.xdata,event.ydata
    global coords
    coords.append((ix,iy))
    



#File names
stdspecfile = 'wnb.GD50_930_blue.ms.fits'
stdfile = 'mgd50.dat'
specfile = 'wnb.WD0122p0030_930_blue.ms.fits'

#Read in the observed spectrum of the standard star
obs_spectra,airmass,exptime,dispersion = st.readspectrum(stdspecfile) #This is an object containing var_farr,farr,sky,sigma,warr
#plt.clf()
#plt.plot(obs_spectra.warr,obs_spectra.farr)
#plt.show()
#read in the sandard file
std_spectra = st.readstandard(stdfile)

#plt.clf()
#plt.plot(std_spectra.warr,std_spectra.magarr)
#plt.show()
#Only keep the part of the standard file that overlaps with observation.
lowwv = np.where(std_spectra.warr >= np.min(obs_spectra.warr))
lowwv = np.asarray(lowwv)
highwv = np.where(std_spectra.warr <= np.max(obs_spectra.warr))
highwv = np.asarray(highwv)
index = np.intersect1d(lowwv,highwv)

std_spectra.warr = std_spectra.warr[index]
std_spectra.magarr = std_spectra.magarr[index]
std_spectra.wbin = std_spectra.wbin[index]

#Convert from AB mag to fnu, then to fwave (ergs/s/cm2/A)
stdzp = 3.68e-20 #The absolute flux per unit frequency at an AB mag of zero
std_spectra.magarr = st.magtoflux(std_spectra.magarr,stdzp)
std_spectra.magarr = st.fnutofwave(std_spectra.warr, std_spectra.magarr)

#plt.clf()
#plt.plot(std_spectra.warr,std_spectra.magarr)
#plt.show()
#exit()

#We want to rebin the observed spectrum to match with the bins in the standard file. This makes summing up counts significantly easier.
#Set the new binning here.
low = np.rint(np.min(obs_spectra.warr)) #Rounds to nearest integer
high = np.rint(np.max(obs_spectra.warr))
size = 0.05 #size in Angstroms you want each bin

num = (high - low) / size + 1. #number of bins. Must add one to get correct number.
wavenew = np.linspace(low,high,num=num) #wavelength of each new bin

#Now do the rebinning using Ian Crossfield's rebinning package
binflux = st.resamplespec(wavenew,obs_spectra.warr,obs_spectra.farr,200.) #200 is the oversampling factor

#plt.clf()
#plt.plot(obs_spectra.warr,obs_spectra.farr)
#plt.plot(wavenew,binflux)
#plt.show()

#Now sum the rebinned spectra into the same bins as the standard star file
counts = st.sum_std(std_spectra.warr,std_spectra.wbin,wavenew,binflux)
#plt.clf()
#plt.plot(std_spectra.warr,std_spectra.magarr)
#plt.plot(obs_spectra.warr,obs_spectra.farr,'b')
#plt.plot(std_spectra.warr,counts,'g+')
#plt.plot(irafwv,irafcounts,'r+')
#plt.show()

#Calculate the sensitivity function
sens_function = st.sensfunc(counts,std_spectra.magarr,exptime,std_spectra.wbin,airmass)
#plt.clf()
#plt.plot(std_spectra.warr,sens_function)
#plt.show()

#Fit a low order polynomial to this function so that it is smooth.
#The sensitivity function is in units of 2.5 * log10[counts/sec/Ang / ergs/cm2/sec/Ang]
#Choose regions to not include in fit
coords = []
fig = plt.figure(1)
ax = fig.add_subplot(111)
ax.plot(std_spectra.warr,sens_function)
cid = fig.canvas.mpl_connect('button_press_event',onclick)
print 'Please click on both sides of regions you want to exclude. Then close the plot.'
plt.show(1)

#Mask our the regions you don't want to git
mask = np.ones(len(std_spectra.warr))
n = 0
if len(coords) > 0:
    while n < len(coords):
        x1 = np.where(std_spectra.warr == (find_nearest(std_spectra.warr,coords[n][0])))
        n += 1
        x2 = np.where(std_spectra.warr == (find_nearest(std_spectra.warr,coords[n][0])))
        mask[x1[0][0]:x2[0][0]] = 0
        n += 1

indices = np.where(mask !=0.)
lambdasfit = std_spectra.warr[indices]
fluxesfit = sens_function[indices]

#Make sure they are finite
ind1 = np.isfinite(lambdasfit) & np.isfinite(fluxesfit)
lambdasfit = lambdasfit[ind1]
fluxesfit = fluxesfit[ind1]


print 'Fitting the sensitivity funtion now. Close the plot to continue.'
order = 4.
repeat = 'yes'
while repeat == 'yes':
    p = np.polyfit(lambdasfit,fluxesfit,order)
    f = np.poly1d(p)
    smooth_sens = f(lambdasfit)
    plt.clf()
    plt.plot(lambdasfit,fluxesfit,'b+')
    plt.plot(lambdasfit,smooth_sens,'r',linewidth=2.0)
    plt.show()
    repeat = raw_input('Do you want to try again (yes/no)? ')
    if repeat == 'yes':
        order = raw_input('New order for polynomial: ')


#Read in the spectrum we want to flux calibrate
WD_spectra,airmass,exptime,dispersion = st.readspectrum(specfile)

#To get the flux calibration, perform the following
#Flux = counts / (Exptime * dispersion * 10**(sens/2.5))

#Get the sensitivity function at the correct wavelength spacing
sens_wave = f(WD_spectra.warr)

#Perform the flux calibration. We do this on the non-variance weighted aperture, the sky spectrum, and the sigma spectrum.

star_opflux = st.cal_spec(WD_spectra.opfarr,sens_wave,exptime,dispersion)
star_flux = st.cal_spec(WD_spectra.farr,sens_wave,exptime,dispersion)
sky_flux = st.cal_spec(WD_spectra.sky,sens_wave,exptime,dispersion)
sigma_flux = st.cal_spec(WD_spectra.sigma,sens_wave,exptime,dispersion)

#plt.clf()
#plt.plot(WD_spectra.warr,star_opflux)
#plt.show()


#Save the wavelenghts, counts, and fluxes
#np.savetxt('python_counts.txt',np.transpose([std_spectra.warr,counts]))
#np.savetxt('python_sens.txt',np.transpose([obs_spectra.warr,sens_wave]))
#np.savetxt('python_flux.txt',np.transpose([WD_spectra.warr,flux]))

#Save the flux-calibrated spectrum and update the header
header = st.readheader(specfile)
header.set('EX-FLAG',-1) #Extiction correction? 0=yes, -1=no
header.set('CA-FLAG',0) #Calibrated to flux scale? 0=yes, -1=no
header.set('BUNIT','erg/cm2/s/A') #physical units of the array value


#Set up size of new fits image
Ni = 3. #Number of extensions
Nx = len(star_flux)
Ny = 1. #All 1D spectra

data = np.empty(shape = (Ni,Ny,Nx))
data[0,:,:] = star_flux
data[1,:,:] = sky_flux
data[2,:,:] = sigma_flux

#newim = pf.PrimaryHDU(data=data,header=header)
#newim.writeto('blah.fits')



