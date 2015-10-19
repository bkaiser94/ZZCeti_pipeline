"""
Written by JT Fuchs in July 2015
Based off pySALT redution routine specsens.py by S. Crawford
And reading the darned IRAF documentation

spec_sens.p calculates the calibration curve given an
observation and a standard star.


"""

import os
import sys
import time
import numpy as np
import pyfits as pf
import spectools as st

import matplotlib.pyplot as plt


#File names
stdspecfile = 'wnb.GD50_930_red.ms.fits'
stdfile = 'mgd50.dat'
specfile = 'wnb.WD0122p0030_930_blue.ms.fits'

#Read in the observed spectrum of the standard star
obs_spectra,airmass,exptime,dispersion = st.readspectrum(stdspecfile) #This is an object containing var_farr,farr,sky,sigma,warr
#plt.clf()
#plt.plot(obs_spectra.warr,obs_spectra.farr)

#read in the sandard file
std_spectra = st.readstandard(stdfile)

#plt.clf()
#plt.plot(std_spectra.warr,std_spectra.magarr)

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

#Now do the rebinning using Ian Crossfield rebinning package
binflux = st.resamplespec(wavenew,obs_spectra.warr,obs_spectra.farr,200.) #200 is the oversampling factor

plt.clf()
plt.plot(obs_spectra.warr,obs_spectra.farr)
plt.plot(wavenew,binflux)
#plt.show()

irafwv,irafflux,irafbin,irafcounts = np.genfromtxt('wnb.GD50_930_red_std',unpack=True,skip_header=1)

#Now sum the rebinned spectra into the same bins as the standard star file
counts = st.sum_std(std_spectra.warr,std_spectra.wbin,wavenew,binflux)
plt.clf()
#plt.plot(std_spectra.warr,std_spectra.magarr)
plt.plot(obs_spectra.warr,obs_spectra.farr,'b')
plt.plot(std_spectra.warr,counts,'g+')
plt.plot(irafwv,irafcounts,'r+')
plt.show()

#Calculate the sensitivity function
sens_function = st.sensfunc(counts,std_spectra.magarr,exptime,std_spectra.wbin,airmass)
#plt.clf()
#plt.plot(std_spectra.warr,sens_function)
#plt.show()

#Fit a low order polynomial to this function so that it is smooth. Want to ignore the first few and last few pixels.
#The sensitivity function is in units of 2.5 * log10[counts/sec/Ang / ergs/cm2/sec/Ang]
new_wave = std_spectra.warr[25:-10]
new_sens = sens_function[25:-10]
p = np.polyfit(new_wave,new_sens,4)
f = np.poly1d(p)
#smooth_sens = f(new_wave)
#plt.clf()
#plt.plot(std_spectra.warr[25:-10],sens_function[25:-10])
#plt.plot(new_wave,smooth_sens)
#plt.show()
#sys.exit()

#Read in the spectrum we want to flux calibrate
WD_spectra,airmass,exptime,dispersion = st.readspectrum(specfile)

#To get the flux calibration, perform the following
#Flux = counts / (Exptime * dispersion * 10**(sens/2.5))

#Get the sensitivity function at the correct wavelength spacing
sens_wave = f(WD_spectra.warr)

#Perform the flux calibration. We do this on the non-variance weighted aperture, the sky spectrum, and the sigma spectrum.
star_flux = st.cal_spec(WD_spectra.farr,sens_wave,exptime,dispersion)
sky_flux = st.cal_spec(WD_spectra.sky,sens_wave,exptime,dispersion)
sigma_flux = st.cal_spec(WD_spectra.sigma,sens_wave,exptime,dispersion)

#plt.clf()
#plt.plot(WD_spectra.warr,smooth_sens)
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



