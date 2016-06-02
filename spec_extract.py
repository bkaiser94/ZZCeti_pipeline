'''
Spectral Extraction from a 2D image.
Uses superextract written by Ian Crossfield and modified by JTF.
superextract is based on optimal spectral extraction as detailed by Marsh (1989) and Horne (1986)
Crossfield's version can be found at www.lpl.arizona.edu/!ianc/python/index.html
Dependencies: superextract.py and superextrac_tools.py

For best results, first bias-subract, flat-field, and trim the 2D image before running this. It is also best to set the extract_radius to be approximately the FWHM. This maximizes the S/N.

To run:
python spec_extract.py filename_of_2d_spectrum
python spec_extract.py tnb.0526.WD1422p095_930_blue.fits

Inputs:
spectral filename, lamp file (?)
extract_radius and  bkg_radii computed automatically. FWHM used for extract_radius
output file name and lamp output filename done automatically

To Do:
- 

'''
import sys
import os
import numpy as np
import pyfits as fits
import matplotlib.pyplot as plt
import datetime
import mpfit
from scipy.interpolate import InterpolatedUnivariateSpline

import spectools as st
import superextract
from superextract_tools import lampextract
from pylab import *

def gauss(x,p): #single gaussian
    return p[0] +  p[1]*np.exp(-(((x-p[2])/(np.sqrt(2)*p[3])))**2.)

def fitgauss(p,fjac=None,x=None,y=None,err=None):
    #Parameter values are passed in p
    #fjac = None just means partial derivatives will not be computed
    model = gauss(x,p)
    status = 0
    return([status,(y-model)/err])

#Read in file from command line
if len(sys.argv) == 3:
    script, specfile, lamp = sys.argv
if len(sys.argv) == 2:
    script, specfile = sys.argv
    lampcheck = raw_input('Do you want to extract a lamp too? (yes/no) ')
    if lampcheck == 'yes':
        lamp = raw_input('Enter lamp name: ')
    else:
        lamp = 'no'


#script, specfile = sys.argv
#specfile = 'tnb.0526.WD1422p095_930_blue.fits'


#Open file and read gain and readnoise
datalist = fits.open(specfile)
data = datalist[0].data
data = data[0,:,:]
data = np.transpose(data)

gain = datalist[0].header['GAIN']
rdnoise = datalist[0].header['RDNOISE']

#Calculate the variance of each pixel in ADU
varmodel = (rdnoise**2. + np.absolute(data)*gain)/gain

#Fit a column of the 2D image to determine the FWHM
if 'blue' in specfile.lower():
    forfit = data[1200,:]
elif 'red' in specfile.lower():
    forfit = data[1000,:]

guess = np.zeros(4)
guess[0] = np.mean(forfit)
guess[1] = np.amax(forfit)
guess[2] = np.argmax(forfit)
guess[3] = 3.

error_fit = np.ones(len(forfit))
xes = np.linspace(0,len(forfit)-1,num=len(forfit))
fa = {'x':xes,'y':forfit,'err':error_fit}
fitparams = mpfit.mpfit(fitgauss,guess,functkw=fa,quiet=True)

fwhm = 2.*np.sqrt(2.*np.log(2.))*fitparams.params[3]
#extraction_rad = 2.*np.round(fwhm,decimals=1)
extraction_rad = 5. * np.round(fwhm,decimals=1) #Extract up to 5 times FWHM


#Check to make sure background region does not go within 10 pixels of edge
background_radii = [35,60]
#First check this against the bottom
if fitparams.params[2] - background_radii[1] < 10.:
    background_radii[1] = fitparams.params[2] - 10.
    background_radii[0] -= 60 - background_radii[1]
#Then check against the top
hold = background_radii[1]
if fitparams.params[2] + background_radii[1] > 190.:
    background_radii[1] = 190. - fitparams.params[2]
    background_radii[0] -= hold - background_radii[1]
#Ensure that the closest point is at least 20 pixels away.
if background_radii[0] < 20.:
    background_radii[0] = 20.
background_radii[0] = np.round(background_radii[0],decimals=1)
background_radii[1] = np.round(background_radii[1],decimals=1)
#plt.plot(data[1200,:])
#plt.plot(xes,gauss(xes,fitparams.params))
#plt.show()
#extraction_rad = 10.
#background_radii = [40,60]


output_spec = superextract.superExtract(data,varmodel,gain,rdnoise,pord=2,tord=2,bord=2,bkg_radii=background_radii,bsigma=3,extract_radius=extraction_rad,dispaxis=1,verbose=False,csigma=5.,polyspacing=1)
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
#Verbose=True if you want lots of output

###########
# In superextract, to plot a 2D frame at any point, use the following
#   plt.clf()
#   plt.imshow(np.transpose(frame),aspect='auto',interpolation='nearest')
#   plt.show()
##########

print 'Done extracting!'

sigSpectrum = np.sqrt(output_spec.varSpectrum)
#plt.clf()
#plt.imshow(data)
#plt.plot(output_spec.spectrum,'b')
#plt.plot(output_spec.raw,'g')
#plt.plot(output_spec.varSpectrum,'r')
#plt.plot(sigSpectrum,'r')
#plt.plot(output_spec.trace,'m')
#plt.plot(output_spec.background,'k')
#plt.show()

#Get the image header and add keywords
header = st.readheader(specfile)
header.set('BANDID1','Optimally Extracted Spectrum')
header.set('BANDID2','Raw Extracted Spectrum')
header.set('BANDID3','Mean Background')
header.set('BANDID4','Sigma Spectrum')
header.set('DISPCOR',0) #Dispersion axis of image
fwhmsave = np.round(fwhm,decimals=4)
header.set('SPECFWHM',fwhmsave,'FWHM of spectrum in pixels') #FWHM of spectrum in pixels

#Save the extracted image
Ni = 4. #Number of extensions
Nx = 1. #All 1D spectra
Ny = len(output_spec.spectrum[:,0])

spectrum = np.empty(shape = (Ni,Nx,Ny))
spectrum[0,:,:] = output_spec.spectrum[:,0]
spectrum[1,:,:] = output_spec.raw[:,0]
spectrum[2,:,:] = output_spec.background
spectrum[3,:,:] = sigSpectrum[:,0]

#Save the extracted spectra with .ms.fits in the filename
#Ask to overwrite if file already exists or provide new name
loc = specfile.find('.fits')
newname = specfile[0:loc] + '.ms.fits'
clob = False

mylist = [True for f in os.listdir('.') if f == newname]
exists = bool(mylist)

if exists:
    print 'File %s already exists.' % newname
    nextstep = raw_input('Do you want to overwrite or designate a new name (overwrite/new)? ')
    if nextstep == 'overwrite':
        clob = True
        exists = False
    elif nextstep == 'new':
        newname = raw_input('New file name: ')
        exists = False
    else:
        exists = False


newim = fits.PrimaryHDU(data=spectrum,header=header)
newim.writeto(newname,clobber=clob)
print 'Wrote %s to file.' % newname

#Save parameters to a file for future reference. 
# specfile,date of extraction, extration_rad,background_radii,newname,newname2
f = open('extraction_params.txt','a')
now = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M")
newinfo = specfile + '\t' + now + '\t' + str(extraction_rad) + '\t' + str(background_radii) + '\t' + newname
f.write(newinfo + "\n")
f.close()

###########################
#Extract a lamp spectrum using the trace from above
##########################

#lamp = 't0025.Fe_ZZCeti_930_blue_long.fits'
if lamp != 'no':
    lamplist = fits.open(lamp)
    lampdata = lamplist[0].data
    lampdata = lampdata[0,:,:]
    lampdata = np.transpose(lampdata)

    #extraction radius will be the same as the star
    #But since the Fe lamps are taken binned by 1 in spectral direction, we need to adjust the trace to match.
    #We do that by expanding out the trace to match the length of the Fe lamps, then interpolate and read off values for every pixel.
    bin2size = np.arange(1,len(output_spec.trace)+1)
    bin1size = np.arange(1,len(lampdata)+1)
    ratio = float(len(bin1size)) / float(len(bin2size))
    interpolates = InterpolatedUnivariateSpline(ratio*bin2size,output_spec.trace,k=1)
    bin1size = np.arange(1,len(lampdata)+1)
    newtrace = interpolates(bin1size)

    #Do the extraction here.
    lampspec = lampextract(lampdata,newtrace,extraction_rad)


    #Save the 1D lamp
    lampheader = st.readheader(lamp)
    lampheader.set('BANDID2','Raw Extracted Spectrum')
    lampheader.set('REF',newname,'Reference Star used for trace')

    Ni = 1. #We are writing just 1 1D spectrum
    Ny = len(lampspec[:,0])
    lampspectrum = np.empty(shape = (Ni,Ny))
    lampspectrum[0,:] = lampspec[:,0]

    #plt.clf()
    #plt.plot(lampspec)
    #plt.show()

    #Save the extracted spectra with .ms.fits in the filename
    #Ask to overwrite if file already exists or provide new name
    loc2 = lamp.find('.fits')
    newname2 = lamp[0:loc2] + '.ms.fits'
    clob = False

    mylist = [True for f in os.listdir('.') if f == newname2]
    exists = bool(mylist)

    if exists:
        print 'File %s already exists.' % newname2
        nextstep = raw_input('Do you want to overwrite or designate a new name (overwrite/new)? ')
        if nextstep == 'overwrite':
            clob = True
            exists = False
        elif nextstep == 'new':
            newname2 = raw_input('New file name: ')
            exists = False
        else:
            exists = False


    lampim = fits.PrimaryHDU(data=lampspectrum,header=lampheader)
    lampim.writeto(newname2,clobber=clob)
    print 'Wrote %s to file.' % newname2

    #Save parameters to a file for future reference. 
    # specfile,date of extraction, extration_rad,background_radii,newname,newname2
    background_radii2 = [0,0] #We do not extract a background for the lamp
    f = open('extraction_params.txt','a')
    now = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M")
    newinfo2 = lamp + '\t' + now + '\t' + str(extraction_rad) + '\t' + str(background_radii2) + '\t' + newname2
    f.write(newinfo2 + "\n")
    f.close()

#######################
# End lamp extraction
#######################


#To unpack these values, use the following
#arr = np.genfromtxt('extraction_params.txt',dtype=None,delimiter='\t')
#names, date, fwhm, back, newname = [], [],np.zeros(len(arr)),[],[]
#for m in np.arange(len(arr)):
#    names.append(arr[m][0])
#    date.append(arr[m][1])
#    fwhm[m] = arr[m][2]    
#    back.append(arr[m][3])
#    newname.append(arr[m][4])
    
