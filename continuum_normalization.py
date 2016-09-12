'''
Initially written by E. Dennihy, May 2015. Modified by J. Fuchs June 2016.

Continuum normalizes a ZZ Ceti spectrum to match a DA model

To run: (Red filename is optional)
python model_calibration.py bluefilename redfilename
python model_calibration.py wtfb.wd1401-147_930_blue.ms.fits wtfb.wd1401-147_930_red.ms.fits

:INPUTS:
     bluefilename: string, filename of wavelength calibrated ZZ Ceti blue spectrum

:OPTIONAL:
     redfilename: string, filename of wavelength calibrated ZZ Ceti red spectrum


'''


import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as inter
import pyfits as pf
import scipy.signal as signal
import spectools as st
import os
import sys
import datetime


if len(sys.argv) == 3:
    script, filenameblue, filenamered = sys.argv
    redfile = True
elif len(sys.argv) == 2:
    script, filenameblue = sys.argv
    redfile = False
else:
    print '\n Incorrect number of arguments. \n'

#Read in the observed spectrum
obs_spectrablue,airmass,exptime,dispersion = st.readspectrum(filenameblue)
datalistblue = pf.open(filenameblue)
if redfile:
    obs_spectrared, airmassred,exptimered,dispersionred = st.readspectrum(filenamered)

FWHMpix = datalistblue[0].header['specfwhm']
FWHM = FWHMpix * (obs_spectrablue.warr[-1] - obs_spectrablue.warr[0])/len(obs_spectrablue.warr)



#Read in DA model
cwd = os.getcwd()
os.chdir('/afs/cas.unc.edu/depts/physics_astronomy/clemens/students/group/modelfitting/Koester_08')
dafile = 'da12500_800.dk'
mod_wav, mod_spec = np.genfromtxt(dafile,unpack=True,skip_header=33)
os.chdir(cwd) #Move back to directory with observed spectra


#Convolve the model to match the seeing of the spectrum
intlambda = np.divide(range(31000),10.) + 3660.0
lowlambda = np.min(np.where(mod_wav > 3600.))
highlambda = np.min(np.where(mod_wav > 6800.))
shortlambdas = mod_wav[lowlambda:highlambda]
shortinten = mod_spec[lowlambda:highlambda]
interp = inter.InterpolatedUnivariateSpline(shortlambdas,shortinten,k=1)
intflux = interp(intlambda)
sig = FWHM / (2. * np.sqrt(2.*np.log(2.)))
gx = np.divide(range(360),10.)
gauss = (1./(sig * np.sqrt(2. * np.pi))) * np.exp(-(gx-18.)**2./(2.*sig**2.))
gf = np.divide(np.outer(intflux,gauss),10.)
length = len(intflux) - 360.
cflux = np.zeros(length)
clambda = intlambda[180:len(intlambda)-180]
x  = 0
while x < length:
    cflux[x] = np.sum(np.diagonal(gf,x,axis1=1,axis2=0),dtype='d')
    x += 1
interp2 = inter.InterpolatedUnivariateSpline(clambda,cflux,k=1)
cflux2blue = interp2(obs_spectrablue.warr)
if redfile:
    cflux2red = interp2(obs_spectrared.warr)

plt.clf()
plt.plot(obs_spectrablue.warr,obs_spectrablue.opfarr,'b')
plt.plot(obs_spectrablue.warr,cflux2blue*8000./cflux2blue[0],'r')
if redfile:
    plt.plot(obs_spectrared.warr,obs_spectrared.opfarr,'b')
    plt.plot(obs_spectrared.warr,cflux2red*8000./cflux2red[0],'r')
#plt.show()


balmer_features_blue = [[3500,3650],[3745,3757],[3760,3780],[3784,3812],[3816,3856],[3865,3921],[3935,4021],[4040,4191],[4223,4460],[4691,5019]]
balmer_features_red = [[6350,6780],[6835,6970]]

balmer_mask_blue = obs_spectrablue.warr == obs_spectrablue.warr
for wavrange in balmer_features_blue:
    inds = np.where((obs_spectrablue.warr > wavrange[0]) & (obs_spectrablue.warr < wavrange[1]))
    balmer_mask_blue[inds] = False

if redfile:
    balmer_mask_red = obs_spectrared.warr == obs_spectrared.warr
    for wavrange in balmer_features_red:
        indxs = np.where((obs_spectrared.warr > wavrange[0]) * (obs_spectrared.warr < wavrange[1]))
        balmer_mask_red[indxs] = False

mod_wav_masked_blue = obs_spectrablue.warr[balmer_mask_blue]
mod_spec_masked_blue = cflux2blue[balmer_mask_blue]/(10**13.6)

spec_wav_masked_blue = obs_spectrablue.warr[balmer_mask_blue]
spec_flux_masked_blue = obs_spectrablue.opfarr[balmer_mask_blue]
weights_masked_blue = obs_spectrablue.sigma[balmer_mask_blue]

if redfile:
    mod_wav_masked_red = obs_spectrared.warr[balmer_mask_red]
    mod_spec_masked_red = cflux2red[balmer_mask_red]/(10**13.6)
    
    spec_wav_masked_red = obs_spectrared.warr[balmer_mask_red]
    spec_flux_masked_red = obs_spectrared.opfarr[balmer_mask_red]
    weights_masked_red = obs_spectrared.sigma[balmer_mask_red]

'''
#This section uses a spline to fit the spectra
def moving_average(series,pts):
    b = signal.gaussian(len(series),pts)
    average = signal.convolve(series,b,mode='same')
    var = signal.convolve(np.power(series-average,2),b,mode='same')
    return average, var

def spline_fit(wav,data,scale,pts,smooth=False,k=3):
    if smooth:
        _ , var = moving_average(data,pts)
        spline = inter.UnivariateSpline(wav,data,w=scale/np.sqrt(var),k=k)
    else:
        spline = inter.UnivariateSpline(wav,data)
    nor_data = data/spline(wav)
    return nor_data,spline

mod_nor_blue, mod_fit_blue = spline_fit(mod_wav_masked_blue,mod_spec_masked_blue,1.0,5,smooth=True,k=5)
spec_nor_blue, spec_fit_blue = spline_fit(spec_wav_masked_blue,spec_flux_masked_blue,1.0,11,smooth=True,k=5)

if redfile:
    mod_nor_red, mod_fit_red = spline_fit(mod_wav_masked_red,mod_spec_masked_red,1.0,5,smooth=True,k=5)
    spec_nor_red, spec_fit_red = spline_fit(spec_wav_masked_red,spec_flux_masked_red,1.0,11,smooth=True,k=5)
'''

#This section uses a polynomial to fit the spectra
spec_poly_order_blue = 5.
mod_poly_order_blue = 6.
spec_fit_blue_poly = np.polyfit(spec_wav_masked_blue,spec_flux_masked_blue,spec_poly_order_blue,w=weights_masked_blue)
spec_fit_blue = np.poly1d(spec_fit_blue_poly)
mod_fit_blue_poly = np.polyfit(mod_wav_masked_blue,mod_spec_masked_blue,mod_poly_order_blue)
mod_fit_blue = np.poly1d(mod_fit_blue_poly)

if redfile:
    spec_poly_order_red = 3.
    mod_poly_order_red = 3.
    spec_fit_red_poly = np.polyfit(spec_wav_masked_red,spec_flux_masked_red,spec_poly_order_red,w=weights_masked_red)
    spec_fit_red = np.poly1d(spec_fit_red_poly)
    mod_fit_red_poly = np.polyfit(mod_wav_masked_red,mod_spec_masked_red,mod_poly_order_red)
    mod_fit_red = np.poly1d(mod_fit_red_poly)


plt.clf()
plt.plot(obs_spectrablue.warr,cflux2blue/(10**13.6),'r')
plt.plot(mod_wav_masked_blue,mod_spec_masked_blue,'b.')
plt.plot(obs_spectrablue.warr,mod_fit_blue(obs_spectrablue.warr),'k--')
plt.plot(obs_spectrablue.warr,obs_spectrablue.opfarr,'r')
plt.plot(spec_wav_masked_blue,spec_flux_masked_blue,'g.')
plt.plot(obs_spectrablue.warr,spec_fit_blue(obs_spectrablue.warr),'k--')



#plt.clf()
if redfile:
    plt.plot(obs_spectrared.warr,cflux2red/(10**13.6),'r')
    plt.plot(mod_wav_masked_red,mod_spec_masked_red,'b.')
    plt.plot(obs_spectrared.warr,mod_fit_red(obs_spectrared.warr),'k--')
    plt.plot(obs_spectrared.warr,obs_spectrared.opfarr,'r')
    plt.plot(spec_wav_masked_red,spec_flux_masked_red,'g.')
    plt.plot(obs_spectrared.warr,spec_fit_red(obs_spectrared.warr),'k--')
plt.show()


wd_response_blue = mod_fit_blue(obs_spectrablue.warr)/spec_fit_blue(obs_spectrablue.warr)
if redfile:
    wd_response_red = mod_fit_red(obs_spectrared.warr)/spec_fit_red(obs_spectrared.warr)

plt.clf()
plt.plot(obs_spectrablue.warr,wd_response_blue,'b')
if redfile:
    plt.plot(obs_spectrared.warr,wd_response_red,'r')
plt.show()

fcorr_wd_blue_opfarr = obs_spectrablue.opfarr * wd_response_blue
fcorr_wd_blue_farr = obs_spectrablue.farr * wd_response_blue
fcorr_wd_blue_sky = obs_spectrablue.sky * wd_response_blue
fcorr_wd_blue_sigma = obs_spectrablue.sigma * wd_response_blue


if redfile:
    fcorr_wd_red_opfarr = obs_spectrared.opfarr * wd_response_red
    fcorr_wd_red_farr = obs_spectrared.farr * wd_response_red
    fcorr_wd_red_sky = obs_spectrared.sky * wd_response_red
    fcorr_wd_red_sigma = obs_spectrared.sigma * wd_response_red

plt.clf()
plt.plot(obs_spectrablue.warr,fcorr_wd_blue_opfarr,'b')
if redfile:
    plt.plot(obs_spectrared.warr,fcorr_wd_red_opfarr,'r')
plt.show()
#exit()

#Save parameters for diagnostics
if redfile:
    bigarray = np.zeros([len(obs_spectrablue.warr),14])
    bigarray[:,0] = obs_spectrablue.warr
    bigarray[:,1] = obs_spectrablue.opfarr
    bigarray[:,2] = spec_fit_blue(obs_spectrablue.warr)
    bigarray[:,3] = cflux2blue/(10**13.6)
    bigarray[:,4] = mod_fit_blue(obs_spectrablue.warr)
    bigarray[:,5] = wd_response_blue
    bigarray[:,6] = fcorr_wd_blue_opfarr
    bigarray[:,7] = obs_spectrared.warr
    bigarray[:,8] = obs_spectrared.opfarr
    bigarray[:,9] = spec_fit_red(obs_spectrared.warr)
    bigarray[:,10] = cflux2red/(10**13.6)
    bigarray[:,11] = mod_fit_red(obs_spectrared.warr)
    bigarray[:,12] = wd_response_red
    bigarray[:,13] = fcorr_wd_red_opfarr
    now = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M")
    endpoint = '_930'
    with open('continuum_normalization_' + filenameblue[5:filenameblue.find(endpoint)] + '_' + now + '.txt','a') as handle:
        header = str(filenameblue) + ',' + str(filenamered) + ',' + dafile + '\n Columns structured as blue then red. If no red file, only blue data given. Columns are: wavelengths, optimized spectra, polynomial fit to data, \n model fluxes, polynomial fit to model, response, continuum normalized spectrum.'
        np.savetxt(handle,bigarray,fmt='%f',header=header)
if not redfile:
    bigarray = np.zeros([len(obs_spectrablue.warr),7])
    bigarray[:,0] = obs_spectrablue.warr
    bigarray[:,1] = obs_spectrablue.opfarr
    bigarray[:,2] = spec_fit_blue(obs_spectrablue.warr)
    bigarray[:,3] = cflux2blue/(10**13.6)
    bigarray[:,4] = mod_fit_blue(obs_spectrablue.warr)
    bigarray[:,5] = wd_response_blue
    bigarray[:,6] = fcorr_wd_blue_opfarr
    now = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M")
    with open('continuum_normalization_' + filenameblue[5:filenameblue.find(endpoint)] + '_' + now + '.txt','a') as handle:
        header = str(filenameblue) + ',' + ',' + dafile + '\n Columns structured as blue then red. If no red file, only blue data given. Columns are: wavelengths, optimized spectra, polynomial fit to data, \n model fluxes, polynomial fit to model, response, continuum normalized spectrum.'
        np.savetxt(handle,bigarray,fmt='%f',header=header)
    
    

Ni = 4. #Number of extensions
Nx1 = len(fcorr_wd_blue_opfarr)
if redfile:
    Nx2 = len(fcorr_wd_red_opfarr)
Ny = 1. #All 1D spectra


header1 = st.readheader(filenameblue)
header1.set('STANDARD',dafile,'DA Model for Continuum Calibration')
header1.set('SPECPOLY',spec_poly_order_blue,'Spectrum polynomial order for Continuum Calib.')
header1.set('MODPOLY',mod_poly_order_blue,'Model polynomial order for Continuum Calib.')

data1 = np.empty(shape = (Ni,Ny,Nx1))
data1[0,:,:] = fcorr_wd_blue_opfarr 
data1[1,:,:] = fcorr_wd_blue_farr
data1[2,:,:] = fcorr_wd_blue_sky
data1[3,:,:] = fcorr_wd_blue_sigma 



loc1 = filenameblue.find('.ms.fits')
newname1 = filenameblue[0:loc1] + '_flux_model.ms.fits'
clob = False
mylist = [True for f in os.listdir('.') if f == newname1]
exists = bool(mylist)

if exists:
    print 'File %s already exists.' % newname1
    nextstep = raw_input('Do you want to overwrite or designate a new name (overwrite/new)? ')
    if nextstep == 'overwrite':
        clob = True
        exists = False
    elif nextstep == 'new':
        newname1 = raw_input('New file name: ')
        exists = False
    else:
        exists = False

newim1 = pf.PrimaryHDU(data=data1,header=header1)
newim1.writeto(newname1,clobber=clob)

if redfile:
    header2 = st.readheader(filenamered)
    header2.set('STANDARD',dafile,'DA Model for Continuum Calibration')
    header1.set('SPECPOLY',spec_poly_order_red,'Spectrum polynomial order for Continuum Calibration')
    header1.set('MODPOLY',mod_poly_order_red,'Model polynomial order for Continuum Calibration')
    data2 = np.empty(shape = (Ni,Ny,Nx2))
    data2[0,:,:] = fcorr_wd_red_opfarr 
    data2[1,:,:] = fcorr_wd_red_farr
    data2[2,:,:] = fcorr_wd_red_sky
    data2[3,:,:] = fcorr_wd_red_sigma 
    
    
    loc2 = filenamered.find('.ms.fits')
    newname2 = filenamered[0:loc2] + '_flux_model.ms.fits'
    clob = False
    mylist = [True for f in os.listdir('.') if f == newname2]
    exists = bool(mylist)

    if exists:
        print 'File %s already exists.' % newname1
        nextstep = raw_input('Do you want to overwrite or designate a new name (overwrite/new)? ')
        if nextstep == 'overwrite':
            clob = True
            exists = False
        elif nextstep == 'new':
            newname2 = raw_input('New file name: ')
            exists = False
        else:
            exists = False
    
    newim2 = pf.PrimaryHDU(data=data2,header=header2)
    newim2.writeto(newname2,clobber=clob)





