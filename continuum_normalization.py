'''
Written by J. T. Fuchs, UNC, 2016. Some initial work done by E. Dennihy, UNC.

Continuum normalizes a ZZ Ceti spectrum to match a DA model. The response function is determined by dividing the observed spectrum by the model spectrum. The response function is fitted with a polynomial. The observed spectrum is then divided by this fit to deliver the continuum normalized spectrum.

To run: (Red filename is optional)
python model_calibration.py bluefilename redfilename
python model_calibration.py wtfb.wd1401-147_930_blue.ms.fits wtfb.wd1401-147_930_red.ms.fits

:INPUTS:
     bluefilename: string, filename of wavelength calibrated ZZ Ceti blue spectrum

:OPTIONAL:
     redfilename: string, filename of wavelength calibrated ZZ Ceti red spectrum

:OUTPUTS:
     continuum normalized spectrum: '_flux_model' added to filename. Name of model used written to header. 
     
     normalization_ZZCETINAME_DATE.txt: File for diagnostics. ZZCETINAME is name of the ZZ Ceti spectrum supplied. DATE is the current date and time. Columns are: blue wavelengths, blue response all data, blue masked wavelengths, blue masked response data, blue response fit, red wavelengths, red response all data, red masked wavelengths, red masked response data, red response fit

'''


import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as inter
#import pyfits as fits
import astropy.io.fits as fits
import spectools as st
import os
import sys
import datetime
from scipy.interpolate import UnivariateSpline


#=============================================
def extinction_correction(lams, flux, airmass):
    # Function inputs are wavelengths and flux values for the spectrum as well 
    # as the airmass the spectrum was measured at
    
    # wavelength-dependent extinction coefficients from CTIO
    # Strizinger et. al. 2005
    ctio_lams = [3050.0, 3084.6500000000001, 3119.3099999999999, 3153.96, 3188.6100000000001, 3223.27, 3257.9200000000001, 3292.5700000000002, 3327.23, 3361.8800000000001, 3396.54, 3431.1900000000001, 3465.8400000000001, 3500.5, 3535.1500000000001, 3569.8000000000002, 3604.46, 3639.1100000000001, 3673.7600000000002, 3708.4200000000001, 3743.0700000000002, 3777.7199999999998, 3812.3800000000001, 3847.0300000000002, 3881.6900000000001, 3916.3400000000001, 3950.9899999999998, 3985.6500000000001, 4020.3000000000002, 4054.9499999999998, 4089.6100000000001, 4124.2600000000002, 4158.9099999999999, 4193.5699999999997, 4228.2200000000003, 4262.8699999999999, 4297.5299999999997, 4332.1800000000003, 4366.8299999999999, 4401.4899999999998, 4436.1400000000003, 4470.79, 4505.4499999999998, 4540.1000000000004, 4574.7600000000002, 4609.4099999999999, 4644.0600000000004, 4678.7200000000003, 4713.3699999999999, 4748.0200000000004, 4782.6800000000003, 4817.3299999999999, 4851.9799999999996, 4886.6400000000003, 4921.29, 4955.9399999999996, 4990.6000000000004, 5025.25, 5059.9099999999999, 5094.5600000000004, 5129.21, 5163.8699999999999, 5198.5200000000004, 5233.1700000000001, 5267.8299999999999, 5302.4799999999996, 5337.1300000000001, 5371.79, 5406.4399999999996, 5441.0900000000001, 5475.75, 5510.3999999999996, 5545.0500000000002, 5579.71, 5614.3599999999997, 5649.0200000000004, 5683.6700000000001, 5718.3199999999997, 5752.9799999999996, 5787.6300000000001, 5822.2799999999997, 5856.9399999999996, 5891.5900000000001, 5926.2399999999998, 5960.8999999999996, 5995.5500000000002, 6030.1999999999998, 6064.8599999999997, 6099.5100000000002, 6134.1700000000001, 6168.8199999999997, 6203.4700000000003, 6238.1300000000001, 6272.7799999999997, 6307.4300000000003, 6342.0900000000001, 6376.7399999999998, 6411.3900000000003, 6446.0500000000002, 6480.6999999999998, 6482.8500000000004, 6535.3800000000001, 6587.9099999999999, 6640.4399999999996, 6692.96, 6745.4899999999998, 6798.0200000000004, 6850.5500000000002, 6903.0699999999997, 6955.6000000000004, 7008.1300000000001, 7060.6499999999996, 7113.1800000000003, 7165.71, 7218.2399999999998, 7270.7600000000002, 7323.29, 7375.8199999999997, 7428.3500000000004, 7480.8699999999999, 7533.3999999999996, 7585.9300000000003, 7638.4499999999998, 7690.9799999999996, 7743.5100000000002, 7796.04, 7848.5600000000004, 7901.0900000000001, 7953.6199999999999, 8006.1499999999996, 8058.6700000000001, 8111.1999999999998, 8163.7299999999996, 8216.25, 8268.7800000000007, 8321.3099999999995, 8373.8400000000001, 8426.3600000000006, 8478.8899999999994, 8531.4200000000001, 8583.9500000000007, 8636.4699999999993, 8689.0, 8741.5300000000007, 8794.0499999999993, 8846.5799999999999, 8899.1100000000006, 8951.6399999999994, 9004.1599999999999, 9056.6900000000005, 9109.2199999999993, 9161.75, 9214.2700000000004, 9266.7999999999993, 9319.3299999999999, 9371.8500000000004, 9424.3799999999992, 9476.9099999999999, 9529.4400000000005, 9581.9599999999991, 9634.4899999999998, 9687.0200000000004, 9739.5499999999993, 9792.0699999999997, 9844.6000000000004, 9897.1299999999992, 9949.6499999999996, 10002.200000000001, 10054.700000000001, 10107.200000000001, 10159.799999999999, 10212.299999999999, 10264.799999999999, 10317.299999999999, 10369.9, 10422.4, 10474.9, 10527.5, 10580.0, 10632.5, 10685.0, 10737.6, 10790.1, 10842.6, 10895.1, 10947.700000000001, 11000.200000000001]
    ctio_ext = [1.395, 1.2830000000000001, 1.181, 1.0880000000000001, 1.004, 0.92900000000000005, 0.86099999999999999, 0.80099999999999993, 0.748, 0.69999999999999996, 0.65900000000000003, 0.623, 0.59099999999999997, 0.56399999999999995, 0.54000000000000004, 0.52000000000000002, 0.502, 0.48700000000000004, 0.47299999999999998, 0.46000000000000002, 0.44799999999999995, 0.436, 0.42499999999999999, 0.41399999999999998, 0.40200000000000002, 0.39100000000000001, 0.38100000000000001, 0.37, 0.35999999999999999, 0.34899999999999998, 0.33899999999999997, 0.33000000000000002, 0.32100000000000001, 0.313, 0.30399999999999999, 0.29600000000000004, 0.28899999999999998, 0.28100000000000003, 0.27399999999999997, 0.26700000000000002, 0.26000000000000001, 0.254, 0.247, 0.24100000000000002, 0.23600000000000002, 0.23000000000000001, 0.22500000000000001, 0.22, 0.215, 0.20999999999999999, 0.20600000000000002, 0.20199999999999999, 0.19800000000000001, 0.19399999999999998, 0.19, 0.187, 0.184, 0.18100000000000002, 0.17800000000000002, 0.17600000000000002, 0.17300000000000001, 0.17100000000000001, 0.16899999999999998, 0.16699999999999998, 0.16600000000000001, 0.16399999999999998, 0.16300000000000001, 0.16200000000000001, 0.16, 0.159, 0.158, 0.158, 0.157, 0.156, 0.155, 0.155, 0.154, 0.153, 0.153, 0.152, 0.151, 0.151, 0.14999999999999999, 0.14899999999999999, 0.14899999999999999, 0.14800000000000002, 0.14699999999999999, 0.14599999999999999, 0.14400000000000002, 0.14300000000000002, 0.14199999999999999, 0.14000000000000001, 0.13800000000000001, 0.13600000000000001, 0.13400000000000001, 0.13200000000000001, 0.129, 0.126, 0.12300000000000001, 0.12, 0.12, 0.115, 0.111, 0.107, 0.10300000000000001, 0.099000000000000005, 0.096000000000000002, 0.091999999999999998, 0.088000000000000009, 0.085000000000000006, 0.08199999999999999, 0.078, 0.074999999999999997, 0.072000000000000008, 0.069000000000000006, 0.066000000000000003, 0.064000000000000001, 0.060999999999999999, 0.057999999999999996, 0.055999999999999994, 0.052999999999999999, 0.050999999999999997, 0.049000000000000002, 0.047, 0.044999999999999998, 0.042999999999999997, 0.040999999999999995, 0.039, 0.037000000000000005, 0.035000000000000003, 0.034000000000000002, 0.032000000000000001, 0.029999999999999999, 0.028999999999999998, 0.027999999999999997, 0.026000000000000002, 0.025000000000000001, 0.024, 0.023, 0.022000000000000002, 0.02, 0.019, 0.019, 0.018000000000000002, 0.017000000000000001, 0.016, 0.014999999999999999, 0.014999999999999999, 0.013999999999999999, 0.013000000000000001, 0.013000000000000001, 0.012, 0.011000000000000001, 0.011000000000000001, 0.011000000000000001, 0.01, 0.01, 0.0090000000000000011, 0.0090000000000000011, 0.0090000000000000011, 0.0080000000000000002, 0.0080000000000000002, 0.0080000000000000002, 0.0069999999999999993, 0.0069999999999999993, 0.0069999999999999993, 0.0069999999999999993, 0.0069999999999999993, 0.0060000000000000001, 0.0060000000000000001, 0.0060000000000000001, 0.0060000000000000001, 0.0060000000000000001, 0.0060000000000000001, 0.0050000000000000001, 0.0050000000000000001, 0.0050000000000000001, 0.0050000000000000001, 0.0050000000000000001, 0.0050000000000000001, 0.0040000000000000001, 0.0040000000000000001, 0.0040000000000000001, 0.0040000000000000001, 0.0030000000000000001, 0.0030000000000000001, 0.0030000000000000001]

    smooth_param = 0.001
    spline_fit = UnivariateSpline(ctio_lams, ctio_ext, s=smooth_param, k=3)

    a_lambda = spline_fit(lams)

    corrected_flux = flux*(10.0**(.4*a_lambda*(1.0+airmass)))    
    
    xx = np.linspace(np.min(ctio_lams), np.max(ctio_lams), 1000)
    yy = spline_fit(xx)
    '''
    plt.figure()
    plt.scatter(ctio_lams, ctio_ext, label=smooth_param)
    plt.axvline(np.min(lams), color='g')
    plt.axvline(np.max(lams), color='g')
    plt.plot(xx, yy)
    plt.xlabel('Wavelength')
    plt.ylabel('Extinction Coefficient')
    plt.title('Gemini Extinction Coefficient Fit')
    '''
    '''
    plt.figure()
    plt.plot(lams,flux)
    plt.plot(lams,corrected_flux)
    plt.show()
    '''
    return corrected_flux


#=============================================


def normalize_now(filenameblue,filenamered,redfile,plotall=True,extinct_correct=False):
    #Read in the observed spectrum
    obs_spectrablue,airmass,exptime,dispersion = st.readspectrum(filenameblue)
    datalistblue = fits.open(filenameblue)

    if redfile:
        obs_spectrared, airmassred,exptimered,dispersionred = st.readspectrum(filenamered)
    
    #Extinction correction
    if extinct_correct:
        print 'Extinction correcting spectra.'
        plt.clf()
        plt.plot(obs_spectrablue.warr,obs_spectrablue.opfarr)
        obs_spectrablue.opfarr = extinction_correction(obs_spectrablue.warr,obs_spectrablue.opfarr,airmass)
        plt.plot(obs_spectrablue.warr,obs_spectrablue.opfarr)
        plt.show()
        
        if redfile:
            plt.clf()
            plt.plot(obs_spectrared.warr,obs_spectrared.opfarr)
            obs_spectrared.opfarr = extinction_correction(obs_spectrared.warr,obs_spectrared.opfarr,airmassred)
            plt.plot(obs_spectrared.warr,obs_spectrared.opfarr)
            plt.show()
    

    #Read in measured FWHM from header. This is used to convolve the model spectrum.
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
    cflux2blue /= 10**13. #Divide by 10**13 to scale
    if redfile:
        cflux2red = interp2(obs_spectrared.warr)
        cflux2red /= 10**13. #Divide by 10**13 to scale


    #plt.clf()
    #plt.plot(obs_spectrablue.warr,obs_spectrablue.opfarr,'b')
    #plt.plot(obs_spectrablue.warr,cflux2blue,'r')
    #if redfile:
    #    plt.plot(obs_spectrared.warr,obs_spectrared.opfarr,'b')
    #    plt.plot(obs_spectrared.warr,cflux2red,'r')
    #plt.show()

    #The response function is the observed spectrum divided by the model spectrum.
    response_blue = obs_spectrablue.opfarr/cflux2blue
    if redfile:
        response_red = obs_spectrared.opfarr/cflux2red

    '''
    plt.clf()
    plt.plot(obs_spectrablue.warr,response_blue,'k')
    if redfile:
        plt.plot(obs_spectrared.warr,response_red,'k')
    plt.show()
    '''

    #We want to mask out the Balmer line features, and the telluric line in the red spectrum. Set up the wavelength ranges to mask here.
    #balmer_features_blue = [[3745,3757],[3760,3780],[3784,3812],[3816,3856],[3865,3921],[3935,4021],[4040,4191],[4223,4460],[4691,5010]] #Keeping ends
    balmer_features_blue = [[3400,3700],[3745,3757],[3760,3780],[3784,3812],[3816,3856],[3865,3921],[3935,4021],[4040,4191],[4223,4460],[4691,5010],[5140,5500]] #Discarding ends
    balmer_features_red = [[6350,6780],[6835,6970]]

    balmer_mask_blue = obs_spectrablue.warr == obs_spectrablue.warr
    for wavrange in balmer_features_blue:
        inds = np.where((obs_spectrablue.warr > wavrange[0]) & (obs_spectrablue.warr < wavrange[1]))
        balmer_mask_blue[inds] = False

    if redfile:
        balmer_mask_red = obs_spectrared.warr == obs_spectrared.warr
        for wavrange in balmer_features_red:
            indxs = np.where((obs_spectrared.warr > wavrange[0]) & (obs_spectrared.warr < wavrange[1]))
            balmer_mask_red[indxs] = False

    spec_wav_masked_blue = obs_spectrablue.warr[balmer_mask_blue]
    response_masked_blue = response_blue[balmer_mask_blue]

    if redfile:    
        spec_wav_masked_red = obs_spectrared.warr[balmer_mask_red]
        response_masked_red = response_red[balmer_mask_red]


    #Fit the response function with a polynomial. The order of polynomial is specified first. 
    response_poly_order_blue = 7.
    response_fit_blue_poly = np.polyfit(spec_wav_masked_blue,response_masked_blue,response_poly_order_blue)
    response_fit_blue = np.poly1d(response_fit_blue_poly)


    if redfile:
        response_poly_order_red = 3.
        response_fit_red_poly = np.polyfit(spec_wav_masked_red,response_masked_red,response_poly_order_red)
        response_fit_red = np.poly1d(response_fit_red_poly)

    #Save response function
    #np.savetxt('response_model_no_extinction.txt',np.transpose([obs_spectrablue.warr,response_fit_blue(obs_spectrablue.warr),obs_spectrared.warr,response_fit_red(obs_spectrared.warr)]))
    #plt.clf()
    #plt.plot(obs_spectrablue.warr,response_fit_blue(obs_spectrablue.warr)/response_fit_blue(obs_spectrablue.warr)[1000])
    #plt.show()
    #exit()
    if plotall:
        plt.clf()
        plt.plot(obs_spectrablue.warr,response_blue,'r')
        plt.plot(spec_wav_masked_blue,response_masked_blue,'g.')
        plt.plot(obs_spectrablue.warr,response_fit_blue(obs_spectrablue.warr),'k--')
        #plt.show()
        
        #plt.clf()
        if redfile:
            plt.plot(obs_spectrared.warr,response_red,'r')
            plt.plot(spec_wav_masked_red,response_masked_red,'g.')
            plt.plot(obs_spectrared.warr,response_fit_red(obs_spectrared.warr),'k--')
        plt.show()

    #Divide by the fit to the response function to get the continuum normalized spectra. Divide every extension by the same polynomial
    fcorr_wd_blue_opfarr = obs_spectrablue.opfarr / response_fit_blue(obs_spectrablue.warr)
    fcorr_wd_blue_farr = obs_spectrablue.farr / response_fit_blue(obs_spectrablue.warr)
    fcorr_wd_blue_sky = obs_spectrablue.sky / response_fit_blue(obs_spectrablue.warr)
    fcorr_wd_blue_sigma = obs_spectrablue.sigma / response_fit_blue(obs_spectrablue.warr)


    if redfile:
        fcorr_wd_red_opfarr = obs_spectrared.opfarr / response_fit_red(obs_spectrared.warr)
        fcorr_wd_red_farr = obs_spectrared.farr / response_fit_red(obs_spectrared.warr)
        fcorr_wd_red_sky = obs_spectrared.sky / response_fit_red(obs_spectrared.warr)
        fcorr_wd_red_sigma = obs_spectrared.sigma / response_fit_red(obs_spectrared.warr)
    
    if plotall:
        plt.clf()
        plt.plot(obs_spectrablue.warr,fcorr_wd_blue_opfarr,'b')
        if redfile:
            plt.plot(obs_spectrared.warr,fcorr_wd_red_opfarr,'r')
        plt.show()
    #exit()

    #Save parameters for diagnostics
    if redfile:
        bigarray = np.zeros([len(obs_spectrablue.warr),12])
        bigarray[0:len(obs_spectrablue.warr),0] = obs_spectrablue.warr
        bigarray[0:len(response_blue),1] = response_blue
        bigarray[0:len(spec_wav_masked_blue),2] = spec_wav_masked_blue
        bigarray[0:len(response_masked_blue),3] = response_masked_blue
        bigarray[0:len(response_fit_blue(obs_spectrablue.warr)),4] = response_fit_blue(obs_spectrablue.warr)
        bigarray[0:len(fcorr_wd_blue_opfarr),5] = fcorr_wd_blue_opfarr
        bigarray[0:len(obs_spectrared.warr),6] = obs_spectrared.warr
        bigarray[0:len(response_red),7] = response_red
        bigarray[0:len(spec_wav_masked_red),8] = spec_wav_masked_red
        bigarray[0:len(response_masked_red),9] = response_masked_red
        bigarray[0:len(response_fit_red(obs_spectrared.warr)),10] = response_fit_red(obs_spectrared.warr)
        bigarray[0:len(fcorr_wd_red_opfarr),11] = fcorr_wd_red_opfarr
        now = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M")
        endpoint = '_930'
        with open('continuum_normalization_' + filenameblue[5:filenameblue.find(endpoint)] + '_' + now + '.txt','a') as handle:
            header = str(filenameblue) + ',' + str(filenamered) + ',' + dafile + '\n Columns structured as blue then red. If no red file, only blue data given. Columns are: blue wavelengths, blue response all data, blue masked wavelengths, blue masked response data, blue response fit, blue continuum-normalize flux, red wavelengths, red response all data, red masked wavelengths, red masked response data, red response fit, red continuum-normalized flux'
            np.savetxt(handle,bigarray,fmt='%f',header=header)
    if not redfile:
        bigarray = np.zeros([len(obs_spectrablue.warr),6])
        bigarray[0:len(obs_spectrablue.warr),0] = obs_spectrablue.warr
        bigarray[0:len(response_blue),1] = response_blue
        bigarray[0:len(spec_wav_masked_blue),2] = spec_wav_masked_blue
        bigarray[0:len(response_masked_blue),3] = response_masked_blue
        bigarray[0:len(response_fit_blue(obs_spectrablue.warr)),4] = response_fit_blue(obs_spectrablue.warr)
        bigarray[0:len(fcorr_wd_blue_opfarr),5] = fcorr_wd_blue_opfarr
        now = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M")
        endpoint = '_930'
        with open('continuum_normalization_' + filenameblue[5:filenameblue.find(endpoint)] + '_' + now + '.txt','a') as handle:
            header = str(filenameblue) + ',' + ',' + dafile + '\n Columns structured as blue then red. If no red file, only blue data given. Columns are: blue wavelengths, blue response all data, blue masked wavelengths, blue masked response data, blue response fit, blue continuum-normalized flux'
            np.savetxt(handle,bigarray,fmt='%f',header=header)
    
    
    #Save the continuum normalized spectra here.
    Ni = 4. #Number of extensions
    Nx1 = len(fcorr_wd_blue_opfarr)
    if redfile:
        Nx2 = len(fcorr_wd_red_opfarr)
    Ny = 1. #All 1D spectra

    #Update header
    header1 = st.readheader(filenameblue)
    header1.set('STANDARD',dafile,'DA Model for Continuum Calibration')
    header1.set('RESPPOLY',response_poly_order_blue,'Polynomial order for Response Function')
    header1.set('DATENORM',datetime.datetime.now().strftime("%Y-%m-%d"),'Date of Continuum Normalization')


    data1 = np.empty(shape = (Ni,Ny,Nx1))
    data1[0,:,:] = fcorr_wd_blue_opfarr 
    data1[1,:,:] = fcorr_wd_blue_farr
    data1[2,:,:] = fcorr_wd_blue_sky
    data1[3,:,:] = fcorr_wd_blue_sigma 

    #Check that filename does not already exist. Prompt user for input if it does.
    loc1 = filenameblue.find('.ms.fits')
    newname1 = filenameblue[0:loc1] + '_flux_model_short.ms.fits'
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
    print 'Writing ', newname1
    newim1 = fits.PrimaryHDU(data=data1,header=header1)
    newim1.writeto(newname1,clobber=clob)


    #Save the red file if it exists.
    if redfile:
        header2 = st.readheader(filenamered)
        header2.set('STANDARD',dafile,'DA Model for Continuum Calibration')
        header2.set('RESPPOLY',response_poly_order_red,'Polynomial order for Response Function')
        header2.set('DATENORM',datetime.datetime.now().strftime("%Y-%m-%d"),'Date of Continuum Normalization')
        data2 = np.empty(shape = (Ni,Ny,Nx2))
        data2[0,:,:] = fcorr_wd_red_opfarr 
        data2[1,:,:] = fcorr_wd_red_farr
        data2[2,:,:] = fcorr_wd_red_sky
        data2[3,:,:] = fcorr_wd_red_sigma 
    
    
        loc2 = filenamered.find('.ms.fits')
        newname2 = filenamered[0:loc2] + '_flux_model_short.ms.fits'
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
        print 'Writing ', newname2
        newim2 = fits.PrimaryHDU(data=data2,header=header2)
        newim2.writeto(newname2,clobber=clob)


#################################
#############################

if __name__ == '__main__':
    if len(sys.argv) == 3:
        script, filenameblue, filenamered = sys.argv
        redfile = True
    elif len(sys.argv) == 2:
        script, filenameblue = sys.argv
        filenamered = None
        redfile = False
    else:
        print '\n Incorrect number of arguments. \n'
    
    normalize_now(filenameblue,filenamered,redfile,extinct_correct=False)

        
