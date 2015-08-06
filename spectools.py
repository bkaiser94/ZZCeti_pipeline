"""
This program contains various definitions and commonly done things
for spectra for the ZZ CETI pipeline.

Written primarily by JT Fuchs
Based on pySALT

"""

import pyfits as fits
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as interpo


class spectrum(object):

    def __init__(self,var_farr,farr,sky,sigma,warr):
        self.var_farr = var_farr
        self.farr = farr
        self.sky = sky
        self.sigma = sigma
        self.warr = warr

class standard(object):

    def __init__(self,warr,magarr,wbin):
        self.warr = warr
        self.magarr = magarr
        self.wbin = wbin

def readspectrum(specfile):
    """ Given a specfile, read in the spectra and return
    a spectrum object
    warr, farr, farr_err
    """

    spec = fits.open(specfile)
    var_farr = spec[0].data[0,0,:]
    farr = spec[0].data[1,0,:]
    sky = spec[0].data[2,0,:]
    sigma = spec[0].data[3,0,:]

    #Read in header info
    airmass = spec[0].header['airmass']
    exptime = spec[0].header['exptime']

    #Set up wavelengths
    specwav0 = spec[0].header['crval1'] #Grab the leftmost wavelength coordinate
    specdeltawav = spec[0].header['cd1_1'] #Grab the delta coordinate
    warr = np.zeros(len(farr)) #Fill an array with appropriate wavelength values
    warr[0] = specwav0
    ival = np.arange(1,len(farr))
    for i in ival:
        warr[i] = warr[i-1] + specdeltawav

    result = spectrum(var_farr,farr,sky,sigma,warr)
    return result,airmass,exptime,specdeltawav

def readheader(specfile):
    spec = fits.open(specfile)
    #Delete the parts of the header that are not uniform in Goodman. These are primarily the parts that contain degree symbols.
    header = spec[0].header
    del header['param0']
    del header['param61']
    del header['param62']
    del header['param63']
    return header

def readstandard(stdfile):
    warr,magarr,wbin = np.genfromtxt(stdfile,unpack=True)
    result = standard(warr,magarr,wbin)
    return result

def magtoflux(marr, fzero):
    """Convert from magnitude to flux
    marr - input array in mags
    fzero - zero point for convertion
    """
    return fzero * 10. ** (-0.4 * marr)

def fnutofwave(warr,farr):
    """Converts farr in ergs/s/cm2/Hz to ergs/s/cm2/A"""
    c = 2.99792458e18 #speed of light in Angstroms/s
    return farr * c / warr**2.

def sum_std(std_warr,wbin,spec_warr,spec_farr,dw,binsize):
    #Sum the standard star spectrum into the same bins as the flux
    #calibration file.
    #numbins is a scaling factor for how the flux should change based on the relative sizes of the old and new bins. This assumes the slope does not change significantly over the size of either bin. We have to rebin because it makes the summing better as we are including the whole size of the bin instead of fractional parts.
    n = 0
    for lambdas in std_warr:
        low = lambdas - wbin[n]/2.
        high = lambdas + wbin[n]/2.
        #print low,high
        c = np.where(spec_warr >= low)
        d = np.where(spec_warr <= high)
        lowflux = np.asarray(c)
        highflux = np.asarray(d)
        index = np.intersect1d(lowflux,highflux)
        fluxnew = spec_farr[index]
        wavenew = spec_warr[index]
        #print wavenew[0],wavenew[-1]
        numbins = binsize/dw
        #print numbins
        total = numbins * np.sum(fluxnew) #numbins is a scaling factor 
        if n == 0:
            result = [total]
        if n > 0:
            result.append(total)
        #blah = np.asarray(result)
        #print blah[n]
        n += 1.
    return np.asarray(result)

def sensfunc(obs_counts,std_flux,exptime,bins,airmass):
    #This function calculates the sensitivity curve for the spectrum
    #It is calculated by:
    #C = 2.5 * log(obs_counts/ (exptime * bin * std_flux)) + airmass * extinction
    n = 0
    for counts in obs_counts:
        cal = 2.5 * np.log10(counts/ (exptime * bins[n] * std_flux[n]))
        if n == 0:
            sens = [cal]
        if n > 0:
            sens.append(cal)
        n += 1
    sensitivity = np.asarray(sens)
    return sensitivity

def cal_spec(counts,sens,exptime,disp):
    #Calibrates a observed star using a sensitivity function
    flux = (counts) / (exptime * disp * 10.**(sens/2.5))
    return flux
