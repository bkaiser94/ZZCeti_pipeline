"""
This program contains various definitions and commonly done things
for spectra for the ZZ CETI pipeline.

Written primarily by JT Fuchs
Based on pySALT

"""

import pyfits as fits
import numpy as np
import scipy
from scipy.interpolate import InterpolatedUnivariateSpline as interpo


class spectrum(object):

    def __init__(self,opfarr,farr,sky,sigma,warr):
        self.opfarr = opfarr
        self.farr = farr
        self.sky = sky
        self.sigma = sigma
        self.warr = warr

# ===========================================================================

class standard(object):

    def __init__(self,warr,magarr,wbin):
        self.warr = warr
        self.magarr = magarr
        self.wbin = wbin

# ===========================================================================


def readspectrum(specfile):
    """ Given a specfile, read in the spectra and return
    a spectrum object consisting of
    opfar(optimally extracted spectrum),farr(raw extracted spectrum),sky(background),sigma(sigma spectrum)
    """

    spec = fits.open(specfile)
    opfarr = spec[0].data[0,0,:]
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

    result = spectrum(opfarr,farr,sky,sigma,warr)
    return result,airmass,exptime,specdeltawav

# ===========================================================================


def readheader(specfile):
    spec = fits.open(specfile)
    #Delete the parts of the header that are not uniform in Goodman. These are primarily the parts that contain degree symbols.
    header = spec[0].header
    del header['param0']
    del header['param61']
    del header['param62']
    del header['param63']
    return header

# ===========================================================================


def readstandard(stdfile):
    warr,magarr,wbin = np.genfromtxt(stdfile,unpack=True)
    result = standard(warr,magarr,wbin)
    return result

# ===========================================================================


def magtoflux(marr, fzero):
    """Convert from magnitude to flux
    marr - input array in mags
    fzero - zero point for convertion
    """
    return fzero * 10. ** (-0.4 * marr)

# ===========================================================================


def fnutofwave(warr,farr):
    """Converts farr in ergs/s/cm2/Hz to ergs/s/cm2/A"""
    c = 2.99792458e18 #speed of light in Angstroms/s
    return farr * c / warr**2.

# ===========================================================================


def sum_std(std_warr,wbin,spec_warr,spec_farr):
    #Sum the standard star spectrum into the same bins as the flux
    #calibration file.
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
        total = np.sum(fluxnew)
        if n == 0:
            result = [total]
        if n > 0:
            result.append(total)
        #blah = np.asarray(result)
        #print blah[n]
        n += 1.
    return np.asarray(result)

# ===========================================================================


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

# ===========================================================================


def cal_spec(counts,sens,exptime,disp):
    #Calibrates a observed star using a sensitivity function
    flux = (counts) / (exptime * disp * 10.**(sens/2.5))
    return flux

# ===========================================================================


def resamplespec(w1, w0, spec0, oversamp=100):
    """
    Resample a spectrum while conserving flux density.
    Written by Ian Crossfield: www.lpl.arizona.edu/~ianc/python/index.html

    :INPUTS:
      w1 : sequence
        new wavelength grid (i.e., center wavelength of each pixel)

      w0 : sequence
        old wavelength grid (i.e., center wavelength of each pixel)

      spec0 : sequence
        old spectrum (e.g., flux density or photon counts)

      oversamp : int
        factor by which to oversample input spectrum prior to
        rebinning.  The worst fractional precision you achieve is
        roughly 1./oversamp.

    :NOTE: 
      Format is the same as :func:`numpy.interp`
      
    :REQUIREMENTS:
      :doc:`tools`

    """
    #from tools import errxy

    # 2012-04-25 18:40 IJMC: Created
    nlam = len(w0)
    x0 = np.arange(nlam, dtype=float)
    x0int = np.arange((nlam-1.)*oversamp + 1., dtype=float)/oversamp
    w0int = np.interp(x0int, x0, w0)
    spec0int = np.interp(w0int, w0, spec0)/oversamp

    # Set up the bin edges for down-binning
    maxdiffw1 = np.diff(w1).max()
    w1bins = np.concatenate(([w1[0] - maxdiffw1], 
                             .5*(w1[1::] + w1[0:-1]), \
                                 [w1[-1] + maxdiffw1]))
    # Bin down the interpolated spectrum:
    junk, spec1, junk2, junk3 = errxy(w0int, spec0int, w1bins, xmode=None, ymode='sum', xerr=None, yerr=None)

    return spec1

# ===========================================================================



def errxy(x,y,xbins, xmode='mean', ymode='mean', xerr='minmax', yerr='sdom', clean=None, binfactor=None, verbose=False,returnstats=False, timing=False):
    """Bin down datasets in X and Y for errorbar plotting
    Written by Ian Crossfield: www.lpl.arizona.edu/~ianc/python/index.html
    
    :INPUTS:
       x -- (array) independent variable data

       y -- (array) dependent variable data

       xbins -- (array) edges of bins, in x-space.  Only x-data
                between two bin edges will be used.  Thus if M bin
                edges are entered, (M-1) datapoints will be returned.
                If xbins==None, then no binning is done.

    :OPTIONAL INPUT:
       xmode/ymode -- (str) method to aggregate x/y data into datapoints:
              'mean' -- use numpy.mean
              'median' -- use numpy.median
              'sum' -- use numpy.sum
              None -- don't compute; return the empty list []

       xerr/yerr -- (str) method to aggregate x/y data into errorbars
              'std' -- sample standard deviation (numpy.std)
              'sdom' -- standard deviation on the mean; i.e., std/sqrt(N)
              'minmax' -- use full range of data in the bin
              None -- don't compute; return the empty list []

       binfactor -- (int) If not None, average over this many
              consecutive values instead of binning explicitly by
              time-based bins.  Can also be a sequence, telling the
              number of values over which to average.  E.g.,
              binfactor=[10,10,20] will bin over the first 10 points,
              the second 10 points, and the next 20 points.

       clean -- (dict) keyword options to clean y-data ONLY, via
                analysis.removeoutliers, with an additional "nsigma"
                keyword.  See removeoutliers for more information.
                E.g.:  clean=dict(nsigma=5,remove='both',niter=1)

    :OUTPUTS:     a tuple of four arrays to be passed to matplotlib.pyplot.errorbar:
       xx -- locations of the aggregated x-datapoint in each bin

       yy -- locations of the aggregated y-datapoint in each bin

       xerr -- x-errorbars

       yerr -- y-errorbars

    :EXAMPLE:
      ::

          x = hstack((arange(10), arange(20)+40))
          y = randn(len(x))
          xbins = [-1,15,70]
          xx,yy,xerr,yerr = errxy(x,y,xbins)
          plot(x,y, '.b')
          errorbar(xx,yy,xerr=xerr,yerr=yerr, fmt='or')
      
    :NOTES:

       To just bin down uncleaned data (i.e., no 'error' terms
          returned), set clean, xerr, yerr to None.  However, when
          computing all values (xerr and yerr not None) it is faster
          to set clean to some rediculous value, i.e.,
          clean=dict(niter=0, nsigma=9e99).  This probably means more
          optimization could be done.

       Be sure you call the errorbar function using the keywords xerr
          and yerr, since otherwise the default order of inputs to the
          function is (x,y,yerr,xerr).

       Data 'x' are determined to be in a bin with sides (L, R) when
          satisfying the condition (x>L) and (x<=R)

    :SEE ALSO:  matplotlib.pyplot.errorbar, :func:`analysis.removeoutliers`

    :REQUIREMENTS:  :doc:`numpy`, :doc:`analysis`
    """
    # 2009-09-29 20:07 IJC: Created w/mean-median and std-sdom-minmax.
    # 2009-12-14 16:01 IJC: xbins can be 'None' for no binning.
    # 2009-12-15 10:09 IJC: Added "binfactor" option.
    # 2009-12-22 09:56 IJC: "binfactor" can now be a sequence.
    # 2009-12-29 01:16 IJC: Fixed a bug with binfactor sequences.
    # 2010-04-29 09:59 IJC: Added 'returnstats' feature
    # 2010-10-19 16:25 IJC: Added 'sum' option for x-data
    # 2011-03-22 12:57 IJC: Added 'none' option for data and errors
    # 2012-03-20 16:33 IJMC: Fixed bug; xmode=='none' now works.
    # 2012-03-27 14:00 IJMC: Now using np.digitize -- speed boost.
    #                        Rewrote code to optimize (somewhat),
    #                        cleaned up 'import' statements.
    # 2012-04-08 15:57 IJMC: New speed boost from adopting
    #                        numpy.histogram-like implementation:
    #                        numpy.searchsorted, etc.

    
    #from analysis import removeoutliers
    import numpy as np

    if timing:
        import time
        tic = time.time()

    def sdom(data):
        """Return standard deviation of the mean."""
        return np.std(data)/np.sqrt(data.size)

    def getcenter(data, cmode):
        """Get data center based on mode.  Helper function."""
        if cmode is None:
            ret = 0
        elif cmode=='mean':
            ret = np.mean(data)
        elif cmode=='median':
            ret = np.median(data)
        elif cmode=='sum':
            ret = np.sum(data)
        return ret

    def geterr(data, emode, cmode):  
        """Get errorbar. Helper function."""
        if emode is None:
            ret = []
        elif emode=='std':
            ret = np.std(data)
        elif emode=='sdom':
            ret = sdom(data)
        elif emode=='minmax':
            if len(data)==0:
                ret = [np.nan, np.nan]
            else:
                center = getcenter(data,cmode)
                ret = [center-min(data), max(data)-center]
        return ret

    def cleandata(data, clean, returnstats=False):
        """Clean data using removeoutliers. Helper function."""
        init_count = np.array(data).size

        if clean==None: # Don't clean at all!
            #clean = dict(nsigma=1000, niter=0)
            if returnstats:
                ret = data, (init_count, init_count)
            else:
                ret = data

        else:  # Clean the data somehow ('clean' must be a dict)
            if not clean.has_key('nsigma'):
                clean.update(dict(nsigma=99999))
            data = removeoutliers(data, **clean)
            if returnstats:
                ret = data, (init_count, np.array(data).size)
            else:
                ret = data

        return ret

    if timing:
        print "%1.3f sec since starting function; helpers defined" % (time.time() - tic)

    ####### Begin main function ##########
    sorted_index = np.argsort(x)
    x = np.array(x, copy=False)[sorted_index]
    y = np.array(y, copy=False)[sorted_index]
    #x = np.array(x,copy=True).ravel()
    #y = np.array(y,copy=True).ravel()
    xbins = np.array(xbins,copy=True).ravel()
    if xbins[0]==None and binfactor==None:
        if returnstats ==False:
            ret = x, y, np.ones(x.shape)*np.nan, np.ones(y.shape)*np.nan
        else:
            ret = x, y, np.ones(x.shape)*np.nan, np.ones(y.shape)*np.nan, (x.size, x.size)
        return ret

    if binfactor==None:  # used passed-in 'xbins'
        xbins = np.sort(xbins)
    elif hasattr(binfactor,'__iter__'): # use variable-sized bins
        binfactor = np.array(binfactor).copy()
        sortedx = np.sort(x)
        betweens = np.hstack((x.min()-1, 0.5*(sortedx[1::]+sortedx[0:len(x)-1]), x.max()+1))
        xbins = []
        counter = 0
        for ii in range(len(binfactor)):
            thisbin = betweens[counter]
            xbins.append(thisbin)
            counter += binfactor[ii]
        xbins.append(x.max() + 1)
    else: # bin down by the same factor throughout
        binfactor = int(binfactor)
        sortedx = np.sort(x)
        betweens = np.hstack((x.min()-1, 0.5*(sortedx[1::]+sortedx[0:len(x)-1]), x.max()+1))
        xbins = betweens[::binfactor]

    if timing:
        print "%1.3f sec since starting function; bins defined" % (time.time() - tic)

    nbins = len(xbins)-1

    arraynan = np.array([np.nan])

    exx = []
    eyy = []
    xx = np.zeros(nbins)
    yy = np.zeros(nbins)
    yy2 = np.zeros(nbins)

    init_count, final_count = y.size, 0
    if timing:
        setuptime = 0
        xdatatime = 0
        ydatatime = 0
        statstime = 0

    #import pylab as py
    #xxx = np.sort(x)

    if timing: tic1 = time.time()
    #inds = np.digitize(x, xbins)
    inds2 = [[x.searchsorted(xbins[ii], side='left'), \
                  x.searchsorted(xbins[ii+1], side='left')] for ii in range(nbins)]
    if timing: setuptime += (time.time() - tic1)
    #pdb.set_trace()
    #bin_means = [data[digitized == i].mean() for i in range(1, len(bins))]



    dox = xmode is not None 
    doy = ymode is not None 
    doex = xerr is not None 
    doey = yerr is not None 

    if clean is None:
        if timing: tic3 = time.time()
        if dox: exec ('xfunc = np.%s' % xmode) in locals()
        if doy: exec ('yfunc = np.%s' % ymode) in locals()
        for ii in range(nbins):
            #index = inds==(ii+1)
            if dox:
                #xx[ii] = xfunc(x[index])
                xx[ii] = xfunc(x[inds2[ii][0]:inds2[ii][1]])
            if doy:
                #yy[ii] = yfunc(y[index])
                yy[ii] = yfunc(y[inds2[ii][0]:inds2[ii][1]])
            if doex:
                #exx.append(geterr(x[index], xerr, xmode))
                exx.append(geterr(x[inds2[ii][0]:inds2[ii][1]], xerr, xmode))
            if doey:
                #eyy.append(geterr(y[index], yerr, ymode))
                eyy.append(geterr(y[inds2[ii][0]:inds2[ii][1]], yerr, ymode))

        if timing: statstime += (time.time() - tic3)
        #pdb.set_trace()
    else:
        for ii in range(nbins):
            if timing: tic1 = time.time()
            #index = inds==(ii+1)
            if timing: setuptime += (time.time() - tic1)

            if timing: tic2 = time.time()
            xdata = x[inds2[ii][0]:inds2[ii][1]]
            if timing: xdatatime += (time.time() - tic2)

            if timing: tic25 = time.time()
            if ymode is None and yerr is None:  # We're free to ignore the y-data:
                ydata = arraynan
            else:  # We have to compute something with the y-data:
                if clean is not None:
                    ydata, retstats = cleandata(y[inds2[ii][0]:inds2[ii][1]], clean, returnstats=True)
                    if returnstats:
                        final_count += retstats[1]
                else:  # We don't have to clean the data
                    ydata = y[inds2[ii][0]:inds2[ii][1]]
                    if returnstats:
                        final_count += ydata.size
            if timing: ydatatime += (time.time() - tic25)

            if timing: tic3 = time.time()
            xx[ii] = getcenter(xdata,xmode)
            if timing: tic4 = time.time()
            yy[ii] = getcenter(ydata,ymode)
            if timing: tic5 = time.time()
            exx.append(geterr(  xdata,xerr,xmode))
            if timing: tic6 = time.time()
            eyy.append(geterr(  ydata,yerr,ymode))
            if timing: tic7 = time.time()
            if timing: statstime += (time.time() - tic3)
            #exx[ii] = geterr(  xdata,xerr,xmode)
            #eyy[ii] = geterr(  ydata,yerr,ymode)

    if timing:
        print "%1.3f sec for setting up bins & indices..." % setuptime
        print "%1.3f sec for getting x data clean and ready." % xdatatime
        print "%1.3f sec for getting y data clean and ready." % ydatatime
        #print "%1.3f sec for computing x-data statistics." % (tic4-tic3)
        #print "%1.3f sec for computing y-data statistics." % (tic5-tic4)
        #print "%1.3f sec for computing x-error statistics." % (tic6-tic5)
        #print "%1.3f sec for computing y-error statistics." % (tic7-tic6)


        print "%1.3f sec for computing statistics........." % statstime

    if timing:
        print "%1.3f sec since starting function; uncertainties defined" % (time.time() - tic)



    #xx = array(xx)
    #yy = array(yy)
    exx = np.array(exx).transpose()  # b/c 2D if minmax option used
    eyy = np.array(eyy).transpose()  # b/c 2D if minmax option used

    #pdb.set_trace()
 
    if returnstats:
        ret= xx,yy,exx,eyy,(init_count, final_count)
    else:
        ret = xx,yy,exx,eyy

    #print 'tools: returnstats, len(ret)>>', returnstats, len(ret)
    if timing:
        print "%1.3f sec since starting function; returning" % (time.time() - tic)

    return ret 

# ===========================================================================


def removeoutliers(data, nsigma, remove='both', center='mean', niter=500, retind=False, verbose=False):
    """Strip outliers from a dataset, iterating until converged.
    Written by Ian Crossfield: www.lpl.arizona.edu/~ianc/python/index.html
    
    :INPUT:
      data -- 1D numpy array.  data from which to remove outliers.

      nsigma -- positive number.  limit defining outliers: number of
                standard deviations from center of data.

    :OPTIONAL INPUTS:               
      remove -- ('min'|'max'|'both') respectively removes outliers
                 below, above, or on both sides of the limits set by
                 nsigma.

      center -- ('mean'|'median'|value) -- set central value, or
                 method to compute it.

      niter -- number of iterations before exit; defaults to Inf,
               which can occasionally result in empty arrays returned

      retind -- (bool) whether to return index of good values as
                second part of a 2-tuple.

    :EXAMPLE: 
       ::

           from numpy import hist, linspace, randn
           from analysis import removeoutliers
           data = randn(1000)
           hbins = linspace(-5,5,50)
           d2 = removeoutliers(data, 1.5, niter=1)
           hist(data, hbins)
           hist(d2, hbins)

       """
    # 2009-09-04 13:24 IJC: Created
    # 2009-09-24 17:34 IJC: Added 'retind' feature.  Tricky, but nice!
    # 2009-10-01 10:40 IJC: Added check for stdev==0
    # 2009-12-08 15:42 IJC: Added check for isfinite

    from numpy import median, ones, isfinite

    def getcen(data, method):
        "Get central value of a 1D array (helper function)"
        if method.__class__==str:
            if method=='median':
                cen = median(data)
            else:
                cen = data.mean()
        else:
            cen = method
        return cen

    def getgoodindex(data, nsigma, center, stdev, remove):
        "Get number of outliers (helper function!)"
        if stdev==0:
            distance = data*0.0
        else:
            distance = (data-center)/stdev
        if remove=='min':
            goodind = distance>-nsigma
        elif remove=='max':
            goodind = distance<nsigma
        else:
            goodind = abs(distance)<=nsigma
        return goodind

    data = data.ravel().copy()

    ndat0 = len(data)
    ndat = len(data)
    iter=0
    goodind = ones(data.shape,bool)
    goodind *= isfinite(data)
    while ((ndat0<>ndat) or (iter==0)) and (iter<niter) and (ndat>0) :
        ndat0 = len(data[goodind])
        cen = getcen(data[goodind], center)
        stdev = data[goodind].std()
        thisgoodind = getgoodindex(data[goodind], nsigma, cen, stdev, remove)
        goodind[find(goodind)] = thisgoodind
        if verbose:
            print "cen>>",cen
            print "std>>",stdev
        ndat = len(data[goodind])
        iter +=1
        if verbose:
            print ndat0, ndat
    if retind:
        ret = data[goodind], goodind
    else:
        ret = data[goodind]
    return ret


# ===========================================================================


def resample(old_dispersion, new_dispersion):
    """
    Written by Andy Casey: https://github.com/andycasey
    This program is found under sick/specutils.py

    The output is a compressed sparse column matrix that tells the ratio between the old and new binning. To get the fluxes of the new binning: newflux = np.dot(oldflux,output.toarray())

    Resample a spectrum to a new dispersion map while conserving total flux.
    :param old_dispersion:
        The original dispersion array.
    :type old_dispersion:
        :class:`numpy.array`
    :param new_dispersion:
        The new dispersion array to resample onto.
    :type new_dispersion:
        :class:`numpy.array`
    """

    data = []
    old_px_indices = []
    new_px_indices = []
    for i, new_wl_i in enumerate(new_dispersion):

        # These indices should span just over the new wavelength pixel.
        indices = np.unique(np.clip(
            old_dispersion.searchsorted(new_dispersion[i:i + 2], side="left") \
                + [-1, +1], 0, old_dispersion.size - 1))
        N = np.ptp(indices)

        if N == 0:
            # 'Fake' pixel.
            data.append(np.nan)
            new_px_indices.append(i)
            old_px_indices.extend(indices)
            continue

        # Sanity checks.
        assert (old_dispersion[indices[0]] <= new_wl_i \
            or indices[0] == 0)
        assert (new_wl_i <= old_dispersion[indices[1]] \
            or indices[1] == old_dispersion.size - 1)

        fractions = np.ones(N)

        # Edges are handled as fractions between rebinned pixels.
        _ = np.clip(i + 1, 0, new_dispersion.size - 1)
        lhs = old_dispersion[indices[0]:indices[0] + 2]
        rhs = old_dispersion[indices[-1] - 1:indices[-1] + 1]
        fractions[0]  = (lhs[1] - new_dispersion[i])/np.ptp(lhs)
        fractions[-1] = (new_dispersion[_] - rhs[0])/np.ptp(rhs)

        # Being binned to a single pixel. Prevent overflow from fringe cases.
        fractions = np.clip(fractions, 0, 1)
        fractions /= fractions.sum()

        data.extend(fractions) 
        new_px_indices.extend([i] * N) # Mark the new pixel indices affected.
        old_px_indices.extend(np.arange(*indices)) # And the old pixel indices.

    return scipy.sparse.csc_matrix((data, (old_px_indices, new_px_indices)),
        shape=(old_dispersion.size, new_dispersion.size))
