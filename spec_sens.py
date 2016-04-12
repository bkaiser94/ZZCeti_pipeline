"""
Written by JT Fuchs in July 2015
Based off pySALT redution routine specsens.py by S. Crawford
And reading the darned IRAF documentation

spec_sens.py calculates the calibration curve given an
observation and a standard star.

To run file:
python spec_sens.py liststandard listflux liststar
######
Each list should have the names of the stars, with blue and red exposures next to each other.
The ordering of the standard star flux files should match the order of the standard star list.
Example:

liststandard:
wtfb.LTT3218_930_blue.ms.fits
wtfb.LTT3218_930_red.ms.fits
wnb.GD50_930_blue.ms.fits
wnb.GD50_930_red.ms.fits

listflux:
mltt3218.dat
mgd50.dat

liststar
wnb.WD0122p0030_930_blue.ms.fits
wnb.WD0122p0030_930_red.ms.fits
wnb.WD0235p069_930_blue.ms.fits
wnb.WD0235p069_930_red.ms.fits

#####

Counting variables are fruits and vegetables.

Inputs: list of 1D standard star spectra,list of filename of standard star fluxes, list of  1D spectrum of observed stars

Outputs: flux calibrated file (_flux is added to the filename), sensitivity_params.txt is updated to list all the parameters used in the flux calibration.

To do:
- Read in inputs

"""

import os
import sys
import time
import numpy as np
import pyfits as pf
import spectools as st
import datetime

import matplotlib.pyplot as plt

#These functions are to help with excluding regions from the sensitivity function
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def onclick(event):
    global ix,iy
    ix, iy = event.xdata,event.ydata
    global coords
    ax.axvline(x=ix,color='k',linewidth='3')
    fig.canvas.draw()
    coords.append((ix,iy))


#Read in lists from command line
script, stdlist, fluxlist, speclist = sys.argv
#script, stdspecfile, stdfile, specfile = sys.argv
#stdspecfile = 'wnb.GD50_930_blue.ms.fits'
#stdfile = 'mgd50.dat'
#specfile = 'wnb.WD0122p0030_930_blue.ms.fits'

#Read in each standard star spectrum 
standards = np.genfromtxt(stdlist,dtype=str)
stdflux = np.genfromtxt(fluxlist,dtype=str)

#Check that the files are set up correctly to avoid mixing standards.
onion = 0
for stanspec in standards:
    quickcheck = stdflux[onion//2].lower()[1:-4] in stanspec.lower()
    if not quickcheck:
        print 'Check your standard star and flux files. They are mixed up.'
        sys.exit()
    onion += 1

orderused = np.zeros([len(standards)])
senspolys = []
airstd = np.zeros([len(standards)])
allexcluded = [[None] for i in range(len(standards))]

#Calculating the sensitivity function of each standard star
cucumber = 0
for stdspecfile in standards:
    #Read in the observed spectrum of the standard star
    obs_spectra,airmass,exptime,dispersion = st.readspectrum(stdspecfile) #This is an object containing var_farr,farr,sky,sigma,warr
    airstd[cucumber] = airmass
    #plt.clf()
    #plt.plot(obs_spectra.warr,obs_spectra.farr)
    #plt.show()
    #read in the sandard file
    placeholder = cucumber // 2
    stdfile = stdflux[placeholder]
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
    print 'Starting to rebin: ',stdspecfile 
    low = np.rint(np.min(obs_spectra.warr)) #Rounds to nearest integer
    high = np.rint(np.max(obs_spectra.warr))
    size = 0.05 #size in Angstroms you want each bin

    num = (high - low) / size + 1. #number of bins. Must add one to get correct number.
    wavenew = np.linspace(low,high,num=num) #wavelength of each new bin

    #Now do the rebinning using Ian Crossfield's rebinning package
    binflux = st.resamplespec(wavenew,obs_spectra.warr,obs_spectra.farr,200.) #200 is the oversampling factor
    print 'Done rebinning. Now summing the spectrum into new bins to match', stdfile
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
    plt.clf()
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.plot(std_spectra.warr,sens_function)
    cid = fig.canvas.mpl_connect('button_press_event',onclick)
    print 'Please click on both sides of regions you want to exclude. Then close the plot.'
    plt.title('Click both sides of regions you want to exclude. Then close the plot.')
    plt.show(1)


    #Mask our the regions you don't want to fit
    #We need make sure left to right clicking and right to left clicking both work.
    mask = np.ones(len(std_spectra.warr))
    excluded = np.zeros(len(coords))
    lettuce = 0
    if len(coords) > 0:
        while lettuce < len(coords):
            x1 = np.where(std_spectra.warr == (find_nearest(std_spectra.warr,coords[lettuce][0])))
            excluded[lettuce] = np.asarray(x1)
            lettuce += 1
            x2 = np.where(std_spectra.warr == (find_nearest(std_spectra.warr,coords[lettuce][0])))
            if x2 < x1:
                x1,x2 = x2,x1
            mask[x1[0][0]:x2[0][0]+1] = 0 #have to add 1 here to the second index so that we exclude through that index. Most important for when we need to exclude the last point of the array.
            excluded[lettuce-1] = np.asarray(x1)
            excluded[lettuce] = np.asarray(x2)
            lettuce += 1

    excluded =  np.array(excluded).tolist()
    allexcluded[cucumber] = excluded
    indices = np.where(mask !=0.)
    lambdasfit = std_spectra.warr[indices]
    fluxesfit = sens_function[indices]

    #Make sure they are finite
    ind1 = np.isfinite(lambdasfit) & np.isfinite(fluxesfit)
    lambdasfit = lambdasfit[ind1]
    fluxesfit = fluxesfit[ind1]


    print 'Fitting the sensitivity funtion now.'
    order = 4
    repeat = 'yes'
    while repeat == 'yes':
        p = np.polyfit(lambdasfit,fluxesfit,order)
        f = np.poly1d(p)
        smooth_sens = f(lambdasfit)
        residual = fluxesfit - smooth_sens
        plt.close()
        plt.ion()
        #ax1 = plt.subplot(211)
        g, (ax1,ax2) = plt.subplots(2,sharex=True)
        ax1.plot(lambdasfit,fluxesfit,'b+')
        ax1.plot(lambdasfit,smooth_sens,'r',linewidth=2.0)
        ax1.set_ylabel('Sensitivity Function')
        #ax2 = plt.subplot(212,sharex=ax1)
        ax2.plot(lambdasfit,residual,'k+')
        ax2.set_ylabel('Residuals')
        ax1.set_title('Current polynomial order: %s' % order)
        g.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in g.axes[:-1]],visible=False)
        plt.show()
        repeat = raw_input('Do you want to try again (yes/no)? ')
        if repeat == 'yes':
            order = raw_input('New order for polynomial: ')

    orderused[cucumber] = order
    senspolys.append(f)
    cucumber += 1


#Outline for next steps:
#Read in both red and blue files
#compute airmass and compare to airstd
#choose best standard and flux calibrate both blue and red
#save files and write to sensitivity_params.txt

specfile = np.genfromtxt(speclist,dtype=str)
length = len(specfile)
airwd = np.zeros([length])

avocado = 0
while avocado < length:
    #Read in the blue and red spectra we want to flux calibrate. Save the airmass
    WD_spectra1,airmass1,exptime1,dispersion1 = st.readspectrum(specfile[avocado])
    WD_spectra2,airmass2,exptime2,dispersion2 = st.readspectrum(specfile[avocado+1])

    airwd[avocado] = airmass1
    airwd[avocado+1] = airmass2
    #Compare the airmasses to determine the best standard star
    tomato = 0
    while tomato < len(airstd):
        diff = np.absolute(np.mean([airwd[avocado],airwd[avocado+1]]) - np.mean([airstd[tomato],airstd[tomato+1]]))
        if tomato == 0:
            difference = diff
            choice = tomato
        if diff < difference:
            difference = diff
            choice = tomato
        tomato += 2
    
    #To get the flux calibration, perform the following
    #Flux = counts / (Exptime * dispersion * 10**(sens/2.5))
    #Get the sensitivity function at the correct wavelength spacing
    sens_wave1 = senspolys[choice](WD_spectra1.warr)
    sens_wave2 = senspolys[choice+1](WD_spectra2.warr)

    #Perform the flux calibration. We do this on the optimal extraction, non-variance weighted aperture, the sky spectrum, and the sigma spectrum.
    print 'Doing the final flux calibration.'
    star_opflux1 = st.cal_spec(WD_spectra1.opfarr,sens_wave1,exptime1,dispersion1)
    star_flux1 = st.cal_spec(WD_spectra1.farr,sens_wave1,exptime1,dispersion1)
    sky_flux1 = st.cal_spec(WD_spectra1.sky,sens_wave1,exptime1,dispersion1)
    sigma_flux1 = st.cal_spec(WD_spectra1.sigma,sens_wave1,exptime1,dispersion1)

    star_opflux2 = st.cal_spec(WD_spectra2.opfarr,sens_wave2,exptime2,dispersion2)
    star_flux2 = st.cal_spec(WD_spectra2.farr,sens_wave2,exptime2,dispersion2)
    sky_flux2 = st.cal_spec(WD_spectra2.sky,sens_wave2,exptime2,dispersion2)
    sigma_flux2 = st.cal_spec(WD_spectra2.sigma,sens_wave2,exptime2,dispersion2)

    #plt.clf()
    #plt.plot(WD_spectra.warr,star_opflux)
    #plt.show()

    #Save the wavelenghts, counts, and fluxes
    #np.savetxt('python_counts.txt',np.transpose([std_spectra.warr,counts]))
    #np.savetxt('python_sens.txt',np.transpose([obs_spectra.warr,sens_wave]))
    #np.savetxt('python_flux.txt',np.transpose([WD_spectra.warr,flux]))

    print 'Saving the final spectrum.'

    #Save the flux-calibrated spectrum and update the header
    header1 = st.readheader(specfile[avocado])
    header1.set('EX-FLAG',-1) #Extiction correction? 0=yes, -1=no
    header1.set('CA-FLAG',0) #Calibrated to flux scale? 0=yes, -1=no
    header1.set('BUNIT','erg/cm2/s/A') #physical units of the array value
    header1.set('STANDARD',str(standards[choice]),'Flux standard used') #flux standard used for flux-calibration

    header2 = st.readheader(specfile[avocado+1])
    header2.set('EX-FLAG',-1) #Extiction correction? 0=yes, -1=no
    header2.set('CA-FLAG',0) #Calibrated to flux scale? 0=yes, -1=no
    header2.set('BUNIT','erg/cm2/s/A') #physical units of the array value
    header2.set('STANDARD',str(standards[choice+1]),'Flux standard used') #flux standard used for flux-calibration

    #Set up size of new fits image
    Ni = 4. #Number of extensions
    Nx1 = len(star_flux1)
    Nx2 = len(star_flux2)
    Ny = 1. #All 1D spectra

    data1 = np.empty(shape = (Ni,Ny,Nx1))
    data1[0,:,:] = star_opflux1
    data1[1,:,:] = star_flux1
    data1[2,:,:] = sky_flux1
    data1[3,:,:] = sigma_flux1
    
    data2 = np.empty(shape = (Ni,Ny,Nx2))
    data2[0,:,:] = star_opflux2
    data2[1,:,:] = star_flux2
    data2[2,:,:] = sky_flux2
    data2[3,:,:] = sigma_flux2

    #Add '_flux' to the end of the filename
    loc1 = specfile[avocado].find('.ms.fits')
    newname1 = specfile[avocado][0:loc1] + '_flux.ms.fits'
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

    loc2 = specfile[avocado+1].find('.ms.fits')
    newname2 = specfile[avocado+1][0:loc2] + '_flux.ms.fits'
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

    newim2 = pf.PrimaryHDU(data=data2,header=header2)
    newim2.writeto(newname2,clobber=clob)

    #Finally, save all the used parameters into a file for future reference.
    # specfile,current date, stdspecfile,stdfile,order,size,newname
    f = open('sensitivity_params.txt','a')
    now = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M")
    newinfo1 = specfile[avocado] + '\t' + now + '\t' + standards[choice] + '\t' + stdflux[choice//2] + '\t' + str(allexcluded[choice]) + '\t' + str(orderused[choice]) + '\t' + str(size) + '\t' + newname1
    newinfo2 = specfile[avocado+1] + '\t' + now + '\t' + standards[choice+1] + '\t' + stdflux[choice//2] + '\t' + str(allexcluded[choice+1]) + '\t' + str(orderused[choice+1]) + '\t' + str(size) + '\t' + newname2
    f.write(newinfo1 + "\n" + newinfo2 + "\n")
    f.close()


    avocado += 2

print 'Done flux calibrating the spectra.'




