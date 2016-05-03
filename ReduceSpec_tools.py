# -*- coding: utf-8 -*-
"""
Created on Sun Aug 23 20:48:10 2015

@author: jmeza
"""


# ===========================================================================
# Packages ==================================================================
# ===========================================================================

import numpy as np
import pyfits as pf
import mpfit
import os

# ===========================================================================
# Lesser Functions Used by Main Functions ===================================
# ===========================================================================

def gauss(x,p): #single gaussian
    return p[0] +  p[1]*np.exp(-(((x-p[2])/(np.sqrt(2)*p[3])))**2.)

def fitgauss(p,fjac=None,x=None,y=None,err=None):
    #Parameter values are passed in p
    #fjac = None just means partial derivatives will not be computed
    model = gauss(x,p)
    status = 0
    return([status,(y-model)/err])

def checkspec(listcheck):
    #Calculates the FWHM and profile postion for two points on each spectrum
    #If these values deviate by more than given values, prints warning.
    #Saves all values in a text file.
    listcheck = np.genfromtxt(listcheck,dtype=str)
    print 'Now checking FWHM and center of spectral profile for stability.'
    #Max values acceptable
    maxcendev = 2. #Deviation from center of gaussian
    maxfwhmdev = 0.5 #deviation of fwhm

    fwhm1 = np.zeros(len(listcheck))
    fwhm2 = np.zeros(len(listcheck))
    center1 = np.zeros(len(listcheck))
    center2 = np.zeros(len(listcheck))

    f = open('FWHM_records.txt','a')
    n = 0.
    for specfile in listcheck:
        datalist = pf.open(specfile)
        data = datalist[0].data
        data = data[0,:,:]
        data = np.transpose(data)

        #Fit a column of the 2D image to determine the center and FWHM 
        forfit1 = data[550,2:] #column 550 and 1750 are good for both setups
        guess1 = np.zeros(4)
        guess1[0] = np.mean(forfit1)
        guess1[1] = np.amax(forfit1)
        guess1[2] = np.argmax(forfit1)
        guess1[3] = 3.
        error_fit1 = np.ones(len(forfit1))
        xes1 = np.linspace(0,len(forfit1)-1,num=len(forfit1))
        fa1 = {'x':xes1,'y':forfit1,'err':error_fit1}
        fitparams1 = mpfit.mpfit(fitgauss,guess1,functkw=fa1,quiet=True)
        fwhm1[n] = 2.*np.sqrt(2.*np.log(2.))*fitparams1.params[3]
        center1[n] = fitparams1.params[2]
        #print np.round(fwhm1[n],decimals=1),np.round(center1[n],decimals=1)
        #plt.clf()
        #plt.plot(forfit1)
        #plt.plot(xes1,gauss(xes1,fitparams1.params))
        #plt.show()

        forfit2 = data[1750,2:] #column 550 and 1750 are good for both setups
        guess2 = np.zeros(4)
        guess2[0] = np.mean(forfit2)
        guess2[1] = np.amax(forfit2)
        guess2[2] = np.argmax(forfit2)
        guess2[3] = 3.

        error_fit2 = np.ones(len(forfit2))
        xes2 = np.linspace(0,len(forfit2)-1,num=len(forfit2))
        fa2 = {'x':xes2,'y':forfit2,'err':error_fit2}
        fitparams2 = mpfit.mpfit(fitgauss,guess2,functkw=fa2,quiet=True)

        fwhm2[n] = 2.*np.sqrt(2.*np.log(2.))*fitparams2.params[3]
        center2[n] = fitparams2.params[2]
        #print np.round(fwhm2[n],decimals=1),np.round(center2[n],decimals=1)
        #plt.clf()
        #plt.plot(forfit2)
        #plt.plot(xes2,gauss(xes2,fitparams2.params))
        #plt.show()

        info = specfile + '\t' + '550' + '\t' + str(np.round(fwhm1[n],decimals=2)) + '\t' + str(np.round(center1[n],decimals=2)) + '\t' + '1750' + '\t' + str(np.round(fwhm2[n],decimals=2)) + '\t' + str(np.round(center2[n],decimals=2))
        f.write(info+ "\n")

        n += 1
    f.close()

    #Check if values deviate by more than a certain amount

    if (np.max(fwhm1) - np.min(fwhm1)) > maxfwhmdev:
        print 'WARNING!!! Left FWHM varying significantly. Values are %s' % fwhm1
    elif (np.max(fwhm2) - np.min(fwhm2)) > maxfwhmdev:
        print 'WARNING!!! Right FWHM varying significantly. Values are %s' % fwhm2
    else:
        print 'FWHM is stable.'

    if (np.max(center1) - np.min(center1)) > maxcendev:
        print 'WARNING!!! Left profile center varying significantly. Values are %s' % center1
    elif (np.max(center2) - np.min(center2)) > maxcendev:
        print 'WARNING!!! Right profile center varying significantly. Values are %s' % center2
    else:
        print 'Profile center is stable.'
        
# ============================================================================

def Read_List( lst ):
    # This function reads a list of images and decomposes them into a python
    # list of image names. 
    list_file = open(lst,'r')
    im_list = list_file.read()
    list_file.close()
    im_list = im_list.split()
    return im_list
    
def List_Combe(img_list):
    # This is meant to combe trough list names to identify seperate sublist of
    # stars / flats / standars 
    sub_lists= [] # list of sub_list of images 
    sl= [] # sub_list of images
    sl.append(img_list[0]) # place first image in sublist
    i= 0; # image counter  
    while i < len(img_list)-1: # run trough all images 
        if img_list[i+1].__contains__(img_list[i][4:]) == True:
            sl.append(img_list[i+1]) # place it in the sub_list 
        else:
            # if the images dont match: 
            sub_lists.append(sl) # write the sublist to the list of sublist 
            sl= [] # clear the sublist
            sl.append(img_list[i+1]) # append the image to the new list 
        i= i+1 # image counter
    sub_lists.append(sl) # append the last sublist to the list of sublist 
    return sub_lists # return the list of sub_list of images
    
def check_file_exist(name):
    # This function is to be called before wirting a file. 
    # This function checks if the file name already exist.
    # If it does it appends a number to the begining until 
    # the name no longer matches the files in the directory. 
    
    # List of files in directory
    listDirFiles = [f for f in os.listdir('.') if f.endswith('.fits')]
    # If "name" is in the derectory append a number i until it doent match 
    # If name is not in directory then we simply return name
    if listDirFiles.__contains__(name):
        i= 2
        while listDirFiles.__contains__(name):
            name= str(i) + name
            i= i+1
    return name

def Fix_Header( header ):
    # This function deletes the header cards that contain the badly coded 
    # degree symbol '\xb0'. If they are not deleted pyfits won't write the
    # headers. 
    bad_key = ['param0', 'param61', 'param62', 'param63']
    for p in bad_key:
        if p in header: 
            bad_str = header.comments[p]
            if '\xb0' in bad_str:
                del header[p]      

def decimal_dec(hdu_str):
    # Read header strings in "hh:mm:ss" or "dd:mm:ss" fromat 
    # and outputs the value as a decimal. 
    val_list = [float(n) for n in hdu_str.split(':')]
    if val_list[0] < 0 :
        sng = -1
        val_list[0] = sng*val_list[0]
    else:
        sng = 1
    val_deci =  sng*(val_list[0]+((val_list[1]+(val_list[2]/60.0))/60.0))
    return val_deci

def decimal_ra(hdu_str):
    # Read header strings in "hh:mm:ss" or "dd:mm:ss" fromat 
    # and outputs the value as a decimal. 
    val_list = [float(n) for n in hdu_str.split(':')]
    if val_list[0] < 0 :
        sng = -1
        val_list[0] = sng*val_list[0]
    else:
        sng = 1
    val_deci =  15*sng*(val_list[0]+((val_list[1]+(val_list[2]/60.0))/60.0))

    return val_deci

def SigClip(data_set, lo_sig, hi_sig):
    # Sigma Cliping Function  #
    # Input is set of counts for a particular pixel, 
    # along with low and high sigma factors. 
    # Output is a list containg only the data that is with the sigma factors.
    # Only a single rejection iteration is made. 
    Avg = np.mean(data_set)
    St_Div = np.std(data_set)
    min_val = Avg-lo_sig*St_Div
    max_val = Avg+hi_sig*St_Div
    cliped_data = []
    for val in data_set:
        if min_val <= val <= max_val:
            cliped_data.append( val )
    return cliped_data  
        
def RaDec2AltAz(ra, dec, lat, lst ):
    # Input: RA in decimal hours; DEC in decimal deg; 
    # LAT in decimal deg; LST in decimal hours; 
    # Output: ALT, AZ, HA in decimal deg. 
     
    # Compute Hour Angle
        ha = lst-ra # hour angle in deg
        if ha < 0 :
            ha = ha+360
        if ha > 360:
            ha = ha-360
    # Convert Qunataties to Radians 
        ra = ra*(np.pi/180.0) 
        dec = dec*(np.pi/180.0) 
        lat = lat*(np.pi/180.0) 
        ha = ha*(np.pi/180.0)
    # Calculate Altitiude 
        a =  np.sin(dec)*np.sin(lat)
        b = np.cos(dec)*np.cos(lat)*np.cos(ha)
        alt = np.arcsin( a+b ) # altitude in radians 
    # Calculate Azimuth 
        a = np.sin(dec)-np.sin(lat)*np.sin(alt)
        b = np.cos(lat)*np.cos(alt)
        az = np.arccos( a/b ) # azumuth in radians
        if np.sin(ha) > 0:
            az = (2*np.pi) - az 
    # Convert Alt, Az, and Ha to decimal deg
        alt = alt*(180.0/np.pi)
        az = az*(180.0/np.pi)
        ha = ha*(180.0/np.pi)
        return alt, az, ha
        
def AirMass(alt, scale):
    # Calculates instantaneus airmass to be called by SetAirMass() #
    # Input: 
    #   scale = atmospheric scale factor (defalut 750)
    #   alt = altitude of star in degrees.  
    # Output: 
    #   AM = airmass from given altitude and scale factor
    x = scale*np.sin(np.pi*alt/180.)
    AM = np.sqrt( x**2 + 2*scale + 1 ) - x
    return AM
        
def EffectiveAirMass(AM_st, AM_mid, AM_end):
    # Calculate effective airmass to be called by SetAirMass() and Imcombine() 
    # Input: airmass at start, middel, and end of an exposure. 
    # Output: Effective Airmass 
    AM_eff = (AM_st + 4*AM_mid + AM_end)/6.  
    return AM_eff

def Trim_Spec(img):
    # Trims Overscan region and final row of of image #
    # The limits of the 2x2 binned trim are: [:, 1:199, 9:2054]
    # The limits of the 1x2 trim are: [:, 1:199, 19:4111]
    print "\n====================\n"  
    print 'Triming Image: %s\n' % img
    img_head= pf.getheader(img) 
    img_data= pf.getdata(img)    
    Fix_Header(img_head)
    length = float(img_head['PARAM17'])
    if length == 2071.:
        img_head.append( ('CCDSEC', '[9:2055,1:200]' ,'Original Pixel Indices'),
                   useblanks= True, bottom= True )
        NewHdu = pf.PrimaryHDU(data= img_data[:, 1:200, 9:2055], header= img_head)
        new_file_name= check_file_exist('t'+img)
        NewHdu.writeto(new_file_name, output_verify='warn', clobber= True )
        return (new_file_name)
    elif length == 4142.:
        img_head.append( ('CCDSEC', '[19:4111,1:200]' ,'Original Pixel Indices'),
                   useblanks= True, bottom= True )
        NewHdu = pf.PrimaryHDU(data= img_data[:, 1:200, 19:4111], header= img_head)
        new_file_name= check_file_exist('t'+img)
        NewHdu.writeto(new_file_name, output_verify='warn', clobber= True )
        return (new_file_name)
    else:
        print 'WARNING. Image not trimmed. \n'

def Add_Scale (img_block):
    # Function to be called by Imcombine. 
    # The function is meant to additively sclae a set of images, (zeros in particular). 
    # The input is a numpy block of pixel values (see imcombine). 
    # The function calculates the average number of 
    # counts of the region [25:75, 1700:1800] of the first image. 
    # Then scales the rest of the images by adding the diffrence between the 
    # average counts of the first image and its own.
    # Returns a scaled image block, and a list of scale values. 
    print("Scaling Counts Additively.\n")
    ni, ny, nx = np.shape(img_block)
    Cavg= [] # Average Counts 
    Sval= []  # Scale Values
    for i in range(0,ni):
        Cavg.append( np.mean(img_block[i, 25:75, 1700:1800]) )
        Sval.append( Cavg[0]-Cavg[i] )
        img_block[i]= img_block[i] + Sval[i]     
    return img_block, Sval
    
def Mult_Scale (img_block):
    # Function to be called by Imcombine. 
    # The function is meant to multiplicative sclae a set of images, (flats in particular). 
    # The input is a numpy block of pixel values (see imcombine). 
    # The function calculates the average number of 
    # counts of the region [25:75, 1700:1800] of the first image. 
    # Then scales the rest of the images by multiplying by the ratio between the 
    # average counts of the first image and its own.
    # Returns a scaled image block, and a list of scale values. 
    print("Scaling Counts Multiplicatively.\n")
    ni, ny, nx = np.shape(img_block)
    Cavg= [] # Average Counts 
    Sval= []  # Scale Values
    for i in range(0,ni):
        Cavg.append( np.mean(img_block[i, 25:75, 1700:1800]) )
        Sval.append( Cavg[0]/Cavg[i] )
        img_block[i]= img_block[i]*Sval[i]     
    return img_block, Sval
        
    
# ===========================================================================
# Main Functions ============================================================
# ===========================================================================

def Bias_Subtract( img_list, zero_img ):
    # This function takes in a list of images and a bias image 'zero_img'
    # and performs a pixel by pixel subtration using numpy.
    # The function writes the bias subtracted images as 'b.Img_Name.fits'.
    # The output is a list of names for the bias subtrated images. 
    print "\n====================\n"  
    print 'Bias Subtracting Images: \n' 
        
    zero_data = pf.getdata(zero_img)
    bias_sub_list = []
    for img in img_list:
        hdu = pf.getheader(img)
        Fix_Header(hdu) 
        img_data = pf.getdata(img)
        img_data[ np.isnan(img_data) ] = 0
        b_img_data = np.subtract(img_data, zero_data)
        print 'b.'+"%s Mean: %.3f StDiv: %.3f" % (img, np.mean(b_img_data), np.std(img_data))
        hdu.append( ('BIASSUB', zero_img ,'Image Used to Bias Subtract.'),
                   useblanks= True, bottom= True )
        NewHdu = pf.PrimaryHDU(b_img_data, hdu)
        bias_sub_name= check_file_exist('b.'+img)
        NewHdu.writeto(bias_sub_name, output_verify='warn', clobber= True)
        bias_sub_list.append( bias_sub_name )
    return bias_sub_list

# ===========================================================================

def Norm_Flat_Avg( flat ):
    # Takes average value of all the pixels and devides the entier flat by 
    # that value using numpy. 
    print "\n====================\n" 
    print 'Normalizing %s By Dividing Each Pixel By Average Value:' % ( flat )
    # Read Data, take average, and divide # 
    flat_data = pf.getdata(flat)
    flat_data[ np.isnan(flat_data) ] = 0
    # Calculate Average of the flat excluding bottom row and overscan regions # 
    avg_flat = np.average( flat_data[:, 1:200, 9:2055] )
    norm_flat_data = np.divide( flat_data, float(avg_flat) )
    print 'Average Value: %s\n' % avg_flat
    # Copy Header, write changes, and write file #
    hdu = pf.getheader(flat)
    Fix_Header(hdu)
    hdu.append( ('NORMFLAT', avg_flat,'Average Used to Normalize the Flat.'), 
               useblanks= True, bottom= True )
    NewHdu = pf.PrimaryHDU(data= norm_flat_data, header= hdu)
    norm_flat_name= check_file_exist('n'+flat)
    NewHdu.writeto(norm_flat_name, output_verify='warn', clobber= True )
    
    print 'Flat: %s Mean: %.3f StDiv: %.3f' % (norm_flat_name, np.mean(norm_flat_data), np.std(norm_flat_data)) 
    return (norm_flat_name)

# ============================================================================    

def Norm_Flat_Poly( flat ):
    print "\n====================\n" 
    print 'Normalizing %s By Fitting Polynomial to center rows [95:105]:' % ( flat )
    # Decide Order # 
    if flat.__contains__("blue")== True:
        order= 4;
    elif flat.__contains__("red")== True:
        order= 4; 
    else:
        print ("Could not identifiy blue or red flat")
        order= raw_input("Fit Order?>>>")
    print "Fit Order: %s" % order
    # Read Flat and Average Center Rows # 
    flat_data = pf.getdata(flat)
    flat_data[ np.isnan(flat_data) ] = 0
    fit_data= np.median(flat_data[0][95:105], axis=0) # Median of center Rows
    X= range(0,len(fit_data)) # Column Numbers 
    # Fit the data removeing the limits of the overscan regions. #
    lo= 10;
    hi= 2055;
    # Calculate Fit # 
    coeff= np.polyfit(X[lo:hi], fit_data[lo:hi], order ) # coefficents of polynomial fit # 
    profile= np.poly1d(coeff)(X) # Profile Along Dispersion axis # 
    # Divide each Row by the Profile # 
    for row in flat_data[0]:
        i= 0; 
        while i < len(row): 
            row[i]= row[i]/profile[i]
            i= i+1   
            
    # Copy Header, write changes, and write file #
    hdu = pf.getheader(flat)
    Fix_Header(hdu)
    hdu.append( ('NORMFLAT ', order,'Flat Polynomial Fit Order'), 
               useblanks= True, bottom= True )
    for i in range(0,len(coeff)):
        coeff_str=  "{0:.5e}".format(coeff[i])
        coeff_order = str(len(coeff)-i-1)
        coeff_title = 'NCOEF%s' %coeff_order
        coeff_expla = 'Flat Polynomial Coefficient - Term %s' %coeff_order
        hdu.append((coeff_title,coeff_str,coeff_expla),
                   useblanks= True, bottom= True )
    NewHdu = pf.PrimaryHDU(data= flat_data, header= hdu)
    norm_flat_name= check_file_exist('n'+flat)
    NewHdu.writeto(norm_flat_name, output_verify='warn', clobber= True )
    
    print '\nFlat: %s Mean: %.3f StDiv: %.3f' % (norm_flat_name, np.mean(flat_data), np.std(flat_data))
    return (norm_flat_name)
    
# ===========================================================================    
    
def Flat_Field( spec_list, flat ):
    # This Function divides each spectrum in spec_list by the flat and writes
    # The new images as fits files. The output is a list of file names of 
    # the flat fielded images. 
    print "\n====================\n" 
    print 'Flat Fielding Images by Dividing by %s\n' % (flat) 
        
    np.seterr(divide= 'warn')
    flat_data = pf.getdata(flat)
    f_spec_list = []
    for spec in spec_list:
        spec_data = pf.getdata(spec)
        f_spec_data = np.divide(spec_data, flat_data)
        f_spec_data[ np.isnan(f_spec_data) ] = 0
        print "f"+"%s Mean: %.3f StDiv: %.3f" % (spec, np.mean(f_spec_data), np.std(f_spec_data) ) 
        hdu = pf.getheader(spec)
        Fix_Header(hdu)
        hdu.append( ('FLATFLD', flat,'Image used to Flat Field.'), 
               useblanks= True, bottom= True )    
        NewHdu = pf.PrimaryHDU(data= f_spec_data, header= hdu)
        new_file_name= check_file_exist('f'+spec)
        NewHdu.writeto(new_file_name, output_verify='warn', clobber= True)
        f_spec_list.append(new_file_name)
    return f_spec_list

# ===========================================================================

def SetAirMass(img, lat= -30.238, scale= 750):
    # This Function Calculates The Effective Airmass of a single image  
    # Inputs:
    #   img = image name
    #   lat = latitude of observer in decimal degrees. 
    #       (Default Soar lat: '-30:14:16.8' = -30.238 deg)
    #   scale = atmospheric scale factor 750
    # Output: 
    #   AMeff = effective airmass for single image
     
    # Image Info #    
    hdulist = pf.open(img, 'update')
    hdu = hdulist[0]
    
    Fix_Header(hdu.header)
            
    ra = decimal_ra( hdu.header['RA'] ) # hours
    dec = decimal_dec( hdu.header['DEC'] ) # deg
    lst_st = decimal_ra( hdu.header['LST'] ) # start exposure LST in hours
    exp = hdu.header['EXPTIME']  # sec
    lst_mid = lst_st + (exp/2.)/3600. # mid exposure LST in hours
    lst_end = lst_st + (exp)/3600. # end exposure LST in hours

    # Air Mass Calculations # 
    times = [lst_st, lst_mid, lst_end]
    AM = []
    for t in times:
        alt, az, ha = RaDec2AltAz(ra, dec, lat, t )
        airmass = AirMass(alt, scale)
        AM.append( airmass )
    AMeff = EffectiveAirMass(AM[0], AM[1], AM[2])
    
    # Print and write to header # 
    print '\nImage:', img
    print 'Observatory Latitude: %s' % lat
    print 'AM_st   AM_mid  AM_end  AM_eff'
    print '%5.4f  %5.4f  %5.4f  %5.4f' % (AM[0], AM[1], AM[2], AMeff)
    hdu.header.set( 'AIRMASS', np.round(AMeff,6) , 
                   'Calculated Effective Airmass' )
    hdulist.close()
    return AMeff    
 
# =========================================================================== 
 
def imcombine(im_list, output_name, method,  
              lo_sig = 10, hi_sig = 3, overwrite= False):
# Image Combination Script # 
# Inputs:
#   im_list = mist be a python list of images or "@listfile"
#   output_name =  name of combined fits image 
#   method = The method to use for combining (median, average, sum)
#   lo_sig = low sigma cliping factor (default = 3 sigma) 
#   hi_sig = high sigma cliping factor (default = 3 sigma)
#   overwrite = if true go ahead and re write existing file 'output_name'
#               if false it will warn you and ask for new output_name. 
#               (default false)
# Output:
#   After succefully combining, calculateing airmass, and writing to fits file, 
#   The return of this function is the name of the combined 
#   image (Output_name).
              
    print "\n====================\n" 
    print "Combining Images:"
    print "Using %s of count values." % method 
    print "Sigma Cliping Factors (low, high): (%s, %s)\n" % (lo_sig, hi_sig)
    
    # Read image data and put it in a numpy block # 
    Ni = len(im_list)
    for i in range(0, Ni):
        # First size the array to contain the data based on 1st image #
        # Create block with 3 axis:
        #   axis[0] has length of number of images.
        #   axis[1] is the vertical axis of the chip.
        #   axis[2] is the horizontal axis of the chip.
        if i == 0:  
            img_data = pf.getdata(im_list[i])
            n,Ny,Nx = np.shape(img_data)
            img_block = np.ndarray( shape= (Ni,Ny,Nx) )
            img_block[i,:,:] = img_data
        # Then go ahead and read the rest of the images into the block #   
        else: 
            img_block[i,:,:] = pf.getdata(im_list[i])
        # set nan values to zero # 
        img_block[ np.isnan(img_block) ] = 0
        
    # If Zero Additive Scale Images # 
    if im_list[0].lower().__contains__("zero"):
        img_block, Scale= Add_Scale(img_block)
    # If Flats Multiplicative Scale Images # 
    elif im_list[0].lower().__contains__("flat"):
        img_block, Scale= Mult_Scale(img_block)
    # If Not, Dont Scale # 
    else: 
        print "Did Not Scale Images.\n" 
        Scale= np.empty(Ni)
        Scale[:]= np.NaN
    
    # Print Name and Statistics of Each image % 
    for i in range(0,Ni):
        Avg= np.mean(img_block[i,25:75,1700:1800])
        Std= np.std(img_block[i,25:75,1700:1800])
        print ( "%02d: %s ScaleValue:% .3f Mean: %.3f StDiv: %.3f" 
                % (i, im_list[i], Scale[i], Avg, Std) )
    
    ## Combine the images acording to input "method" using SigmaClip() above ## 
    comb_img = np.ndarray( shape= (1,Ny,Nx), dtype='float32')
    while True: # Contunualy askes for method if input is wierd # 
        
        if method == 'median':
            for y in range(0,Ny):
                for x in range(0,Nx):
                    counts = img_block[:,y,x]
                    val = np.median( SigClip(counts, lo_sig, hi_sig) )
                    comb_img[0,y,x] = np.float32(val)
            break # exit while loop 
    
        elif method == 'average':
            for y in range(0,Ny):
                for x in range(0,Nx):
                    counts = img_block[:,y,x]
                    val = np.average( SigClip(counts, lo_sig, hi_sig) )
                    comb_img[0,y,x] = np.float32(val)
            break # exit while loop
        
        elif method == 'sum':
            for y in range(0,Ny):
                for x in range(0,Nx):
                    counts = img_block[:,y,x]
                    val = np.sum( SigClip(counts, lo_sig, hi_sig) )
                    comb_img[0,y,x] = np.float32(val)
            #print img_block[:,100,50]
            #print comb_img[:,100,50]
            break # exit while loop
        
        else:
            # if 'method' input is wanky, ask for method again. 
            print "\nError: Method NOT AVALABLE." 
            print "Available Methods: ('median', 'average', 'sum')"
            print "Enter Valid Method"
            method = raw_input('>>>')
    
    # Set NAN values to zero 
    comb_img[ np.isnan(comb_img) ] = np.float32(0)
    
    ###### Calculate Effetive Airmass for combined image ######
    # The EffAM value is writen into the header in the next section #
    print '\nCalculating Effective Airmass:'    
    
    # if were just combining 2 images #
    if Ni == 2:
        AM0 = SetAirMass(im_list[0])
        AM2 = SetAirMass(im_list[1])
        AM1 = (AM0+AM2)/2
        EffAM = EffectiveAirMass(AM0, AM1, AM2)
        print '\nEffective Airmass of combined image: %5.4f' % EffAM
    # if were combining an odd number of images # 
    elif Ni%2 == 1: 
        images = [ im_list[0], im_list[Ni//2], im_list[-1] ] 
        AM = [ SetAirMass(img) for img in images ] 
        EffAM = EffectiveAirMass( AM[0], AM[1], AM[2] )
        print '\nEffective Airmass of combined image: %5.4f' % EffAM
    # if were combing an even number of images #  
    elif Ni%2 == 0:
        images = [im_list[0], im_list[(Ni//2)-1], im_list[Ni//2], im_list[-1]]
        AM = [ SetAirMass(img) for img in images ]
        EffAM = EffectiveAirMass( AM[0], (AM[1]+AM[2])/2, AM[3])
        print '\nEffective Airmass of combined image: %5.4f' % (EffAM)
    # Otherwise we fail # 
    else:
        print "Eff AirMass calculation failed? This never happens!"
    
    ###### Overwrite Protection loop, just in case ######
    if overwrite == False:
        from os.path import isfile
    
    while overwrite == False: # Outer Loop #  
    # Breaks if file name doesnot exist or overwrite == true # 
        exist = isfile(output_name) # Asks computer if file name exist # 
        
        if exist == False: 
            print "\nWriting combined image to fits file",output_name,"..." 
            break # Break out of outter loop and continue writing # 
        elif exist == True:
            while True: # Inner Loop # 
            # Breaks if user wishes to overwite, abort, or gives new name.
            # loop also checks new names for existance.  
                print "\nFile name",output_name,
                print "already exist do you wish to overwrite?"
                yes_no = raw_input('yes or no ?>>>')
            
                if yes_no == 'no':
                    # If overwrite no: prompt new name or abort # 
                    print"\nEnter new file name or Ctrl-c to Abort "
                    output_name = raw_input('>>>')
                    print "\nNew File Name:= ", output_name
                    break # breaks out of Inner loop only.
                          # Code proceeds to Outer Loop to 
                          # ask computer if new file name exist.     
                
                elif yes_no == 'yes':
                    # If overwrite yes: Break Inner Loop and Outer loop # 
                    overwrite = True
                    print "\nOverwriting Image:", output_name
                    break
                
                else: 
                    # If yes_no input is wierd return to Inner loop
                    # to ask question again. 
                    print "\nInput Not Recognized."
                    
    ###### The following part only runs if above while loop is satisfied ######
    
    # Copy header of first image in im_list and fix degree symbol issue. 
    hdulist = pf.open(im_list[0])
    hdu = hdulist[0]
    # This checks the string and deletes the bad keywords from header. 
    Fix_Header(hdu.header)
    
    # Write Effective Airmass into header # 
    hdu.header.set('AIRMASS',np.round(EffAM,6),'Calculated Effective Airmass')
    
    # Write the imcombine information into header #
    N = len(im_list)
    for i in range(0,N):
        num = str(i+1).zfill(3)
        key = 'IMCMB'+num
        hdu.header.append( (key, im_list[i]), useblanks= True, bottom= True )
    hdu.header.append( ('NCOMBINE', N), useblanks= True, bottom = True )
    hdu.header.append( ('COMBTYPE', method,'Operation Used to Combine'),
                      useblanks= True, bottom= True )
    
    # Make sure header BITPIX reflects data encodeing as float 32 ie: -32 
    hdu.header['BITPIX'] = -32
    
    # Write header to new fits file  
    new_file_name= check_file_exist(output_name)
    hdu.writeto(new_file_name, output_verify='warn', clobber= True)
    
    # write combined data to new fits file  # 
    pf.update(output_name, data= comb_img, header= hdu.header, 
                output_verify='warn')
                         
    print ( "\nCombined Image: %s Mean: %.3f StDiv: %.3f" 
            % (new_file_name, np.mean(comb_img), np.std(comb_img)) ) 
    return new_file_name           

# ===========================================================================
# ===========================================================================            
# ===========================================================================
