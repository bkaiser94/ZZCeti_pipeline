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

# ===========================================================================
# Lesser Functions Used by Main Functions ===================================
# ===========================================================================

def Read_List( lst ):
    # This function reads a list of images and decomposes them into a python
    # list of image names. 
    list_file = open(lst,'r')
    im_list = list_file.read()
    list_file.close()
    im_list = im_list.split()
    return im_list

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
    # The limits of the trim are: [:, 1:199, 9:2054]
    img_head= pf.getheader(img) 
    img_data= pf.getdata(img)    
    Fix_Header(img_head)
    img_head.append( ('TRIM', '[:, 1:200, 9:2055]' ,'Original Pixel Indices'),
                   useblanks= True, bottom= True )
    NewHdu = pf.PrimaryHDU(data= img_data[:, 1:200, 9:2055], header= img_head)
    NewHdu.writeto('t'+img, output_verify='warn', clobber= True )
    return ('t'+img)
    

# ===========================================================================
# Main Functions ============================================================
# ===========================================================================

def Bias_Subtract( img_list, zero_img ):
    # This function takes in a list of images and a bias image 'zero_img'
    # and performs a pixel by pixel subtration using numpy.
    # The function writes the bias subtracted images as 'bImg_Name.fits'.
    # The output is a list of names for the bias subtrated images. 
    print "\n====================\n"  
    print 'Bias Subtracting Images.\n' 
    zero_data = pf.getdata(zero_img)
    bias_sub_list = []
    for img in img_list:
        hdu = pf.getheader(img)
        Fix_Header(hdu) 
        img_data = pf.getdata(img)
        b_img_data = np.subtract(img_data, zero_data)
        hdu.append( ('BIASSUB', zero_img ,'Image Used to Bias Subtract.'),
                   useblanks= True, bottom= True )
        NewHdu = pf.PrimaryHDU(b_img_data, hdu)
        NewHdu.writeto('b.'+img, output_verify='warn', clobber= True)
        bias_sub_list.append( 'b.'+img )
    return bias_sub_list

# ===========================================================================

def Norm_Flat( flat ):
    # Takes average value of all the pixels and devides the entier flat by 
    # that value using numpy. 
    print "\n====================\n" 
    print 'Normalizing %s By Dividing Each Pixel By Average Value:' % ( flat )
    # Read Data, take average, and divide # 
    flat_data = pf.getdata(flat)
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
    NewHdu.writeto('n'+flat, output_verify='warn', clobber= True )
    return ('n'+flat)
    
# ===========================================================================    
    
def Flat_Field( spec_list, flat ):
    # This Function divides each spectrum in spec_list by the flat and writes
    # The new images as fits files. The output is a list of file names of 
    # the flat fielded images. 
    print "\n====================\n" 
    print 'Flat Fielding %s by Dividing by %s\n' % (spec_list, flat) 
    np.seterr(divide= 'warn')
    flat_data = pf.getdata(flat)
    f_spec_list = []
    for spec in spec_list:
        spec_data = pf.getdata(spec)
        f_spec_data = np.divide(spec_data, flat_data)
        hdu = pf.getheader(spec)
        Fix_Header(hdu)
        hdu.append( ('FLTFIELD', flat,'Image used to Flat Field.'), 
               useblanks= True, bottom= True )    
        NewHdu = pf.PrimaryHDU(data= f_spec_data, header= hdu)
        NewHdu.writeto('f'+spec, output_verify='warn', clobber= True)
        f_spec_list.append('f'+spec)
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
    lst_st = decimal_dec( hdu.header['LST'] ) # start exposure LST in hours
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
    # print images being combined
    for img in im_list:
        print "%s" % img
    
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
    
    # Make sure header BITPIX reflects data encodeing as float 32 ie: -32 
    hdu.header['BITPIX'] = -32
    
    # Write header to new fits file  
    hdu.writeto(output_name, output_verify='warn', clobber= True)
    
    # write combined data to new fits file  # 
    pf.update(output_name, data= comb_img, header= hdu.header, 
                output_verify='warn')
                        
    print '\nDone!' 
    print 'File', output_name, 'was succesfuly writen.' 
    return output_name             

# ===========================================================================
# ===========================================================================            
# ===========================================================================