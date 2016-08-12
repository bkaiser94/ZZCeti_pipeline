# -*- coding: utf-8 -*-

# This is my first attempt at combining images. 
# it works on laptop but not in infierno or cielo. 

# ===========================================================================
# Packages ==================================================================
# ===========================================================================

import numpy as np
import pyfits as pf
import ReduceSpec_tools as rt
import warnings

# ===========================================================================
# Code to Reduce Spectrum ===================================================
# ===========================================================================
# This peice takes in arguments from the command line to 
# All the functions required are called from ReduceSpec_tools.py
            
if __name__ == "__main__":
    from sys import argv
    args = argv # arguments from comand line #
    nargs = len(args) # number of arguments # 
    
    if (nargs < 5):
        print "\n====================\n"
        print "\nNot Enough Inputs." 
        print "Need at least 4 inputs: listZero, listFlat, listSpec, listFe"
        print "Optional inputs: overwrite= , low_sig= , high_sig=  "
        print "Example:"
        print "\n>>> python imcombine.py listZero listFlat listSpec listFe \n"
        print "\n====================\n"
    
    # Unpack list from command line and combe trough them for diffrent observations # 
    scriptname = args[0]
    zero_lists = rt.List_Combe( rt.Read_List( args[1] ) )
    flat_lists = rt.List_Combe( rt.Read_List( args[2] ) )
    spec_lists = rt.List_Combe( rt.Read_List( args[3] ) )
    fe_lists = rt.List_Combe( rt.Read_List( args[4] ) )
     
    # Select names from the first image of each observation # 
    zero_names= []
    for zero in zero_lists:
        zero_names.append(zero[0][5:])
    flat_names= []
    for flat in flat_lists:
        flat_names.append(flat[0][5:])
    spec_names= []
    for spec in spec_lists:
        spec_names.append(spec[0][5:])
    fe_names = []
    for lamp in fe_lists:
        fe_names.append(lamp[0][5:])

    # Default values for special commands if none are given these dont change #   
    overwrite = False # dont give imcombine permision to overwrite files # 
    lo_sig = 10
    hi_sig = 3
    method = 'median' # method used to combine images 
    
    # If overwrite special comand is given # 
    if nargs >= 6:
        overwrite = args[5]
        warnings.filterwarnings('ignore', category=UserWarning, append=True)
    # If low_sigma and high_sigma values are given # 
    if nargs >= 8: 
        lo_sig = float(args[6])
        hi_sig = float(args[7]) 
    # If method is given #  
    if nargs >= 9:
        method = args[8]
        
    #Set up array to save for diagnostics. This is defined in rt.init()
    rt.init()
    
    # The rest of the code runs the reduction procces up to apall #  =========
    # Combine Zeros # 
    comb_zero = rt.imcombine(zero_lists[0], zero_names[0], 'average', lo_sig= lo_sig, 
                        hi_sig= hi_sig, overwrite= overwrite)
    
    # Bias Subtract Flats # 
    nf= len(flat_lists) # number of flats
    b_flat_lists= []
    i= 0
    while i < nf:
        b_flat_lists.append( rt.Bias_Subtract(flat_lists[i], comb_zero ) )
        i= i+1
    
    # Combine Bias Subtracted Flats # 
    i= 0
    comb_flat= []
    while i < nf:
        comb_flat.append( rt.imcombine(b_flat_lists[i], 'b.'+flat_names[i], 'median', 
                        lo_sig= lo_sig, hi_sig= hi_sig, overwrite= overwrite) )
        i= i+1
    
    # Normalize Flat (divide by average of counts) # 
    i= 0
    nb_flat= []
    while i < nf:
        nb_flat.append( rt.Norm_Flat_Poly(comb_flat[i]) )
        i= i+1
                        
    # Bias Subtract Spec # 
    i= 0
    b_spec_list= []
    nsp= len(spec_lists); # number of spectra
    while i < nsp:
        b_spec_list.append( rt.Bias_Subtract(spec_lists[i], comb_zero) )
        i= i+1
    
    # Flat Field Individual Spectra #
    blueindex = [i for i, s in enumerate(nb_flat) if 'blue' in s.lower()]
    nbflatblue = nb_flat[blueindex[0]]
    redindex = [i for i, s in enumerate(nb_flat) if 'red' in s.lower()]
    if len(redindex) > 0:
        nbflatred = nb_flat[redindex[0]]
    i= 0
    fb_spec_list = []
    while i < nsp:
        if b_spec_list[i][0].__contains__('blue') == True:
            fb_spec_list.append( rt.Flat_Field(b_spec_list[i], nbflatblue) )
        elif b_spec_list[i][0].__contains__('red') == True:
            fb_spec_list.append( rt.Flat_Field(b_spec_list[i], nbflatred) )
        else: 
            print ("Problem applying the Flats." )
            print ("Could not identify blue or red setup.")
        i= i+1
        
    # Combine Spectra # 
    i= 0 
    comb_fb_spec = []
    while i < nsp:
        rt.checkspec(fb_spec_list[i])
        comb_fb_spec.append ( rt.imcombine(fb_spec_list[i], 'fb.'+spec_names[i], 'average', 
                        lo_sig= lo_sig, hi_sig= hi_sig, overwrite= overwrite) )
        i= i+1
                        
    # Trim Spectra # 
    i= 0
    while i < nsp:
        rt.Trim_Spec(comb_fb_spec[i]); 
        i= i+1
                        
    print "\n====================\n"

    #########################################
    # Combine Fe lamps # 
    print "Combining and trimming Fe lamps."
    nf = len(fe_lists) #number of fe lamps
    i = 0
    comb_lamp = []
    while i < nf:
        comb_lamp.append( rt.imcombine(fe_lists[i], fe_names[i], 'average', lo_sig= lo_sig, 
                        hi_sig= hi_sig, overwrite= overwrite) )
        i = i+1

    # Trim lamps # 
    i= 0
    while i < nf:
        rt.Trim_Spec(comb_lamp[i]); 
        i= i+1
 
    ########################################

    print "Done. Ready for Apeture Extraction.\n"
    
# ===========================================================================
# ===========================================================================
# ===========================================================================
