# -*- coding: utf-8 -*-

# This is my first attempt at combining images. 
# it works on laptop but not in infierno or cielo. 

# ===========================================================================
# Packages ==================================================================
# ===========================================================================

import numpy as np
import pyfits as pf
import ReduceSpec_tools as rt

# ===========================================================================
# Code to Reduce Spectrum ===================================================
# ===========================================================================
# This peice takes in arguments from the command line to 
# All the functions required are called from ReduceSpec_tools.py
            
if __name__ == "__main__":
    from sys import argv
    args = argv # arguments from comand line #
    print args
    nargs = len(args) # number of arguments # 
    
    if (nargs < 4):
        print "\n====================\n"
        print "\nNot Enough Inputs." 
        print "Need at least 3 inputs: listZero, listFlat, listSpec"
        print "Optional inputs: overwrite= , low_sig= , high_sig=  "
        print "Example:"
        print "\n>>> python imcombine.py listZero listFlat listSpec \n"
        print "\n====================\n"
    
    # Unpack list from command line and take their names # 
    scriptname = args[0]
    zero_list = rt.Read_List( args[1] )
    name_zero = args[1].lstrip("list") + '.fits'
    flat_list = rt.Read_List( args[2] )
    name_flat = args[2].lstrip("list") + '.fits'
    spec_list = rt.Read_List( args[3] ) 
    name_spec = args[3].lstrip("list") + '.fits'    
    
    # Default values for special commands if None are Given these dont change #   
    overwrite = False # dont give imcombine permision to overwrite files # 
    lo_sig = 10
    hi_sig = 3
    method = 'median' # method used to combine images 
    
    # If overwrite special comand is given # 
    if nargs >= 5:
        overwrite = args[4]
    # If low_sigma and high_sigma values are given # 
    if nargs >= 7: 
        lo_sig = float(args[5])
        hi_sig = float(args[6]) 
    # If method is given #  
    if nargs >= 8:
        method = args[7]
    
    # The rest of the code runs the reduction procces up to apall # 
    
    # Combine Zeros # 
    comb_zero = rt.imcombine(zero_list, name_zero, method, lo_sig= lo_sig, 
                        hi_sig= hi_sig, overwrite= overwrite)
    # Bias Subtract Flats # 
    b_flat_list = rt.Bias_Subtract(flat_list, comb_zero )
    # Combine Bias Subtracted Flats # 
    comb_flat = rt.imcombine(b_flat_list, 'b.'+name_flat, method, 
                        lo_sig= lo_sig, hi_sig= hi_sig, overwrite= overwrite)
    # Normalize Flat (divide by average of counts) # 
    nb_flat = rt.Norm_Flat(comb_flat)
    # Bias Subtract Spectra # 
    b_spec_list = rt.Bias_Subtract(spec_list, comb_zero)
    # Flat Field Individual Spectra # 
    fb_spec_list = rt.Flat_Field(b_spec_list, nb_flat )
    # Combine Spectra # 
    comb_fb_spec = rt.imcombine(fb_spec_list, 'fb.'+name_spec, method, 
                        lo_sig= lo_sig, hi_sig= hi_sig, overwrite= overwrite)
    # Trim Spectra # 
    rt.Trim_Spec(comb_fb_spec); 
                        
    print "\n====================\n"
    print "Done. Ready for Apeture Extraction.\n"
    
# ===========================================================================
# ===========================================================================
# ===========================================================================