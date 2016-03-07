# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 21:53:13 2016

@author: jmeza
"""



# Imports # ==================================================================

import numpy as np
import pyfits as pf
import matplotlib.pyplot as plt
   
# ============================================================================
    
# Combined Flat # 
# flat= "b.QuartzFlat_ZZCeti_930_blue.fits"; order= 4; 
flat= "b.DomeFlat_ZZCeti_930_red.fits"; order= 4; 

# Read Flat # 
flat_data = pf.getdata(flat)
#plt.imshow(flat_data[0,:,:], origin= 'lower' )#cmap='bone')
#plt.show()

# Gather the rows from 95 to 105 and combine them % 
fit_data= np.median(flat_data[0][95:105], axis=0)
plt.plot(fit_data)

# trim the overscan from the fi data # 
lo= 10;
hi= 2055; 
 
fit_data= fit_data[lo:hi]

# Fit polynomial to fit data # 
X= range(lo,hi)
coeff= np.polyfit(X, fit_data, order )
line= np.poly1d(coeff)


# Calculate Residuals # 
Res= fit_data-line(X)
Per= Res/fit_data

# Plot Fit # 
plt.plot(X, line(X),'r')
plt.title(flat+"      Fit Order: "+str(order))
plt.show()
 
# Plot Fit # 
plt.plot(X, fit_data/line(X),'r')
plt.title("Data Divide Fit"+"      Fit Order: "+str(order))
plt.show()

# Plot Res and Diff # 
plt.subplot(2,1,1)
plt.plot(X,Res)
plt.title('Residuals in Counts')
plt.xlabel('Pixels')
plt.ylabel('Residuals')
plt.subplot(2,1,2)
plt.plot(X,Per)
plt.title('Percent Difference')
plt.xlabel('Pixels')
plt.ylabel('Percent')
plt.show()


# Divide Flat by Polynomial # 




