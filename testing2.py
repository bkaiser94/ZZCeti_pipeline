'''
From S. Crawford on github

st.
'''
import numpy as np

def boxcar_smooth(spec, smoothwidth):
    # get the average wavelength separation for the observed spectrum
    # This will work best if the spectrum has equal linear wavelength spacings
    wavespace = np.diff(spec.warr).mean()
    # kw
    kw = int(smoothwidth / wavespace)
    # make sure the kernel width is odd
    if kw % 2 == 0:
        kw += 1
    kernel = np.ones(kw)
    # Conserve flux
    kernel /= kernel.sum()
    smoothed = spec.farr.copy()
    smoothed[(kw / 2):-(kw / 2)] = np.convolve(spec.farr, kernel, mode='valid')
    return smoothed
