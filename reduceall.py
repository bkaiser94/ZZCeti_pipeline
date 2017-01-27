import ReduceSpec
import spectral_extraction
import Wavelength_Calibration
import continuum_normalization
from glob import glob



#=========================
#Begin Fits Reduction
#=========================

ReduceSpec.reduce_now(['script_name','listZero','listFlat','listSpec','listFe'])


#========================
#Begin Spectral Extraction
#========================
print 'Beginning spectral extraction.'
spec_files = sorted(glob('ctfb*fits'))
single_spec_list = []
for x in spec_files:
    if ('ctfb.0' or 'ctfb.1') in x:
        single_spec_list.append(x)
for x in single_spec_list:
    spec_files.remove(x)
spec_files = sorted(spec_files)

lamp_file_blue = sorted(glob('tFe*blue*fits'))
lamp_file_red = sorted(glob('tFe*red*fits'))


#Get the trace and FWHM files if they exist. 
'''
trace_files = sorted(glob('*trace*npy'))
if len(trace_files) != len(spec_files):
    trace_files = [None] * len(spec_files)
FWHM_files = sorted(glob('*poly.npy'))
if len(FWHM_files) != len(spec_files):
    FWHM_files = [None] * len(spec_files)
'''

#Search for FWHM and trace file for each spectrum. If it does not exist, these go to None and will be fit and saved during the extraction.
trace_files = []
FWHM_files = []
for x in spec_files:
    trace_name = '*' + x[5:-5] + '*trace.npy'
    new_trace = glob(trace_name)
    if len(new_trace) == 0:
        trace_files.append(None)
    else:
        trace_files.append(new_trace[0])
    fwhm_name = '*' + x[5:-5] + '*poly.npy'
    new_fwhm = glob(fwhm_name)
    if len(new_fwhm) == 0:
        FWHM_files.append(None)
    else:
        FWHM_files.append(new_fwhm[0])


for x in spec_files:
    if 'blue' in x.lower():
        lamp_file = lamp_file_blue[0]
    elif 'red' in x.lower():
        lamp_file = lamp_file_red[0]
    FWHM_thisfile = FWHM_files[spec_files.index(x)]
    trace_thisfile = trace_files[spec_files.index(x)]
    if trace_thisfile != None:
        trace_exist_file = True
    else:
        trace_exist_file = False
    print ''
    print x, lamp_file,trace_thisfile, FWHM_thisfile
    #Must add in option of not have trace file or FWHM file
    #if no FWHMfile, FWHMfile=None
    spectral_extraction.extract_now(x,lamp_file,FWHMfile=FWHM_thisfile,tracefile=trace_thisfile,trace_exist=trace_exist_file)
    #spectral_extraction.extract_now('ctfb.wd1307-017_930_blue.fits','tFe_ZZCeti_930_blue_long.fits',FWHMfile='ctfb.wd1307-017_930_blue_poly.npy',tracefile='tfb.wd1307-017_930_blue_trace.npy',trace_exist=True)


#=========================
# Begin Wavelength Calibration
#=========================
print '\n Beginning Wavelength Calibration'
spec_files = sorted(glob('ctfb*ms.fits'))
lamp_files = sorted(glob('tFe*ms.fits'))
offset_file = glob('offsets.txt') #Offset file must be structured as blue, then red
if len(offset_file) == 0:
    offset_file = None
else:
    offset_file = offset_file[0]

#print spec_files
#print lamp_files
for x in lamp_files:
    #print ''
    #print x
    for y in spec_files:
        #print y[5:y.find('_930')], y[y.find('930_'):y.find('.ms')]
        #print y[5:y.find('_930')], y[y.find('_930'):y.find('_930')+8]
        if (y[5:y.find('_930')] in x) and (y[y.find('_930'):y.find('_930')+8] in x):
            print x, y, offset_file
            #Wavelength_Calibration.calibrate_now('tFe_ZZCeti_930_blue_long_.wd1307-017.ms.fits','ctfb.wd1307-017_930_blue.ms.fits','yes','yes')
            Wavelength_Calibration.calibrate_now(x,y,'yes','yes',offset_file,plotall=False)


#=========================
#Begin Continuum Normalization
#=========================
print '\n Begin continuum normalization.'
continuum_files = sorted(glob('wctfb*ms.fits'))
print continuum_files
x = 0
while x < len(continuum_files):
    if x == len(continuum_files)-1:
        print continuum_files[x]
        continuum_normalization.normalize_now(continuum_files[x],None,False,plotall=False)
        x += 1
    elif continuum_files[x][0:continuum_files[x].find('930')] == continuum_files[x+1][0:continuum_files[x].find('930')]:
        print continuum_files[x],continuum_files[x+1]
        continuum_normalization.normalize_now(continuum_files[x],continuum_files[x+1],True,plotall=False)
        x += 2
    else:
        print continuum_files[x]
        continuum_normalization.normalize_now(continuum_files[x],None,False,plotall=False)
        x += 1

#continuum_normalization.normalize_now('wctfb.wd1307-017_930_blue.ms.fits',None,False)
#exit()
