"""Settings file for astrometry.py.

Author Lukas Wenzl
written in python 3

"""

#Author Lukas Wenzl
#written in python 3

###########################
#settings
#calibration upper and lower limit
# FOV_MAX= 10 #arcmin (Field of View)
# FOV_MIN= 2 #arcmin (Field of View)
#
# #WCS initial guess if no info in header or wrong info in header (currently hardcoded)
# CDELT1_GUESS = -8.0000000000000E-5#-5!!!! #about 0.288 arcsec per pixel
# CDELT2_GUESS = 8.0000000000003E-5#-5!!!
# PC1_1_GUESS =1. #no rotation or scaling
# PC1_2_GUESS =0.
# PC2_1_GUESS=0.
# PC2_2_GUESS =1.


#source detection
#FWHM = 4. #pixels, seeing in pixel: standard:4
N_BRIGHTEST_SOURCES = 200 # only use the XXX brightest sources in the image
DETECTION_ABSOLUTE_THRESHOLD = None #set to replace sigma threshold by absolute threshold


#RMS calculation
RMS_PX_THRESHOLD = 10 #threshold for rms calculation



#transformation
USE_N_SOURCES = 30 #number of sources to be used in fast mode
FASTMODE_THRESHOLD = 0.5 #half the sources in fast mode have to be detected otherwise I try again without fast mode
OFFSET_BINWIDTH = 1 #binning for peak finding to determine x y offset, default: 1px

# #Hubbe Deep Field:
# FWHM = 7. #pixels, seeing in pixel
# N_BRIGHTEST_SOURCES = 300 # only use the XXX brightest sources in the image
# DETECTION_SIGMA_THRESHOLD = 5 #threshold for detection. standard only use detections above 5 sigma
# DETECTION_ABSOLUTE_THRESHOLD = 0.05
# OFFSET_BINWIDTH = 3
# RMS_PX_THRESHOLD = 20

#query catalog_data
#MAG_PS = "brightest" #magnitue of PANSTARRS to be used: either: 'gmag', "zmag", .., "brightest"
#MAG_GAIA = "phot_g_mean_mag" #magnitude of GAIA to be used
