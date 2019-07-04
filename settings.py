"""Settings file for astrometry.py.

Author Lukas Wenzl
written in python 3

"""

#Author Lukas Wenzl
#written in python 3

###########################
#settings
#calibration upper and lower limit
FOV_MAX= 10 #arcmin (Field of View)
FOV_MIN= 2 #arcmin (Field of View)

#WCS initial guess if no info in header or wrong info in header
CDELT1_GUESS = -8.0000000000000E-5#-5!!!! #about 0.288 arcsec per pixel
CDELT2_GUESS = 8.0000000000003E-5#-5!!!
PC1_1_GUESS =1. #no rotation or scaling
PC1_2_GUESS =0.
PC2_1_GUESS=0.
PC2_2_GUESS =1.


#source detection


#peak detection


#transformation
USE_N_SOURCES = 30 #number of sources to be used in fast mode
FASTMODE_THRESHOLD = 0.5 #half the sources in fast mode have to be detected otherwise I try again without fast mode

#query catalog_data
#MAG_PS = "brightest" #magnitue of PANSTARRS to be used: either: 'gmag', "zmag", .., "brightest"
#MAG_GAIA = "phot_g_mean_mag" #magnitude of GAIA to be used
