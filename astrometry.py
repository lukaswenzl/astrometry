"""This aims to do very fast and simple astrometric calibration of images.

Author Lukas Wenzl
written in python 3

"""

#Author Lukas Wenzl
#written in python 3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# from sklearn.externals import joblib ##load pkl (pickle) file
#
from datetime import datetime
#
# #for parsing the arguments for the file
from argparse import ArgumentParser

import get_catalog_data as query
import get_transformation as register
import settings as s
#
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord

#import photutils
from photutils import DAOStarFinder
from photutils import aperture_photometry, CircularAperture
#from astropy.stats import mad_std
from astropy.stats import sigma_clipped_stats
from matplotlib.colors import LogNorm
from astropy.wcs import WCS
from astropy.wcs import Wcsprm
from astropy.table import Table
import warnings
import os
import copy


def determine_if_fit_converged(dic_rms, catalog, observation, wcsprm, NAXIS1, NAXIS2, silent=False):
    """Cut off criteria for deciding if the astrometry fit worked.

    Parameters
    ----------
    dic_rm : dictionary
        Info about the rms
    Returns
    -------
    converegd : boolean
        Boolean about wether it converged

    """
    converged = False

    N_obs = observation.shape[0]

    catalog_on_sensor = wcsprm.s2p(catalog[["ra", "dec"]], 1)['pixcrd']
    cat_x = np.array( [catalog_on_sensor[:,0] ])
    cat_y = np.array( [catalog_on_sensor[:,1] ])
    mask = (cat_x >= 0) & (cat_x <= NAXIS1) & (cat_y >= 0) & (cat_y <= NAXIS2)
    N_cat = mask.sum()

    N = min([N_obs, N_cat])
    if(not silent):
        print("Conditions for converegence of fit: at least 5 matches within {}px and at least 0.2 * minimum ( observed sources, catalog sources in fov) matches with the same radius".format(dic_rms["radius_px"]))
    if (dic_rms["matches"] > 5 and  dic_rms["matches"] > 0.1 * N ):
        converged = True

    return converged

def find_sources(image, vignette=3,vignette_rectangular=1., cutouts=None,sigma_threshold_for_source_detection=5, only_rectangle=None, FWHM=4):
    """Find surces in the image. Uses DAOStarFinder with symmetric gaussian kernels. Only uses 5 sigma detections. It only gives the 200 brightest sources or less.
    This has to work well for the later calculations to work. Possible issues: low signal to noise image, gradient in the background

    Parameters
    ----------
    image
        Observed image (without background)
    vignette : float
        Cut off courners with a vignette. Default: nothing cut off

    Returns
    -------
    observaion : dataframe
        Pandas dataframe with the sources

    """
    #find sources
    #bkg_sigma = mad_std(image)
    mean, median, std = sigma_clipped_stats(image, sigma=3.0)
    image = copy.copy(image)
    # #only search sources in a circle with radius <vignette>
    if(vignette < 3):
        sidelength = np.max(image.shape)
        x = np.arange(0, image.shape[1])
        y = np.arange(0, image.shape[0])
        vignette = vignette * sidelength/2
        mask = (x[np.newaxis,:]-sidelength/2)**2 + (y[:,np.newaxis]-sidelength/2)**2 < vignette**2
        image[~mask] = median

    #ignore a fraction of the image at the corner
    if(vignette_rectangular < 1.):
        sidelength_x = image.shape[1]
        sidelength_y = image.shape[0]
        cutoff_left = (1.-vignette_rectangular)*sidelength_x
        cutoff_right = vignette_rectangular*sidelength_x
        cutoff_bottom = (1.-vignette_rectangular)*sidelength_y
        cutoff_top = vignette_rectangular*sidelength_y
        x = np.arange(0, image.shape[1])
        y = np.arange(0, image.shape[0])
        left =    x[np.newaxis,:] > cutoff_left 
        right =   x[np.newaxis,:] < cutoff_right
        bottom =  y[:,np.newaxis] > cutoff_bottom
        top =     y[:,np.newaxis] < cutoff_top
        mask = (left*bottom)*(right*top)
        image[~mask] = median

    
    #cut out rectangular regions of the image, [(xstart, xend, ystart, yend)]
    if(cutouts != None):
        sidelength = np.max(image.shape)
        x = np.arange(0, image.shape[1])
        y = np.arange(0, image.shape[0])
        for cutout in cutouts:
            print("cutting out {}".format(cutout))
            left =    x[np.newaxis,:] > cutout[0] 
            right =   x[np.newaxis,:] < cutout[1] 
            bottom =  y[:,np.newaxis] > cutout[2]
            top =     y[:,np.newaxis] < cutout[3]
            mask = (left*bottom)*(right*top)
            image[mask] = median

    #use only_rectangle within image format: (xstart, xend, ystart, yend)
    if(only_rectangle != None):
        sidelength = np.max(image.shape)
        x = np.arange(0, image.shape[1])
        y = np.arange(0, image.shape[0])
        print("only using {}".format(only_rectangle))
        left =    x[np.newaxis,:] > only_rectangle[0] 
        right =   x[np.newaxis,:] < only_rectangle[1] 
        bottom =  y[:,np.newaxis] > only_rectangle[2]
        top =     y[:,np.newaxis] < only_rectangle[3]
        mask = (left*bottom)*(right*top)
        image[~mask] = median

    #daofind = DAOStarFinder(fwhm=4., threshold=5.*std, brightest=200)
    #daofind = DAOStarFinder(fwhm=7., threshold=0.6, brightest=400 )
    daofind = DAOStarFinder(fwhm=FWHM, threshold=sigma_threshold_for_source_detection *std, brightest=s.N_BRIGHTEST_SOURCES )
    if(s.DETECTION_ABSOLUTE_THRESHOLD is not None):
        daofind = DAOStarFinder(fwhm=FWHM, threshold=s.DETECTION_ABSOLUTE_THRESHOLD, brightest=s.N_BRIGHTEST_SOURCES )




    sources = daofind(image)
    for col in sources.colnames:
        sources[col].info.format = '%.8g'  # for consistent table output

    #changed order of positions to [(x,y), (x,y),...] for compatibility with photutils 1.4
    xcenters = np.array(sources['xcentroid'])
    ycenters = np.array(sources['ycentroid'])
    print("Changed convention to start with idx 1 for pixels")
    xcenters = xcenters +1.
    ycenters = ycenters +1.
    positions = [(xcenters[i], ycenters[i]) for i in range(len(xcenters))]
    apertures = CircularAperture(positions, r=4.)
    phot_table = aperture_photometry(image, apertures)
    for col in phot_table.colnames:
        phot_table[col].info.format = '%.8g'  # for consistent table output

    observation = Table(phot_table).to_pandas()

    #through out candidates where the star finder messed up
    observation = observation.query("aperture_sum > "+str(5*std))
    return  observation


def write_wcs_to_hdr(original_filename, wcsprm, report, hdul_idx=0):
    """Update the header of the fits file itself.

    Parameters
    ----------
    original_filename : str
        Original filename of the fits file
    wcsprm : astropy.wcs.wcsprm
        World coordinate system object decsribing translation between image and skycoord

    """
    with fits.open(original_filename) as hdul:

        hdu = hdul[hdul_idx]
        hdr_file = hdu.header

        #new_header_info = wcs.to_header()

        wcs =WCS(wcsprm.to_header())

        #I will through out CD which contains the scaling and separate into pc and Cdelt
        for old_parameter in ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', "PC1_1", "PC1_2", "PC2_1", "PC2_2"]:
            if (old_parameter in hdr_file):
                del hdr_file[old_parameter]

        hdr_file.update(wcs.to_header())
        repr(hdr_file)

        #adding report
        hdr_file['AST_SCR'] = ("Astrometry", 'Astrometry by Lukas Wenzl')
        hdr_file['AST_CAT'] = (report["catalog"], 'Catalog used')
        hdr_file['AST_MAT'] = (report["matches"], 'Number of catalog matches')
        hdr_file['AST_RAD'] = (report["match_radius"], 'match radius in pixel')
        hdr_file['AST_CONV'] = (report["converged"], "T or F for astroemtry converged or not")



        hdu.header = hdr_file
        hdul[hdul_idx] = hdu
        #removing fits ending
        name_parts = original_filename.rsplit('.', 1)
        hdul.writeto(name_parts[0]+'_astro.fits', overwrite=True)
        print("file written.")



def read_additional_info_from_header(wcsprm, hdr, RA_input=None, DEC_input=None, projection_ra=None, projection_dec=None, ignore_header_rot=False, radius = -1., silent=False):
    """Tries to handle additional or missing data from the header.
    If your picture does not conform with fits standards this might be the section to edit to make the code work

    Parameters
    ----------
    wcsprm : astropy.wcs.wcsprm
        World coordinate system object decsribing translation between image and skycoord
    hdr : header
    ...
    """
    fov_radius = 4 #arcmin radius to include field of view
    if(radius > 0):
        fov_radius = radius
    INCREASE_FOV_FLAG = False # increase the field to view by 50% to search in catalog if position on sky is inaccurate
    PIXSCALE_UNCLEAR = False

    keywords_check = ["PIXSCALE", "NAXIS1", "NAXIS2", "RA", "DEC"] #list of possible keywords the scs parser might miss
    keywords_present = [] # list of keywords that are actually present
    for i in keywords_check:
        if(i in hdr.keys()):
            keywords_present.append(i)

    if("NAXIS1" not in keywords_present or "NAXIS2" not in keywords_present ):
        print("ERROR: NAXIS1 or NAXIS2 missing in file. Please add!")
    else:
        axis1 = hdr["NAXIS1"]
        axis2 = hdr["NAXIS2"]

    pc = wcsprm.get_pc()
    cdelt = wcsprm.get_cdelt()
    wcs_pixscale = (pc @ cdelt )
    if((np.abs(wcs_pixscale[0])) < 1e-7 or (np.abs(wcs_pixscale[1])) < 1e-7 or
        (np.abs(wcs_pixscale[0])) > 5e-3 or (np.abs(wcs_pixscale[1])) > 5e-3):
        if(not silent):
            print("pixelscale is completely unrealistic. Will guess")
        print(wcs_pixscale)
        guess = 8.43785734e-05
        #guess = 6.94444461259988e-05
        wcsprm.pc = [[1,0],[0,1]]
        wcsprm.cdelt = [guess, guess]
        if(not silent):
            print("Changed pixelscale to {:.3g} deg/arcsec".format(guess))
        PIXSCALE_UNCLEAR = True
    if(ignore_header_rot):
        wcsprm.pc = [[1,0],[0,1]]
        #wcsprm.cdelt = [8.0006871225376e-05, 8.0006871225376e-05]
    if("PIXSCALE" in keywords_present):
        #normal around 0.450000 / arcsec/pixel, for now i assume arcsec per pixel
        pixscale = hdr["PIXSCALE"]
        if("deg" in hdr.comments['PIXSCALE']): #correction if in deg/pixel
            pixscale = pixscale *60*60
        x_size = axis1 * pixscale /60# arcmin
        y_size = axis2 * pixscale /60# arcmin

        if 20 > x_size > 0.5 and 20 > y_size> 0.5 :
            #pixscale is sensical
            #Now: is the pixscale of the current wcs realistic?
            pc = wcsprm.get_pc()
            cdelt = wcsprm.get_cdelt()
            wcs_pixscale = (pc @ cdelt )
            pixscale = pixscale /60 /60 #pixelscale now in deg / pixel
            if( wcs_pixscale[0]/pixscale < 0.1 or wcs_pixscale[0]/pixscale > 10 or wcs_pixscale[1]/pixscale < 0.1 or wcs_pixscale[1]/pixscale > 10):
                #check if there is a huge difference in the scales
                #if yes then replace the wcs scale with the pixelscale info
                wcsprm.pc = [[1,0],[0,1]]

                wcsprm.cdelt = [pixscale, pixscale]
                if(not silent):
                    print("changed pixelscale to {:.3g} deg/arcsec".format(pixscale))
                fov_radius = (x_size/2+y_size/2)/np.sqrt(2) #try to get corners
                PIXSCALE_UNCLEAR=True


    if(np.array_equal(wcsprm.crpix, [0,0])):
        #centrl pixel seems to not be in header, better set in middle
        wcsprm.crpix = [axis1/2, axis2/2]

    if(np.array_equal(wcsprm.crval, [0,0] )):
        ###sky position not found. Maybe there is some RA and DEC info in the header:
        INCREASE_FOV_FLAG = True
        if ("RA" in keywords_present and "DEC" in keywords_present): ##carefull degree and hourangle!!!
            wcsprm.crval = [hdr["RA"], hdr["DEC"]]
            if(not silent):
                print("Found ra and dec information in the header")
            print(wcsprm.crval)
            if(not silent):
                print("Is this position within the field of view in degrees? otherwise it will not work. In that case give a more accurate position as an argument: -ra XX -dec XX both in degrees")

    if (RA_input is not None): #use user input if provided
        wcsprm.crval = [RA_input, wcsprm.crval[1]]
        wcsprm.crpix = [axis1/2, wcsprm.crpix[1]]

    if (DEC_input is not None):
        wcsprm.crval = [wcsprm.crval[0], DEC_input]
        wcsprm.crpix = [wcsprm.crpix[0], axis2/2, ]


    if(np.array_equal(wcsprm.crval, [0,0] )):
        print(">>>>>>>>>WARNING")
        print("No rough sky position was found for this object. Please add as -ra XX -dex XX both in degress. Adding the position as keywords in the fits file header will also work. The keywords are RA and DEC. The program expects the values in degrees. ")

    if(np.array_equal(wcsprm.ctype, ["",""])):
        INCREASE_FOV_FLAG = True
        if(projection_ra is not None and projection_dec is not None):
            wcsprm.ctype = [ projection_ra,  projection_dec]
        else:
            wcsprm.ctype = [ 'RA---TAN',  'DEC--TAN'] #this is a guess
            print(">>>>>>>>>WARNING")
            print("The wcs in the header has no projection specified. Will guess 'RA---TAN', 'DEC--TAN' (gnomonic projection) if this is incorrect the fit will fail. You can specify the projection via -projection_ra XX -projection_dec XX")
            print("make sure you do not use quotations, example: -proj1 RA---TAN -proj2 DEC--TAN")
    if(INCREASE_FOV_FLAG):
        fov_radius = fov_radius*2.5
    return wcsprm, fov_radius, INCREASE_FOV_FLAG, PIXSCALE_UNCLEAR



def parseArguments():
    """Parse the given Arguments when calling the file from the command line.

    Returns
    -------
    arg
        The result from parsing.

    """
    # Create argument parser
    #parser = argparse.ArgumentParser()
    parser = ArgumentParser()


    # Positional mandatory arguments
    parser.add_argument("input", nargs='+',help="Input Image with .fits ending. A folder or multiple Images also work", type=str)
    parser.add_argument("-c", "--catalog", help="Catalog to use for position reference ('PS', '2MASS' or 'GAIA') or name of local file", type=str, default="PS")

    parser.add_argument("-s", "--save_images", help="Set True to create _image_before.pdf and _image_after.pdf", type=bool, default=False)
    parser.add_argument("-p", "--images", help="Set 0 to not create the plots.", type=int, default=1)

    parser.add_argument("-v", "--verbose", help="More console output about what is happening. Helpfull for debugging.", type=bool, default=False)
    parser.add_argument("-silent", "--silent", help="Almost no console output. Use when running many files. turns verbose off, turns plots off, turns ignore warning on. ", type=bool, default=False)



    parser.add_argument("-w", "--ignore_warnings", help="Set False to see all Warnings about the header if there is problems. Default is to ignore most warnings.", type=bool, default=True)

    parser.add_argument("-ra", "--ra", help="Set ra by hand in degrees. Should not be necessary if the info is in the fits header.", type=float, default=None)
    parser.add_argument("-dec", "--dec", help="Set dec by hand in degrees. Should not be necessary if the info is in the fits header.", type=float, default=None)

    parser.add_argument("-proj1", "--projection_ra", help="Set the projection by hand. Should not be necessary if the info is in the fits header. example: -proj1 RA---TAN -proj2 DEC--TAN", type=str, default=None)
    parser.add_argument("-proj2", "--projection_dec", help="Set the projection by hand. Should not be necessary if the info is in the fits header. example: -proj1 RA---TAN -proj2 DEC--TAN", type=str, default=None)


    parser.add_argument("-rot_scale", "--rotation_scaling", help="By default rotation and scaling is determined. If wcs already contains this info and the fit fails you can try deactivating this part by setting it to 0", type=int, default=1)
    parser.add_argument("-xy_trafo", "--xy_transformation", help="By default the x and y offset is determined. If wcs already contains this info and the fit fails you can try deactivating this part by setting it to 0", type=int, default=1)
    parser.add_argument("-fine", "--fine_transformation", help="By default a fine transformation is applied in the end. You can try deactivating this part by setting it to 0", type=int, default=1)

    parser.add_argument("-vignette", "--vignette", help="Do not use corner of the image. Only use the data in a circle around the center with certain radius. Default: not used. Set to 1 for circle that touches the sides. Less to cut off more", type=float, default=3)
    parser.add_argument("-vignette_rec", "--vignette_rectangular", help="Do not use corner of the image. Cutoff 1- <value> on each side of the image. Default: not used. 1 uses the full image, 0.9 cuts 10%% on each side", type=float, default=1.)

    parser.add_argument("-cutout", "--cutout", help="Cutout bad recangle from Image. Specify corners in pixel as -cutout xstart xend ystart yend. You can give multiple cutouts ", nargs='+',action="append", type=int)

    parser.add_argument("-ignore_header_rot", "--ignore_header_rot", help="Set to 0 to ignore rotation information contained in header", type=int, default=0)

    parser.add_argument("-r", "--radius", help="Change the download radius for how many catalog objects should be downloaded in arcmin. ", type=float, default=-1)
    parser.add_argument("-save_bad_result", "--save_bad_result", help="Save result even when it is considered bad by setting this to 1. By default bad results will not be saved", type=int, default=0)

    parser.add_argument("-sigma_threshold_for_source_detection", "--sigma_threshold_for_source_detection", help="Sigma threshold for source detection on the image", type=float, default=5)
    parser.add_argument("-seeing", "--seeing", help="Average FWHM of sources. Default 4px", type=float, default=4)

    parser.add_argument("-high_res", "--high_resolution", help="Set true to indicate higher resolution image like HST to allow larger steps for fine transformation. Default False for smaller telescopes", type=bool, default=False)
    parser.add_argument("-hdul_idx", "--hdul_idx", help="Index of the image data in the fits file. Default 0", type=int, default=0)
    parser.add_argument("-filename_for_sources", "--filename_for_sources", help="Save the sky positions of the sources to file to calibrate another image with it. Set to the filename without extension. Default: not used", type=str, default=None)

    # Print version
    parser.add_argument("--version", action="version", version='%%(prog)s - Version 1.5') #
    #changelog
    #version 0.0 proof of concept
    #version 0.1 alpha version
    #version 1.0 first public version
    #version 1.2 included option to not show images
    #version 1.3 improved support for high resolution images
    #version 1.4 bugfixes and ability to calibrate based on another image
    #version 1.5 added compatibility with photutils v1.4

    # Parse arguments
    args = parser.parse_args()

    return args


def main():
    """Perform astrometry for the given file."""
    print("Program version: 1.5")
    StartTime = datetime.now()
    args = parseArguments()

    verbose = args.verbose
    images = args.images
    ignore_warnings = args.ignore_warnings
    if(args.silent):
        verbose = False
        images = False
        ignore_warnings = True

    if(args.images):
        plt.ioff()

    if(args.ignore_warnings):
        warnings.simplefilter('ignore', UserWarning)

    #sample header keywords
    #     OBJECT  = 'P016+03_P1_JKdeep'  / Original target
    # RA      = ' 01:06:37.759'               / 01:06:37.7 RA (J2000) pointing
    # DEC     = ' 03:32:36.096'               / 03:32:36.0  DEC (J2000) pointing
    # EQUINOX =                2000.          / Standard FK5 (years)
    # RADECSYS= 'FK5     '                    / Coordinate reference frame
    # CRVAL1  =             16.65733          / 01:06:37.7, RA at ref pixel
    # CRVAL2  =              3.54336          / 03:32:36.0, DEC at ref pixel
    # CRPIX1  =                 447. /Ref pixel in X
    # CRPIX2  =                 452. / Ref pixel in Y
    # CDELT1  =  -8.0000000000000E-5 / SS arcsec per pixel in RA
    # CDELT2  =  8.00000000000003E-5 / SS arcsec per pixel in DEC
    # CTYPE1  = 'RA---TAN'                    / pixel coordinate system
    # CTYPE2  = 'DEC--TAN'                    / pixel coordinate system
    # PC1_1   =             0.000000          / Translation matrix element
    # PC1_2   =             1.000000          / Translation matrix element
    # PC2_1   =            -1.000000          / Translation matrix element
    # PC2_2   =             0.000000          / Translation matrix element

    fits_image_filenames = args.input

    #if directory given search for appropriate fits files

    if(os.path.isdir(fits_image_filenames[0])):
        print("detected a directory. Will search for fits files in it")
        path = fits_image_filenames[0]
        fits_image_filenames = []
        for file in os.listdir(path):
            if file.endswith(".fits") and "_astro" not in file:
                fits_image_filenames.append(path+"/"+file)
        print(fits_image_filenames)

    multiple = False
    if(len(fits_image_filenames)>1):
        multiple = True
    not_converged = []
    converged_counter = 0
    for fits_image_filename in fits_image_filenames:

        result,_  = astrometry_script(fits_image_filename, catalog=args.catalog, rotation_scaling=0, xy_transformation=args.xy_transformation, fine_transformation=args.fine_transformation,
        images=images, vignette=args.vignette,vignette_rectangular=args.vignette_rectangular, cutouts=args.cutout, ra=args.ra, dec=args.dec, projection_ra=args.projection_ra, projection_dec=args.projection_dec, verbose=verbose, save_images=args.save_images, ignore_header_rot=args.ignore_header_rot, radius = args.radius, save_bad_result=args.save_bad_result, silent =args.silent, sigma_threshold_for_source_detection= args.sigma_threshold_for_source_detection, high_res=args.high_resolution, hdul_idx=args.hdul_idx, filename_for_sources=args.filename_for_sources, FWHM=args.seeing)

        if((not result) and args.rotation_scaling):
            print("Did not converge. Will try again with full rotation and scaling")
            result, _ = astrometry_script(fits_image_filename, catalog=args.catalog, rotation_scaling=args.rotation_scaling, xy_transformation=args.xy_transformation, fine_transformation=args.fine_transformation,
            images=images, vignette=args.vignette,vignette_rectangular=args.vignette_rectangular, cutouts=args.cutout, ra=args.ra, dec=args.dec, projection_ra=args.projection_ra, projection_dec=args.projection_dec, verbose=verbose, save_images=args.save_images, ignore_header_rot=args.ignore_header_rot, radius = args.radius, save_bad_result=args.save_bad_result, silent=args.silent, sigma_threshold_for_source_detection=args.sigma_threshold_for_source_detection, high_res=args.high_resolution, hdul_idx=args.hdul_idx, filename_for_sources=args.filename_for_sources, FWHM=args.seeing)

        if(result):
            print("Astrometry was determined to be good.")
            converged_counter = converged_counter+1
        else:
            print("Astrometry was determined to be bad.")
            not_converged.append(fits_image_filename)
            if(args.save_bad_result):
                print("Result was saved anyway")
            else:
                print("Result was not saved.")
        # print("")
        # print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        # print("> Astrometry for {} ".format(fits_image_filename))
        #
        # with fits.open(fits_image_filename) as hdul:
        #     #print(hdul.info())
        #     if(args.verbose):
        #         print("if image is not at first position in the fits file the program will break later on")
        #     #print(hdul[0].header)
        #
        #     hdu = hdul[0]
        #     #hdu.verify('fix')
        #     hdr = hdu.header
        #
        #
        #     image_or = hdul[0].data.astype(float)
        #     median = np.nanmedian(image_or)
        #     image_or[np.isnan(image_or)]=median
        #     image = image_or - median
        #
        # observation = find_sources(image, args.vignette)
        # #print(observation)
        #
        # positions = (observation['xcenter'], observation['ycenter'])
        # apertures = CircularAperture(positions, r=4.)
        #
        #
        # #world coordinates
        # print(">Info found in the file -- (CRVAl: position of central pixel (CRPIX) on the sky)")
        # print(WCS(hdr))
        #
        # hdr["NAXIS1"] = image.shape[0]
        # hdr["NAXIS2"] = image.shape[1]
        #
        # #wcsprm = Wcsprm(hdr.tostring().encode('utf-8')) #everything else gave me errors with python 3, seemed to make problems with pc conversios, so i wwitched to the form below
        # wcsprm = WCS(hdr).wcs
        # wcsprm_original = WCS(hdr).wcs
        # if(args.verbose):
        #     print(WCS(wcsprm.to_header()))
        # wcsprm, fov_radius, INCREASE_FOV_FLAG, PIXSCALE_UNCLEAR = read_additional_info_from_header(wcsprm, hdr, args.ra, args.dec, args.projection_ra, args.projection_dec)
        # if(args.verbose):
        #     print(WCS(wcsprm.to_header()))
        #
        # #print(wcsprm)
        # #wcsprm.pc = [[2, 0],[0,1]]
        #
        #
        # #Possibly usefull examples of how to use wcsprm:
        # #print(wcsprm.set())
        # #print(wcsprm.get_pc())
        # #pc = wcsprm.get_pc()
        # #print(np.linalg.det(pc))
        # #print(wcsprm.get_cdelt())
        # #wcs.fix()
        # #print(wcsprm.print_contents())
        # #print(repr(hdr.update(wcsprm.to_header().encode('utf-8')))) #not working
        #
        # #hdu.verify("fix")
        # #print(repr(hdr))
        # #wcs.wcs_pix2world(pixcrd, 1)
        # #wcs.wcs_world2pix(world, 1)
        # #wcs.wcs.crpix = [-234.75, 8.3393]
        # # wcs.wcs.cdelt = np.array([-0.066667, 0.066667])
        # # wcs.wcs.crval = [0, -90]
        # # wcs.wcs.ctype = ["RA---AIR", "DEC--AIR"]
        # # wcs.wcs.set_pv([(2, 1, 45.0)])
        # # For historical compatibility, three alternate specifications of the linear transformations
        # # are available in wcslib. The canonical PCi_ja with CDELTia, CDi_ja, and the deprecated CROTAia
        # # keywords. Although the latter may not formally co-exist with PCi_ja,
        # # the approach here is simply to ignore them if given in conjunction with PCi_ja.
        # # has_pc, has_cd and has_crota can be used to determine which of these alternatives are present in the header.
        # # These alternate specifications of the linear transformation matrix are translated immediately to PCi_ja by set
        # # and are nowhere visible to the lower-level routines. In particular, set resets cdelt to unity if CDi_ja is present
        # # (and no PCi_ja). If no CROTAia is associated with the latitude axis, set reverts to a unity PCi_ja matrix.
        #
        #
        #
        #
        #
        # #get rough coordinates
        # #print(hdr["RA"])
        # #coord = SkyCoord(hdr["RA"], hdr["DEC"], unit=(u.hourangle, u.deg), frame="icrs")
        # coord = SkyCoord(wcsprm.crval[0], wcsprm.crval[1], unit=(u.deg, u.deg), frame="icrs")
        # if(not PIXSCALE_UNCLEAR):
        #     if(wcsprm.crpix[0] < 0 or wcsprm.crpix[1] < 0 or wcsprm.crpix[0] > image.shape[0] or wcsprm.crpix[1] > image.shape[1] ):
        #         print("central value outside of the image, moving it to the center")
        #         coord_radec = wcsprm.p2s([[image.shape[0]/2, image.shape[1]/2]], 1)["world"][0]
        #         coord = SkyCoord(coord_radec[0], coord_radec[1], unit=(u.deg, u.deg), frame="icrs")
        #         #print(wcsprm)
        #
        #
        #
        # #better: put in nice wrapper! with repeated tries and maybe try synchron!
        # print(">Dowloading catalog data")
        # radius = u.Quantity(fov_radius, u.arcmin)#will prob need more
        # catalog_data = query.get_data(coord, radius, args.catalog)
        # #reference = reference.query("mag <20")
        # max_sources = 500
        # if(INCREASE_FOV_FLAG):
        #     max_sources= max_sources*2.25 #1.5 times the radius, so 2.25 the area
        # if(catalog_data.shape[0]>max_sources):
        #     catalog_data = catalog_data.nsmallest(400, "mag")
        #
        # if(args.catalog == "GAIA" and catalog_data.shape[0] < 5):
        #     print("GAIA seems to not have enough objects, will enhance with PS1")
        #     catalog_data2 = query.get_data(coord, radius, "PS")
        #     catalog_data = pd.concat([catalog_data, catalog_data2])
        #     #apertures_catalog = CircularAperture(wcs.wcs_world2pix(catalog_data[["ra", "dec"]], 1), r=5.)
        #     print("Now we have a total of {} sources. Keep in mind that there might be duplicates now  since we combined 2 catalogs".format(catalog_data.shape[0]))
        # elif(args.catalog == "PS" and (catalog_data is None or catalog_data.shape[0] < 5)):
        #     print("We seem to be outside the PS footprint, enhance with GAIA data")
        #     catalog_data2 = query.get_data(coord, radius, "GAIA")
        #     catalog_data = pd.concat([catalog_data, catalog_data2])
        #     #apertures_catalog = CircularAperture(wcs.wcs_world2pix(catalog_data[["ra", "dec"]], 1), r=5.)
        #     print("Now we have a total of {} sources. Keep in mind that there might be duplicates now since we combined 2 catalogs".format(catalog_data.shape[0]))
        #
        # #remove duplicates in catalog?
        #
        # apertures_catalog = CircularAperture(wcsprm.s2p(catalog_data[["ra", "dec"]], 1)['pixcrd'], r=5.)
        #
        #
        # #plotting what we have, I keep it in the detector field, world coordinates are more painfull to plot
        # if(args.images):
        #     fig = plt.figure()
        #     fig.canvas.set_window_title('Input for {}'.format(fits_image_filename))
        #     plt.xlabel("pixel x direction")
        #     plt.ylabel("pixel y direction")
        #     plt.title("Input - red: catalog sources, blue: detected sources in img")
        #     plt.imshow(image,cmap='Greys', origin='lower', norm=LogNorm())
        #     apertures.plot(color='blue', lw=1.5, alpha=0.5)
        #     apertures_catalog.plot(color='red', lw=1.5, alpha=0.5)
        #
        #     plt.xlim(-200,image.shape[0]+200)
        #     plt.ylim(-200,image.shape[1]+200)
        #     if(args.save_images):
        #         name_parts = fits_image_filename.rsplit('.', 1)
        #         plt.savefig(name_parts[0]+"_image_before.pdf")
        #
        # ###tranforming to match the sources
        # print("---------------------------------")
        # print(">Finding the transformation")
        # if(args.rotation_scaling):
        #     print("Finding scaling and rotation")
        #     wcsprm = register.get_scaling_and_rotation(observation, catalog_data, wcsprm, scale_guessed=PIXSCALE_UNCLEAR, verbose=args.verbose)
        # if(args.xy_transformation):
        #     print("Finding offset")
        #     wcsprm,_,_ = register.offset_with_orientation(observation, catalog_data, wcsprm, fast=False , INCREASE_FOV_FLAG=INCREASE_FOV_FLAG, verbose= args.verbose)
        #
        # #correct subpixel error
        # obs_x, obs_y, cat_x, cat_y, distances = register.find_matches(observation, catalog_data, wcsprm, threshold=3)
        # rms = np.sqrt(np.mean(np.square(distances)))
        # best_score = len(obs_x)/(rms+10) #start with current best score
        # fine_transformation = False
        # if(args.fine_transformation):
        #     for i in [2,3,5,8,10,6,4, 20,2,1,0.5]:
        #         wcsprm_new, score = register.fine_transformation(observation, catalog_data, wcsprm, threshold=i)
        #         if(score> best_score):
        #             wcsprm = wcsprm_new
        #             best_score = score
        #             fine_transformation = True
        #     if not fine_transformation:
        #         print("Fine transformation did not improve result so will be discarded.")
        #     else:
        #         print("Fine transformation applied to improve result")
        # #register.calculate_rms(observation, catalog_data,wcs)
        #
        # #make wcsprim more physical by moving scaling to cdelt, out of the pc matrix
        # wcs =WCS(wcsprm.to_header())
        # if(args.verbose):
        #     print(wcs)
        #
        # from astropy.wcs import utils
        # scales = utils.proj_plane_pixel_scales(wcs)
        # print(scales)
        # cdelt = wcsprm.get_cdelt()
        # print(cdelt)
        # scale_ratio = scales/cdelt
        # #print(scale_ratio)
        # pc = np.array(wcsprm.get_pc())
        # pc[0,0] = pc[0,0]/scale_ratio[0]
        # pc[1,0] = pc[1,0]/scale_ratio[1]
        # pc[0,1] = pc[0,1]/scale_ratio[0]
        # pc[1,1] = pc[1,1]/scale_ratio[1]
        # wcsprm.pc = pc
        # wcsprm.cdelt = scales
        # if(args.verbose):
        #     print("moved scaling info to CDelt")
        #     print(WCS(wcsprm.to_header()))
        #
        # #WCS difference before and after
        # print("> Compared to the input the Wcs was changed by: ")
        # scales_original = utils.proj_plane_pixel_scales(WCS(hdr))
        # print("WCS got scaled by {} in x direction and {} in y direction".format(scales[0]/scales_original[0], scales[1]/scales_original[1]))
        # #sources:
        # #https://math.stackexchange.com/questions/2113634/comparing-two-rotation-matrices
        # #https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python/13849249#13849249
        # def unit_vector(vector):
        #     """ Returns the unit vector of the vector.  """
        #     return vector / max(np.linalg.norm(vector), 1e-10)
        # def matrix_angle( B, A ):
        #     """ comment cos between vectors or matrices """
        #     Aflat = A.reshape(-1)
        #     Aflat = unit_vector(Aflat)
        #     Bflat = B.reshape(-1)
        #     Bflat = unit_vector(Bflat)
        #     #return np.arccos((np.dot( Aflat, Bflat ) / max( np.linalg.norm(Aflat) * np.linalg.norm(Bflat), 1e-10 )))
        #     return np.arccos(np.clip(np.dot(Aflat, Bflat), -1.0, 1.0))
        # #print(matrix_angle(wcsprm.get_pc(), wcsprm_original.get_pc()) /2/np.pi*360)
        # rotation_angle = matrix_angle(wcsprm.get_pc(), wcsprm_original.get_pc()) /2/np.pi*360
        # if((wcsprm.get_pc() @ wcsprm_original.get_pc() )[0,1] > 0):
        #     text = "counterclockwise"
        # else:
        #     text = "clockwise"
        # print("Rotation of WCS by an angle of {} deg ".format(rotation_angle)+text)
        # old_central_pixel = wcsprm_original.s2p([wcsprm.crval], 1)["pixcrd"][0]
        # print("x offset: {} px, y offset: {} px ".format(wcsprm.crpix[0]- old_central_pixel[0], wcsprm.crpix[1]- old_central_pixel[1]))
        #
        #
        # #check final figure
        # if(args.images):
        #     fig = plt.figure()
        #     fig.canvas.set_window_title('Result for {}'.format(fits_image_filename))
        #     plt.xlabel("pixel x direction")
        #     plt.ylabel("pixel y direction")
        #     plt.title("Result - red: catalog sources, blue: detected sources in img")
        #     plt.imshow(image,cmap='Greys', origin='lower', norm=LogNorm())
        #     apertures.plot(color='blue', lw=1.5, alpha=0.5)
        #     #apertures_catalog = CircularAperture(wcs.wcs_world2pix(catalog_data[["ra", "dec"]], 1), r=5.)
        #     apertures_catalog = CircularAperture(wcsprm.s2p(catalog_data[["ra", "dec"]], 1)['pixcrd'], r=5.)
        #
        #     apertures_catalog.plot(color='red', lw=1.5, alpha=0.5)
        #     if(args.save_images):
        #         name_parts = fits_image_filename.rsplit('.', 1)
        #         plt.savefig(name_parts[0]+"_image_after.pdf")
        #
        # print("--- Evaluate how good the transformation is ----")
        # register.calculate_rms(observation, catalog_data,wcsprm)
        #
        #
        # #updating file
        # write_wcs_to_hdr(fits_image_filename, wcsprm)
        #
        #
        # print("overall time taken")
        # print(datetime.now()-StartTime)
        # if(args.images):
        #     plt.show()
    if(multiple):
        print(">> Final report:")
        print("Processed {} files, {} of them did converge. The following files failed:".format(len(fits_image_filenames), converged_counter))
        print(not_converged)
    print("-- finished --")

def astrometry_script(filename, catalog="PS", rotation_scaling=True, xy_transformation=True, fine_transformation=True, images=False, vignette=3,vignette_rectangular=1., cutouts=None, ra=None, dec=None, projection_ra=None, projection_dec=None, verbose=False, save_images=False, ignore_header_rot=False, radius=-1., save_bad_result=False, silent=False, sigma_threshold_for_source_detection=5, high_res = False, hdul_idx=0, filename_for_sources=None, FWHM=4):
    """Perform astrometry for the given file. Scripable version"""
    #print("Program version: 1.2")

    report = {}
    if(images):
        plt.ioff()
    warnings.simplefilter('ignore', UserWarning)
    fits_image_filename = filename

    print("> Astrometry for {} ".format(fits_image_filename))

    with fits.open(fits_image_filename) as hdul:
        #print(hdul.info())
        #print(hdul[0].header)

        hdu = hdul[hdul_idx]
        #hdu.verify('fix')
        hdr = hdu.header


        image_or = hdul[hdul_idx].data.astype(float)
        median = np.nanmedian(image_or)
        image_or[np.isnan(image_or)]=median
        image = image_or - median

    observation = find_sources(image, vignette,vignette_rectangular,cutouts, sigma_threshold_for_source_detection, FWHM=FWHM)
    #print(observation)

    #changed order of positions to [(x,y), (x,y),...] for compatibility with photutils 1.4
    xcenters = np.array(observation['xcenter'])
    ycenters = np.array(observation['ycenter'])
    positions = [(xcenters[i], ycenters[i]) for i in range(len(xcenters))]
    apertures = CircularAperture(positions, r=4.)


    #world coordinates
    if(not silent):
        print(">Info found in the file -- (CRVAl: position of central pixel (CRPIX) on the sky)")
        print(WCS(hdr))

    hdr["NAXIS1"] = image.shape[0]
    hdr["NAXIS2"] = image.shape[1]

    #wcsprm = Wcsprm(hdr.tostring().encode('utf-8')) #everything else gave me errors with python 3, seemed to make problems with pc conversios, so i wwitched to the form below
    wcsprm = WCS(hdr).wcs
    wcsprm_original = WCS(hdr).wcs
    wcsprm, fov_radius, INCREASE_FOV_FLAG, PIXSCALE_UNCLEAR = read_additional_info_from_header(wcsprm, hdr, ra, dec,projection_ra, projection_dec, ignore_header_rot, radius)
    if(verbose):
        print(WCS(wcsprm.to_header()))
    coord = SkyCoord(wcsprm.crval[0], wcsprm.crval[1], unit=(u.deg, u.deg), frame="icrs")
    if(not PIXSCALE_UNCLEAR):
        if(wcsprm.crpix[0] < 0 or wcsprm.crpix[1] < 0 or wcsprm.crpix[0] > image.shape[0] or wcsprm.crpix[1] > image.shape[1] ):
            if(not silent):
                print("central value outside of the image, moving it to the center")
            coord_radec = wcsprm.p2s([[image.shape[0]/2, image.shape[1]/2]], 1)["world"][0]
            coord = SkyCoord(coord_radec[0], coord_radec[1], unit=(u.deg, u.deg), frame="icrs")
            #print(wcsprm)



    #better: put in nice wrapper! with repeated tries and maybe try synchron!
    if(not silent):
        print(">Dowloading catalog data")
    radius = u.Quantity(fov_radius, u.arcmin)#will prob need more
    catalog_data = query.get_data(coord, radius, catalog)
    report["catalog"] = catalog
    #reference = reference.query("mag <20")
    

    if(catalog == "GAIA" and catalog_data.shape[0] < 5):
        if(not silent):
            print("GAIA seems to not have enough objects, will enhance with PS1")
        catalog_data2 = query.get_data(coord, radius, "PS")
        report["catalog"] = "PS"
        catalog_data = pd.concat([catalog_data, catalog_data2])
        #apertures_catalog = CircularAperture(wcs.wcs_world2pix(catalog_data[["ra", "dec"]], 1), r=5.)
        if(not silent):
            print("Now we have a total of {} sources. Keep in mind that there might be duplicates now  since we combined 2 catalogs".format(catalog_data.shape[0]))
    elif(catalog == "PS" and (catalog_data is None or catalog_data.shape[0] < 5)):
        if(not silent):
            print("We seem to be outside the PS footprint, enhance with GAIA data")
        catalog_data2 = query.get_data(coord, radius, "GAIA")
        report["catalog"] = "GAIA"
        catalog_data = pd.concat([catalog_data, catalog_data2])
        #apertures_catalog = CircularAperture(wcs.wcs_world2pix(catalog_data[["ra", "dec"]], 1), r=5.)
        if(not silent):
            print("Now we have a total of {} sources. Keep in mind that there might be duplicates now since we combined 2 catalogs".format(catalog_data.shape[0]))

    max_sources = 400
    if(INCREASE_FOV_FLAG):
        max_sources= max_sources*2.25 #1.5 times the radius, so 2.25 the area
    if(catalog_data.shape[0]>max_sources):
        catalog_data = catalog_data.nsmallest(400, "mag")
    #remove duplicates in catalog?

    apertures_catalog = CircularAperture(wcsprm.s2p(catalog_data[["ra", "dec"]], 1)['pixcrd'], r=5.)
    #plotting what we have, I keep it in the detector field, world coordinates are more painfull to plot
    if(images):
        fig = plt.figure()
        fig.canvas.manager.set_window_title('Input for {}'.format(fits_image_filename))
        plt.xlabel("pixel x direction")
        plt.ylabel("pixel y direction")
        plt.title("Input - red: catalog sources, blue: detected sources in img")
        #plt.imshow(image,cmap='Greys', origin='lower', norm=LogNorm())
        plt.imshow(image,cmap='Greys', origin="lower", norm=LogNorm() , extent=[0.5, image.shape[1]+0.5, 0.5, image.shape[0]+0.5])
        apertures.plot(color='blue', lw=1.5, alpha=0.5)
        apertures_catalog.plot(color='red', lw=1.5, alpha=0.5)

        plt.xlim(-200,image.shape[0]+200)
        plt.ylim(-200,image.shape[1]+200)
        if(save_images):
            name_parts = fits_image_filename.rsplit('.', 1)
            plt.savefig(name_parts[0]+"_image_before.pdf")

    ###tranforming to match the sources
    if(not silent):
        print("---------------------------------")
        print(">Finding the transformation")
    if(rotation_scaling):
        if(not silent):
            print("Finding scaling and rotation")
        wcsprm = register.get_scaling_and_rotation(observation, catalog_data, wcsprm, scale_guessed=PIXSCALE_UNCLEAR, verbose=verbose)
    if(xy_transformation):
        if(not silent):
            print("Finding offset")
        wcsprm,_,_ = register.offset_with_orientation(observation, catalog_data, wcsprm, fast=False , INCREASE_FOV_FLAG=INCREASE_FOV_FLAG, verbose= verbose, silent=silent)

    #correct subpixel error
    compare_threshold = 3
    if(high_res):
        compare_threshold = 100
    obs_x, obs_y, cat_x, cat_y, distances = register.find_matches(observation, catalog_data, wcsprm, threshold=compare_threshold)#3
    if (len(distances) == 0): #meaning the list is empty
        best_score = 0
    else:
        rms = np.sqrt(np.mean(np.square(distances)))
        best_score = len(obs_x)/(rms+10) #start with current best score
    fine_transformation_success = False
    if(fine_transformation):
        print("Finding scaling and rotation")
        lis = [2,3,5,8,10,6,4, 20,2,1,0.5]
        if(high_res):
            lis = [200,300,100,150,80,40,70, 20, 100, 30,9,5]
        skip_rot_scale = True
        for i in lis:
            wcsprm_new, score = register.fine_transformation(observation, catalog_data, wcsprm, threshold=i, compare_threshold=compare_threshold, skip_rot_scale=skip_rot_scale)
            if(i == 20):
                #only allow rot and scaling for the last few tries
                skip_rot_scale = False
            if(score> best_score):
                print("fine transformation successfull - new best fit")
                wcsprm = wcsprm_new
                best_score = score
                fine_transformation_success = True
        if not fine_transformation_success:
            if(not silent):
                print("Fine transformation did not improve result so will be discarded.")
        else:
            if(not silent):
                print("Fine transformation applied to improve result")
    #register.calculate_rms(observation, catalog_data,wcs)

    #make wcsprim more physical by moving scaling to cdelt, out of the pc matrix
    wcs =WCS(wcsprm.to_header())
    if(verbose):
        print(wcs)
    from astropy.wcs import utils
    scales = utils.proj_plane_pixel_scales(wcs)
    #print(scales)
    cdelt = wcsprm.get_cdelt()
    #print(cdelt)
    scale_ratio = scales/cdelt
    #print(scale_ratio)
    pc = np.array(wcsprm.get_pc())
    pc[0,0] = pc[0,0]/scale_ratio[0]
    pc[1,0] = pc[1,0]/scale_ratio[1]
    pc[0,1] = pc[0,1]/scale_ratio[0]
    pc[1,1] = pc[1,1]/scale_ratio[1]
    wcsprm.pc = pc
    wcsprm.cdelt = scales

    #WCS difference before and after
    if(not silent):
        print("> Compared to the input the Wcs was changed by: ")
    scales_original = utils.proj_plane_pixel_scales(WCS(hdr))
    if(not silent):
        print("WCS got scaled by {} in x direction and {} in y direction".format(scales[0]/scales_original[0], scales[1]/scales_original[1]))
    #sources:
    #https://math.stackexchange.com/questions/2113634/comparing-two-rotation-matrices
    #https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python/13849249#13849249
    def unit_vector(vector):
        """ Returns the unit vector of the vector.  """
        return vector / max(np.linalg.norm(vector), 1e-10)
    def matrix_angle( B, A ):
        """ comment cos between vectors or matrices """
        Aflat = A.reshape(-1)
        Aflat = unit_vector(Aflat)
        Bflat = B.reshape(-1)
        Bflat = unit_vector(Bflat)
        #return np.arccos((np.dot( Aflat, Bflat ) / max( np.linalg.norm(Aflat) * np.linalg.norm(Bflat), 1e-10 )))
        return np.arccos(np.clip(np.dot(Aflat, Bflat), -1.0, 1.0))
    #print(matrix_angle(wcsprm.get_pc(), wcsprm_original.get_pc()) /2/np.pi*360)
    #bugfix: multiplying by cdelt otherwise the calculated angle is off by a tiny bit
    rotation_angle = matrix_angle(wcsprm.get_pc()@wcsprm.get_cdelt(), wcsprm_original.get_pc()@wcsprm_original.get_cdelt()) /2./np.pi*360.
    if((wcsprm.get_pc() @ wcsprm_original.get_pc() )[0,1] > 0):
        text = "counterclockwise"
    else:
        text = "clockwise"
    if(not silent):
        print("Rotation of WCS by an angle of {} deg ".format(rotation_angle)+text)
    old_central_pixel = wcsprm_original.s2p([wcsprm.crval], 1)["pixcrd"][0]
    if(not silent):
        print("x offset: {} px, y offset: {} px ".format(wcsprm.crpix[0]- old_central_pixel[0], wcsprm.crpix[1]- old_central_pixel[1]))

    if(not silent):
        print("new WCS:")
        print(WCS(wcsprm.to_header()))

    #check final figure
    if(images):
        print("Plotting final figue. Forcing origin to be (1,1) to be consistent with ds9 and fits file conventions.")
        fig = plt.figure()
        fig.canvas.manager.set_window_title('Result for {}'.format(fits_image_filename))
        plt.xlabel("pixel x direction")
        plt.ylabel("pixel y direction")
        plt.title("Result - red: catalog sources, blue: detected sources in img \n (Note origin = (1,1))")
        print(image.shape)
        plt.imshow(image,cmap='Greys', origin="lower", norm=LogNorm() , extent=[0.5, image.shape[1]+0.5, 0.5, image.shape[0]+0.5])
        apertures.plot(color='blue', lw=1.5, alpha=0.5)
        #apertures_catalog = CircularAperture(wcs.wcs_world2pix(catalog_data[["ra", "dec"]], 1), r=5.)
        apertures_catalog = CircularAperture(wcsprm.s2p(catalog_data[["ra", "dec"]], 1)['pixcrd'], r=5.)

        apertures_catalog.plot(color='red', lw=1.5, alpha=0.5)
        if(save_images):
            name_parts = fits_image_filename.rsplit('.', 1)
            plt.savefig(name_parts[0]+"_image_after.pdf")
    if(not silent):
        print("--- Evaluate how good the transformation is ----")
    dic_rms = register.calculate_rms(observation, catalog_data,wcsprm)
    #updating file
    converged = determine_if_fit_converged(dic_rms, catalog_data, observation, wcsprm, image.shape[0], image.shape[1], silent)
    report["converged"] = converged
    report["matches"] = dic_rms["matches"]
    report["match_radius"] = dic_rms["radius_px"]
    if(converged or save_bad_result):
        write_wcs_to_hdr(fits_image_filename, wcsprm, report, hdul_idx=hdul_idx)
        if(filename_for_sources != None):
            wcs =WCS(wcsprm.to_header())
            observation_on_sky = wcs.wcs_pix2world(observation[["xcenter","ycenter"]], 1)
            #catalog_from_obs = np.zeros(observation_on_sky.shape[0], dtype={'names':('ra', 'dec', 'aperture_sum'),'formats':('f8', 'f8', 'f8')})
            catalog_from_obs = pd.DataFrame()
            catalog_from_obs["ra"]= observation_on_sky[:,0]
            catalog_from_obs["dec"]= observation_on_sky[:,1]
            catalog_from_obs["aperture_sum"]= observation["aperture_sum"]
            catalog_from_obs["mag"]= -1.* observation["aperture_sum"]#this is fine since we only use the mag to order the sources!
            catalog_from_obs.to_csv(filename_for_sources+".csv")
    if(images):
        plt.show()

    return converged, dic_rms #dictionary with short info about fit, "matches" gives a number of objects matched within certain radius


if __name__ == '__main__':
    main()
