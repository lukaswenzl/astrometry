"""This aims to do very fast and simple astrometric calibration of images.

Author Lukas Wenzl
written in python 3

NOTES:
-not sure if reflected images will work
"""



#
#
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



def find_sources(image):
    """Find surces in the image. Uses DAOStarFinder with symmetric gaussian kernels. Only uses 5 sigma detections. It only gives the 200 brightest sources or less.

    Parameters
    ----------
    image
        Observed image (without background)

    Returns
    -------
    observaion : dataframe
        Pandas dataframe with the sources

    """
    #find sources
    #bkg_sigma = mad_std(image)
    mean, median, std = sigma_clipped_stats(image, sigma=3.0)
    #sigma = np.std(image)
    #print(bkg_sigma)
    #print(std)
    daofind = DAOStarFinder(fwhm=4., threshold=5.*std, brightest=200)
    sources = daofind(image)
    for col in sources.colnames:
        sources[col].info.format = '%.8g'  # for consistent table output
    #print(sources)

    positions = (sources['xcentroid'], sources['ycentroid'])
    apertures = CircularAperture(positions, r=4.)
    phot_table = aperture_photometry(image, apertures)
    for col in phot_table.colnames:
        phot_table[col].info.format = '%.8g'  # for consistent table output
    #print(phot_table)

    observation = Table(phot_table).to_pandas()

    #through out candidates where the star finder messed up
    observation = observation.query("aperture_sum > "+str(5*std))
    return  observation


def write_wcs_to_hdr(original_filename, wcsprm):
    """Update the header of the fits fileself.

    Parameters
    ----------
    original_filename : str
        Original filename of the fits file
    wcsprm : astropy.wcs.wcsprm
        World coordinate system object decsribing translation between image and skycoord

    """
    with fits.open(original_filename) as hdul:

        hdu = hdul[0]
        hdr_file = hdu.header

        #new_header_info = wcs.to_header()

        wcs =WCS(wcsprm.to_header())

        #I will through out CD which contains the scaling and separate into pc and Cdelt
        for old_parameter in ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', "PC1_1", "PC1_2", "PC2_1", "PC2_2"]:
            if (old_parameter in hdr_file):
                del hdr_file[old_parameter]

        #if new file doe not have off diagonal px values or unity delete the old ones
        #pc = wcsprm.get_pc()
        #if (pc[0, 0] ==1 and "PC1_1" in hdr_file):
        #    del hdr_file["PC1_1"]
        #if (pc[1, 1] ==1 and "PC2_2" in hdr_file):
        #    del hdr_file["PC1_1"]
        #if ( (pc[0, 1] ==0 or pc[1,0] == 0):
        #    if "PC1_1" in hdr_file):
        #    del hdr_file["PC1_1"]
        #if (pc[0, 0] ==1 and "PC1_1" in hdr_file):
        #    del hdr_file["PC1_1"]
        hdr_file.update(wcs.to_header())

        #print(repr(wcs.to_header()))
        #print("-----------------------")
        #hdr_old.update(wcs.to_header())#.keys(),wcs.to_header().values())
        repr(hdr_file)

        hdu.header = hdr_file
        hdul[0] = hdu
        #removing fits ending
        name_parts = original_filename.rsplit('.', 1)
        hdul.writeto(name_parts[0]+'_astro.fits', overwrite=True)
        print("file written.")



def read_additional_info_from_header(wcsprm, hdr, RA_input=None, DEC_input=None, projection_ra=None, projection_dec=None):
    fov_radius = 4 #arcmin radius to include field of view
    INCREASE_FOV_FLAG = False # increase the field to view by 50% to search in catalog if position on sky is inaccurate
    PIXSCALE_UNCLEAR = False

    keywords_check = ["PIXSCALE", "NAXIS1", "NAXIS2", "RA", "DEC"] #list of possible keywords the scs parser might miss
    keywords_present = [] # list of keywords that are actually present
    for i in keywords_check:
        if(i in hdr.keys()):
            keywords_present.append(i)
    # print("Additional keywords found by header:")
    # print(keywords_present)

    if("NAXIS1" not in keywords_present or "NAXIS2" not in keywords_present ):
        print("ERROR: NAXIS1 or NAXIS2 missing in file. Please add!")
    else:
        axis1 = hdr["NAXIS1"]
        axis2 = hdr["NAXIS2"]
    #if(wcsprm.isunity()):

    pc = wcsprm.get_pc()
    cdelt = wcsprm.get_cdelt()
    wcs_pixscale = (pc @ cdelt )
    if((np.abs(wcs_pixscale[0])) < 1e-7 or (np.abs(wcs_pixscale[1])) < 1e-7 or
        (np.abs(wcs_pixscale[0])) > 5e-3 or (np.abs(wcs_pixscale[1])) > 5e-3):
        print("pixelscale is completely unrealistic. Will guess")
        print(wcs_pixscale)
        guess = 8.43785734e-05
        wcsprm.pc = [[1,0],[0,1]]
        wcsprm.cdelt = [guess, guess]
        print("Changed pixelscale to {:.3g} deg/arcsec".format(guess))
        PIXSCALE_UNCLEAR = True
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
                print("changed pixelscale to {:.3g} deg/arcsec".format(pixscale))
                fov_radius = (x_size/2+y_size/2)/np.sqrt(2) #try to get corners
                PIXSCALE_UNCLEAR=True


    if(np.array_equal(wcsprm.crpix, [0,0])):
        #seems to not be in header, better set in middle
        wcsprm.crpix = [axis1/2, axis2/2]

    if(np.array_equal(wcsprm.crval, [0,0] )): ###what to do??
        INCREASE_FOV_FLAG = True
        if ("RA" in keywords_present and "DEC" in keywords_present): ##carefull degree and hourangle!!!
            wcsprm.crval = [hdr["RA"], hdr["DEC"]]
            print("Found ra and dec information in the header")
            print(wcsprm.crval)
            print("Is this position within the field of view in degrees? otherwise it will not work. In that case give a more accurate position as an argument: -ra XX -dec XX both in degrees")

    if (RA_input is not None): #use user input if provided
        wcsprm.crval = [RA_input, wcsprm.crval[1]]
    if (DEC_input is not None):
        wcsprm.crval = [wcsprm.crval[0], DEC_input]

    if(np.array_equal(wcsprm.crval, [0,0] )):
        print(">>>>>>>>>WARNING")
        print("No rough sky position was found for this object. Please add as -ra XX -dex XX both in degress. Adding the position as keywords in the fits file header will also work. The keywords are RA and DEC. The program expects the values in degrees. ")
    #wcsprm.crval = [172.3556987, 18.77342494] #deg  , deg
    #wcsprm.pc = [[-1,0],[0,1]]
    #wcsprm.crpix = [wcsprm.crpix[0]+22, wcsprm.crpix[1]+70]
    #
    if(np.array_equal(wcsprm.ctype, ["",""])):
        INCREASE_FOV_FLAG = True
        if(projection_ra is not None and projection_dec is not None):
            wcsprm.ctype = [ projection_ra,  projection_dec]
            print("reached")
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
    parser.add_argument("input", nargs='+',help="Input Image with .fits ending. A folder or multiple Images also work", type=str)#, default="sample_images/pso016p03_Jrot.fits")  #pso016p03_Jrot.fits")
    #parser.add_argument('-l','--list', nargs='+', help='<Required> Set flag', required=True)
    # Use like:
    # python arg.py -l 1234 2345 3456 4567
    parser.add_argument("-c", "--catalog", help="Catalog to use for position reference ('PS' or 'GAIA')", type=str, default="PS")

    parser.add_argument("-s", "--save_images", help="Set True to create _image_before.pdf and _image_after.pdf", type=bool, default=False)
    parser.add_argument("-p", "--show_images", help="Set False to not have the plots pop up.", type=bool, default=True)

    parser.add_argument("-v", "--verbose", help="More console output about what is happening. Helpfull for debugging.", type=bool, default=False)


    parser.add_argument("-w", "--ignore_warnings", help="Set False to see all Warnings about the header if there is problems. Default is to ignore most warnings.", type=bool, default=True)

    parser.add_argument("-ra", "--ra", help="Set ra by hand in degrees. Should not be necessary if the info is in the fits header.", type=float, default=None)
    parser.add_argument("-dec", "--dec", help="Set ra by hand in degrees. Should not be necessary if the info is in the fits header.", type=float, default=None)

    parser.add_argument("-proj1", "--projection_ra", help="Set the projection by hand. Should not be necessary if the info is in the fits header. example: -proj1 RA---TAN -proj2 DEC--TAN", type=str, default=None)
    parser.add_argument("-proj2", "--projection_dec", help="Set the projection by hand. Should not be necessary if the info is in the fits header. example: -proj1 RA---TAN -proj2 DEC--TAN", type=str, default=None)







    # Print version
    parser.add_argument("--version", action="version", version='%(prog)s - Version 0.1') #
    #changelog
    #version 0.0 proof of concept
    #version 0.1 alpha version


    # Parse arguments
    args = parser.parse_args()

    return args


def main():
    """Perform astrometry for the given file."""
    print("Program version: 0.1")
    StartTime = datetime.now()
    args = parseArguments()

    if(args.show_images):
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
    #print(fits_image_filenames)

    #for directories search for appropriate fits files
    if(os.path.isdir(fits_image_filenames[0])):
        print("detected a directory. Will search for fits files in it")
        path = fits_image_filenames[0]
        fits_image_filenames = []
        for file in os.listdir(path):
            if file.endswith(".fits") and "_astro" not in file:
                fits_image_filenames.append(path+"/"+file)
        print(fits_image_filenames)

    for fits_image_filename in fits_image_filenames:
        print("")
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        print("> Astrometry for {} ".format(fits_image_filename))

        with fits.open(fits_image_filename) as hdul:
            #print(hdul.info())
            if(args.verbose):
                print("if image is not at first position in the fits file the program will break later on")
            #print(hdul[0].header)

            hdu = hdul[0]
            #hdu.verify('fix')
            hdr = hdu.header


            image_or = hdul[0].data.astype(float)
            image = image_or - np.median(image_or)

        #pixel scale by hand for testing
        # hdr["CDELT1"] = -6.0000000000000E-5#-5!!!!
        # hdr["CDELT2"] = 6.0000000000003E-5#-5!!!
        # hdr["PC1_1"] =0.98006658
        # hdr["PC1_2"] =0.19866933
        # hdr["PC2_1"]=-0.19866933
        # hdr["PC2_2"] =0.98006658


        observation = find_sources(image)
        #print(observation)

        positions = (observation['xcenter'], observation['ycenter'])
        apertures = CircularAperture(positions, r=4.)


        #world coordinates
        #wcs = WCS(hdr)
        print(">Info found in the file -- (CRVAl: position of central pixel (CRPIX) on the sky)")
        print(WCS(hdr))
        #print(wcs)
        #




        #header = hdr)
        wcsprm = Wcsprm(hdr.tostring().encode('utf-8')) #everything else gave me errors with python 3
        #print(wcsprm.get_pc())
        wcsprm, fov_radius, INCREASE_FOV_FLAG, PIXSCALE_UNCLEAR = read_additional_info_from_header(wcsprm, hdr, args.ra, args.dec, args.projection_ra, args.projection_dec)
        if(args.verbose):
            print(WCS(wcsprm.to_header()))
        #wcsprm.pc = [[2, 0],[0,1]]

        #print(wcsprm.set())
        #print(wcsprm.get_pc())
        #pc = wcsprm.get_pc()
        #print(np.linalg.det(pc))
        #print(wcsprm.get_cdelt())
        #wcs.fix()
        #print(wcsprm.print_contents())
        #print(repr(hdr.update(wcsprm.to_header().encode('utf-8')))) #not working

        #hdu.verify("fix")
        #print(repr(hdr))
        #wcs.wcs_pix2world(pixcrd, 1)
        #wcs.wcs_world2pix(world, 1)
        #wcs.wcs.crpix = [-234.75, 8.3393]
        # wcs.wcs.cdelt = np.array([-0.066667, 0.066667])
        # wcs.wcs.crval = [0, -90]
        # wcs.wcs.ctype = ["RA---AIR", "DEC--AIR"]
        # wcs.wcs.set_pv([(2, 1, 45.0)])
        # For historical compatibility, three alternate specifications of the linear transformations
        # are available in wcslib. The canonical PCi_ja with CDELTia, CDi_ja, and the deprecated CROTAia
        # keywords. Although the latter may not formally co-exist with PCi_ja,
        # the approach here is simply to ignore them if given in conjunction with PCi_ja.
        # has_pc, has_cd and has_crota can be used to determine which of these alternatives are present in the header.
        # These alternate specifications of the linear transformation matrix are translated immediately to PCi_ja by set
        # and are nowhere visible to the lower-level routines. In particular, set resets cdelt to unity if CDi_ja is present
        # (and no PCi_ja). If no CROTAia is associated with the latitude axis, set reverts to a unity PCi_ja matrix.

        #ax = plt.subplot(projection=wcs)
        #fig = plt.figure()
        #ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs)

        # ax = plt.subplot(projection=wcs)
        # plt.imshow(image, cmap='Greys', origin='lower', norm=LogNorm())
        # apertures.plot(color='blue', lw=1.5, alpha=0.5)
        # ra = ax.coords[0]
        # dec = ax.coords[1]




        #get rough coordinates
        #print(hdr["RA"])
        #coord = SkyCoord(hdr["RA"], hdr["DEC"], unit=(u.hourangle, u.deg), frame="icrs")
        coord = SkyCoord(wcsprm.crval[0], wcsprm.crval[1], unit=(u.deg, u.deg), frame="icrs")


        #put in nice wrapper! with repeated tries and maybe try synchron!
        print(">Dowloading catalog data")
        radius = u.Quantity(fov_radius, u.arcmin)#will prob need more
        catalog_data = query.get_data(coord, radius, args.catalog)
        #reference = reference.query("mag <20")
        max_sources = 500
        if(INCREASE_FOV_FLAG):
            max_sources= max_sources*2.25 #1.5 times the radius, so 2.25 the area
        if(catalog_data.shape[0]>max_sources):
            catalog_data = catalog_data.nsmallest(400, "mag")

        if(args.catalog == "GAIA" and catalog_data.shape[0] < 5):
            print("GAIA seems to not have enough objects, will enhance with PS1")
            catalog_data2 = query.get_data(coord, radius, "PS")
            catalog_data = pd.concat([catalog_data, catalog_data2])
            #apertures_catalog = CircularAperture(wcs.wcs_world2pix(catalog_data[["ra", "dec"]], 1), r=5.)
            print("Now we have a total of {} sources. Keep in mind that there might be duplicates now  since we combined 2 catalogs".format(catalog_data.shape[0]))
        elif(args.catalog == "PS" and (catalog_data is None or catalog_data.shape[0] < 5)):
            print("We seem to be outside the PS footprint, enhance with GAIA data")
            catalog_data2 = query.get_data(coord, radius, "GAIA")
            catalog_data = pd.concat([catalog_data, catalog_data2])
            #apertures_catalog = CircularAperture(wcs.wcs_world2pix(catalog_data[["ra", "dec"]], 1), r=5.)
            print("Now we have a total of {} sources. Keep in mind that there might be duplicates now since we combined 2 catalogs".format(catalog_data.shape[0]))

        #remove duplicates in catalog


        #apertures_catalog = CircularAperture(wcs.wcs_world2pix(catalog_data[["ra", "dec"]], 1), r=5.)
        apertures_catalog = CircularAperture(wcsprm.s2p(catalog_data[["ra", "dec"]], 1)['pixcrd'], r=5.)


        #plotting what we have, I keep it in the detector field, world coordinates are a pain to plot
        fig = plt.figure()
        fig.canvas.set_window_title('Input for {}'.format(fits_image_filename))
        plt.xlabel("pixel x direction")
        plt.ylabel("pixel y direction")
        plt.title("Input - red: catalog sources, blue: detected sources in img")
        #ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs)
        #plt.plot(r["RAJ2000"], r["DEJ2000"], "x")
        #plt.xlabel("X")
        plt.imshow(image,cmap='Greys', origin='lower', norm=LogNorm())
        apertures.plot(color='blue', lw=1.5, alpha=0.5)
        apertures_catalog.plot(color='red', lw=1.5, alpha=0.5)
        #apertures_GAIA.plot(color='black', lw=1.5, alpha=0.5)

        plt.xlim(-200,image.shape[0]+200)
        plt.ylim(-200,image.shape[1]+200)
        if(args.save_images):
            name_parts = fits_image_filename.rsplit('.', 1)
            plt.savefig(name_parts[0]+"_image_before.pdf")


        ###tranforming to match the sources
        print("---------------------------------")
        print(">Finding the transformation")
        print("Finding scaling and rotation")
        wcsprm = register.get_scaling_and_rotation(observation, catalog_data, wcsprm, scale_guessed=PIXSCALE_UNCLEAR, verbose=args.verbose)
        print("Finding offset")
        wcsprm,_,_ = register.offset_with_orientation(observation, catalog_data, wcsprm, fast=False , INCREASE_FOV_FLAG=INCREASE_FOV_FLAG, verbose= args.verbose)
        #wcs,_,_ = register.simple_offset(observation, catalog_data, wcs)
        #register.general_transformation(observation, catalog_data, wcs)

        #correct subpixel error
        #register.calculate_rms(observation, catalog_data,wcsprm)
        obs_x, obs_y, cat_x, cat_y, distances = register.find_matches(observation, catalog_data, wcsprm, threshold=3)
        rms = np.sqrt(np.mean(np.square(distances)))
        best_score = len(obs_x)/(rms+1) #start with current best score
        fine_transformation = False
        for i in [2,3,5,8,10,6,4,2,1,0.5]:
            wcsprm_new, score = register.fine_transformation(observation, catalog_data, wcsprm, threshold=i)
            if(score> best_score):
                wcsprm = wcsprm_new
                best_score = score
                fine_transformation = True
        if not fine_transformation:
            print("Fine transformation did not imporve result so will be discarded.")
        #register.calculate_rms(observation, catalog_data,wcs)


        #check final figure
        fig = plt.figure()
        fig.canvas.set_window_title('Result for {}'.format(fits_image_filename))
        plt.xlabel("pixel x direction")
        plt.ylabel("pixel y direction")
        plt.title("Result - red: catalog sources, blue: detected sources in img")
        plt.imshow(image,cmap='Greys', origin='lower', norm=LogNorm())
        apertures.plot(color='blue', lw=1.5, alpha=0.5)
        #apertures_catalog = CircularAperture(wcs.wcs_world2pix(catalog_data[["ra", "dec"]], 1), r=5.)
        apertures_catalog = CircularAperture(wcsprm.s2p(catalog_data[["ra", "dec"]], 1)['pixcrd'], r=5.)

        apertures_catalog.plot(color='red', lw=1.5, alpha=0.5)
        if(args.save_images):
            name_parts = fits_image_filename.rsplit('.', 1)
            plt.savefig(name_parts[0]+"_image_after.pdf")

        print("--- Evaluate how good the transformation is ----")
        register.calculate_rms(observation, catalog_data,wcsprm)


        #updating file
        write_wcs_to_hdr(fits_image_filename, wcsprm)


        print("overall time taken")
        print(datetime.now()-StartTime)
        if(args.show_images):
            plt.show()
    print("-- finished --")





if __name__ == '__main__':
    main()
