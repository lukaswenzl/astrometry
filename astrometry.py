"""This aims to do very fast and simple astrometric calibration of images.

Author Lukas Wenzl
written in python 3

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


def write_wcs_to_hdr(original_filename, wcs):
    """Update the header of the fits fileself.

    Parameters
    ----------
    original_filename : str
        Original filename of the fits file
    wcs : astropy.wcs
        World coordinate system object decsribing translation between image and skycoord

    """
    with fits.open(original_filename) as hdul:

        hdu = hdul[0]
        hdr_old = hdu.header

        #new_header_info = wcs.to_header()
        pc = wcs.wcs.pc
        hdr_old["PC1_1"] =pc[0][0]
        hdr_old["PC1_2"] = pc[0][1]
        hdr_old["PC2_1"]= pc[1][0]
        hdr_old["PC2_2"] = pc[1][1]
        #print(repr(wcs.to_header()))
        #print("-----------------------")
        hdr_old.update(wcs.to_header())#.keys(),wcs.to_header().values())
        repr(hdr_old)

        hdu.header = hdr_old
        hdul[0] = hdu
        #removing fits ending
        name_parts = original_filename.rsplit('.', 1)
        hdul.writeto(name_parts[0]+'_astro.fits', overwrite=True)
        print("file written.")








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
    parser.add_argument("input", nargs='+',help="Input Image with .fits ending.", type=str)#, default="sample_images/pso016p03_Jrot.fits")  #pso016p03_Jrot.fits")
    #parser.add_argument('-l','--list', nargs='+', help='<Required> Set flag', required=True)
    # Use like:
    # python arg.py -l 1234 2345 3456 4567
    parser.add_argument("-c", "--catalog", help="Catalog to use for position reference ('PS' or 'GAIA')", type=str, default="PS")

    parser.add_argument("-s", "--save_images", help="Set True to create _image_before.pdf and _image_after.pdf", type=bool, default=False)
    parser.add_argument("-p", "--show_images", help="Set False to not have the plots pop up.", type=bool, default=True)

    parser.add_argument("-w", "--ignore_warnings", help="Set False to see all Warnings about the header if there is problems. Default is to ignore most warnings.", type=bool, default=True)





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
                fits_image_filenames.append(file)
        print(fits_image_filenames)

    for fits_image_filename in fits_image_filenames:

        print("")
        print("--------------------------------------------------------------------------------------")
        print("---- Astrometry for {} ----".format(fits_image_filename))

        #header_keywords = [""]
        with fits.open(fits_image_filename) as hdul:
            #print(hdul.info())
            #print("if image is not at 0 the program will break later on")

            #print(hdul[0].header)

            hdu = hdul[0]
            #hdu.verify('fix')
            hdr = hdu.header


            image_or = hdul[0].data.astype(float)  #hdu.data[500:700, 500:700].astype(float)
            # plt.imshow(image)
            # plt.show()
            image = image_or - np.median(image_or)

        #pixel scale by hand for testing
        #hdr["CDELT1"] = -8.0000000000000E-5#-5!!!!
        #hdr["CDELT2"] = 8.0000000000003E-5#-5!!!
        # hdr["PC1_1"] =1.
        # hdr["PC1_2"] =0.
        # hdr["PC2_1"]=0.
        # hdr["PC2_2"] =1.



        observation = find_sources(image)

        positions = (observation['xcenter'], observation['ycenter'])
        apertures = CircularAperture(positions, r=4.)


        #world coordinates
        wcs = WCS(hdr)
        print("--- Info found in the file -- (CRVAl: position of central pixel (CRPIX) on the sky)")
        print(wcs)
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
        coord = SkyCoord(wcs.wcs.crval[0], wcs.wcs.crval[1], unit=(u.deg, u.deg), frame="icrs")

        #print(coord)


        #put in nice wrapper! with repeated tries and maybe try synchron!
        print("--- Dowloading catalog data -- ")
        radius = u.Quantity(4, u.arcmin)#will prob need more
        catalog_data = query.get_data(coord, radius, args.catalog)
        #reference = reference.query("mag <20")

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

        apertures_catalog = CircularAperture(wcs.wcs_world2pix(catalog_data[["ra", "dec"]], 1), r=5.)

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
        print("--- Finding the tranformation -- ")
        ###
        wcs,_,_ = register.offset_with_orientation(observation, catalog_data, hdr, fast=True)
        #wcs,_,_ = register.simple_offset(observation, catalog_data, wcs)
        #register.general_transformation(observation, catalog_data, wcs)

        #correct subpixel error
        register.fine_transformation(observation, catalog_data, wcs)
        #register.calculate_rms(observation, catalog_data,wcs)


        #check final figure
        fig = plt.figure()
        fig.canvas.set_window_title('Result for {}'.format(fits_image_filename))
        plt.xlabel("pixel x direction")
        plt.ylabel("pixel y direction")
        plt.title("Result - red: catalog sources, blue: detected sources in img")
        plt.imshow(image,cmap='Greys', origin='lower', norm=LogNorm())
        apertures.plot(color='blue', lw=1.5, alpha=0.5)
        apertures_catalog = CircularAperture(wcs.wcs_world2pix(catalog_data[["ra", "dec"]], 1), r=5.)
        apertures_catalog.plot(color='red', lw=1.5, alpha=0.5)
        if(args.save_images):
            name_parts = fits_image_filename.rsplit('.', 1)
            plt.savefig(name_parts[0]+"_image_after.pdf")

        print("--- Evaluate how goot the transformation is ----")
        register.calculate_rms(observation, catalog_data,wcs)


        #updating file
        write_wcs_to_hdr(fits_image_filename, wcs)


        print("overall time taken")
        print(datetime.now()-StartTime)
        if(args.show_images):
            plt.show()
    print("-- finished --")





if __name__ == '__main__':
    main()
