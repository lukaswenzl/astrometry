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
from astropy.wcs import Wcsprm


from astropy.table import Table

import warnings
import os

"""
EXAMPLES:
pso200m10_i_astro.fits -ra 200.38292412 -dec 10.30514953 -b i

photometry r_eduardo_3_pisco_J_astro.fits -b J -name eduardo_3
"""



def find_sources(image, aperture):
    """Find surces in the image. Uses DAOStarFinder with symmetric gaussian kernels. Only uses 5 sigma detections. It only gives the 200 brightest sources or less.

    Parameters
    ----------
    image
        Observed image (without background)
    aperture : float
        aperture in pixel


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

    #changed order of positions to [(x,y), (x,y),...] for compatibility with photutils 1.4
    xcenters = np.array(sources['xcentroid'])
    ycenters = np.array(sources['ycentroid'])
    positions = [(xcenters[i], ycenters[i]) for i in range(len(xcenters))]
    apertures = CircularAperture(positions, r=aperture)
    phot_table = aperture_photometry(image, apertures)
    for col in phot_table.colnames:
        phot_table[col].info.format = '%.8g'  # for consistent table output
    #print(phot_table)

    observation = Table(phot_table).to_pandas()

    #through out candidates where the star finder messed up
    observation = observation.query("aperture_sum > "+str(5*std))
    return  observation




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
    parser.add_argument("-b", "--band", help="Photometric band to use. Options are g,r,i,z,y,J,K,H", type=str, default="z")
    parser.add_argument("-c", "--catalog", help="Catalog to use for photometric reference ('PS'), this can be left to  the standard 'auto' to let the program choose automaticaly", type=str, default="auto")

    parser.add_argument("-ra", "--ra", help="Set ra of the object targeted for observation", type=float, default=None)
    parser.add_argument("-dec", "--dec", help="Set dec of the object targeted for observation", type=float, default=None)
    parser.add_argument("-name", "--name", help="Provide name of object with coordinates given in targets.csv (comma separated, columns NAME, RA, DEC, both in degrees)", type=str, default=None)


    parser.add_argument("-v", "--verbose", help="More console output about what is happening. Helpfull for debugging.", type=bool, default=False)


    parser.add_argument("-w", "--ignore_warnings", help="Set False to see all Warnings about the header if there is problems. Default is to ignore most warnings.", type=bool, default=True)


    parser.add_argument("-a", "--aperture", help="Aperture in arcsec to use for extraction. Default 2 arcseconds.", type=float, default="2")





    # Print version
    parser.add_argument("--version", action="version", version='%(prog)s - Version 1.5') #
    #changelog
    #version 0.0 proof of concept
    #version 0.1 alpha version
    #version 0.2 bugfixes, 5sigma magnitue now more consistent and panstarrs data handling is more robust
    #version 1.5 added compatibility with photutils v1.4 and bumped version number to match other files


    # Parse arguments
    args = parser.parse_args()

    return args


def main():
    """Perform photometry for the given file."""
    print("Program version: 1.5")
    StartTime = datetime.now()
    args = parseArguments()


    if(args.ignore_warnings):
        warnings.simplefilter('ignore', UserWarning)


    fits_image_filenames = args.input
    #print(fits_image_filenames)

    #for directories search for appropriate fits files
    if(os.path.isdir(fits_image_filenames[0])):
        print("detected a directory. Will search for fits files in it that already have astrometry calibration and therefor contain _astro")
        path = fits_image_filenames[0]
        fits_image_filenames = []
        for file in os.listdir(path):
            if file.endswith(".fits") and "_astro" in file:
                fits_image_filenames.append(path+"/"+file)
        print(fits_image_filenames)

    for fits_image_filename in fits_image_filenames:
        print("")
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
        print("> Photometry for {} ".format(fits_image_filename))

        with fits.open(fits_image_filename) as hdul:
            #print(hdul.info())
            if(args.verbose):
                print("if image is not at first position in the fits file the program will break later on")
            #print(hdul[0].header)

            hdu = hdul[0]
            hdr = hdu.header


            image_or = hdul[0].data.astype(float)
            image = image_or - np.median(image_or)



        wcsprm = Wcsprm(hdr.tostring().encode('utf-8')) #everything else gave me errors with python 3

        #tranlating aperture into pixel:
        on_sky = wcsprm.p2s([[0,0],[1,1]], 1)["world"]
        px_scale = np.sqrt((on_sky[0,0]-on_sky[1,0])**2+(on_sky[0,1]-on_sky[1,1])**2)
        px_scale = px_scale*60*60 #in arcsec
        aperture = args.aperture / px_scale  #aperture n pixel
        print("aperture")
        observation = find_sources(image, aperture)
        #print(observation)

        #changed order of positions to [(x,y), (x,y),...] for compatibility with photutils 1.4
        xcenters = np.array(observation['xcenter'])
        ycenters = np.array(observation['xcenter'])
        positions = [(xcenters[i], ycenters[i]) for i in range(len(xcenters))]
        apertures = CircularAperture(positions, r=4.)


        #get rough coordinates
        coord = SkyCoord(wcsprm.crval[0], wcsprm.crval[1], unit=(u.deg, u.deg), frame="icrs")


        #put in nice wrapper! with repeated tries and maybe try synchron!
        print(">Dowloading catalog data")
        #WCS.calc_footprint(header=None, undistort=True, axes=None, center=True)
        radius = u.Quantity(5, u.arcmin)#should be enough for all images
        catalog_data, band_name, catalog_name, mag_sys = query.get_photometry_data(coord, radius, args.band, args.catalog)

        #throwing out blended sources (should be improved, TODO)


        apertures_catalog = CircularAperture(wcsprm.s2p(catalog_data[["ra", "dec"]], 1)['pixcrd'], r=5.)

        obs_matched, cat_matched, distances = register.find_matches_keep_catalog_info(observation, catalog_data, wcsprm, threshold=3)
        print("Found {} matches".format(obs_matched.shape[0]))
        MAG_CALC = True
        if(obs_matched.shape[0] == 0):
            MAG_CALC = False

        obs_matched["aperture_sum"]
        cat_matched[band_name]

        #mag = -2.5 log10(cts) + ZP
        ZP = -1*(-2.5 *np.log10(obs_matched["aperture_sum"].values) - cat_matched[band_name].values)

        #FIGURE OUT LINEAR RANGE TOTO
        ZP_median = np.median(ZP[~np.isnan(ZP)])
        new_magnitudes = -2.5 *np.log10(observation["aperture_sum"].values) + ZP_median




        ra = args.ra
        dec = args.dec
        if(args.name):
            targets = pd.read_csv("targets.csv")

            if ((targets["NAME"] == args.name).sum()):
                ra = targets.loc[targets["NAME"]==args.name, "RA"].values[0]
                dec = targets.loc[targets["NAME"]==args.name, "DEC"].values[0]
                print("For {} the following coordinates where found".format(args.name))
                print("ra: {} dec: {}".format(ra,dec))
            else:
                print("{} not found".format(args.name))
        if(ra and dec):
            c_unknown = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame="icrs")
            #print(observation[["xcenter", "ycenter"]].values)
            coordinats_obs = wcsprm.p2s(observation[["xcenter", "ycenter"]].values,1)["world"]
            #print(coordinats_obs[17])
            c_obs = SkyCoord(coordinats_obs[:,0], coordinats_obs[:,1], unit=(u.deg, u.deg), frame="icrs")
            #print(coordinats_obs)
            from astropy.coordinates import match_coordinates_sky
            idx, d2d, d3d = match_coordinates_sky(c_unknown, c_obs)
            cand_dups = d2d < 5*u.arcsec

            mean, median, std = sigma_clipped_stats(image, sigma=3.0)
            ap_area= CircularAperture((2,2), r=aperture)
            sig3_limiting_mag = -2.5*np.log10(3*std*np.sqrt(ap_area.area)) + ZP_median
            sig5_limiting_mag = -2.5*np.log10(5*std*np.sqrt(ap_area.area)) + ZP_median

            #######################
            #estimating the error via random apertures
            # print(std*np.sqrt(aperture_obj.area()))
            ran_x = 2*aperture + np.random.random(1000) * (image.shape[0]-4*aperture)
            ran_y = 2*aperture + np.random.random(1000) * (image.shape[0]-4*aperture)
            ran_pos = [(ran_x[i], ran_y[i]) for i in range(len(ran_x))]
            apertures_error = CircularAperture(ran_pos, r=aperture)
            phot_table = aperture_photometry(image, apertures_error)
            phot_table['aperture_sum'].info.format = '%.8g'  # for consistent table output
            random_extractions = Table(phot_table).to_pandas()
            #print(np.mean(np.abs(random_extractions["aperture_sum"])))
            # print(np.median(np.abs(random_extractions["aperture_sum"])))
            # sig5_limiting_mag_apertures = -2.5*np.log10(5*np.median(np.abs(random_extractions["aperture_sum"]))) + ZP_median
            # print(sig5_limiting_mag_apertures)
            m,_, std_apertures = sigma_clipped_stats(random_extractions["aperture_sum"], sigma=3)
            print(std_apertures)
            sig5_limiting_mag_apertures = -2.5*np.log10(5*std_apertures) + ZP_median
            # print(sig5_limiting_mag_apertures)
            # print(np.mean(random_extractions.loc[random_extractions["aperture_sum"]<0,"aperture_sum"]))
            ###########################
            mag = 0
            if(cand_dups.sum() > 0):
                print("Position was given. Found {} sources in a 5 arcsec radius. Here is the magnitude:".format(cand_dups.sum()))
                print(new_magnitudes[idx])
                print("----")


                aperture_obj = CircularAperture(((observation["xcenter"].values)[idx], (observation["ycenter"].values)[idx]), r=aperture)
                #noise = (observation["aperture_sum"].values)[idx] /(std*np.sqrt(aperture_obj.area()))
                #https://en.wikipedia.org/wiki/Sum_of_normally_distributed_random_variables
                #the std**2 is added, so we get a square root for the area!
                noise = (observation["aperture_sum"].values)[idx] /std_apertures


                #print(noise_prob_wrong)
                #print(noise_prob_wrong2)

                print("We get a signal to noise of {} for the fixed aperture".format(noise))

                mag = new_magnitudes[idx]
                mag_err = 2.5 * np.log10(1+ 1/noise)
                text = "found detection in {} band,\n {:.4g} +- {:.2g} {}mag, {:.3g} S/N \n 5 sig limiting mag is {:.4g}".format(args.band, mag, mag_err, mag_sys, noise, sig5_limiting_mag_apertures)

            else:
                print("Position was given. No object was found at that position.")
                pix_unknown = wcsprm.s2p([[ra, dec]], 1)
                pix_unknown  = pix_unknown['pixcrd']

                aperture_obj = CircularAperture(pix_unknown, r=aperture)
                phot_table = aperture_photometry(image, aperture_obj)
                phot_table['aperture_sum'].info.format = '%.8g'  # for consistent table output
                #print(phot_table)
                forced_phot = Table(phot_table).to_pandas()

                forced_mag = -2.5 *np.log10(forced_phot["aperture_sum"].values) + ZP_median
                forced_mag = forced_mag[0]

                #noise = (forced_phot["aperture_sum"].values)[0] /(std*np.sqrt(aperture_obj.area()))
                noise = (forced_phot["aperture_sum"].values)[0] /std_apertures


                print("Forced photometry gives a magnitude of {}".format(forced_mag))


                #local error
                # pix_unknown = pix_unknown[0]
                # mean_loc, median_loc, std_loc = sigma_clipped_stats(image[int(pix_unknown[0])-50:int(pix_unknown[0])+50,int(pix_unknown[1])-50:int(pix_unknown[1])+50], sigma=3.0)
                # ap_area= CircularAperture((2,2), r=aperture)
                # sig5_limiting_mag_local = -2.5*np.log10(5*std_loc*np.sqrt(ap_area.area())) + ZP_median
                # print((std_loc*np.sqrt(aperture_obj.area())))
                # print(sig5_limiting_mag_local)


                print("We get a signal to noise of {} for the fixed aperture".format(noise))
                mag = forced_mag
                mag_err = 2.5 * np.log10(1+ 1/noise)
                text = "forced photometry in {} band,\n {:.4g} +- {:.2g} {}mag, {:.3g} S/N \n 5 sig limiting mag is {:.4g} ".format(args.band, mag, mag_err,mag_sys, noise, sig5_limiting_mag_apertures)


        else:
            print(new_magnitudes)
            mag = 0

        obs_with_photometry = observation.copy()
        obs_with_photometry["mag"] = new_magnitudes

        #search for targeted object and print out its magnitude TODO

        if(ra and dec):
            plt.figure(figsize=(15,8))
            plt.subplot(1,2,1)
            pix_unknown = wcsprm.s2p([[ra, dec]], 1)
            pix_unknown  = pix_unknown['pixcrd']
            pix_unknown = pix_unknown[0]

            plt.xlabel("pixel x direction")
            plt.ylabel("pixel y direction")
            plt.title(text, size=20)
            #plt.imshow(image[int( pix_unknown[0])-30: int(pix_unknown[0])+30, int(pix_unknown[1])-30:int(pix_unknown[1])+30],cmap='Greys', origin='lower', norm=LogNorm())
            #plt.plot(30, 30, "+", color="red", markersize=20)
            plt.imshow(image,cmap='Greys', origin='lower', norm=LogNorm())
            plt.xlim(int( pix_unknown[0])-20, int( pix_unknown[0])+20)
            plt.ylim(int( pix_unknown[1])-20, int( pix_unknown[1])+20)
            aperture_obj.plot(color='blue', lw=1.5, alpha=0.5, label=str(args.aperture)+" arcsec aperture")
            plt.plot(pix_unknown[0],pix_unknown[1], "+", color="red", markersize=20, linewidth=10, label="predicted target position")
            plt.legend(bbox_to_anchor=(1.1, -0.1), ncol=2)

            plt.subplot(1,2,2)
            plt.title("Photometric calibration with {} sources from {}".format(obs_matched.shape[0], catalog_name))
            plt.plot(cat_matched[band_name].values,ZP, "o", markersize=20, label="catalog objects")
            if(MAG_CALC):
                plt.ylim(np.nanmin(ZP)-1, np.nanmax(ZP)+1)
            else:
                plt.ylim(10,30)
            if(mag):
                plt.axvline(x=mag, linewidth=4, color="green", label="target")
            plt.xlabel("catalog "+ band_name+" magnitude")
            plt.ylabel("Zeropoint (should be constant in linear part of the detector)")
            plt.legend()
            outputname = fits_image_filename.replace('.fits','')
            plt.savefig(outputname+"_ra{}dec{}.pdf".format(ra, dec))
            plt.show()

        else:

            plt.figure()
            plt.plot(cat_matched[band_name].values,ZP, "o", markersize=20, label="catalog objects")
            plt.ylim(np.nanmin(ZP)-1, np.nanmax(ZP)+1)
            if(mag):
                plt.axvline(x=mag, linewidth=4, color="green", label="target")
            plt.xlabel("catalog "+ band_name+" magnitude")
            plt.ylabel("Zeropoint (should be constant in linear part of the detector)")
            plt.legend()
            plt.show()


        print("overall time taken")
        print(datetime.now()-StartTime)
        # if(args.show_images):
        #     plt.show()
    print("-- finished --")


if __name__ == '__main__':
    main()
