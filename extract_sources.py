"""This aims extract the sources from one image to calibrate another (possible smaller one with no catalog sources)

Author Lukas Wenzl
written in python 3

"""

from astropy.wcs import WCS
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from datetime import datetime
# #for parsing the arguments for the file
from argparse import ArgumentParser
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
from matplotlib.colors import LogNorm

import astrometry

#import photutils
from photutils import DAOStarFinder
from photutils import aperture_photometry, CircularAperture

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
    parser.add_argument("source_image", help="Source image for the calibration. The program extracts and saves sources from this image. ", type=str)
    parser.add_argument("filename_for_sources", help="Filename for the sources to save to (without ending)", type=str)

    parser.add_argument("-t", "--target_image", help="Target image for the calibration (so only roughly the area of the image is used for source extraction)", type=str, default = None)

    parser.add_argument("-target_hdul_idx", "--target_hdul_idx", help="Index of the target image data in the fits file. Default 0", type=int, default=0)
    parser.add_argument("-source_hdul_idx", "--source_hdul_idx", help="Index of the source image data in the fits file. Default 0", type=int, default=0)

    parser.add_argument("-det_thres", "--detection_threshold", help="Detection threshold for sources in sigma above mean. Default 5 sigma", type=float, default=5)
    parser.add_argument("-padding", "--padding", help="Padding around target image for extraction. Default 50px", type=int, default=50)
    parser.add_argument("-cutout", "--cutout", help="Cutout bad recangle from Image. Specify corners in pixel as -cutout xstart xend ystart yend. You can give multiple cutouts ", nargs='+',action="append", type=int)

    parser.add_argument("-v", "--verbose", help="More console output about what is happening. Helpfull for debugging.", type=bool, default=False)

    # Print version
    parser.add_argument("--version", action="version", version='%(prog)s - Version 1.0') #
    #changelog
    #version 1.0


    # Parse arguments
    args = parser.parse_args()

    return args


def main():
    """Perform photometry for the given file."""
    print("Program version: 1.5")
    StartTime = datetime.now()
    args = parseArguments()


    with fits.open(args.source_image) as hdul:
        hdu = hdul[args.source_hdul_idx]
        hdr = hdu.header

        image_or = hdul[args.source_hdul_idx].data.astype(float)
        median = np.nanmedian(image_or)
        image_or[np.isnan(image_or)]=median
        image = image_or - median

    wcs = WCS(hdr)


    #if sources specified ... laod images
    only_rectangle = None
    if(args.target_image != None):
        with fits.open(args.target_image) as hdul:
            target_hdu = hdul[args.target_hdul_idx]
            target_hdr = target_hdu.header
            target_image = hdul[args.target_hdul_idx].data.astype(float)
        target_wcs = WCS(target_hdr)
        corners_sky = target_wcs.wcs_pix2world( [[0.,0.], [target_image.shape[0], target_image.shape[1]]] , 1)
        corners = wcs.wcs_world2pix(corners_sky, 1)
        corners = np.rint(corners)

        #test
        # corners_all_sky = target_wcs.wcs_pix2world( [[0.,0.], [0, target_image.shape[1]], [target_image.shape[0], 0], [target_image.shape[0], target_image.shape[1]]] , 1)
        # corners_all = wcs.wcs_world2pix(corners_all_sky, 1)
        # tmp =wcs.wcs_pix2world(corners_all,1)
        # tmp2 = target_wcs.wcs_world2pix(tmp,1)
        # print(corners_all)
        #import pdb; pdb.set_trace()

        left = np.min(corners[:,0])   - args.padding
        right = np.max(corners[:,0])  + args.padding
        bottom = np.min(corners[:,1]) - args.padding
        top = np.max(corners[:,1])    + args.padding
        only_rectangle = [left, right, bottom, top]



    sources_in_target_image = astrometry.find_sources(image, only_rectangle=only_rectangle, cutouts=args.cutout , sigma_threshold_for_source_detection=args.detection_threshold)
    print("Found {} sources".format(sources_in_target_image.shape[0]))

    observation_on_sky = wcs.wcs_pix2world(sources_in_target_image[["xcenter","ycenter"]], 1)
    #catalog_from_obs = np.zeros(observation_on_sky.shape[0], dtype={'names':('ra', 'dec', 'aperture_sum'),'formats':('f8', 'f8', 'f8')})
    catalog_from_obs = pd.DataFrame()
    catalog_from_obs["ra"]= observation_on_sky[:,0]
    catalog_from_obs["dec"]= observation_on_sky[:,1]
    catalog_from_obs["aperture_sum"]= sources_in_target_image["aperture_sum"]
    catalog_from_obs["mag"]= -1.* sources_in_target_image["aperture_sum"]#this is fine since we only use the mag to order the sources!
    catalog_from_obs.to_csv(args.filename_for_sources+".csv")

    #show as image
    fig = plt.figure()
    fig.canvas.set_window_title('Sources extracted from {}'.format(args.source_image))
    plt.xlabel("pixel x direction")
    plt.ylabel("pixel y direction")
    plt.imshow(image,cmap='Greys', origin='lower', norm=LogNorm())

    if(args.target_image != None):
        print(corners)
        ax = plt.gca()
        rect = Rectangle(corners[0,:],corners[1,0]-corners[0,0],corners[1,1]-corners[0,1],linewidth=1,edgecolor='r',facecolor='none')
        ax.add_patch(rect)

    positions = (sources_in_target_image['xcenter'], sources_in_target_image['ycenter'])
    apertures = CircularAperture(positions, r=4.)
    apertures.plot(color='blue', lw=1.5, alpha=0.5)
    plt.show()


if __name__ == '__main__':
    main()
