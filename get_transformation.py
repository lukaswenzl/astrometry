"""Calculating the tranformation to make astrometry (wcs) consistent with the literature.

Author Lukas Wenzl
written in python 3

"""

#import matplotlib.pyplot as plt
import numpy as np

from astropy.wcs import WCS


import copy

import settings as s

def simple_offset(observation, catalog, wcs, report=""):
    """Get best offset in x, y direction.

    Parameters
    ----------
    observation : dataframe
        pandas dataframe with sources on the observation
    catalog : dataframe
        pandas dataframe with nearby sources from online catalogs with accurate astrometric information
    wcs
        Wold coordinates file
    report : str
        Previous part of the final report that will be extended by the method.

    Returns
    -------
    wcs, signal, report

    """
    report = report+"simple_offset aproach via a histogram \n"
    # distances_x = []
    # distances_y = []
    # for _, obs in observation.iterrows():
    #     for _, cat in catalog.iterrows():
    #         #using PS!
    #         cat_on_sensor = wcs.wcs_world2pix([[cat["ra"],cat["dec"]]],1)[0]
    #
    #         dist_x = (obs['xcenter']-cat_on_sensor[0])
    #         dist_y = (obs['ycenter']-cat_on_sensor[1])
    #
    #         distances_x.append(dist_x)
    #         distances_y.append(dist_y)

    #vectorized distances, better by a factor. Went from a second for this method to well below a second,
    #complete runtime for sample image went from 5.3 to 2.9 seconds (mostly overhead from download)
    catalog_on_sensor = wcs.wcs_world2pix(catalog[["ra", "dec"]], 1)
    #catalog_on_sensor[:,1]
    obs = [observation["xcenter"].values]
    cat = np.array( [catalog_on_sensor[:,0] ])
    distances_x =  (obs - cat.T).flatten()

    obs = [observation["ycenter"].values]
    cat = np.array( [catalog_on_sensor[:,1] ])
    distances_y =  (obs - cat.T).flatten()

    # plt.figure()
    # plt.hist(distances_x, bins=1500)
    # plt.xlim(-100,100)
    #plt.show()

    #to find the statistics i have to use the center! otherwise I get selection effects because i only downloaded a small radius

    #hist = plt.figure()
    binwidth= 1 #would there be a reason to make it bigger?
    bins = [np.arange(min(distances_x), max(distances_x) + binwidth, binwidth), np.arange(min(distances_y), max(distances_y) + binwidth, binwidth)]
    #bins = [range(-100,100, 4), range(-100,100, 4)]


    #H, x_edges, y_edges,tmp = plt.hist2d(distances_x, distances_y, bins=bins)
    H, x_edges, y_edges = np.histogram2d(distances_x, distances_y, bins=bins)

    #weights with brightness?
    # plt.xlim(-100,100)
    # plt.ylim(-100,100)
    #plt.show()
    #
    #find background IRRELEVANT
    # H_center,_,_ = np.histogram2d(distances_x, distances_y, bins=bins, range=[[-100,100], [-100,100]]) #same calibration as above just without image output and only for central range
    # bkg = np.mean(H_center) #median would just be 0 or 1, not usefull, maybe I need to exclude the peak?
    # bkg_sigma = np.std(H_center)

    peak = np.argwhere(H == H.max())[0] #take fisrt peak
    signal = np.sum(H[peak[0]-1:peak[0]+2, peak[1]-1:peak[1]+2])   #sum up signal in fixed aperture 1 pixel in each direction around the peak, so a 3x3 array, total 9 pixel
    signal_wide = np.sum(H[peak[0]-4:peak[0]+5, peak[1]-4:peak[1]+5])
    report = report+"signal wide (64pixel) - signal (9pixel)  = {}. If this value is large then there might be rotation or scaling issues. \n".format(signal_wide-signal)
    #aperture = 9 #pixel
    ##signal wide: 64 pixel total
    # print(signal, signal_wide)

    x_shift = (x_edges[peak[0]] + x_edges[peak[0]+1])/2
    y_shift = (y_edges[peak[1]] + y_edges[peak[1]+1])/2

    report = report+"We find an offset of {} in the x direction and {} in the y direction \n".format(x_shift, y_shift)
    report = report+"{} sources are fitting well with this offset. \n".format(signal)

    current_central_pixel = wcs.wcs.crpix
    new_central_pixel = [current_central_pixel[0] + x_shift, current_central_pixel[1] +y_shift]
    wcs.wcs.crpix = new_central_pixel
    #wcs.wcs
    return wcs, signal, report



def rotate(hdr, rot):
    """Help method for offset_with_orientation. Set the different rotations in the header."""
    hdr["PC1_1"] =rot[0][0]
    hdr["PC1_2"] =rot[1][0]
    hdr["PC2_1"]=rot[0][1]
    hdr["PC2_2"] =rot[1][1]
    return hdr
def offset_with_orientation(observation, catalog, hdr_file, verbose=True, fast=False, report_global=""):
    """Use simple_offset(...) but tries 0,90,180,270 rotation and reflections.

    Parameters
    ----------
    observation : dataframe
        pandas dataframe with sources on the observation
    catalog : dataframe
        pandas dataframe with nearby sources from online catalogs with accurate astrometric information
    hdr_file
        header from fits files
    verbose : boolean
        Set to False to supress output to the console
    fast : boolean
        If true will run with subset of the sources to increase speed.


    Returns
    -------
    wcs, signal, report

    """
    observation_all = copy.copy(observation)
    N_SOURCES = observation.shape[0]
    if(fast):
        if(N_SOURCES > s.USE_N_SOURCES):
            N_SOURCES = s.USE_N_SOURCES
        observation = observation.nlargest(N_SOURCES, 'aperture_sum')
        catalog = catalog.nsmallest(N_SOURCES*4, 'mag')
    if(verbose):
        #print("--------------------------------------")
        print("offset_with_orientation, seaching for offset while considering reflections and 0,90,180,270 rotations")
        if(fast):
            print("running in fast mode")
    rotations = [ [[1,0],[0,1]], [[-1,0],[0,-1]],
                  [[-1,0],[0,1]], [[1,0],[0,-1]],
                  [[0,1],[1,0]], [[0,-1],[-1,0]],
                  [[0,-1],[1,0]], [[0,1],[-1,0]],
                ]
    hdr_local = copy.copy(hdr_file)
    results = []
    for rot in rotations:
        if(verbose):
            print("Trying rotation {}".format(rot))
            #print(report)
        hdr_local = rotate(hdr_local,rot)
        wcs = WCS(hdr_local)
        report = report_global +  "---- Report for rotation {} ---- \n".format(rot)
        wcs, signal, report = simple_offset(observation, catalog, wcs, report)
        results.append([wcs,signal,report])


    signals = [i[1] for i in results]
    median = np.median(signals)
    i = np.argmax(signals)
    wcs = results[i][0]
    #hist = results[i][3]
    report = results[i][2]
    report = report + "A total of {} sources from the fits file where used. \n".format(N_SOURCES)
    report = report + "The signal (#stars) is {} times higher than noise outlierers for other directions. (more than 2 would be nice, typical: 8 for PS)\n".format(signals[i]/median)

    if(fast & (signals[i] < N_SOURCES * s.FASTMODE_THRESHOLD)):
        print("Not enough sources were matched. trying again without fast mode")
        report_global = report_global + "Turned fast mode off because not enough sources were detected"
        wcs, signal, report = offset_with_orientation(observation_all, catalog, hdr_file, report_global=report_global)
    elif(verbose):
        print("We found the following world coordinates: ")
        print(wcs)
        print("And here is the report:")
        print(report)
        print("-----------------------------")

    return wcs, signal, report


# def general_transformation(observation, catalog, wcs, verbose=True):
#     """Find general transformation. assumes square pixel.
#
#     Parameters
#     ----------
#     observation : dataframe
#         pandas dataframe with sources on the observation
#     catalog : dataframe
#         pandas dataframe with nearby sources from online catalogs with accurate astrometric information
#     wcs
#         Wold coordinates file
#
#     Returns
#     -------
#
#     """
#     N_SOURCES = observation.shape[0]
#     if(N_SOURCES > s.USE_N_SOURCES):
#         N_SOURCES = s.USE_N_SOURCES
#     observation = observation.nlargest(N_SOURCES, 'aperture_sum')
#     catalog = catalog.nsmallest(N_SOURCES*4, 'mag')
#     if(verbose):
#         print("--------------------------------------")
#         print("Searching for a general transformation")
#
#     #I need one match. I will guess that the two brightest objects align and will go down from there
#     match_guess = []
#     for i in range(3):
#         match_guess.append([i,i])
#         guesses = [[i,j] for j in range(i)]
#         match_guess.append(guesses)
#         guesses = [[j,i] for j in range(i)]
#         match_guess.append(guesses)
#
#     print(match_guess)


def calculate_rms(observation, catalog, wcs):
    """Calculate the rms of the final transformation."""
    catalog_on_sensor = wcs.wcs_world2pix(catalog[["ra", "dec"]], 1)
    obs = [observation["xcenter"].values]
    cat = np.array( [catalog_on_sensor[:,0] ])
    distances_x =  (obs - cat.T)#.flatten()

    obs = [observation["ycenter"].values]
    cat = np.array( [catalog_on_sensor[:,1] ])
    distances_y =  (obs - cat.T)#.flatten()

    #find closest points to each source is observaion
    distances = np.min(distances_x**2 + distances_y**2, axis=0)
    #print(distances)
    #indices for the 10 brightest sources

    N_OBS = observation.shape[0]
    rms = np.sqrt(np.mean(np.square(distances)))
    print("rms {:.3g} for all {} sources detected. (likely high because not all will be in the catalog)".format(rms, N_OBS))
    if(N_OBS > 50):
        N_OBS = 50
        ind = np.argpartition(observation["aperture_sum"].values, -N_OBS)[-N_OBS:]
        rms = np.sqrt(np.mean(np.square(distances[ind])))
        print("rms {:.3g} for the {} brightest sources detected. (likely high because not all will be in the catalog)".format(rms, N_OBS))
    if(N_OBS > 10):
        N_OBS = 10
        ind = np.argpartition(observation["aperture_sum"].values, -N_OBS)[-N_OBS:]
        rms = np.sqrt(np.mean(np.square(distances[ind])))
        print("rms {:.3g} for the {} brightest sources detected".format(rms, N_OBS))

    N_OBS = 50
    if(len(distances) < 50):
        N_OBS = len(distances)
    closest = np.partition(distances, N_OBS)[:N_OBS]
    rms = np.sqrt(np.mean(np.square(closest)))
    print("rms {:.3g} for the {} closest sources detected. (If below one the tranforation likely worked)".format(rms, N_OBS))

    if(N_OBS > 10):
        N_OBS = 10
    closest = np.partition(distances, N_OBS)[:N_OBS]
    rms = np.sqrt(np.mean(np.square(closest)))
    print("rms {:.3g} for the {} closest sources detected".format(rms, N_OBS))
    print("-------------------")




def fine_transformation(observation, catalog, wcs, verbose=True):
    """TODO Final improvement of registration. This requires that the wcs is already accurate to a few pixels.

    Parameters
    ----------
    observation : dataframe
        pandas dataframe with sources on the observation
    catalog : dataframe
        pandas dataframe with nearby sources from online catalogs with accurate astrometric information

    Returns
    -------
    wcs

    """
    if(verbose):
        print("---- Final fine subpixel improvement of the registration. ----")
    catalog_on_sensor = wcs.wcs_world2pix(catalog[["ra", "dec"]], 1)
    obs = [observation["xcenter"].values]
    cat = np.array( [catalog_on_sensor[:,0] ])
    distances_x =  (obs - cat.T)#.flatten()

    obs = [observation["ycenter"].values]
    cat = np.array( [catalog_on_sensor[:,1] ])
    distances_y =  (obs - cat.T)#.flatten()

    #find closest points to each source is observaion
    distances = np.min(distances_x**2 + distances_y**2, axis=0)
    distances_index = np.argmin(distances_x**2 + distances_y**2, axis=0)
    index1 = distances<=2.12 #1,5 pixel diagonally so sqrt(2)*3/2
    if(verbose):
        print("{} sources are matched within 2.12 pixels. Will register them.".format(index1.sum()))
    index2 =distances_index[index1] #reconstruct the corresponding observation and catalog items

    x_shift = np.mean(distances_x[index2,index1])
    y_shift = np.mean(distances_y[index2,index1])
    total_offset = (x_shift**2 + y_shift**2)**(0.5)


    if(total_offset > 5):
        print("Fine transformation failed, offset would be more than 5 pixel. Will not use fine transformation and quit")
        return wcs

    print("We find a subpixel offset of {:.3g} in the x direction and {:.3g} in the y direction \n".format(x_shift, y_shift))

    current_central_pixel = wcs.wcs.crpix
    new_central_pixel = [current_central_pixel[0] + x_shift, current_central_pixel[1] +y_shift]
    wcs.wcs.crpix = new_central_pixel
    #wcs.wcs
    return wcs
