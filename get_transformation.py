"""Calculating the tranformation to make astrometry (wcs) consistent with the literature.

Author Lukas Wenzl
written in python 3

"""

#import matplotlib.pyplot as plt
import numpy as np

from astropy.wcs import WCS
from astropy.wcs import Wcsprm

import matplotlib.pyplot as plt

#from scipy.signal import correlate2d


import copy

import settings as s



def simple_offset(observation, catalog, wcsprm, report=""):
    """Get best offset in x, y direction.

    Parameters
    ----------
    observation : dataframe
        pandas dataframe with sources on the observation
    catalog : dataframe
        pandas dataframe with nearby sources from online catalogs with accurate astrometric information
    wcsprm
        Wold coordinates file
    report : str
        Previous part of the final report that will be extended by the method.

    Returns
    -------
    wcsprm, signal, report

    """
    report = report+"simple_offset aproach via a histogram \n"

    #catalog_on_sensor = wcsprm.wcs_world2pix(catalog[["ra", "dec"]], 1) #now using wcsprm
    catalog_on_sensor = wcsprm.s2p(catalog[["ra", "dec"]], 1)
    catalog_on_sensor  = catalog_on_sensor['pixcrd']
    #catalog_on_sensor[:,1]
    obs = [observation["xcenter"].values]
    cat = np.array( [catalog_on_sensor[:,0] ])
    distances_x =  (obs - cat.T).flatten()

    obs = [observation["ycenter"].values]
    cat = np.array( [catalog_on_sensor[:,1] ])
    distances_y =  (obs - cat.T).flatten()
    #vectorized distances, better by a factor. Went from a second for this method to well below a second,

    #plot to visualize if something is not working
    # plt.figure()
    # plt.hist(distances_x, bins=1500)
    # plt.xlim(-100,100)
    #plt.show()

    #to find the statistics i have to use the center! otherwise I get selection effects because i only downloaded a small radius

    #hist = plt.figure()
    binwidth= s.OFFSET_BINWIDTH #would there be a reason to make it bigger than 1?
    bins = [np.arange(min(distances_x), max(distances_x) + binwidth, binwidth), np.arange(min(distances_y), max(distances_y) + binwidth, binwidth)]

    #H, x_edges, y_edges,tmp = plt.hist2d(distances_x, distances_y, bins=bins) #to visualize for testing
    H, x_edges, y_edges = np.histogram2d(distances_x, distances_y, bins=bins)

    #finding the peak for the x and y distance where the two sets overlap
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

    current_central_pixel = wcsprm.crpix
    new_central_pixel = [current_central_pixel[0] + x_shift, current_central_pixel[1] +y_shift]
    wcsprm.crpix = new_central_pixel
    return wcsprm, signal, report


def rotation_matrix(angle):
    rot = [[np.cos(angle), np.sin(angle)], [-np.sin(angle), np.cos(angle)]]
    return rot

def rotate(wcsprm, rot):
    """Help method for offset_with_orientation. Set the different rotations in the header."""
    #hdr["PC1_1"] =rot[0][0]
    #hdr["PC1_2"] =rot[1][0]
    #hdr["PC2_1"]=rot[0][1]
    #hdr["PC2_2"] =rot[1][1]
    pc = wcsprm.get_pc()
    pc_rotated =rot@ pc
    wcsprm.pc = pc_rotated
    return wcsprm

def offset_with_orientation(observation, catalog, wcsprm, verbose=True, fast=False, report_global="", INCREASE_FOV_FLAG=False, silent=False):
    """Use simple_offset(...) but with trying 0,90,180,270 rotation.

    Parameters
    ----------
    observation : dataframe
        pandas dataframe with sources on the observation
    catalog : dataframe
        pandas dataframe with nearby sources from online catalogs with accurate astrometric information
    wcsprm
        Wcsprm file
    verbose : boolean
        Set to False to supress output to the console
    fast : boolean
        If true will run with subset of the sources to increase speed.


    Returns
    -------
    wcs, signal, report

    """
    observation = copy.copy(observation)
    N_SOURCES = observation.shape[0]
    if(fast):
        if(N_SOURCES > s.USE_N_SOURCES):
            N_SOURCES = s.USE_N_SOURCES
        observation = observation.nlargest(N_SOURCES, 'aperture_sum')
        if(INCREASE_FOV_FLAG):
            N_CATALOG = N_SOURCES *12
        else:
            N_CATALOG = N_SOURCES * 4
        catalog = catalog.nsmallest(N_CATALOG, 'mag')
    if(verbose):
        #print("--------------------------------------")
        #print("offset_with_orientation, seaching for offset while considering reflections and 0,90,180,270 rotations")
        print("offset_with_orientation, seaching for offset while considering 0,90,180,270 rotations")
        if(fast):
            print("running in fast mode")
    rotations = [ [[1,0],[0,1]], [[-1,0],[0,-1]], #already checking for reflections with the scaling and general rotations, but in case rot_scale is of it is nice to have
                  [[-1,0],[0,1]], [[1,0],[0,-1]],
                  [[0,1],[1,0]], [[0,-1],[-1,0]],
                  [[0,-1],[1,0]], [[0,1],[-1,0]],
                ]
    # rotations = [ [[1,0],[0,1]], [[-1,0],[0,-1]],
    #               [[0,-1],[1,0]], [[0,1],[-1,0]],
    #             ]
    wcsprm_global = copy.copy(wcsprm)
    results = []
    for rot in rotations:
        if(verbose):
            print("Trying rotation {}".format(rot))
            #print(report)
        wcsprm = rotate(copy.copy(wcsprm_global),rot)
        #wcs = WCS(...)
        report = report_global +  "---- Report for rotation {} ---- \n".format(rot)
        wcsprm, signal, report = simple_offset(observation, catalog, wcsprm, report)
        results.append([copy.copy(wcsprm),signal,report])


    signals = [i[1] for i in results]
    median = np.median(signals)
    i = np.argmax(signals)
    wcsprm = results[i][0]
    signal = signals[i]
    #hist = results[i][3]
    report = results[i][2]
    report = report + "A total of {} sources from the fits file where used. \n".format(N_SOURCES)
    report = report + "The signal (#stars) is {} times higher than noise outlierers for other directions. (more than 2 would be nice, typical: 8 for PS)\n".format(signals[i]/median)


    if(verbose):
        print("We found the following world coordinates: ")
        print(WCS(wcsprm.to_header()))
        print("And here is the report:")
        print(report)
        print("-----------------------------")
    off = wcsprm.crpix - wcsprm_global.crpix
    if(not silent):
        print("Found offset {:.3g} in x direction and {:.3g} in y direction".format(off[0], off[1]))

    return wcsprm, signal, report


def calculate_dist(data_x, data_y):
    data_x = np.array(data_x)
    data_y = np.array(data_y)

    distances_x =  (data_x - data_x.T)
    distances_y =  (data_y - data_y.T)

    #only use off diagonal elements
    distances_x = distances_x[np.where(~np.eye(distances_x.shape[0],dtype=bool))]
    distances_y = distances_y[np.where(~np.eye(distances_y.shape[0],dtype=bool))]

    distances = np.sqrt(distances_x**2 + distances_y**2)

    return distances


def calculate_log_dist(data_x, data_y):
    log_distances = np.log(calculate_dist(data_x, data_y)+np.finfo(float).eps)
    return log_distances

def calculate_angles(data_x, data_y):
    #plan: i need all pairs: vector differences, then angle with x axis.
    #so
    data_x = np.array(data_x)
    data_y = np.array(data_y)
    vec_x = data_x - data_x.T
    vec_y = data_y - data_y.T
    vec_x = vec_x[np.where(~np.eye(vec_x.shape[0],dtype=bool))]
    vec_y = vec_y[np.where(~np.eye(vec_y.shape[0],dtype=bool))]
    angles = np.arctan2(vec_x,vec_y)
    angles = angles % (2*np.pi) #make sure angles are between 0 and 2 pi
    angles[np.where(angles > np.pi)] = -1*(2*np.pi - angles[np.where(angles > np.pi)]) #shift to -pi to pi
    return angles

def peak_with_histogram(obs_x, obs_y, cat_x, cat_y):
    """Find the relation between the two sets. Either the positional offset (not used for that at the moment) or the scale+angle between them.

    This should be replaced with a method using convolution instead. Also currently the bandwith is choosen quite random

    Parameters
    ----------
    obs_x: array
        first axis to consider of observations (x pos or log distance)
    obs_y: array
        second axis to consider of observations (x pos or log distance)
    cat_x: array
        first axis to consider of catalog data (x pos or log distance)
    cat_x: array
        first axis to consider of catalog data (x pos or log distance)

    Returns
    -------
    x_shift, y_shift

    """
    obs_x = obs_x[:, np.newaxis]#to make transpose work!
    cat_x = cat_x[:, np.newaxis]#to make transpose work!
    distances_x =  (obs_x - cat_x.T).flatten()
    obs_y = obs_y[:, np.newaxis]#to make transpose work!
    cat_y = cat_y[:, np.newaxis]#to make transpose work!
    distances_y =  (obs_y - cat_y.T).flatten()
    #to find the statistics i have to use the center! otherwise I get selection effects because i only downloaded a small radius
    binwidth= 0.001 #would there be a reason to make it bigger?

    bins = [np.arange(min(distances_x), max(distances_x) + binwidth, binwidth), np.arange(min(distances_y), max(distances_y) + binwidth*10, binwidth*10)]

    plt.figure()
    H, x_edges, y_edges,tmp = plt.hist2d(distances_x, distances_y, bins=bins)

    plt.show()

    H, x_edges, y_edges = np.histogram2d(distances_x, distances_y, bins=bins)

    peak = np.argwhere(H == H.max())[0] #take fisrt peak


    x_shift = (x_edges[peak[0]] + x_edges[peak[0]+1])/2
    y_shift = (y_edges[peak[1]] + y_edges[peak[1]+1])/2

    return x_shift,y_shift

def cross_corr_to_fourier_space(a):
    "Tranform 2D array into fourier space. Uses padding and normalization."
    aa = (a - np.mean(a))/np.std(a)
    #aaa = np.pad(aa, (aa.shape[0]//2,aa.shape[0]//2), 'constant') #wraps around so half the size should be fine, padds 2D array with zeros
    aaa = np.pad(aa, (2,2), 'constant')
    #I think for the angle the padding is not relevant since it is cyclic, for the scale i just don't care about the extremes
    ff_a = np.fft.fft2(aaa)
    return ff_a

def peak_with_cross_correlation(log_distance_obs, angle_obs, log_distance_cat, angle_cat, scale_guessed=False):
    """Find the relation between the two sets. Either the positional offset (not used for that at the moment) or the scale+angle between them.

    This is using cross correlation

    Parameters
    ----------
    log_distance_obs: array
        first axis to consider of observations ( log distance)
    angle_obs: array
        second axis to consider of observations ( angle)
    log_distance_cat: array
        first axis to consider of catalog data (log distance)
    angle_cat: array
        first axis to consider of catalog data (angle)

    Returns
    -------
    x_shift, y_shift

    """
    if(scale_guessed==False):
        minimum_distance = np.log(8)#minimum pixel distance
        maximum_distance = max(log_distance_obs)
    else:
        #broader distance range if the scale is just a guess so there is a higher chance to find the correct one
        minimum_distance = min([min(log_distance_cat), min(log_distance_obs)])
        maximum_distance = max([max(log_distance_cat), max(log_distance_obs)])

    bins_dist, binwidth_dist = np.linspace(minimum_distance, maximum_distance, 3000, retstep=True)
    # print(binwidth_dist)
    # print(np.e**(binwidth_dist))
    bins_ang, binwidth_ang = np.linspace(min([min(angle_cat), min(angle_obs)]), max([max(angle_cat),max(angle_obs)]), 360*3, retstep=True)
    #print(binwidth_ang/2/np.pi*360)
    bins = [bins_dist, bins_ang]#min max of both
    H_obs, x_edges_obs, y_edges_obs = np.histogram2d(log_distance_obs, angle_obs, bins=bins)
    H_cat, x_edges_cat, y_edges_cat = np.histogram2d(log_distance_cat,angle_cat, bins=bins)

    #print(bins)
    H_obs = (H_obs- np.mean(H_obs))/ np.std(H_obs)
    H_cat = (H_cat- np.mean(H_cat))/ np.std(H_cat)

    ff_obs = cross_corr_to_fourier_space(H_obs)
    ff_cat = cross_corr_to_fourier_space(H_cat)


    cross_corr = ff_obs*np.conj(ff_cat)

    #frequency cut off
    step = 1 # maybe in arcsec??, this is usually the timestep to get a frequency
    frequ = np.fft.fftfreq(ff_obs.size, d=step).reshape(ff_obs.shape)
    max_frequ = np.max(frequ)#frequ are symmatric - to +
    threshold = 0.02* max_frequ #todo make threshold changable
    #print(threshold)
    cross_corr[(frequ<threshold)&(frequ>-threshold)] = 0 ##how to choose the frequency cut off?


    cross_corr = np.real(np.fft.ifft2(cross_corr))
    cross_corr = np.fft.fftshift(cross_corr) #zero shift is at (0,0), this move it to the middle


    peak = np.argwhere(cross_corr == cross_corr.max())[0] #take fisrt peak

    around_peak = cross_corr[peak[0]-1:peak[0]+2, peak[1]-1:peak[1]+2]

    #finding the sub pixel shift of the true peak
    # print("sub pixel offset:")
    peak_x_subpixel = np.sum(np.sum(around_peak, axis=1)*(np.arange(around_peak.shape[0])+1))/np.sum(around_peak)-2#should be peak[0] offset
    peak_y_subpixel = np.sum(np.sum(around_peak, axis=0)*(np.arange(around_peak.shape[1])+1))/np.sum(around_peak)-2#should be peak[1] offset

    signal = np.sum(cross_corr[peak[0]-1:peak[0]+2, peak[1]-1:peak[1]+2])   #sum up signal in fixed aperture 1 pixel in each direction around the peak, so a 3x3 array, total 9 pixel
    #signal_wide = np.sum(cross_corr[peak[0]-4:peak[0]+5, peak[1]-4:peak[1]+5])
    # #report = report+"signal wide (64pixel) - signal (9pixel)  = {}. If this value is large then there might be rotation or scaling issues. \n".format(signal_wide-signal)
    # #aperture = 9 #pixel
    # ##signal wide: 64 pixel total
    # # print(signal, signal_wide)
    middle_x = cross_corr.shape[0]/2#is that corroct? yes I think so, shape is uneven number and index counting starts at 0
    middle_y = cross_corr.shape[1]/2

    x_shift = (peak[0]+peak_x_subpixel-middle_x)*binwidth_dist
    y_shift = (peak[1]+peak_y_subpixel-middle_y)*binwidth_ang

    scaling= np.e**(-x_shift)
    rotation = y_shift#/2/np.pi*360 maybe easier in rad

    # plt.figure()
    # #plt.imshow(correlate2d(H_obs,H_cat), interpolation="nearest")
    # plt.imshow(cross_corr, interpolation="nearest")
    # plt.colorbar()
    # plt.show()

    return scaling, rotation, signal

def scale(wcsprm, scale_factor):
    pc = wcsprm.get_pc()
    pc_scaled =scale_factor* pc
    wcsprm.pc = pc_scaled
    return wcsprm


def get_scaling_and_rotation(observation, catalog, wcsprm, scale_guessed, verbose=True, report_global={}): #INCREASE_FOV_FLAG=False
    """Calculate the scaling and rotation compared to the catalog based on the method of Kaiser et al. (1999).

    This should be quite similar to the approach by SCAMP.

    Parameters
    ----------
    observation : dataframe
        pandas dataframe with sources on the observation
    catalog : dataframe
        pandas dataframe with nearby sources from online catalogs with accurate astrometric information
    wcsprm
        Wcsprm file
    verbose : boolean
        Set to False to supress output to the console


    Returns
    -------
    wcs, signal, report

    """
    catalog_on_sensor = wcsprm.s2p(catalog[["ra", "dec"]], 1)
    catalog_on_sensor  = catalog_on_sensor['pixcrd']
    obs_x = [observation["xcenter"].values]
    cat_x = np.array( [catalog_on_sensor[:,0] ])
    obs_y = [observation["ycenter"].values]
    cat_y = np.array( [catalog_on_sensor[:,1] ])

    log_distances_obs = calculate_log_dist(obs_x, obs_y)
    log_distances_cat = calculate_log_dist(cat_x, cat_y)

    angles_obs = calculate_angles(obs_x,obs_y)
    angles_cat = calculate_angles(cat_x, cat_y)


    scaling, rotation, signal = peak_with_cross_correlation(log_distances_obs, angles_obs, log_distances_cat, angles_cat, scale_guessed= scale_guessed)
    scaling_reflected, rotation_reflected, signal_reflected =  peak_with_cross_correlation(log_distances_obs, -angles_obs, log_distances_cat, angles_cat, scale_guessed=scale_guessed)

    if(signal_reflected > signal):
        is_reflected = True
        confidence = signal_reflected/signal
        scaling = scaling_reflected
        rotation = rotation_reflected
    else:
        is_reflected = False
        confidence = signal/signal_reflected

    rot = rotation_matrix(rotation)
    if(is_reflected):
        rot = np.array([[1, 0], [0,-1]])@rot #reflecting, but which direction?? this is a flip of the y axis
    wcsprm_new = rotate(copy.copy(wcsprm), rot)

    wcsprm_new = scale(wcsprm_new, scaling)

    if(is_reflected):
        refl = ""
    else:
        refl = "not "
    print("Found a rotation of {:.3g} deg and the pixelscale was scaled with the factor {:.3g}.".format(rotation/2/np.pi*360,scaling)+"The image was "+refl+"mirrored.")
    if(verbose):
        print("The confidence level is {}. values between 1 and 2 are bad. Much higher values are best.".format(confidence))
        print("Note that there still might be a 180deg rotation. If this is the case it should be correct in the next step")
    #report_global["rotation"] +

    #scaling, rotation = peak_with_histogram(log_distances_obs, angles_obs, log_distances_cat, angles_cat)

    # print(scaling)
    # print(rotation)


    return wcsprm_new

def calculate_rms(observation, catalog, wcsprm):
    """Calculate the root mean square deviation of the astrometry fit"""
    #finding pixel scale
    on_sky = wcsprm.p2s([[0,0],[1,1]], 0)["world"]
    px_scale = np.sqrt((on_sky[0,0]-on_sky[1,0])**2+(on_sky[0,1]-on_sky[1,1])**2)
    px_scale = px_scale*60*60 #in arcsec
    obs_x, obs_y, cat_x, cat_y, distances = find_matches(observation, catalog, wcsprm, threshold=3)
    rms = np.sqrt(np.mean(np.square(distances)))
    print("Within 3  pixel or {:.3g} arcsec {} sources where matched. The rms is {:.3g} pixel or {:.3g} arcsec".format(px_scale*3, len(obs_x), rms, rms*px_scale))
    obs_x, obs_y, cat_x, cat_y, distances = find_matches(observation, catalog, wcsprm, threshold=5)
    rms = np.sqrt(np.mean(np.square(distances)))
    print("Within 5  pixel or {:.3g} arcsec {} sources where matched. The rms is {:.3g} pixel or {:.3g} arcsec".format(px_scale*5,len(obs_x), rms, rms*px_scale))
    obs_x, obs_y, cat_x, cat_y, distances = find_matches(observation, catalog, wcsprm, threshold=s.RMS_PX_THRESHOLD)
    rms = np.sqrt(np.mean(np.square(distances)))
    print("Within {} pixel or {:.3g} arcsec {} sources where matched. The rms is {:.3g} pixel or {:.3g} arcsec".format(s.RMS_PX_THRESHOLD, px_scale*s.RMS_PX_THRESHOLD,len(obs_x), rms, rms*px_scale))
    return {"radius_px":s.RMS_PX_THRESHOLD, "matches":len(obs_x), "rms":rms}

def find_matches_keep_catalog_info(observation, catalog, wcsprm, threshold=5):
    catalog_on_sensor = wcsprm.s2p(catalog[["ra", "dec"]], 1)
    catalog_on_sensor  = catalog_on_sensor['pixcrd']
    obs_x = [observation["xcenter"].values]
    cat_x = np.array( [catalog_on_sensor[:,0] ])
    obs_y = [observation["ycenter"].values]
    cat_y = np.array( [catalog_on_sensor[:,1] ])

    distances_x =  (obs_x - cat_x.T)#.flatten()
    distances_y =  (obs_y - cat_y.T)#.flatten()
    obs_x = obs_x[0]
    obs_y = obs_y[0]
    cat_x = cat_x[0]
    cat_y = cat_y[0]


    #find closest points to each source in observaion
    distances = np.min(distances_x**2 + distances_y**2, axis=0)
    matches = np.argmin(distances_x**2 + distances_y**2, axis=0)

    #search for all matches within threshold
    obs_matched = observation.iloc[distances<threshold]
    #obs_y = obs_y[distances<threshold]
    cat_matched = catalog.iloc[matches[distances<threshold]]
    #cat_y = cat_y[matches[distances<threshold]]
    distances = distances[distances < threshold]
    return obs_matched, cat_matched, distances

def find_matches(observation, catalog, wcsprm, threshold=5):
    catalog_on_sensor = wcsprm.s2p(catalog[["ra", "dec"]], 1)
    catalog_on_sensor  = catalog_on_sensor['pixcrd']
    #catalog_on_sensor[:,1]
    obs_x = [observation["xcenter"].values]
    cat_x = np.array( [catalog_on_sensor[:,0] ])
    obs_y = [observation["ycenter"].values]
    cat_y = np.array( [catalog_on_sensor[:,1] ])

    distances_x =  (obs_x - cat_x.T)#.flatten()
    distances_y =  (obs_y - cat_y.T)#.flatten()
    obs_x = obs_x[0]
    obs_y = obs_y[0]
    cat_x = cat_x[0]
    cat_y = cat_y[0]


    #find closest points to each source in observaion
    distances = np.min(distances_x**2 + distances_y**2, axis=0)
    matches = np.argmin(distances_x**2 + distances_y**2, axis=0)

    #search for all matches within threshold
    obs_x = obs_x[distances<threshold]
    obs_y = obs_y[distances<threshold]
    cat_x = cat_x[matches[distances<threshold]]
    cat_y = cat_y[matches[distances<threshold]]
    distances = distances[distances < threshold]
    return obs_x, obs_y, cat_x, cat_y, distances






def fine_transformation(observation, catalog, wcsprm, threshold=1, verbose=True, compare_threshold=3, skip_rot_scale=False):
    """Final improvement of registration. This requires that the wcs is already accurate to a few pixels.

    Parameters
    ----------
    observation : dataframe
        pandas dataframe with sources on the observation
    catalog : dataframe
        pandas dataframe with nearby sources from online catalogs with accurate astrometric information
    wcsprm
        Wcsprm file
    threshold : float
        maximum separation to consider two sources matches
    verbose : boolean
        print details

    Returns
    -------
    wcsprm

    """
    wcsprm_original = wcsprm
    wcsprm = copy.copy(wcsprm)

    if(threshold == 20 or threshold == 100):
        observation = observation.nlargest(5, "aperture_sum")
        #print("using 5 brightest sources")
    obs_x, obs_y, cat_x, cat_y, _ = find_matches(observation, catalog, wcsprm, threshold=threshold)
    if(len(obs_x)<4):
        return wcsprm_original, 0 #not enough matches
    #seems to work

    #angle:
    angle_offset = -calculate_angles([obs_x],[obs_y])+calculate_angles([cat_x],[cat_y])
    log_distances_obs = calculate_log_dist([obs_x],[obs_y])
    log_distances_cat = calculate_log_dist([cat_x],[cat_y])
    threshold_min = np.log(20) #minimum distance to make usefull scaling or angle estimation
    if(threshold == 10):
        threshold_min = np.log(200)
    #if(threshold == 20):
    #    threshold_min = np.log(250)
    mask = (log_distances_obs>threshold_min)&(log_distances_cat>threshold_min)
    scale_offset = -log_distances_obs+log_distances_cat

    #only take usefull datapoints
    angle_offset = angle_offset[mask]
    scale_offset = scale_offset[mask]

    rotation = np.mean(angle_offset)
    scaling = np.e**(np.mean(scale_offset))

    rot = rotation_matrix(rotation)
    if (not skip_rot_scale):
        wcsprm = rotate(wcsprm, rot)
        if (scaling > 0.9 and scaling < 1.1):
            wcsprm = scale(wcsprm, scaling)
        else:
            #print("fine transformation failed. Scaling calculation failed")
            return wcsprm_original,0

    #need to recalculate positions
    obs_x, obs_y, cat_x, cat_y, _ = find_matches(observation, catalog, wcsprm, threshold=threshold)
    if(len(obs_x)<4):
        return wcsprm_original,0

    #offset:
    x_shift = np.mean(obs_x- cat_x)
    y_shift = np.mean(obs_y- cat_y)

    current_central_pixel = wcsprm.crpix
    new_central_pixel = [current_central_pixel[0] + x_shift, current_central_pixel[1] +y_shift]
    wcsprm.crpix = new_central_pixel

    obs_x, obs_y, cat_x, cat_y, distances = find_matches(observation, catalog, wcsprm, threshold=compare_threshold)
    rms = np.sqrt(np.mean(np.square(distances)))
    score = len(obs_x)/(rms+10) #number of matches within 3 pixel over rms+1 (so its bigger than 0)
    return wcsprm, score
