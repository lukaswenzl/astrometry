"""Methods to get catalog data online.

Author Lukas Wenzl
written in python 3

"""

#
#
#Author Lukas Wenzl
#written in python 3



#from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier

#import pandas as pd

#import settings as s






def get_GAIA_data(coord, radius):
    """Query Gaia database.

    Parameters
    ----------
    coord  : SkyCoord
        Position to search.
    radius : float
        Radius of search cone.

    Returns
    -------
    objects
        table with the objects including the following info: ra (deg), ra_error (milliarcsec), dec (deg), dec_error (milliarcsec), mag

    """
    j = Gaia.cone_search_async(coord, radius)
    r = j.get_results()
    #r.pprint()
    #r.show_in_browser()

    catalog_data = r.to_pandas()
    catalog_data["mag"] = catalog_data["phot_g_mean_mag"]
    catalog_data = catalog_data[["ra", "ra_error", "dec", "dec_error", "mag"]]
    print("Found {} sources in GAIA within a radius of {}".format(catalog_data.shape[0], radius))
    #columns: ra, ra_error, dec, dec_error, mag
    return catalog_data


def get_PS_data(coord, radius):
    """Query Vizier for Panstarrs Dr1 data. PS1 was calibrated with GAIA and goes darker so if it is availabe it is likely the best choice.

    Parameters
    ----------
    coord  : SkyCoord
        Position to search.
    radius : float
        Radius of search cone.

    Returns
    -------
    objects
        table with the objects including the following info: ra (deg), ra_error (milliarcsec), dec (deg), dec_error (milliarcsec), mag

    """
    Vizier.ROW_LIMIT = -1
    j = Vizier.query_region(coord, radius=radius, catalog="II/349/ps1")
    if (j==[]):
        return None
    r = j[0]
    #r.pprint()
    #r.show_in_browser()

    catalog_data = r.to_pandas()
    magnitudes = ["gmag", "rmag", "imag", "zmag", "ymag"]
    catalog_data["mag"] = catalog_data[magnitudes].min(axis=1) #there should not be an issue with -999 values since in Vizier nulls are used fro missing values not -999
    catalog_data.rename(index=str, columns={"RAJ2000": "ra", "e_RAJ2000":"ra_error", "DEJ2000": "dec", "e_DEJ2000":"dec_error"}, inplace=True)
    catalog_data = catalog_data[["ra", "ra_error", "dec", "dec_error", "mag"]]
    #print(catalog_data)
    #columns: ra, ra_error, dec, dec_error, mag
    #
    print("Found {} sources in PS1 within a radius of {}".format(catalog_data.shape[0], radius))
    #
    return catalog_data





def get_data(pos, radius, source):
    """Query databases.

    Parameters
    ----------
    pos : SkyCoord
        Position to search.
    radius : float
        Radius of search cone.
    source : string
        Define which catalog to query. Standard is 'PS' for Panstarrs DR1, alternatives: GAIA

    Returns
    -------
    objects
        table with the objects.

    """
    if(source == "PS" or source == "PANSTARRS" or source == "Panstarrs" or source == "PS1"):
        return get_PS_data(pos,radius)

    if(source == "GAIA" or source == "GAIADR1"):
        return get_GAIA_data(pos,radius)
