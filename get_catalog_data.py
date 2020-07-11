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

from astropy.coordinates import SkyCoord
from astropy import units as u
import pandas as pd

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
    """Query Vizier for Panstarrs1 data. PS1 was calibrated with GAIA and goes darker so if it is availabe it is likely the best choice.

    Parameters
    ----------
    coord  : SkyCoord
        Position to search.
    radius : float*unit
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


def get_PS_photometry_data(ra=0, dec=0, radius=60, coord=None, band="z"):
    """Query Vizier for Panstarrs1 data. This method is meant to get the color info for nearby stars

    Parameters
    ----------
    ra : float
        ra to search
    dec : float
        dec to search
    radius : float
        Radius of search cone in arcsec.
    coord  : SkyCoord
        Position to search. Use instead of ra and dec
    band : str
        photometric band that is needed. PanStarrs has g,r,i,z,y

    Returns
    -------
    objects
        table with the objects including the following info: ra (deg), ra_error (milliarcsec), dec (deg), dec_error (milliarcsec), mag
    band
        photometric band with the info required

    """
    dict = {"g":"gmag", "r":"rmag", "i":"imag","z":"zmag","y":"ymag"}
    band_name = dict[band]
    #coord = SkyCoord(ra=ra*u.deg,dec=dec*u.deg, frame="icrs")
    radius = u.Quantity(radius, u.arcsec)
    if(coord == None):
        coord = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame="icrs")
    Vizier.ROW_LIMIT = -1
    j = Vizier.query_region(coord, radius=radius, catalog="II/349/ps1")
    if (j==[]):
        return None
    r = j[0]


    catalog_data = r.to_pandas()
    catalog_data.rename(index=str, columns={"RAJ2000": "ra", "e_RAJ2000":"ra_error", "DEJ2000": "dec", "e_DEJ2000":"dec_error"}, inplace=True)
    #columns: ra, ra_error, dec, dec_error, mag
    #
    print("Found {} sources in PS1 within a radius of {}".format(catalog_data.shape[0], radius))
    #

    return catalog_data, band_name, "PS1", "AB"



def get_2MASS_data(coord, radius):
    """Query Vizier for 2MASS astronometry data.

    Parameters
    ----------
    coord  : SkyCoord
        Position to search. Use instead of ra and dec
    radius : float*unit
        Radius of search cone in arcsec.

    Returns
    -------
    objects
        table with the objects including the following info: ra (deg), ra_error (milliarcsec), dec (deg), dec_error (milliarcsec), mag
    band
        photometric band with the info required

    """

    #coord = SkyCoord(ra=ra*u.deg,dec=dec*u.deg, frame="icrs")
    Vizier.ROW_LIMIT = -1
    j = Vizier.query_region(coord, radius=radius, catalog="II/246/out")
    if (j==[]):
        return None
    r = j[0]


    catalog_data = r.to_pandas()
    magnitudes = ["Jmag", "Kmag", "Hmag"]
    catalog_data["mag"] = catalog_data[magnitudes].min(axis=1) #there should not be an issue with -999 values since in Vizier nulls are used fro missing values not -999
    catalog_data.rename(index=str, columns={"RAJ2000": "ra", "DEJ2000": "dec"}, inplace=True)
    catalog_data = catalog_data[["ra", "dec", "mag"]]
    print("WARNING: 2MASS has no error for the positions.")
    #columns: ra, ra_error, dec, dec_error, mag
    #
    print("Found {} sources in 2MASS within a radius of {}".format(catalog_data.shape[0], radius))
    #

    return catalog_data


def get_2MASS_photometry_data(ra=0, dec=0, radius=60, coord=None, band="z"):
    """Query Vizier for Panstarrs1 data. This method is meant to get the color info for nearby stars

    Parameters
    ----------
    ra : float
        ra to search
    dec : float
        dec to search
    radius : float
        Radius of search cone in arcsec.
    coord  : SkyCoord
        Position to search. Use instead of ra and dec
    band : str
        photometric band that is needed. PanStarrs has g,r,i,z,y

    Returns
    -------
    objects
        table with the objects including the following info: ra (deg), ra_error (milliarcsec), dec (deg), dec_error (milliarcsec), mag
    band
        photometric band with the info required

    """
    dict = {"J":"Jmag", "H":"Hmag", "K":"Kmag"}
    dict_error = {"J":"e_Jmag", "H":"e_Hmag", "K":"e_Kmag"}
    band_name = dict[band]
    #coord = SkyCoord(ra=ra*u.deg,dec=dec*u.deg, frame="icrs")
    radius = u.Quantity(radius, u.arcsec)
    if(coord == None):
        coord = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame="icrs")
    Vizier.ROW_LIMIT = -1
    j = Vizier.query_region(coord, radius=radius, catalog="II/246")
    if (j==[]):
        return None
    r = j[0]


    catalog_data = r.to_pandas()
    catalog_data.rename(index=str, columns={"RAJ2000": "ra", "DEJ2000": "dec"}, inplace=True)
    #columns: ra, ra_error, dec, dec_error, mag
    #
    print("Found {} sources in 2MASS within a radius of {}".format(catalog_data.shape[0], radius))
    #

    return catalog_data, band_name, "2MASS", "VEGA"

def get_file_data(filename):
    """Read catalog data from local file. At the minimum the columns ra, dec and mag are needed.
    The magnitudes do not need to be calibrated only ordered.
    """
    catalog_data = pd.read_csv(filename)
    print("Found {} sources in file: {}".format(catalog_data.shape[0], filename))
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

    if(source == "PS_photometry"):
        return get_PS_data(pos,radius)

    if(source == "GAIA" or source == "GAIADR1"):
        return get_GAIA_data(pos,radius)

    if(source == "2MASS" or source == "TWOMASS" or source == "2mass" or source =="twomass"):
        return get_2MASS_data(pos, radius)

    print("trying to read local file")
    return get_file_data(source)

def get_photometry_data(pos, radius, band, source="auto"):
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
    if(source == "auto"):
        dict_source = {"g":"PS", "r":"PS", "i":"PS", "z":"PS", "y":"PS", "J":"2MASS", "H":"2MASS", "K":"2MASS"}
        source = dict_source[band]

    if(source == "PS"):
        return get_PS_photometry_data(radius=radius, coord=pos, band=band)

    if(source == "2MASS"):
        return get_2MASS_photometry_data(radius=radius, coord=pos, band = band)
