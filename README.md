# Astrometry

Simple python3 tool to quickly correct the rough astrometry given by a telescope for a fits image. For the calibration of the position Panstarrs Dr1 or GAIA data can be used. The program considers translation within about an arcminute and tries 0, 90, 180, 270 degree rotations as well as reflections. To start with it needs the pixelscale and a rough position from the fits header.


## How to install it

### Prerequisites

Apart from standard packages this package needs astroquery, astropy and photutils

```
pip install astroquery astropy photutils pandas

```

Standard packages that are required: matplotlib pandas numpy

### Installing

To use astrometry just clone this repository to your local machine.

```
git clone https://github.com/lukaswenzl/astrometry.git
```

With this you are all set. If you want to run astroquery from anywhere just install it as a package as follows

```
pip install -e PATH/TO/CLONED/GITHUB
```

## How to use it

Without installing you have to be in the cloned git folder and then you can run

```
python astrometry.py sample_images/PSO016m06_K.fits
```

This will perform an astrometric calibration on the sample file. The result will be stored as FILENAME_astro.fits. 

If you install the package you can use it anywhere like this:

```
astrometry sample_images/PSO016m06_K.fits
```

This installation of the package also takes care of the required packages.


## Further options

Get a list of all options with

```
astrometry --help
```

You can give a list of filenames

```
astrometry sample_images/PSO016m06_K.fits sample_images/PSO016m06_K.fits
```

Or a directory where every fits file without the addition _astro will be used.

```
astrometry sample_images/
```






All of these also work for python astrometry.py ...


## Author

Lukas Wenzl 
