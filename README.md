# Astrometry

Simple python3 tool to quickly correct the rough astrometry given by a telescope for a fits image. For the calibration of the position Panstarrs Dr1 or GAIA data can be used. The program considers scaling, rotation, reflection and translation. To start with it needs a rough position from the fits header or specified as input paramters.

general notes: If already present in the header wcs (World Coordinate System) information is used. If not the wcs is build from scratch. If not specified genomic projection is assumed.

The resulting wcs is saved to the fits header and the data is saved as a new file <filename>_astro.fits.



## How to install it

### Prerequisites

The program was tested for python 3.6 or newer
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

With this you are all set. If you want to run astrometry from anywhere just install it as a package as follows*

```
pip install -e PATH/TO/CLONED/GITHUB
```

## How to use it

Without installing you have to be in the cloned git folder and then you can run

```
python astrometry.py sample_images/sample_file.fits
```

This will perform an astrometric calibration on the sample file. The result will be stored as FILENAME_astro.fits. 

*If you install the package you can use it anywhere like this:

```
astrometry sample_images/sample_file.fits
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

If the header is missing a rough position you can directly input the position as follows. Note that ra and dec have to be in degrees

```
astrometry sample_images/sample_file.fits -ra 16.65733  -dec 3.54336
```

If the header is missing the projection you can specify it directly as follows. Note that None is found it assumes genomic projection (TAN). If this is incorrect the fit will fail. Also do NOT put the RA---TAN within quotation marks.

```
astrometry sample_images/sample_file.fits -proj1 RA---TAN -proj2 DEC--TAN
```

All of these also work for python astrometry.py of course.
The full list of parameters can be accessed with astrometry --help

## Example

Input:

![Input read from file and loaded from online catalog](sample_files/sample_file_input.gif)

Result:

![Result](sample_files/sample_file_result.gif)

## Author

Lukas Wenzl 
