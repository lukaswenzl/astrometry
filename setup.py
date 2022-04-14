from setuptools import setup
setup(
    name='astrometry',
    version='1.5.0',
    author='Lukas Wenzl',
    description='Simple python3 tool to quickly correct the rough astronometry given by a telescope for a fits image.',
    entry_points={
        'console_scripts': [
            'astrometry=astrometry:main',
            'photometry=photometry:main'
        ]
    },
    install_requires = [ 'astroquery',
   'astropy',
   'photutils',
   'matplotlib',
   'pandas',
   'numpy',
   'scipy'
   ]
)
