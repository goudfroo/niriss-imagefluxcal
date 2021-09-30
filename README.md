# NIRISS Alternative Analysis for CAP-020 (NIRISS Imaging Flux Calibration) Workflow
This repository contains scripts that can be used for determining **Encircled Energy** and **Integrated Intensity** measurements from a set of (dithered) JWST images of a flux standard star that were taken with different filters.

## Installation
This package has some external dependencies. Specifically, it currently requires the packages NumPy, SciPy, AstroPy, PySynphot and Photutils. It also uses functions defined in `image_interpolation.py` that was taken from the [niriss-commissioning repo](https://github.com/spacetelescope/niriss-commissioning). Python version 3.6 or higher is required; version 3.8 or later is preferred but not required if one selects the default settings of the scripts in this package. This package will work fine within the conda environment used for the general [niriss-commissioning repo](https://github.com/spacetelescope/niriss-commissioning) or a clean conda environment in which the packages mentioned above have been installed using pip install. The instructions below assume that one of those environments have been installed on your machine and that you have loaded into that environment.

# Usage
Below are brief descriptions on the use of the main scripts.

(1) `getcircEE1.py` - This script performs multi-aperture photometry of the (single) flux standard star in the input image. The measurement radii are currently defined near the bottom of the script in parameter "myradii". It first uses `photutils` functions to find the accurate centroid of the source and lists how many pixels in the innermost 5 pixels of the source area are saturated or flagged "bad" by the pipeline (this is just for info - pixel values are not changed by the script in this central region of the PSF). In the outer regions of the source area (between 5 pixels and the outermost measurement radius), the script allows two ways to identify and/or deal with bad pixels. By default, the script identifies pixels flagged as "bad" by the pipeline and replaces their values by a polynomial surface fitted to the surrounding pixels. The other option, which is useful for cases where the detector noise is significant, is selected by calling the script with "usemedian" as 6th input parameter (see below). In this mode, the script also identifies pixels with values lower than 3 sigma below the clipped mean sky background level, and replaces those by the same pixels in a median-filtered version of the image. The script is run as follows:

```
python getcircEE1.py <input FITS file> <x> <y> <output flux table> <output EE table> <usemedian>
```

where the first three parameters are required: the input FITS file (a _rate.fits or _cal.fits file), and <x> and <y> which are the approximated pixel coordinates of the source on the image. It also takes three optional parameters, namely (1) an output file name for an ASCII table with all the photometric measurements (defaulted to "getcircEE1.fluxes"); (2) an output file name for an ASCII table with encircled energy measurements (defaulted to "getcircEE1.radtab"). It is important to use a common file root name for all exposures taken with a given combination of filter and subarray (an example root name is "f090w_sub128") for the next steps; and (3) when the text "usemedian" is used as 6th input parameter, the script will replace anomalously low pixel values in the outer source area by median-filtered values (see above). A wrapper shell script `do_getcircEE1` is included in this repo to provide an example of running `getcircEE1.py` on a set of data files. 


(2) `meanflux1.py` - This script takes the root name of the EE output tables produced by `getcircEE1.py` as input. It looks for all EE tables with that root name (in the case of the data for this CAP, there two such tables for each passband), and proceeds to calculate weighted average EE measurements and their uncertainties. This script also flags entries for measurement radii at which the EE is smaller than at the previous (smaller) measurement radius, which would indicate a problem with the data. In case all (2 in this case) EE measurements at a given measurement radius are flagged this way, the script calculates non-weighted average fluxes and their uncertainties, along with an "average flag" of 0. This script is run as follows:

```
python meanflux1.py <root name of EE tables> <outtab>
```

where "outtab" is defaulted to "(root name of EE tables).avtab". A wrapper shell script `do_meanflux1` is included in this repo to provide an example of running `meanflux1.py` on a set of data files. 


(3) `getcircSB1.py` - This script is similar to `getcircEE1.py`, except that it calculates surface brightness (SB), rather than encircled energy, in a set of concentric annuli around the source. The radii at the centers of the annuli are currently defined near the bottom of the script in parameter "myradii". By default, the script flags pixels with values lower than 3 sigma below the background level established far away from the source, and does not take those pixels into account in the surface brightness calculation. The script is run as follows:

```
python getcircSB1.py <input FITS file> <x> <y> <output SB table> <noclip>
```

where parameter <output SB table> is defaulted to "getcircSB1.radtab", and optional parameter "noclip", when used, disables the flagging of pixels with values more than 3 sigma below the background value. This script is useful for checking the impact of possibly noisy pixels in the outer regions of the PSF on the EE profiles measured with `getcircEE1.py`. A wrapper shell script `do_SB1` is included in this repo to provide an example of running `getcircSB1.py` on a set of data files. 


(4) `meanSB1.py` - This script is very similar to `meanflux1.py`, except that it now calculates average SB values as a function of measurement radius. The syntax is the same as for `meanflux1.py`:

```
python meanSB1.py <root name of SB tables> <outtab>
```

where "outtab" is defaulted to "<root name of SB tables>.avsbtab". A wrapper shell script `do_meanSB1` is included in this repo to provide an example of running `meanSB1.py` on a set of data files. 


(5) `totfluxes.py` - This script combines all EE results for a given subarray (and a given flux standard star) and calculates a total integrated count rate for all NIRISS filters that were used for the measurements, using aperture corrections that were derived from WebbPSF PSFs (the EE tables for the latter are provided in this repo as files *psf.radtab). In the current version, these aperture corrections were only derived from measurement radii between 5-10 pixels from the center of the PSF. Depending on the results seen during commissioning, the aperture correction part of this script may need to be adjusted. The script is run as follows:

```
python totfluxes.py <subarray name> <outtab>
```

where "outtab" is an output ASCII-format table, defaulted to "totfluxes.out". The resulting count rates and their uncertainties are also displayed on the terminal.


(6) `NIRISScountrates.py` - This pysynphot script calculates expected total NIRISS imaging count rates for a given flux standard in the `calspec` directory (/grp/hst/cdbs/calspec at STScI, also available at [this web address](https://archive.stsci.edu/hlsps/reference-atlases/cdbs/current_calspec/) in all NIRISS filters). For this script to work correctly, the user needs to have the NIS_SYNPHOT environment variable point to a directory that contains all NIRISS-related throughput files for (py)synphot. This directory on STScI Central Store is /ifs/jwst/wit/niriss/pysynphot (i.e., `export NIS_SYNPHOT="/ifs/jwst/wit/niriss/pysynphot/"`) One must also have the PYSYN_CDBS environment variable defined correctly (see [pysynphot docs](https://pysynphot.readthedocs.io)). The syntax of this script is as follows for the example of the current spectrum for WD1057+719 in the `calspec` directory:

```
python NIRISScountrates.py wd1057_719_stisnic_008.fits
```

(see the table on [this web page](https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/calspec.html) for the names of current calspec files of flux standards).


The document "NIS-020_Analysis_PG.pdf", which is included in this repo, further describes the use of these scripts and illustrates the results of running them on simulated data of flux standard star WD1057+719.


## Contributing

Any changes to code should be made on a feature branch. The workflow should be:

1. Create a branch
2. Commit changes on branch
3. Submit pull request for branch to be merged into master

Detailed instructions can be found in the 
[contribution guidelines](https://github.com/spacetelescope/niriss-commissioning/blob/master/CONTRIBUTING.md).

## Code of Conduct

All users are expected to adopt a 
[code of conduct](https://github.com/spacetelescope/niriss-commissioning/blob/master/CODE_OF_CONDUCT.md) 
that ensures a productive, respectful environment for all contributors and 
participants. Any issues or violations should be reported to conduct@stsci.edu. 

## Questions or Issues

Any code problems or questions should be noted by 
[opening an issue](https://github.com/spacetelescope/niriss-commissioning/issues).

If you have git questions, see 
[these resources](https://github.com/spacetelescope/niriss-commissioning/blob/master/CONTRIBUTING.md#Resources).


