#!/usr/bin/env python
#
# Syntax: python getcircSB1.py <input fits file> <x> <y>
#                              <output flux file> <output radial flux file>
#     Calculates surface brightness as function of source radius
# This version is for one source and uses an approximate (x,y) as input.
#
import sys
import os
import math
import numpy as np
from photutils import DAOStarFinder
import photutils.detection as detection
from photutils.morphology import data_properties
from photutils import CircularAperture, CircularAnnulus
from photutils import aperture_photometry as apphot
from photutils.centroids import centroid_sources, centroid_com
from image_interpolation import fix_bad_pixels_surface, pixel_interpolation
from astropy import units
from astropy import stats
from astropy import units as u
from astropy.io import fits
from astropy.table import QTable, Table, Column
from scipy.ndimage.filters import median_filter
from astropy.nddata import Cutout2D

#------------- DEFINITIONS -----------------------------------------

def sb_circ(image, dqarr, x, y, radtab, radii, inttime, dqflag, skip=5.,
              w_ann=2., pixscale=0.0656, gain=1.6, badclip=True):
    # Refine the center coordinates using the 2D image moments
    newx, newy = centroid_sources(image, round(x), round(y), box_size=7,
                                        centroid_func=centroid_com)
    coords = np.column_stack((newx,newy))
    if dqflag:
    # Check for saturated pixels and bad DQ values in the innermost annulus
    # to avoid spurious results. Just to check the source is OK.
        srcaper = CircularAnnulus(coords, r_in=radii[0]-w_ann, r_out=radii[0]+w_ann)
        srcaper_masks = srcaper.to_mask(method='center')
        satflag = np.zeros((len(newx),), dtype=int)
        i = 0
        for mask in srcaper_masks:
            srcaper_dq = mask.multiply(dqarr)
            srcaper_dq_1d = srcaper_dq[mask.data > 0]
            badpix0 = np.logical_and(srcaper_dq_1d >= 2, srcaper_dq_1d <= 3)
            badpix1 = np.logical_and(srcaper_dq_1d >= 5, srcaper_dq_1d < 4095)
            badpix2 = np.where(srcaper_dq_1d >= 8388608)
            donotuse = np.where(srcaper_dq_1d == 1)
            nbadpix0 = len(srcaper_dq_1d[badpix0])
            nbadpix1 = len(srcaper_dq_1d[badpix1])
            nbadpix2 = len(srcaper_dq_1d[badpix2])
            nbadpix = nbadpix0 + nbadpix1 + nbadpix2
            ndonotuse = len(srcaper_dq_1d[donotuse])
            if ((nbadpix > 0) or (len(srcaper_dq_1d[donotuse]) > 0)):
                print('   {} sat, {} dead, {} unreliable, {} DONT USE pixels in source {}'.format(nbadpix0,nbadpix1,nbadpix2,ndonotuse,i))
                satflag[i] = 1
                i = i+1
        goodx = newx[np.where(satflag == 0)]
        goody = newy[np.where(satflag == 0)]
        #print('Number of sources without saturated or bad pixels: ', len(goodx))
        coords = np.column_stack((newx,newy))
    #
    # Get clipped statistics of background region
    # First eliminate bad pixels from the tally, using DQ array
    #
    maxrad = max(radii)
    bckaper = CircularAnnulus(coords, r_in=maxrad+2., r_out=maxrad+2.+skip)
    bckaper_masks = bckaper.to_mask(method='center')
    bck_mean = []
    bck_stdev = []
    bck_npix = []
    badlimits = []
    for mask in bckaper_masks:
        bckaper_data = mask.multiply(image)
        bckaper_data_1d = bckaper_data[mask.data > 0]
        if dqflag:
            bckaper_dq = mask.multiply(dqarr)
            bckaper_dq_1d = bckaper_dq[mask.data > 0]
            okpix1 = np.logical_or(bckaper_dq_1d == 0, bckaper_dq_1d == 4)
            okpix2 = np.logical_and(bckaper_dq_1d >= 4096, bckaper_dq_1d < 4194304)
            goodpix = np.logical_or(okpix1, okpix2)
            thebckaper = bckaper_data_1d[goodpix]
        else:
            thebckaper = bckaper_data_1d
        mean_sigclip, _, sig_sigclip = stats.sigma_clipped_stats(thebckaper, sigma=3.0)
        clipdata = stats.sigma_clip(thebckaper, sigma=3.0, maxiters=10)
        bck_mean.append(mean_sigclip)
        bck_stdev.append(sig_sigclip)
        bck_npix.append(len(np.where(clipdata.mask == False)[0]))
        badlimits.append(mean_sigclip-2.*sig_sigclip)
    bck_mean = np.array(bck_mean)
    bck_stdev = np.array(bck_stdev)
    bck_npix = np.array(bck_npix)
    #
    # Now interpolate over badly negative pixels in the
    # source area, before doing the photometry. Use stats in background area to
    # identify those pixels (i.e., array 'badlimits' above):
    #
    srcaper = CircularAnnulus(coords, r_in=1., r_out=15.)
    srcaper_masks = srcaper.to_mask(method='center')
    #
    i = 0
    for mask in srcaper_masks:
        srcaper_data = mask.multiply(image)
        badpixels = np.where(srcaper_data < badlimits[i])
        if len(badpixels[i]) > 0:
            newimage = fix_bad_pixels_surface(image, badpixels)
            #print(len(newimage))
            #print(len(newimage[0][0]), ' ', len(newimage[0][1]))
            image = newimage[0]
            #print(len(image))
            #print(len(image[0]), ' ', len(image[1]))
        i = i+1
    #
    # Now measure signal level in (by default) 5-pix wide annuli
    # By default, exclude pixel values more than 3 sigma below background value.
    # Selected as such if "noclip" is NOT provided as 5th parameter.
    #
    sb = []
    errsb = []
    for rad in radii:
        srcaper = CircularAnnulus(coords, r_in=rad-w_ann, r_out=rad+w_ann)
        srcaper_masks = srcaper.to_mask(method='center')
        for mask in srcaper_masks:
            srcaper_data = mask.multiply(image)
            srcaper_data_1d = srcaper_data[mask.data > 0]
            if badclip:
                srcaper_goodpix = np.where(srcaper_data_1d > badlimits[0])
                srcaper_thedata = srcaper_data_1d[srcaper_goodpix]
                totsig = np.sum(srcaper_thedata)
            else:
                totsig = np.sum(srcaper_data_1d)
            if totsig > 0.:
                devsig = np.sqrt(totsig)
            else:
                devsig = 0.
            if badclip:
                sigperpix = totsig/len(srcaper_thedata)
                errperpix = devsig/len(srcaper_thedata)
            else:
                sigperpix = totsig/srcaper.area
                errperpix = devsig/srcaper.area
            sb.append(sigperpix)
            errsb.append(errperpix)
    num_sb = np.array(sb)
    num_errsb = np.array(errsb)
    #
    # Now write new table: aperture radius vs. sb and its error
    nradii = np.array(radii)
    rad_tab = Table([nradii, num_sb, num_errsb],
                    names=['Radius', 'SB', 'errSB'])
    rad_tab['Radius'].info.format = '%6.0f'
    rad_tab['SB'].info.format = '11.5f'
    rad_tab['errSB'].info.format = '11.5f'
    rad_tab.write(radtab, format='ascii', overwrite=True)

#--------------------------------------------------------------------
# Main script starts below
#--------------------------------------------------------------------

if len(sys.argv) < 4: 
    print(' Input FITS file name, x, and y are required.')
    sys.exit()

infile = str(sys.argv[1])
x = float(sys.argv[2])
y = float(sys.argv[3])
if len(sys.argv) > 4:
    radtab = str(sys.argv[4])
else:
    radtab = 'getcircSB1.radtab'
if len(sys.argv) > 5:
    clipstr = str(sys.argv[5]).lower()
    noclipstr = "no"
    if clipstr != None and noclipstr in clipstr:
        doclip = False
    else:
        doclip = True
else:
    doclip = True

# Divide science extension (in units of MJy/sr) by the number of MJy/sr producing 1 cps.
# Can only be done right now if the photom step was performed on the input file.
numext = len(fits.open(infile))
if numext > 1:
    image = fits.getdata(infile, ext=1)
    dqarr = fits.getdata(infile, ext=3)
    hdr0 = fits.getheader(infile, ext=0)
    hdr1 = fits.getheader(infile, ext=1)
    do_dqflag = True
    inttime = hdr0["EFFINTTM"]
    ngroups = hdr0["NGROUPS"]
    nints = hdr0["NINTS"]
    pupil = hdr0["PUPIL"]
    filter = hdr0["FILTER"]
    subarray = hdr0["SUBARRAY"]
    # If input image is _cal.fits, then divide by PHOTMJSR value
    if 'PHOTMJSR' in hdr1:
        photmjsr = hdr1['PHOTMJSR']
        image = image/photmjsr
else:
    image = fits.getdata(infile)
    hdr0 = fits.getheader(infile)
    # For now just set it to 1, need to know actual header keyword for a given instrument
    #inttime = hdr0["EXPTIME"]
    inttime = 1.
    do_dqflag = False
    dqarr = np.uint32(np.array(image)*0)
if hdr0 is None:
    print("***** ERROR opening {}, skipping *****".format(infile))
    sys.exit()
# Extract 7x7 subarray around input x,y and get max. count rate and counts per int.
stamp = image[round(y)-4:round(y)+3,round(x)-4:round(x)+3]
maxrate = np.max(stamp)
maxint = round(maxrate*inttime)
if numext > 1:
    print('{:39s}: {:6s}, PW: {:6s}, FW: {:5s}, NGROUPS: {:3d}, NINTS: {:2d}, Max. counts: {:5d}'.format(infile,subarray,pupil,filter,ngroups,nints,maxint))
else: 
    print('  Max. counts: {:5d}'.format(maxint))
# Define multiple measurement radii here
myradii = [5.,9.,13.,17.,21.,25.]
# Define pixel scale in arcsec here
mypixscale = 0.0656
# Define gain factor here
mygain = 1.6

#source_list = find_sources(image, nsigma=nsig, roundness=1.0, sharpness=0.1)
sb_circ(image, dqarr, x, y, radtab, myradii, inttime, do_dqflag, skip=5.,
          w_ann=2., pixscale=mypixscale, gain=mygain, badclip=doclip)
