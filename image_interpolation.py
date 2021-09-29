#! /usr/bin/env python
#
"""
Utilities for interpolation of an image to replace bad pixels.  A few options
are provided for this type of activity.
Code by Kevin Volk along with some changes by PG
"""
import numpy
import scipy.stats


def fix_bad_pixels_surface(image, pixels, box=3, degree=3):
    """
    Interpolate to a set of bad pixels with a polynominal.

    Adapted from IDL code by David Lafreniere.

    Parameters
    ----------

    image:     A two-dimension real numpy array, the image values

    pixels:    A list of integer pixel positions, from nump.where for example,
               length 2 with y positions in 0 and x positions in 1,
               with each list element being a numpy array of integers

    box:       Optional parameter, the size of the polynomial cut-out, an
               integer value.  The cut-out is +/- this number around the
               flagged pixel, so the cut-out is 2*box+1 pixels square

    degree:    Option integer value for the polynomial order for the surface
               fitting

    Returns
    -------

    newimage:  A numpy float image of the same dimensions as image with the
               flagged pixels interpolated if possible; if there is an error
               the original image is returned

    flag:   A boolean value for whether the fitting succeeded or not.

    """
    try:
        workimage = numpy.copy(image)
        sh1 = workimage.shape
        dimx = sh1[1]
        npar = (degree+1)*(degree+1)
        maskimage = numpy.zeros_like(image, dtype=numpy.int16)
        maskimage[pixels] = 1
        goodpixels = numpy.where(maskimage == 0)
        xvalues = numpy.linspace(0, 1, 2*box+1)
        yvalues = numpy.linspace(0, 1, 2*box+1)
        xs, ys = numpy.meshgrid(xvalues, yvalues, copy=False)
        x = xs.flatten()
        y = ys.flatten()
        basis = []
        tdbasis = []
        for ind1 in range(degree+1):
            for ind2 in range(degree+1):
                if ind1 == 0:
                    x1 = x*0+1
                    xs1 = xs*0 + 1
                else:
                    x1 = numpy.power(x, ind1)
                    xs1 = numpy.power(xs, ind1)
                if ind2 == 0:
                    y1 = y*0+1
                    ys1 = ys*0 + 1
                else:
                    y1 = numpy.power(y, ind2)
                    ys1 = numpy.power(ys, ind2)
                basis1 = x1*y1
                basis2 = xs1*ys1
                basis.append(basis1)
                tdbasis.append(basis2)
        basis = numpy.array(basis)
        tdbasis = numpy.array(tdbasis)
        sh2 = basis.shape
        for loop in range(len(pixels[0])):
            xpos = pixels[1][loop]
            ypos = pixels[0][loop]
            subimage = numpy.copy(
                workimage[ypos-box:ypos+box+1, xpos-box:xpos+box+1])
            submask = numpy.copy(
                maskimage[ypos-box:ypos+box+1, xpos-box:xpos+box+1])
            dvector = subimage.flatten()
            mvector = submask.flatten()
            inds = numpy.where(mvector == 0)
            funct = []
            for ind1 in range(sh2[0]):
                funct.append(basis[ind1][inds])
            funct = numpy.asarray(funct)
            funct = numpy.transpose(funct)
            data = subimage.flatten()
            data = data[inds]
            coeff, residuals, rank, singular = numpy.linalg.lstsq(funct, data, rcond=None)
            if numpy.min(coeff) != numpy.max(coeff):
                fitimage = subimage*0.
                for ind1 in range(sh2[0]):
                    fitimage = fitimage + \
                        coeff[ind1]*numpy.squeeze(tdbasis[ind1,:,:])
                newinds = numpy.where(submask == 0)
                fitimage[newinds] = subimage[newinds]
                workimage[
                    ypos-box:ypos+box+1, xpos-box:xpos+box+1] = numpy.copy(
                        fitimage)
        return workimage, True
    except:
        return numpy.copy(image), False


def pixel_interpolation(image, xpixel, ypixel):
    """
    Replace a pixel value by the median of the adjacent 4 values.

    A pixel at the edge of the array is not replaced by the median.

    Parameters
    ----------

    image:   A two-dimensional numpy real or integer array

    xpixel:  An integer value, the x pixel index

    ypixel:  An integer value, the y pixel index

    Returns
    -------

    newimage:   A copy of image with the single pixel value corrected, or the
                original image if the pixel position is at the edge of the
                image and does not have four neighbours

    """
    try:
        values = []
        values.append(image[ypixel-1, xpixel-1])
        values.append(image[ypixel-1, xpixel+1])
        values.append(image[ypixel+1, xpixel+1])
        values.append(image[ypixel+1, xpixel-1])
    except (IndexError, ValueError):
        return image
    values = numpy.asarray(values)
    replace = numpy.median(values)
    newimage = numpy.copy(image)
    newimage[ypixel, xpixel] = replace
    return newimage


def mask_adjacent(maskimage, xpixel, ypixel, npixel=1):
    """
    Given an integer mask return the most common adjacent value.

    If there is an issue with looking at the adjacent values, return the
    current mask value at the pixel of interest.  If the pixel is out of
    bounds, return 0.

    Parameters
    ----------

    maskimage:   A two-dimension numpy array

    xpixel:      An integer value, the x position of the target pixel

    ypixel:      An integer value, the y position of the target pixel

    npixel:      An optional integer parameter, the size of the region
                 to look at to see what the most common value is

    Returns
    -------

    maskvalue:   The most common value in the adjacent pixels

    Note:  This is not likely to work for a flat array and is assumed to
           be applied to an integer array so that there is a clear majority
           value (one hopes) amoung the adjacent pixels.
    """
    try:
        subimage = numpy.copy(maskimage[ypixel-npixel:ypixel+npixel+1,
                                        xpixel-npixel:xpixel+npixel+1])
    except (IndexError, ValueError):
        try:
            return maskimage[ypixel, xpixel]
        except (IndexError, ValueError):
            return 0
    maskvalue, counts = scipy.stats.mode(subimage, axis=None)
    # return the most common mask value in the adjacent pixels
    return maskvalue[0]
