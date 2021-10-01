#!/usr/bin/env python
#
# Syntax: python meansb1a.py <root name of flux tables> <outtab>
#     Specifically for SB tables with one (and the same) source each, 
#     measured using a given set of measurement radii.
#     Goal is to get average SB and error.
#     parameter 'outtab' is defaulted to <input root name>.avsb
# Example: python meansb1.py f115w_sub64 
#
import sys
import os
import math
import glob
import numpy as np
from astropy import stats
from astropy import units as u
from astropy.io import ascii
from astropy.table import QTable, Table, Column

#------------- DEFINITIONS -----------------------------------------
if len(sys.argv) < 2: 
    print(' ERROR: Input root name is required.')
    sys.exit()
inroot = str(sys.argv[1])
infiles = sorted(glob.glob(inroot+'*.sbtab'))
if (len(infiles) < 1) or (infiles is None):
    print(' ERROR: no candidate input files are present in the directory.')
    sys.exit()

if len(sys.argv) > 2: 
    outtab = str(sys.argv[2])
else:
    outtab = inroot+'.avsbtab'

t = Table.read(infiles[0], format='ascii')

sbarr = np.zeros((len(infiles), len(t)))
sberrarr = np.zeros((len(infiles), len(t)))
sbsum = np.zeros(len(t))
wsbsum = np.zeros(len(t))
totweight = np.zeros(len(t))
avsb = np.zeros(len(t))
avsberr = np.zeros(len(t))

for i in range(len(infiles)):
    tabfile = ascii.read(infiles[i],
                         names=['rad','SBMean','SBError'])
    if i == 0:
        pixradii = tabfile['rad']
    sbarr[i,:] = tabfile['SBMean']
    sberrarr[i,:] = tabfile['SBError']
# For entries with SBError = 0, assign them an error of 1000.
thesberrarr = sberrarr
badsbs = np.where(sberrarr <= 0.)
thesberrarr[badsbs] = 1000.
# Now calculate "weighted average" fluxes, by 
# using inverse variances.
for i in range(len(infiles)):
    sbsum = sbsum + sbarr[i,:]
    totweight = totweight + (thesberrarr[i,:])**(-2.)
    wsbsum = wsbsum + (thesberrarr[i,:])**(-2.)*sbarr[i,:]
# Calculate weighted average & median sb values
for j in range(len(t)):
    avsb[j] = wsbsum[j]/totweight[j]
    if totweight[j] > 0.1:
        avsberr[j] = 1./np.sqrt(totweight[j])
    else:
        avsberr[j] = np.std(sbarr[:,j])

# table with radii, fluxes, flux errors, and flags.
# Add data for all dithers as additional columns.
tab = Table([pixradii], names=['Rad_Pix'])
tab['Rad_Pix'].info.format = '%7.0f'
for i in range(len(infiles)):
    tab.add_column(sbarr[i,:], name='sb'+str(i))
    tab.add_column(sberrarr[i,:], name='sberr'+str(i))
    tab['sb'+str(i)].info.format = '10.5f'
    tab['sberr'+str(i)].info.format = '10.5f'
tab.add_column(avsb, name='av_sb')
tab.add_column(avsberr, name='av_sberr')
tab['av_sb'].info.format = '10.5f'
tab['av_sberr'].info.format = '10.5f'

tab.write(outtab, format='ascii.fixed_width', delimiter=None, overwrite=True)
