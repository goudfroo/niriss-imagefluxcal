#!/usr/bin/env python
#
# Syntax: python meanflux1.py <root name of flux tables> <outtab>
#     Specifically for flux tables with one (and the same) source each, 
#     measured using a given set of measurement radii.
#     Goal is to get average flux and error.
#     parameter 'outtab' is defaulted to <input root name>.avtab
# Example: python meanflux1.py f115w_sub64 
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
allfiles = sorted(glob.glob(inroot+'*.radtab'))
infiles = [f for f in allfiles if not 'psf' in f]

if (len(infiles) < 1) or (infiles is None):
    print(' ERROR: no candidate input files are present in the directory.')
    sys.exit()

if len(sys.argv) > 2: 
    outtab = str(sys.argv[2])
else:
    outtab = inroot+'.avtab'

t = Table.read(infiles[0], format='ascii')

fluxarr = np.zeros((len(infiles), len(t)))
fluxerrarr = np.zeros((len(infiles), len(t)))
flagarr = np.ones((len(infiles), len(t)))
avflux = np.zeros(len(t))
avfluxvar = np.zeros(len(t))
avfluxerr = np.zeros(len(t))
avflag = np.ones(len(t))
fluxsum = np.zeros(len(t))
wfluxsum = np.zeros(len(t))
plainavg = np.zeros(len(t))
flagsum = np.zeros(len(t))
totweight = np.zeros(len(t))

for i in range(len(infiles)):
    tabfile = ascii.read(infiles[i], names=['rad','rad_as','flux','fluxerr'])
    if i == 0:
        pixradii = tabfile['rad']
        asradii = tabfile['rad_as']
    fluxarr[i,:] = tabfile['flux']
    fluxerrarr[i,:] = tabfile['fluxerr']
# Now enter zero flags for cases where flux(n) < flux(n-1)
# (i.e., due to pixels with negative values)
for i in range(len(infiles)):
    for j in range(len(t)):
        if j > 0:
            if fluxarr[i,j] <= fluxarr[i,j-1]:
                flagarr[i,j] = 0
# Now calculate "weighted average" fluxes, by ignoring flagged entries
# and using inverse variances.
# Also calculate "plain" average fluxes to cover cases where all 
# individual measurements have flux(n) < flux(n-1).
for i in range(len(infiles)):
    fluxsum = fluxsum + fluxarr[i,:]*flagarr[i,:]
    plainavg = plainavg + fluxarr[i,:]/len(infiles)
    flagsum = flagsum + flagarr[i,:]
    totweight = totweight + (fluxerrarr[i,:])**(-2.)*flagarr[i,:]
    wfluxsum = wfluxsum + (fluxerrarr[i,:])**(-2.)*fluxarr[i,:]*flagarr[i,:]
# Enter avflag = 0 for cases where all measurements have flux(n) < flux(n-1)
for j in range(len(t)):
    if sum(flagarr[:,j]) == 0:
        avflag[j] = sum(flagarr[:,j])
        avflux[j] = plainavg[j]
        avfluxerr[j] = np.sqrt(sum((fluxerrarr[:,j])**2))/len(infiles)
    else:
        avflag[j] = sum(flagarr[:,j])
        avflux[j] = wfluxsum[j]/totweight[j]
        avfluxerr[j] = 1./np.sqrt(totweight[j])

# table with radii, fluxes, flux errors, and flags.
# Add data for all dithers as additional columns.
tab = QTable([pixradii, asradii], 
                   names=('Rad_Pix', 'Arcsec'))
tab['Rad_Pix'].info.format = '%7.0f'
tab['Arcsec'].info.format = '10.4f'
for i in range(len(infiles)):
    tab.add_column(fluxarr[i,:], name='Flux'+str(i))
    tab.add_column(fluxerrarr[i,:], name='Fluxerr'+str(i))
    tab.add_column(flagarr[i,:], name='Flag'+str(i))
    tab['Flux'+str(i)].info.format = '13.4e'
    tab['Fluxerr'+str(i)].info.format = '13.4e'
    tab['Flag'+str(i)].info.format = '5.0f'
tab.add_column(avflux, name='avFlux')
tab.add_column(avfluxerr, name='avFluxerr')
tab.add_column(avflag, name='n_Good')
tab['avFlux'].info.format = '13.4e'
tab['avFluxerr'].info.format = '13.4e'
tab['n_Good'].info.format = '5.0f'

tab.write(outtab, format='ascii', overwrite=True)
