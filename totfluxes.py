import os
import sys
import math
import glob
import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table,QTable

if len(sys.argv) < 2: 
    print(' ERROR: Input subarray name is required.')
    sys.exit()
subarr = str(sys.argv[1])
if len(sys.argv) > 2:
    outtab = str(sys.argv[2])
else:
    outtab = 'totfluxes.out'

nis_subarrays = ['sub64', 'sub128', 'sub256', 'sub80', 'subtasoss', 'subtaami',
                 'subampcal', 'wfss64r', 'wfss64c', 'wfss128r', 'wfss128c',
                 'substrip96', 'substrip256']
if subarr.lower() in nis_subarrays:
    pass
else:
    raise IOError('Error: {} is not a NIRISS subarray.'.format(subarr))

filters = ['f090w', 'f115w', 'f140m', 'f150w', 'f158m', 'f200w', 'f277w', 'f356w',
           'f380m', 'f430m', 'f444w', 'f480m']
infiles = sorted(glob.glob(filters[0]+'_'+subarr+'*a.radtab'))
t = Table.read(infiles[0], format='ascii')
maxcntsarr = np.zeros(len(t))
errcntsarr = np.zeros(len(t))
totcounts = np.zeros(len(filters))
toterr = np.zeros(len(filters))
i = 0
for filter in filters:
    infiles = sorted(glob.glob(filter+'_'+subarr+'*a.radtab'))
    if (len(infiles) < 1) or (infiles is None):
        print('Filter ', filter, ': no candidate input files are present.')
        continue
    psffile = filter + '_psf.radtab'
    if os.path.isfile(psffile):
        pass
    else:
        print('PSF EE file {} not found. Continuing with next filter.'.format(psffile))
        continue
    # For each pair (or set) of input files, determine max(count rate) for each radius
    for j in range(len(infiles)):
        tabfile = ascii.read(infiles[j],
                             names=['rad','rad_as','CountRate','errCountRate'])
        if j == 0:
            pixradii = tabfile['rad']
            maxcntsarr = tabfile['CountRate']
            errcntsarr = tabfile['errCountRate']
        if j > 0:
            for k in range(len(t)):
                if tabfile['CountRate'][k] > maxcntsarr[k]:
                    maxcntsarr[k] = tabfile['CountRate'][k]
                    errcntsarr[k] = tabfile['errCountRate'][k]
    # Read in psf EE table and calculate 'total' count rates
    psftab = ascii.read(psffile, names=['rad', 'rad_as', 'Flux', 'errFlux'])
    maxcntsarr = maxcntsarr/psftab['Flux']
    errcntsarr = errcntsarr/psftab['Flux']
    # Take average of total count rates at radii 5-10 pix
    totcounts[i] = np.mean(maxcntsarr[6:-3])
    toterr[i] = np.max(errcntsarr[6:-3])
    # For filters in the FW: Multiply by 0.84 because the WebbPSFs include the
    # 0.84 reduction of the CLEARP PW element
    if i > 5:
        totcounts[i] = totcounts[i]*0.84
        toterr[i] = toterr[i]*0.84
    i = i+1
# Write results into table 
tab = Table([filters, totcounts, toterr], 
            names=('Filter', 'CountRate', 'CountRateErr'))
tab['CountRate'].info.format = '12.3e'
tab['CountRateErr'].info.format = '12.3e'
tab.write(outtab, format='ascii.fixed_format', delimiter=None, overwrite=True)
# Print out results by filter to Terminal. If results present for all filters,
# then print them after transposing (much more compact).
if len(totcounts) == len(filters):
    print('# WD 1057+719 Measured Count Rates')
    print('#   F090W     F115W     F140M     F150W     F158M     F200W     F277W'
          '     F356W     F380M     F430M     F444W     F480M')
    print('%.3e %.3e %.3e %.3e %.3e'
          ' %.3e %.3e %.3e %.3e %.3e %.3e'
          ' %.3e' % (totcounts[0], totcounts[1], totcounts[2], totcounts[3],
                     totcounts[4], totcounts[5], totcounts[6], totcounts[7],
                     totcounts[8], totcounts[9], totcounts[10], totcounts[11]))
    print('%.3e %.3e %.3e %.3e %.3e'
          ' %.3e %.3e %.3e %.3e %.3e %.3e'
          ' %.3e' % (toterr[0], toterr[1], toterr[2], toterr[3], toterr[4],
                     toterr[5], toterr[6], toterr[7], toterr[8],
                     toterr[9], toterr[10], toterr[11]))
else:
    print(tab)
