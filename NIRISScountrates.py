import os
import sys
import pysynphot as S
from astropy.io import fits

if len(sys.argv) < 2: 
    print(' Input Calspec FITS file name is required.')
    sys.exit()

infile = str(sys.argv[1])

mandatory_env_vars = {"NIS_SYNPHOT", "PYSYN_CDBS"}
envdiff = mandatory_env_vars.difference(os.environ)
if len(envdiff) > 0:
    raise EnvironmentError(f'Failed because variables {envdiff} are not set.')

# Check existence of one file in each of the two directories that are
# environment variables, just making sure
wavefile = os.path.join(
    os.environ['NIS_SYNPHOT'], 'NIRISS_waveset.fits')
if os.path.isfile(wavefile):
    with fits.open(wavefile) as wvf:
        wave = wvf[1].data.field('wavelength')
else:
    raise IOError('Failed because file {} does not exist.'.format(wavefile))
#                  
spfile = os.path.join(
    os.environ['PYSYN_CDBS'], 'calspec', infile)
if os.path.isfile(spfile):
    sp = S.FileSpectrum(spfile)
    hdr = fits.getheader(spfile, ext=0)
    target = hdr['TARGETID']
else:
    raise IOError('Failed because file {} does not exist.'.format(spfile))
#
S.setref(area=254009.0)
#
file090 = os.path.join(os.environ['NIS_SYNPHOT'], 'jwst_nis_f090w_clear.fits') 
file115 = os.path.join(os.environ['NIS_SYNPHOT'], 'jwst_nis_f115w_clear.fits') 
file140 = os.path.join(os.environ['NIS_SYNPHOT'], 'jwst_nis_f140m_clear.fits') 
file150 = os.path.join(os.environ['NIS_SYNPHOT'], 'jwst_nis_f150w_clear.fits') 
file158 = os.path.join(os.environ['NIS_SYNPHOT'], 'jwst_nis_f158m_clear.fits') 
file200 = os.path.join(os.environ['NIS_SYNPHOT'], 'jwst_nis_f200w_clear.fits') 
file277 = os.path.join(os.environ['NIS_SYNPHOT'], 'jwst_nis_clearp_f277w.fits')
file356 = os.path.join(os.environ['NIS_SYNPHOT'], 'jwst_nis_clearp_f356w.fits')
file444 = os.path.join(os.environ['NIS_SYNPHOT'], 'jwst_nis_clearp_f444w.fits')
file380 = os.path.join(os.environ['NIS_SYNPHOT'], 'jwst_nis_clearp_f380m.fits')
file430 = os.path.join(os.environ['NIS_SYNPHOT'], 'jwst_nis_clearp_f430m.fits')
file480 = os.path.join(os.environ['NIS_SYNPHOT'], 'jwst_nis_clearp_f480m.fits')
#
f090w = S.FileBandpass(file090)
f115w = S.FileBandpass(file115)
f140m = S.FileBandpass(file140)
f150w = S.FileBandpass(file150)
f158m = S.FileBandpass(file158)
f200w = S.FileBandpass(file200)
f277w = S.FileBandpass(file277)
f356w = S.FileBandpass(file356)
f444w = S.FileBandpass(file444)
f380m = S.FileBandpass(file380)
f430m = S.FileBandpass(file430)
f480m = S.FileBandpass(file480)
#
obs_f090w = S.Observation(sp, f090w, binset=wave)
obs_f115w = S.Observation(sp, f115w, binset=wave)
obs_f140m = S.Observation(sp, f140m, binset=wave)
obs_f150w = S.Observation(sp, f150w, binset=wave)
obs_f158m = S.Observation(sp, f158m, binset=wave)
obs_f200w = S.Observation(sp, f200w, binset=wave)
obs_f277w = S.Observation(sp, f277w, binset=wave)
obs_f356w = S.Observation(sp, f356w, binset=wave)
obs_f444w = S.Observation(sp, f444w, binset=wave)
obs_f380m = S.Observation(sp, f380m, binset=wave)
obs_f430m = S.Observation(sp, f430m, binset=wave)
obs_f480m = S.Observation(sp, f480m, binset=wave)
#
rate090 = obs_f090w.countrate()
rate115 = obs_f115w.countrate()
rate140 = obs_f140m.countrate()
rate150 = obs_f150w.countrate()
rate158 = obs_f158m.countrate()
rate200 = obs_f200w.countrate()
rate277 = obs_f277w.countrate()
rate356 = obs_f356w.countrate()
rate444 = obs_f444w.countrate()
rate380 = obs_f380m.countrate()
rate430 = obs_f430m.countrate()
rate480 = obs_f480m.countrate()
#
print('# Predicted Count Rates for {}'.format(target))
print('#   F090W     F115W     F140M     F150W     F158M     F200W     F277W'
      '     F356W     F380M     F430M     F444W     F480M')
print('%.3e %.3e %.3e %.3e %.3e'
      ' %.3e %.3e %.3e %.3e %.3e %.3e'
      ' %.3e' % (rate090, rate115, rate140, rate150, rate158, rate200, rate277,
                      rate356, rate380, rate430, rate444, rate480))
