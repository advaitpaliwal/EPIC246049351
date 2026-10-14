"""Same-campaign known-planet reference for the pixel null test: K2-325 b (EPIC 246074965, C12, Kp 16.28, NEA: P 6.93 d, T0 2457743.818 BJD, depth 0.672%,
Castro Gonzalez et al. 2020). Downloads its C12 TPF + EVEREST LC, then refines P, T0, T14 by a box grid on the flattened EVEREST LC
(P 6.90-6.96 step 0.0005, T0 +-0.1 d step 0.002, T14 1-3.5 h) because the NEA P has only 3 decimals."""
import sys, glob, json, numpy as np, lightkurve as lk, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0, 'k2-deep/scripts'); import lcs
from astropy.io import fits
D = 'k2-deep/246049351/pix/ref'
st = lk.search_targetpixelfile('EPIC 246074965', mission='K2', campaign=12); print(st); st[0].download(download_dir=D)
sr = lk.search_lightcurve('EPIC 246074965', mission='K2', campaign=12); print(sr)
for r in sr:
    if str(r.author[0]).lower() == 'everest': r.download(download_dir=D)
open(D + '/query.url', 'w').write('lightkurve.search_targetpixelfile/search_lightcurve EPIC 246074965 K2 campaign 12 -> MAST https://archive.stsci.edu/missions-and-data/k2\n')
p = glob.glob(D + '/mastDownload/**/hlsp_everest*.fits', recursive=True)[0]; d = fits.getdata(p, 1)
m = np.isfinite(d['TIME']) & np.isfinite(d['FLUX']) & ((d['QUALITY'] & ~((1 << 23) | (1 << 24))) == 0); t, f = d['TIME'][m], d['FLUX'][m]/np.nanmedian(d['FLUX'][m])
from wotan import flatten
fl = flatten(t, f, method='biweight', window_length=1.0); ok = np.isfinite(fl); t, fl = t[ok], fl[ok]
best = (-1, 0, 0, 0)
for P in np.arange(6.90, 6.96, 0.0005):
    for T0 in np.arange(2910.818 - 0.1, 2910.818 + 0.1, 0.002):
        for dur in [1.0, 1.5, 2.0, 2.5, 3.0, 3.5]:
            dd, e, n = lcs.boxdepth(t, fl, P, T0, dur/24)
            if n > 5 and dd/e > best[0]: best = (dd/e, P, T0, dur, dd)
print('best', best); json.dump(dict(snr=best[0], P=best[1], T0_BKJD=best[2], T14_h=best[3], depth_ppm=best[4]*1e6), open(D + '/ref_ephem.json', 'w'), indent=1)
