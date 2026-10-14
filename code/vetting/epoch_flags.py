# Per-epoch K2 quality flags, thruster firings and K2SFF arclength near each C12 transit of EPIC 246049351.
# usage: python -I epoch_flags.py K2DEEP_DIR
import sys, glob, json, numpy as np
from astropy.io import fits
from lightkurve.utils import KeplerQualityFlags as Q
D = sys.argv[1]; R = json.load(open(f'{D}/fit/result_prior_everest.json'))
P, T0, T14 = R['P_d'][0], R['T0_BKJD'][0], R['T14_h'][0]/24
q = fits.getdata(glob.glob(f'{D}/lc/mastDownload/K2/ktwo246049351-c12_lc/*llc.fits')[0], 1)
s = fits.getdata(glob.glob(f'{D}/lc/mastDownload/HLSP/*k2sff*c12*/*.fits')[0], 1)
ev = np.load(f'{D}/lc/everest.npz'); print('everest keys', ev.files)
for b in [131072, 262144, 524288, 1048576, 32, 32768]: print(b, Q.decode(b))
t = q['TIME']; ok = np.isfinite(t); t, qq = t[ok], q['SAP_QUALITY'][ok]
fire = (qq & 1048576) > 0
for e in range(-8, 11):
    tm = T0 + e*P; w = np.abs(t - tm) < 0.5*T14; wn = np.abs(t - tm) < T14
    ws = np.abs(s['T'] - tm) < 0.5*T14
    nf = int(fire[wn].sum()); tf = (t[wn][fire[wn]] - tm)*24
    print(e, round(tm, 4), 'cad', int(w.sum()), 'flags', sorted({Q.decode(int(x)).__str__() for x in qq[w]}), 'fire_within_T14', nf, 'at_h', np.round(tf, 2).tolist(),
          'k2sff_arclen', np.round(s['ARCLENGTH'][ws], 2).tolist(), 'moving', s['MOVING'][ws].tolist())
