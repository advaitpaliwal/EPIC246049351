"""Light-curve prep for EPIC 246049351 (C12; C19 K2SFF + TPF also exist). Adapted from ../../scripts/prep.py + lcs.py + ../../248780202/scripts/prep2.py (not modified).
Per reduction: EVEREST v2 C12 (QUALITY==0 except bits 23,24), K2SFF v1 C12 FCOR (BESTAPER ext 1), fresh C12 (threshold-3 SAP from the TPF, quality_bitmask='default',
lightkurve SFFCorrector windows=20 bins=15 with transits masked at |phase|<0.75*T14, restore_trend=True), K2SFF C19 and fresh C19 (same recipe on the C19 TPF).
Then wotan biweight, window 0.75 d, in-transit points (|phase|<0.75*T14) masked from the trend and the trend interpolated across them; out-of-transit clip
+4/-8 MAD-sigma; everything below -20 sigma clipped. Box depth: mean of |phase|<0.3*T14 vs median of T14<|phase|<4*T14, error = std(out)/sqrt(n_in)."""
import sys, json, glob, numpy as np, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0, 'k2-deep/scripts'); import lcs, lightkurve as lk
from astropy.io import fits
epic = 246049351; D = f'k2-deep/{epic}'
P, T0, dur = 4.33316777, 2906.7617, 1.69/24
def k2sff_c(c):
    p = [x for x in lcs._f(epic, 'hlsp_k2sff*.fits') if f'-c{c}_' in x][0]; d = fits.getdata(p, 1)
    t, f = d['T'], d['FCOR']; m = np.isfinite(t) & np.isfinite(f); return t[m], f[m]/np.nanmedian(f[m])
def fresh_c(c):
    tp = lcs.tpf(epic, f'c{c}'); ap = tp.create_threshold_mask(threshold=3, reference_pixel='center')
    lc = tp.to_lightcurve(aperture_mask=ap).remove_nans(); lc = lc[np.isfinite(lc.centroid_col.value) & np.isfinite(lc.centroid_row.value)]
    t = lc.time.value; intr = np.abs(lcs.phase(t, P, T0)) < 0.75*dur
    c_ = lk.SFFCorrector(lc).correct(windows=20, bins=15, cadence_mask=~intr, restore_trend=True)
    f = c_.flux.value; m = np.isfinite(f); print('fresh', c, 'aperture npix', int(ap.sum()), 'n', int(m.sum())); return t[m], f[m]/np.nanmedian(f[m])
out = {}
srcs = {'everest': lambda: lcs.everest(epic), 'k2sff': lambda: k2sff_c(12), 'k2sc': lambda: lcs.k2sc(epic), 'fresh': lambda: fresh_c(12),
        'k2sff_c19': lambda: k2sff_c(19), 'fresh_c19': lambda: fresh_c(19)}
for name, fn in srcs.items():
    try: r = fn()
    except Exception as e: out[name] = f'failed: {e}'; print(name, e); continue
    if r is None: out[name] = 'not available at MAST'; print(name, 'none'); continue
    t, f = r; tf, ff = lcs.flat(t, f, P, T0, dur); np.savez(f'{D}/lc/{name}.npz', t_raw=t, f_raw=f, t=tf, f=ff)
    d, e, n = lcs.boxdepth(tf, ff, P, T0, dur)
    oot = np.abs(lcs.phase(tf, P, T0)) > dur
    out[name] = dict(n=int(len(tf)), t_span=[float(tf.min()), float(tf.max())], rms_ppm=float(np.std(ff[oot])*1e6), box_depth_ppm=float(d*1e6), err_ppm=float(e*1e6), n_in=n, snr=float(d/e))
    print(name, out[name], flush=True)
json.dump(out, open(f'{D}/lc/depths_box.json', 'w'), indent=1)
