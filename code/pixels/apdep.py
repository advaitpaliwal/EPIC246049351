"""Depth vs aperture from the TPF (SAP + lightkurve SFFCorrector, transits masked), box depth at fixed ephemeris.
usage: apdep.py TPF P T0 T14_H OUTJSON"""
import sys, json, numpy as np, warnings; warnings.filterwarnings('ignore')
import lightkurve as lk
sys.path.insert(0, 'k2-deep/scripts'); import lcs
tpp, P, T0, T14, OUT = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])/24, sys.argv[5]
tp = lk.read(tpp, quality_bitmask='default'); med = np.nanmedian(tp.flux.value, 0)
order = np.argsort(np.nan_to_num(med, nan=-1).ravel())[::-1]
def mk(n): m = np.zeros(med.shape, bool); m.flat[order[:n]] = True; return m
aps = {'brightest1': mk(1), 'brightest2': mk(2), 'core4': mk(4), 'core6': mk(6), 'thresh3': tp.create_threshold_mask(3, reference_pixel='center'), 'thresh1': tp.create_threshold_mask(1, reference_pixel='center')}
aps['halo_thresh1_minus_core4'] = aps['thresh1'] & ~aps['core4']
res = {}
for k, ap in aps.items():
    lc = tp.to_lightcurve(aperture_mask=ap).remove_nans(); lc = lc[np.isfinite(lc.centroid_col.value) & np.isfinite(lc.centroid_row.value)]
    t = lc.time.value; intr = np.abs(lcs.phase(t, P, T0)) < 0.75*T14
    try: c = lk.SFFCorrector(lc).correct(windows=20, bins=15, cadence_mask=~intr, restore_trend=True)
    except Exception as e: res[k] = str(e); continue
    f = c.flux.value; m = np.isfinite(f); tt, ff = lcs.flat(t[m], f[m]/np.nanmedian(f[m]), P, T0, T14)
    d, e, n = lcs.boxdepth(tt, ff, P, T0, T14)
    res[k] = dict(npix=int(ap.sum()), depth_ppm=float(d*1e6), err_ppm=float(e*1e6), rms_ppm=float(np.std(ff)*1e6), flux_frac=float(np.nansum(med[ap])/np.nansum(med[aps['thresh1']])))
    print(k, res[k], flush=True)
json.dump(res, open(OUT, 'w'), indent=1)
