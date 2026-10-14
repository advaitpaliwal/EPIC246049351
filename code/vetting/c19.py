"""C19 check (K2SFF v1 C19 FCOR, flattened as in prep3.py). Propagates the C12 prior-fit ephemeris; for shifts of the predicted mid-times within
+-3 sigma_T (step T14/8): box depth (in |phase|<0.5*T14 vs |phase|>T14 within the same data) and SNR; lists predicted transits inside the data span,
and the number of in-transit cadences per transit at the nominal ephemeris."""
import sys, json, numpy as np
sys.path.insert(0, 'k2-deep/scripts'); import lcs
D = 'k2-deep/246049351'; R = json.load(open(f'{D}/fit/result_prior_everest.json'))
P, T0, T14, dep = R['P_d'][0], R['T0_BKJD'][0], R['T14_h'][0]/24, R['depth_ppm'][0]
eP, eT0 = max(R['P_d'][1:]), max(R['T0_BKJD'][1:])
z = np.load(f'{D}/lc/k2sff_c19.npz'); t, f = z['t'], z['f']
n = np.round((np.median(t) - T0)/P); sT = float(np.hypot(eT0, n*eP))
eps = np.arange(np.ceil((t.min() - T0)/P), np.floor((t.max() - T0)/P) + 1).astype(int)
tr = [dict(epoch=int(e), tmid=float(T0 + e*P), n_in=int(np.sum(np.abs(t - (T0 + e*P)) < 0.5*T14)), n_within_3sig=int(np.sum(np.abs(t - (T0 + e*P)) < 0.5*T14 + 3*sT))) for e in eps]
rows = []
for sh in np.arange(-3*sT, 3*sT + 1e-9, T14/8):
    ph = lcs.phase(t, P, T0 + sh); i = np.abs(ph) < 0.5*T14; o = np.abs(ph) > T14
    if i.sum() < 2: rows.append([float(sh*24), None, None, int(i.sum())]); continue
    d = np.median(f[o]) - np.mean(f[i]); e = np.std(f[o])/np.sqrt(i.sum()); rows.append([float(sh*24), float(d*1e6), float(e*1e6), int(i.sum())])
ok = [r for r in rows if r[1] is not None]
res = dict(span=[float(t.min()), float(t.max())], npts=int(len(t)), rms_ppm=float(np.std(f)*1e6), n_epochs_from_ref=int(n), sigma_T_h=sT*24, transits=tr,
           expected_depth_ppm=dep, scan_hshift_depth_err_nin=rows, scan_max_snr=max(r[1]/r[2] for r in ok), scan_min_snr=min(r[1]/r[2] for r in ok),
           frac_shifts_consistent_with_K2_2sig=float(np.mean([abs(r[1]-dep) < 2*r[2] for r in ok])))
json.dump(res, open(f'{D}/fit/c19_check.json', 'w'), indent=1); print(json.dumps({k: v for k, v in res.items() if k != 'scan_hshift_depth_err_nin'})); [print(r) for r in rows]
