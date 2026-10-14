"""Thruster-firing and alias checks, single planet. Adapted from ../../248780202/scripts/thruster2.py (not modified).
Firings = K2 llc SAP_QUALITY bits 20/21 (1-indexed; masks 1<<19, 1<<20). P_thr = median spacing of consecutive firings in 0.1-0.4 d.
Rayleigh coherence |mean exp(2 pi i phase)| of firing times folded on P vs 1/sqrt(N). Box depth (lcs.boxdepth) of the flattened EVEREST LC at P, P/2 (odd phases),
2P (both parities), phase 1/3 and phase 0.5."""
import numpy as np, glob, json, sys
from astropy.io import fits
sys.path.insert(0, 'k2-deep/scripts'); import lcs
D = 'k2-deep/246049351'; R = json.load(open(f'{D}/fit/result_prior_everest.json'))
P, T0, T14 = R['P_d'][0], R['T0_BKJD'][0], R['T14_h'][0]/24
d = fits.getdata(glob.glob(f'{D}/lc/mastDownload/K2/ktwo246049351-c12_lc/*llc.fits')[0], 1)
t, q = d['TIME'], d['SAP_QUALITY']; ok = np.isfinite(t); t, q = t[ok], q[ok]
fire = ((q & (1 << 19)) | (q & (1 << 20))) > 0; tf = t[fire]
dt = np.diff(tf); dt = dt[(dt > 0.1) & (dt < 0.4)]; Pthr = float(np.median(dt))
ph = lcs.phase(t, P, T0); intr = np.abs(ph) < 0.5*T14; phf = lcs.phase(tf, P, T0)
out = dict(n_cad=int(len(t)), n_fire=int(fire.sum()), P_thr_d=Pthr, P_thr_h=Pthr*24, P_over_Pthr=P/Pthr, frac_from_int=float(P/Pthr - round(P/Pthr)),
           frac_fire_in_transit=float(np.mean(np.abs(phf) < 0.5*T14)), expected_frac=float(T14/P), flagged_frac_in=float(np.mean(fire[intr])), flagged_frac_out=float(np.mean(fire[~intr])),
           rayleigh=float(np.abs(np.mean(np.exp(2j*np.pi*phf/P)))), rayleigh_random=float(1/np.sqrt(len(tf))))
for src in ['everest', 'k2sff']:
    z = np.load(f'{D}/lc/{src}.npz'); tt, ff = z['t'], z['f']
    for lab, PP, TT in [('P', P, T0), ('P/2_odd_phase', P/2, T0 + P/2), ('2P_epoch0', 2*P, T0), ('2P_epoch1', 2*P, T0 + P), ('P_phase1/3', P, T0 + P/3), ('P_phase0.5', P, T0 + P/2)]:
        dd, e, n = lcs.boxdepth(tt, ff, PP, TT, T14); out[f'{src}_box_{lab}'] = [float(dd*1e6), float(e*1e6), n]
json.dump(out, open(f'{D}/fit/thruster_alias.json', 'w'), indent=1); print(json.dumps(out, indent=1))
