"""Explain per-epoch K2SFF depths. For each C12 epoch: single-epoch box depth (vet.py estimator: in = |phase|<0.35*T14, out = that epoch's T14<|phase|<0.5 d)
in K2SFF, EVEREST and fresh; K2SFF MOVING flag (thruster), ARCLENGTH, K2 llc SAP_QUALITY bits in-transit; K2SFF segment boundaries (jumps in the fitted
FRAW/FCOR ratio); cadences missing from K2SFF vs EVEREST in-transit. usage: epoch_k2sff.py"""
import sys, glob, json, numpy as np
from astropy.io import fits
sys.path.insert(0, 'k2-deep/scripts'); import lcs
D = 'k2-deep/246049351'; R = json.load(open(f'{D}/fit/result_prior_everest.json'))
P, T0, T14 = R['P_d'][0], R['T0_BKJD'][0], R['T14_h'][0]/24
L = {k: np.load(f'{D}/lc/{k}.npz') for k in ['everest', 'k2sff', 'fresh']}
s = fits.getdata(glob.glob(f'{D}/lc/mastDownload/HLSP/*k2sff*c12*/*.fits')[0], 1)
q = fits.getdata(glob.glob(f'{D}/lc/mastDownload/K2/ktwo246049351-c12_lc/*llc.fits')[0], 1)
corr = s['FCOR']/s['FRAW']
def single(t, f, e):
    ph = lcs.phase(t, P, T0); ep = np.round((t - T0)/P).astype(int); m = (np.abs(ph) < 0.35*T14) & (ep == e); o = (np.abs(ph) > T14) & (np.abs(ph) < 0.5) & (ep == e)
    if m.sum() == 0 or o.sum() < 6: return None
    return [round(float((np.median(f[o]) - np.mean(f[m]))*1e6)), round(float(np.std(f[o])/np.sqrt(m.sum())*1e6)), int(m.sum()), int(o.sum())]
eps = np.unique(np.round((L['everest']['t'] - T0)/P).astype(int)); rows = []
for e in eps:
    tm = T0 + e*P; r = dict(epoch=int(e), tmid=round(float(tm), 4))
    for k, z in L.items(): r[k] = single(z['t'], z['f'], e)
    w = np.abs(s['T'] - tm) < 0.35*T14; r['k2sff_raw_incad'] = int(w.sum()); r['k2sff_moving_in'] = int(s['MOVING'][w].sum())
    r['k2sff_flat_incad'] = int(np.sum(np.abs(L['k2sff']['t'] - tm) < 0.35*T14)); r['everest_flat_incad'] = int(np.sum(np.abs(L['everest']['t'] - tm) < 0.35*T14))
    wq = np.abs(q['TIME'] - tm) < 0.35*T14; r['llc_quality_in'] = [int(x) for x in q['SAP_QUALITY'][wq]]
    w2 = np.abs(s['T'] - tm) < 0.5; r['k2sff_corr_range_ppm_pm0.5d'] = round(float((corr[w2].max() - corr[w2].min())*1e6)) if w2.sum() else None
    r['k2sff_corr_in_minus_out_ppm'] = round(float((np.mean(corr[w]) - np.median(corr[w2 & ~w]))*1e6)) if w.sum() and (w2 & ~w).sum() else None
    r['k2sff_fraw_in_minus_out_ppm'] = round(float((np.median(s['FRAW'][w2 & ~w]) - np.mean(s['FRAW'][w]))/np.median(s['FRAW'])*1e6)) if w.sum() else None
    rows.append(r); print(r)
json.dump(dict(P=P, T0=T0, T14_h=T14*24, columns='single-epoch [depth_ppm, err_ppm, n_in, n_out]', rows=rows), open(f'{D}/fit/epoch_k2sff.json', 'w'), indent=1)
