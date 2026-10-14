"""Adapted from ../../248780202/scripts/pertr_null.py (single planet, no second-planet mask). Empirical per-transit error: the vet2.py single-epoch estimator (in = |phase|<0.35*T14 vs that epoch's T14<|phase|<0.5 d) applied at 300 fake
ephemerides (same P, T0 shifted by U(0.15,0.85)*P, avoiding both planets' transits by 1.5*T14). chi2_emp = sum((d_i-d_w)^2/s_emp^2) with s_emp = std of fake
single-epoch depths (pooled), and the same chi2 computed for each fake ephemeris -> fraction of fakes with chi2 >= real. usage: pertr_null.py SRC PLANET"""
import sys, json, numpy as np
sys.path.insert(0, 'k2-deep/scripts'); import lcs
src, pl = sys.argv[1], 'b'; D = 'k2-deep/246049351'
R = {'b': json.load(open(f'{D}/fit/result_prior_everest.json'))}  # single planet
P, T0, T14 = R[pl]['P_d'][0], R[pl]['T0_BKJD'][0], R[pl]['T14_h'][0]/24
Po, T0o, T14o = 1e9, 0.0, -1.0  # no second planet: mask keeps everything
z = np.load(f'{D}/lc/{src}.npz'); t, f = z['t'], z['f']
keep = np.abs(lcs.phase(t, Po, T0o)) > T14o; t, f = t[keep], f[keep]
def per_epoch(Tx, excl_real):
    ph = lcs.phase(t, P, Tx); ep = np.round((t - Tx)/P).astype(int); inn = np.abs(ph) < 0.35*T14; oot = (np.abs(ph) > T14) & (np.abs(ph) < 0.5)
    if excl_real: oot &= np.abs(lcs.phase(t, P, T0)) > T14
    d = []
    for e in np.unique(ep[inn]):
        m = inn & (ep == e); mo = oot & (ep == e)
        if m.sum() < 2 or mo.sum() < 6: continue
        d.append(np.median(f[mo]) - np.mean(f[m]))
    return np.array(d)*1e6
real = per_epoch(T0, False); rng = np.random.default_rng(3); fakes = []
for _ in range(300):
    Tx = T0 + rng.uniform(0.15, 0.85)*P
    fakes.append(per_epoch(Tx, True))
allf = np.concatenate(fakes); s = 1.4826*np.median(np.abs(allf - np.median(allf))); s_std = np.std(allf)
def chi(d): w = np.mean(d); return float(np.sum((d - w)**2)/s_std**2)
cr = chi(real); cf = np.array([chi(d) for d in fakes if len(d) > 5])
res = dict(source=src, planet=pl, n_epochs=int(len(real)), per_epoch_ppm=real.round(0).tolist(), mean_ppm=float(real.mean()), sem_ppm=float(s_std/np.sqrt(len(real))),
           empirical_single_epoch_sd_ppm=float(s_std), empirical_sd_mad_ppm=float(s), real_sd_ppm=float(np.std(real)), chi2_emp=cr, dof=int(len(real)-1),
           frac_fake_chi2_ge_real=float(np.mean(cf >= cr)), fake_chi2_p50_p95=np.percentile(cf, [50, 95]).tolist(),
           first_half_mean=float(real[:len(real)//2].mean()), second_half_mean=float(real[len(real)//2:].mean()))
json.dump(res, open(f'{D}/fit/pertr_null_{src}_{pl}.json', 'w'), indent=1); print(json.dumps(res))
