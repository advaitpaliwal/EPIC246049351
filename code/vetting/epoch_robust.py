# Robustness of the EVEREST box depth (vet.py estimator) to the shallow epoch -3 and to epochs with a thruster firing within T14 of mid-transit.
# usage: python -I epoch_robust.py K2DEEP_DIR
import sys, json, numpy as np
D = sys.argv[1]; R = json.load(open(f'{D}/fit/result_prior_everest.json'))
P, T0, T14 = R['P_d'][0], R['T0_BKJD'][0], R['T14_h'][0]/24
V = json.load(open(f'{D}/fit/vet_everest.json'))
for src in ['everest', 'k2sff', 'fresh']:
    z = np.load(f'{D}/lc/{src}.npz'); t, f = z['t'], z['f']
    ph = ((t - T0 + 0.5*P) % P) - 0.5*P; ep = np.round((t - T0)/P).astype(int)
    oot = (np.abs(ph) > T14) & (np.abs(ph) < 0.5); inn = np.abs(ph) < 0.35*T14
    def dep(drop):
        k = ~np.isin(ep, drop); d = np.median(f[oot & k]) - np.mean(f[inn & k]); e = np.std(f[oot & k])/np.sqrt((inn & k).sum())
        return round(d*1e6), round(e*1e6), int((inn & k).sum())
    print(src, 'all', dep([]), 'drop -3', dep([-3]), 'drop firing-in-T14 epochs (-8,-3,-2,1,4,6,7,9)', dep([-8, -3, -2, 1, 4, 6, 7, 9]))
pt = V['per_transit']; d = np.array([p['depth_ppm'] for p in pt]); e = np.array([p['err_ppm'] for p in pt])
mu = np.average(d, weights=1/e**2); r = (d - mu)/e
for p, x in zip(pt, r): print('epoch %d depth %.0f err %.0f resid %.2f sigma' % (p['epoch'], p['depth_ppm'], p['err_ppm'], x))
k = np.array([p['epoch'] for p in pt]) != -3; mu2 = np.average(d[k], weights=1/e[k]**2)
print('chi2 all %.2f/%d; without -3 %.2f/%d; weighted mean %.0f -> %.0f' % (np.sum(r**2), len(d) - 1, np.sum(((d[k] - mu2)/e[k])**2), k.sum() - 1, mu, mu2))
