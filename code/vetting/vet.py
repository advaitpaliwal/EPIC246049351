"""Odd/even, secondary eclipse, per-transit depths, negative-dip search. usage: vet.py EPIC SRC P T0 T14_H"""
import sys, json, numpy as np
sys.path.insert(0, 'k2-deep/scripts'); import lcs
epic, src, P, T0, T14 = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])/24
D = f'k2-deep/{epic}'
z = np.load(f'{D}/lc/{src}.npz'); t, f = z['t'], z['f']
ph = lcs.phase(t, P, T0); ep = np.round((t - T0)/P).astype(int)
oot = (np.abs(ph) > T14) & (np.abs(ph) < 0.5)
sd = 1.4826*np.median(np.abs(f[oot] - np.median(f[oot])))
def dep(mask_in, mask_out):
    d = np.median(f[mask_out]) - np.mean(f[mask_in]); e = np.std(f[mask_out])/np.sqrt(mask_in.sum()); return float(d*1e6), float(e*1e6), int(mask_in.sum())
inn = np.abs(ph) < 0.35*T14
res = {}
res['all'] = dep(inn, oot)
res['odd'] = dep(inn & (ep % 2 == 1), oot); res['even'] = dep(inn & (ep % 2 == 0), oot)
res['odd_even_sigma'] = abs(res['odd'][0]-res['even'][0])/np.hypot(res['odd'][1], res['even'][1])
per = []
for e in np.unique(ep[inn]):
    m = inn & (ep == e); mo = oot & (ep == e)
    per.append(dict(epoch=int(e), tmid=float(T0+e*P), n_in=int(m.sum())) | dict(zip(['depth_ppm', 'err_ppm', 'n'], dep(m, mo if mo.sum() > 5 else oot))))
res['per_transit'] = per
ch = [p['depth_ppm'] for p in per]; ce = [p['err_ppm'] for p in per]
mu = np.average(ch, weights=1/np.square(ce)); res['per_transit_chi2'] = float(np.sum(((np.array(ch)-mu)/np.array(ce))**2)); res['per_transit_dof'] = len(ch)-1
# secondary eclipse at phase 0.5 (circular) with same duration
ph2 = lcs.phase(t, P, T0 + P/2); sec_in = np.abs(ph2) < 0.5*T14; sec_out = (np.abs(ph2) > T14) & (np.abs(ph2) < 0.5) & (np.abs(ph) > T14)
res['secondary_phase0.5'] = dep(sec_in, sec_out)
# scan all phases for the deepest box of width T14 (excluding primary) -> secondary at any phase (eccentric)
grid = np.arange(0.05, 0.95, 0.005); sc = []
for g in grid:
    p3 = lcs.phase(t, P, T0 + g*P); m = (np.abs(p3) < 0.5*T14) & (np.abs(ph) > T14)
    if m.sum() < 4: continue
    sc.append((g, *dep(m, oot)))
sc = np.array(sc); snr = sc[:, 1]/sc[:, 2]
res['phase_scan_min_snr'] = float(snr.min()); res['phase_scan_max_snr'] = float(snr.max()); res['phase_scan_argmax'] = float(sc[np.argmax(snr), 0])
res['phase_scan_max_depth_ppm'] = float(sc[np.argmax(snr), 1])
res['phase_scan_std_snr'] = float(np.std(snr))
res['primary_snr'] = res['all'][0]/res['all'][1]
# geometric albedo bound from 3-sigma secondary: depth_sec = Ag*(Rp/a)^2
json.dump(res, open(f'{D}/fit/vet_{src}.json', 'w'), indent=1)
for k, v in res.items():
    if k != 'per_transit': print(k, v)
for p in per: print(p)
