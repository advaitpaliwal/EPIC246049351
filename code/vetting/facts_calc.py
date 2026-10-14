# Recompute the numbers quoted in FACTS.md / the note from the vetting files. usage: python -I facts_calc.py K2DEEP_DIR
import sys, json, glob, csv, numpy as np
D = sys.argv[1]
R = json.load(open(f'{D}/fit/result_prior_everest.json')); F = json.load(open(f'{D}/fit/result_free_everest.json'))
S = json.load(open(f'{D}/fit/result_prior_k2sff.json')); H = json.load(open(f'{D}/fit/result_prior_fresh.json'))
f3 = lambda x: '%.5g (+%.3g/-%.3g)' % (x[0], x[2], x[1])
for k in ['P_d', 'T0_BKJD', 'RpRs', 'b', 'a_Rs', 'inc_deg', 'rho_sun', 'depth_ppm', 'T14_h', 'Rp_earth', 'Teq_K', 'S_earth', 'jitter_ppm', 'q1', 'q2']: print('EVEREST prior', k, f3(R[k]))
print('T0 BJD_TDB %.5f' % (R['T0_BKJD'][0] + 2454833))
print('a AU %.4f' % (R['a_Rs'][0]*0.377432*6.957e8/1.495978707e11))
print('tau max %.0f, post-burn steps %d -> %.0f tau, acc %.3f, npts %d' % (max(R['tau']), R['nstep']*2//3, R['nstep']*2/3/max(R['tau']), R['acc'], R['npts']))
for nm, X in [('free', F), ('k2sff', S), ('fresh', H)]:
    print(nm, 'depth', f3(X['depth_ppm']), 'Rp', f3(X['Rp_earth']), 'b', f3(X['b']), 'rho', f3(X['rho_sun']), 'jitter %.0f' % X['jitter_ppm'][0])
    s = np.hypot(max(X['depth_ppm'][1:]), max(R['depth_ppm'][1:])); print('  depth - EVEREST = %.0f ppm, %.2f sigma' % (X['depth_ppm'][0] - R['depth_ppm'][0], (X['depth_ppm'][0] - R['depth_ppm'][0])/s))
print('free rho vs TIC 6.803: diff %.2f sigma' % ((F['rho_sun'][0] - 6.80341)/F['rho_sun'][1]))
for v in ['gr', 'grs', 'none', 'cc']:
    fs = sorted(glob.glob(f'{D}/fpp_joint/out/res_{v}_*.json')); x = np.array([json.load(open(f))['FPP'] for f in fs]); n = np.array([json.load(open(f))['NFPP'] for f in fs])
    r0 = json.load(open(fs[0]))
    print('FPP %s n=%d mean %.4f sd %.4f min %.4f max %.4f NFPP max %.2g n_stars %d depth %.1f err %.1f nbins %d filt %s' % (v, len(x), x.mean(), x.std(ddof=1) if len(x) > 1 else 0, x.min(), x.max(), n.max(), r0['n_stars'], r0['depth_ppm'], r0['err_ppm'], r0['nbins'], r0['filt']))
g = [r for r in csv.DictReader(open(f'{D}/catalogue/gaia_dr3_60as.csv'))]; t = min(g, key=lambda r: float(r['sep']))
print('Gaia plx %s -> %.1f pc; G %s; BP-RP %s; RUWE %s' % (t['parallax'], 1000/float(t['parallax']), t['phot_g_mean_mag'], t['bp_rp'], t['ruwe']))
for r in sorted(g, key=lambda r: float(r['sep']))[1:]: print('  Gaia neighbour sep %.1f G %.2f plx %.2f' % (float(r['sep']), float(r['phot_g_mean_mag']), float(r['parallax'])))
cc = np.loadtxt(f'{D}/fpp_joint/cc832_arcsec_dmag.csv', delimiter=',')
for s in [0.1, 0.2, 0.5, 1.0, 1.2]: print('832 nm contrast at %.1f arcsec: %.2f mag' % (s, np.interp(s, cc[:, 0], cc[:, 1])))
print('832 curve range %.2f-%.2f arcsec' % (cc[0, 0], cc[-1, 0]))
z = np.load(f'{D}/lc/everest.npz'); print('EVEREST span %.3f-%.3f (%.1f d), raw n %d, flat n %d' % (z['t_raw'].min(), z['t_raw'].max(), np.ptp(z['t_raw']), len(z['t_raw']), len(z['t'])))
ep = np.round((z['t'] - R['T0_BKJD'][0])/R['P_d'][0]).astype(int); ph = z['t'] - R['T0_BKJD'][0] - ep*R['P_d'][0]
print('epochs with >=1 cadence within T14/2:', sorted(set(ep[np.abs(ph) < R['T14_h'][0]/48]))); print('predicted epochs in span:', list(range(ep.min(), ep.max() + 1)))
