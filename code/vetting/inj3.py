"""Copied from ../../scripts/inj.py; only the result path changed (fit/result_prior_everest.json). Injection-recovery using the hunt's own detrend+TLS (k2tls_hunt.py clean/power settings).
usage: inj.py EPIC SRC RPRS ARS_AT_P0 P0 N SEED OUT [mode]
mode=inj (default): inject batman transits at random P in [P0*0.7, P0*1.3] and random T0 into the raw signal-masked light curve.
mode=null: no injection, real signal masked -> TLS false-alarm level.
mode=real: real data unmasked -> SDE of the candidate itself."""
import sys, json, numpy as np, batman, warnings; warnings.filterwarnings('ignore')
from transitleastsquares import transitleastsquares
sys.path.insert(0, 'k2-hunt/scripts'); from k2tls_hunt import clean
epic, src, k, aR0, P0, N, seed, OUT = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]), sys.argv[8]
mode = sys.argv[9] if len(sys.argv) > 9 else 'inj'
T0r, Pr, T14r = float(sys.argv[10]) if len(sys.argv) > 10 else None, P0, None
D = f'k2-deep/{epic}'
z = np.load(f'{D}/lc/{src}.npz'); t0, f0 = z['t_raw'], z['f_raw']
R = json.load(open(f'{D}/fit/result_prior_everest.json')); Tc, Pc, T14 = R['T0_BKJD'][0], R['P_d'][0], R['T14_h'][0]/24
intr = np.abs(((t0 - Tc + 0.5*Pc) % Pc) - 0.5*Pc) < 0.75*T14
rng = np.random.default_rng(seed); out = []
targ = [l for l in open('k2-hunt/targets_all.csv') if l.startswith(epic)][0].split(',')
Rs, Ms = float(targ[7]), float(targ[8])
for i in range(N):
    if mode == 'real': t, f = t0, f0.copy(); Pi = Pc; Ti = Tc
    else:
        t, f = t0[~intr], f0[~intr].copy()
        if mode == 'inj':
            Pi = rng.uniform(0.7*P0, 1.3*P0); Ti = t.min() + rng.uniform(0, Pi)
            pm = batman.TransitParams(); pm.t0=Ti; pm.per=Pi; pm.rp=k; pm.a=aR0*(Pi/P0)**(2/3); pm.inc=np.degrees(np.arccos(0.4/pm.a)); pm.ecc=0; pm.w=90; pm.limb_dark='quadratic'; pm.u=[0.3,0.3]
            f = f*batman.TransitModel(pm, t, supersample_factor=15, exp_time=1765.5/86400).light_curve(pm)
        else: Pi, Ti = np.nan, np.nan
    tc, fc = clean(t, f)
    r = transitleastsquares(tc, fc).power(period_min=0.5, period_max=min(30.0, (tc.max()-tc.min())/2), R_star=Rs, M_star=Ms, R_star_min=0.08, R_star_max=1.0,
        M_star_min=0.08, M_star_max=1.0, oversampling_factor=2, duration_grid_step=1.15, use_threads=1, show_progress_bar=False)
    rec = bool(np.isfinite(Pi) and (abs(r.period - Pi)/Pi < 0.01 or abs(r.period - 2*Pi)/Pi < 0.02 or abs(2*r.period - Pi)/Pi < 0.02))
    out.append(dict(i=i, mode=mode, P_inj=float(Pi), T0_inj=float(Ti), P_rec=float(r.period), SDE=float(r.SDE), snr=float(r.snr), depth=float(1-r.depth), recovered_period=rec))
    print(out[-1], flush=True)
    if mode == 'real': break
json.dump(out, open(OUT, 'w'), indent=1)
