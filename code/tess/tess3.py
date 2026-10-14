"""TESS check at the propagated K2 ephemeris, with dilution, a period-grid search, a look-elsewhere null and injection-recovery.
Adapted from ../../248780202/scripts/tess2.py (not modified). Sectors from tess-point (tess/tesspoint.txt). TESScut 15x15 FFI cutouts (lightkurve).
Photometry: apertures 'pix1' (pixel containing the target, WCS) and 'pix2x2' (that pixel + the 3 others of the 2x2 block nearest the target's sub-pixel position);
background = per-cadence median of the faintest 30% of pixels; quality==0; wotan biweight window 0.75 d; upward 4-sigma clip.
Dilution: Gaussian PRF (sigma = 0.8 px, assumed) for every TIC star within 2.5 arcmin, Tmag from TIC v8.2 -> target fraction of aperture flux.
Expected TESS depth = K2 fitted depth (result_prior_everest.json) x target fraction.
Search: P on a grid P0 +- 4 sigma_P (step chosen so the shift at the last sector is <= T14/4), T0 fixed at the K2 value (sigma_T0 negligible beside n*sigma_P).
For each trial: in-transit = |phase|<0.5*T14; depth = 1 - mean(in), err = rms_sector/sqrt(n_in); sectors combined by inverse-variance weighting; SNR = depth/err.
Null: the same grid at 200 phase-shifted ephemerides T0 + U(0.15,0.85)*P (real-transit windows masked) -> distribution of the max SNR over the grid.
Injection: 100 batman transits at the expected depth (b=0.27, aRs, P0, phase-shifted T0, real windows masked), multiplied into the raw aperture flux before detrending,
then the same grid search; recovered if max SNR over the grid exceeds the 99th percentile of the null max-SNR. usage: tess3.py APERTURE"""
import sys, json, glob, os, numpy as np, warnings; warnings.filterwarnings('ignore')
import lightkurve as lk, batman
from wotan import flatten
D = 'k2-deep/246049351'; RA, DEC = 351.49778, -8.27045; APN = sys.argv[1]
R = json.load(open(f'{D}/fit/result_prior_everest.json'))
T0, eT0 = R['T0_BKJD'][0], max(R['T0_BKJD'][1:]); P0, eP = R['P_d'][0], max(R['P_d'][1:]); T14 = R['T14_h'][0]/24; dep = R['depth_ppm'][0]*1e-6
aRs, b, k = R['a_Rs'][0], R['b'][0], R['RpRs'][0]
secs = [int(l.split('|')[5]) for l in open(f'{D}/tess/tesspoint.txt') if l.startswith('0000')]
tic = json.load(open(f'{D}/catalogue/mast_tic.json'))['data']
out = dict(T0_BKJD=T0, eT0=eT0, P=P0, eP=eP, T14_h=T14*24, K2_depth_ppm=dep*1e6, aperture=APN, sectors={})
LC = {}
for s in secs:
    ex = glob.glob(f'{D}/tess/tesscut/*s{s:04d}*.fits') + glob.glob(f'{D}/tess/tesscut/*-s{s:04d}-*')
    if not ex:
        sr = lk.search_tesscut(f'{RA} {DEC}', sector=s)
        if len(sr) == 0: out['sectors'][s] = 'no TESScut product at MAST'; print(s, 'none', flush=True); continue
        tpf = sr.download(cutout_size=15, download_dir=f'{D}/tess/tesscut')
        open(f'{D}/tess/tesscut_s{s}.url', 'w').write(f'https://mast.stsci.edu/tesscut/api/v0.1/astrocut?ra={RA}&dec={DEC}&y=15&x=15&sector={s}\n')
    else: tpf = lk.read(ex[0])
    cx, cy = tpf.wcs.all_world2pix([[RA, DEC]], 0)[0]; ix, iy = int(np.round(cx)), int(np.round(cy))
    ap = np.zeros(tpf.shape[1:], bool)
    if APN == 'pix1': ap[iy, ix] = True
    else:
        x0 = ix if cx >= ix else ix - 1; y0 = iy if cy >= iy else iy - 1; ap[y0:y0+2, x0:x0+2] = True
    yy, xx = np.mgrid[:ap.shape[0], :ap.shape[1]]; ft = 0.0; fa = 0.0
    for r in tic:
        if r.get('Tmag') is None: continue
        px, py = tpf.wcs.all_world2pix([[r['ra'], r['dec']]], 0)[0]; fl = 10**(-0.4*(r['Tmag'] - 10))
        g = np.exp(-0.5*((xx-px)**2 + (yy-py)**2)/0.8**2); g /= 2*np.pi*0.8**2; fin = fl*g[ap].sum(); fa += fin
        if r['ID'] == 49678675: ft = fin
    frac = ft/fa if fa else np.nan
    med = np.nanmedian(tpf.flux.value, 0); bkm = med < np.nanpercentile(med, 30); bkg = np.nanmedian(tpf.flux.value[:, bkm], axis=1)
    raw = np.nansum(tpf.flux.value[:, ap], axis=1) - bkg*ap.sum()
    t = tpf.time.value + 2167.0; q = tpf.quality == 0; m = q & np.isfinite(raw) & (raw > 0) & np.isfinite(t)
    t, raw = t[m], raw[m]/np.nanmedian(raw[m]); LC[s] = (t, raw)
    n = np.round((np.median(t) - T0)/P0)
    out['sectors'][s] = dict(npts=int(len(t)), cadence_min=float(np.median(np.diff(t))*1440), target_col_row=[float(cx), float(cy)], aperture_npix=int(ap.sum()),
        target_flux_fraction_in_ap=float(frac), expected_depth_ppm=float(dep*frac*1e6), n_epochs_from_K2=int(n), sigma_Tpred_h=float(np.hypot(eT0, n*eP)*24))
def detr(t, f):
    fl = flatten(t, f, method='biweight', window_length=0.75); ok = np.isfinite(fl); t, fl = t[ok], fl[ok]
    md = np.median(fl); sd = 1.4826*np.median(np.abs(fl-md)); ok = fl < md + 4*sd; return t[ok], fl[ok]
FL = {s: detr(*LC[s]) for s in LC}
for s in FL: out['sectors'][s]['rms_ppm'] = float(np.std(FL[s][1])*1e6)
nmax = max(abs(np.round((FL[s][0].max() - T0)/P0)) for s in FL); dPs = (T14/4)/nmax; Pgrid = np.arange(P0 - 4*eP, P0 + 4*eP + dPs/2, dPs)
def search(fl, Tx, excl=None):
    res = []
    for P in Pgrid:
        num = 0.0; den = 0.0; per = {}
        for s, (t, f) in fl.items():
            ph = ((t - Tx + 0.5*P) % P) - 0.5*P; ii = np.abs(ph) < 0.5*T14
            if excl is not None: ii &= excl[s]
            if ii.sum() < 3: continue
            oo = np.abs(ph) > T14
            if excl is not None: oo &= excl[s]
            d = np.median(f[oo]) - np.mean(f[ii]); e = np.std(f[oo])/np.sqrt(ii.sum()); w = 1/e**2/out['sectors'][s]['target_flux_fraction_in_ap']**2
            dc, ec = d/out['sectors'][s]['target_flux_fraction_in_ap'], e/out['sectors'][s]['target_flux_fraction_in_ap']   # undiluted
            num += dc/ec**2; den += 1/ec**2; per[s] = (float(d*1e6), float(e*1e6), int(ii.sum()))
        if den == 0: res.append((P, np.nan, np.nan, per)); continue
        dd = num/den; ee = 1/np.sqrt(den); res.append((P, dd, ee, per))
    return res
real = search(FL, T0); snr = np.array([r[1]/r[2] for r in real])
i0 = int(np.argmin(np.abs(Pgrid - P0))); ib = int(np.nanargmax(snr))
out['grid'] = dict(nP=int(len(Pgrid)), dP=float(dPs), P_range=[float(Pgrid[0]), float(Pgrid[-1])])
out['nominal'] = dict(P=float(Pgrid[i0]), undiluted_depth_ppm=float(real[i0][1]*1e6), err_ppm=float(real[i0][2]*1e6), snr=float(snr[i0]), per_sector=real[i0][3])
out['best'] = dict(P=float(Pgrid[ib]), undiluted_depth_ppm=float(real[ib][1]*1e6), err_ppm=float(real[ib][2]*1e6), snr=float(snr[ib]), per_sector=real[ib][3])
out['marginal'] = dict(min_snr=float(np.nanmin(snr)), median_err_ppm=float(np.nanmedian([r[2] for r in real])*1e6),
    frac_grid_depth_consistent_with_K2_within_2sig=float(np.nanmean([abs(r[1]-dep) < 2*r[2] for r in real])),
    snr_of_K2_depth_per_trial_median=float(np.nanmedian([dep/r[2] for r in real])))
rng = np.random.default_rng(5)
excl = {s: np.ones(len(FL[s][0]), bool) for s in FL}
exclr = {}
for s in LC:
    exclr[s] = np.ones(len(LC[s][0]), bool)
nullmax = []; 
for j in range(200):
    Tx = T0 + rng.uniform(0.15, 0.85)*P0
    rr = search(FL, Tx); nullmax.append(np.nanmax([r[1]/r[2] for r in rr]))
nullmax = np.array(nullmax); thr = float(np.percentile(nullmax, 99))
out['null'] = dict(n=200, max_snr_p50_p95_p99=np.percentile(nullmax, [50, 95, 99]).tolist(), frac_null_ge_real_best=float(np.mean(nullmax >= snr[ib])))
inj = []
for j in range(100):
    Tx = T0 + rng.uniform(0.15, 0.85)*P0; FLi = {}
    for s, (t, f) in LC.items():
        pm = batman.TransitParams(); pm.t0 = Tx; pm.per = P0; pm.rp = k*np.sqrt(out['sectors'][s]['target_flux_fraction_in_ap']); pm.a = aRs; pm.inc = np.degrees(np.arccos(b/aRs))
        pm.ecc = 0; pm.w = 90; pm.limb_dark = 'quadratic'; pm.u = [0.3, 0.3]
        cad = out['sectors'][s]['cadence_min']/1440
        mdl = batman.TransitModel(pm, t, supersample_factor=5, exp_time=cad).light_curve(pm); FLi[s] = detr(t, f*mdl)
    rr = search(FLi, Tx); sn = np.array([r[1]/r[2] for r in rr]); jb = int(np.nanargmax(sn)); jn = int(np.argmin(np.abs(Pgrid - P0)))
    inj.append(dict(T0=float(Tx), max_snr=float(sn[jb]), P_best=float(Pgrid[jb]), snr_at_true_P=float(sn[jn]), depth_at_true_P_ppm=float(rr[jn][1]*1e6), err_at_true_P_ppm=float(rr[jn][2]*1e6)))
ms = np.array([x['max_snr'] for x in inj])
out['injection'] = dict(n=len(inj), threshold_null_p99=thr, frac_recovered=float(np.mean(ms > thr)), max_snr_p16_p50_p84=np.percentile(ms, [16, 50, 84]).tolist(),
    frac_Pbest_within_1step=float(np.mean([abs(x['P_best']-P0) <= dPs for x in inj])), depth_at_true_P_median_ppm=float(np.median([x['depth_at_true_P_ppm'] for x in inj])),
    snr_at_true_P_median=float(np.median([x['snr_at_true_P'] for x in inj])))
out['real_vs_injected'] = dict(real_best_snr=float(snr[ib]), frac_inj_maxsnr_le_real=float(np.mean(ms <= snr[ib])), real_nominal_snr=float(snr[i0]),
    frac_inj_snr_at_true_P_le_real_nominal=float(np.mean([x['snr_at_true_P'] <= snr[i0] for x in inj])))
for s in FL: np.savez(f'{D}/tess/s{s}_{APN}_lc.npz', t=FL[s][0], f=FL[s][1], t_raw=LC[s][0], f_raw=LC[s][1])
np.savez(f'{D}/tess/grid_{APN}.npz', P=Pgrid, depth=np.array([r[1] for r in real]), err=np.array([r[2] for r in real]), nullmax=nullmax, injmax=ms)
json.dump(dict(out, injections=inj), open(f'{D}/tess/tess_check_{APN}.json', 'w'), indent=1, default=float)
print(json.dumps({k: v for k, v in out.items()}, indent=1, default=float))
