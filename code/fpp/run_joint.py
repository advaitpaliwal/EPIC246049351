# Joint TRICERATOPS(K2)+MOLUSC FPP for EPIC 246049351 (C12). Adapted from ../../248780202/fpp_joint/scripts/run_joint.py (single planet).
# Target built with search_radius=10 K2 px (39.8 arcsec), NO stars removed. Light curve: flattened EVEREST (prep3.py), phase-folded on the prior-fit ephemeris and
# binned to 15 min within +-0.3 d; flux_err = std of binned points with |phase|>T14. Depth passed = prior-fit median model mid-transit depth (result_prior_everest.json).
# Aperture = lightkurve threshold-3 mask of the C12 TPF (same as fresh photometry). N=1e6, exptime 29.4 min, nsamples 20, seed 2000+k, parallel=False.
# Variants: gr = MOLUSC RUWE+Gaia companions, no contrast curve, filt='Kepler' (as 248780202);
#           grs = MOLUSC RUWE+Gaia+speckle companions AND the 832 nm Zorro curve passed to TRICERATOPS (filt='TESS', the band TRICERATOPS uses for visible/red curves);
#           none = no MOLUSC, no curve (reference); cc = curve only, no MOLUSC (reference).
# usage: run_joint.py VARIANT K   (VARIANT=init builds the target pkl)
import sys, os, json, pickle, time, glob, types, numpy as np, warnings; warnings.filterwarnings('ignore')
D = 'k2-deep/246049351/fpp_joint'; os.chdir(D)
import triceratops.triceratops as tr, triceratops.marginal_likelihoods as ml, lightkurve as lk
class _F:
    def __init__(s, d): s.d = d
    def __getitem__(s, k):
        if isinstance(k, str): return types.SimpleNamespace(values=s.d[k].to_numpy(copy=True))
        return _F(s.d[np.asarray(k)])
_rc = ml.read_csv
ml.read_csv = lambda f, *a, **kw: _F(_rc(f, *a, **kw)) if 'molusc' in str(f) else _rc(f, *a, **kw)
var, k = sys.argv[1], int(sys.argv[2]); exptime = 29.4/1440; epic = 246049351; pk = f'out/target_{epic}.pkl'
if not os.path.exists(pk):
    if os.environ.get('GAIAPATCH'): sys.path.insert(0, D + '/scripts'); import gaia_patch
    tg = tr.target(ID=epic, sectors=np.array([12]), mission='K2', lightkurve_cache_dir=D + '/lkcache', search_radius=10)
    pickle.dump(tg, open(pk, 'wb')); tg.stars.to_csv(f'out/stars_{epic}.csv', index=False)
if var == 'init': print(pickle.load(open(pk, 'rb')).stars.to_string()); sys.exit()
tg = pickle.load(open(pk, 'rb'))
R = json.load(open('../fit/result_prior_everest.json'))
P, T0, dur, depth = R['P_d'][0], R['T0_BKJD'][0], R['T14_h'][0]/24, R['depth_ppm'][0]*1e-6
tpp = glob.glob(f'k2-deep/{epic}/lc/mastDownload/K2/ktwo{epic}-c12_lc/*lpd-targ.fits.gz')[0]
tp = lk.read(tpp); ap = tp.create_threshold_mask(threshold=3, reference_pixel='center')
rows, cols = np.nonzero(ap); apx = np.column_stack([cols + tp.column, rows + tp.row])
z = np.load('../lc/everest.npz'); t, f = z['t'], z['f']
ph = t - T0 - np.round((t - T0)/P)*P; w = np.abs(ph) < 0.3
edges = np.arange(-0.3, 0.3 + 1e-9, 15/1440); idx = np.digitize(ph[w], edges)
tb = np.array([ph[w][idx == i].mean() for i in np.unique(idx)]); fb = np.array([f[w][idx == i].mean() for i in np.unique(idx)])
err = float(np.std(fb[np.abs(tb) > dur]))
tg.calc_depths(tdepth=depth, all_ap_pixels=[apx])
mf = {'gr': D + '/molusc_out/molusc_gr_for_trc_N1e6.csv', 'grs': D + '/molusc_out/molusc_grs_for_trc_N1e6.csv'}.get(var)
cc = D + '/cc832_arcsec_dmag.csv' if var in ('grs', 'cc') else None; filt = 'TESS' if cc else 'Kepler'
np.random.seed(2000 + k); t0 = time.time()
tg.calc_probs(time=tb, flux_0=fb, flux_err_0=err, P_orb=P, contrast_curve_file=cc, filt=filt, N=1000000, parallel=False,
              exptime=exptime, nsamples=20, verbose=0, molusc_file=mf)
tag = f'{var}_{k}'; tg.probs.to_csv(f'out/probs_{tag}.csv', index=False)
json.dump(dict(epic=epic, variant=var, run=k, seed=2000+k, N=1000000, exptime_min=exptime*1440, P=P, T0=T0, depth_ppm=depth*1e6, err_ppm=err*1e6, nbins=len(tb),
               molusc=mf, contrast_curve=cc, filt=filt, FPP=float(tg.FPP), NFPP=float(tg.NFPP), n_stars=len(tg.stars), sec=time.time() - t0, aperture=apx.tolist()), open(f'out/res_{tag}.json', 'w'), indent=1)
print(tag, tg.FPP, tg.NFPP, len(tg.stars), round(time.time() - t0), flush=True)
