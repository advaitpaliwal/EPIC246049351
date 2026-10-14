"""Per-epoch breakdown of the TESS grid result (tess3.py output): for the nominal and best-SNR periods, each predicted transit inside each sector:
n in-transit cadences, single-epoch box depth (in |phase|<0.5*T14 vs that epoch's T14<|phase|<0.5 d), mean background (faintest-30% pixel median) in vs out.
usage: tess_diag.py APERTURE"""
import sys, json, glob, numpy as np, lightkurve as lk, warnings; warnings.filterwarnings('ignore')
D = 'k2-deep/246049351'; A = sys.argv[1]; R = json.load(open(f'{D}/tess/tess_check_{A}.json'))  # A may carry the _bkgcut suffix; T0 = R['T0_BKJD']; T14 = R['T14_h']/24
out = {}
for lab in ['nominal', 'best']:
    P = R[lab]['P']; out[lab] = dict(P=P, epochs=[])
    for s in R['sectors']:
        if not isinstance(R['sectors'][s], dict): continue
        z = np.load(f'{D}/tess/s{s}_{A}_lc.npz'); t, f = z['t'], z['f']
        tp = lk.read(glob.glob(f'{D}/tess/tesscut/**/*s{int(s):04d}*.fits', recursive=True)[0]); med = np.nanmedian(tp.flux.value, 0)
        bk = np.nanmedian(tp.flux.value[:, med < np.nanpercentile(med, 30)], axis=1); tb = tp.time.value + 2167.0
        ep = np.round((t - T0)/P).astype(int); ph = t - T0 - ep*P
        for e in np.unique(ep[np.abs(ph) < 0.5*T14]):
            i = (ep == e) & (np.abs(ph) < 0.5*T14); o = (ep == e) & (np.abs(ph) > T14) & (np.abs(ph) < 0.5)
            if o.sum() < 5: continue
            tm = T0 + e*P; bi = np.abs(tb - tm) < 0.5*T14; bo = (np.abs(tb - tm) > T14) & (np.abs(tb - tm) < 0.5)
            out[lab]['epochs'].append(dict(sector=s, epoch=int(e), tmid_BTJD=round(float(tm - 2167), 4), n_in=int(i.sum()), depth_ppm=round(float((np.median(f[o]) - np.mean(f[i]))*1e6)),
                err_ppm=round(float(np.std(f[o])/np.sqrt(i.sum())*1e6)), bkg_in_minus_out=round(float(np.nanmean(bk[bi]) - np.nanmedian(bk[bo])), 2), bkg_out_median=round(float(np.nanmedian(bk[bo])), 2)))
    for x in out[lab]['epochs']: print(lab, x)
json.dump(out, open(f'{D}/tess/tess_diag_{A}.json', 'w'), indent=1)
