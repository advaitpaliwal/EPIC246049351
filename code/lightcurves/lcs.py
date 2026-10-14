"""Light-curve loaders + detrending for K2 deep vet. Glue around astropy/wotan/lightkurve."""
import numpy as np, glob, os, warnings; warnings.filterwarnings('ignore')
from astropy.io import fits
from wotan import flatten
B = os.path.expanduser('k2-deep')
def _f(epic, pat): 
    g = sorted(glob.glob(f'{B}/{epic}/lc/mastDownload/**/{pat}', recursive=True)); return g
def everest(epic):
    p = _f(epic, 'hlsp_everest*.fits')[0]; d = fits.getdata(p, 1)
    t, f, q = d['TIME'], d['FLUX'], d['QUALITY']
    m = np.isfinite(t) & np.isfinite(f) & ((q & ~((1 << 23) | (1 << 24))) == 0)
    return t[m], f[m]/np.nanmedian(f[m])
def k2sff(epic):
    p = _f(epic, 'hlsp_k2sff*.fits')[0]; d = fits.getdata(p, 1)  # BESTAPER ext 1
    t, f = d['T'], d['FCOR']; m = np.isfinite(t) & np.isfinite(f)
    return t[m], f[m]/np.nanmedian(f[m])
def k2sc(epic):
    p = _f(epic, 'hlsp_k2sc*.fits')
    if not p: return None
    d = fits.getdata(p[0], 1); t = d['time']; f = d['flux']  # K2SC 'flux' column is already systematics-corrected (1-99% span 0.6%)
    q = d['quality'] if 'quality' in d.names else np.zeros(len(t), int)
    m = np.isfinite(t) & np.isfinite(f) & (q == 0) & (d['mflags'] == 0)
    return t[m], f[m]/np.nanmedian(f[m])
def flat(t, f, P, T0, dur_d, window=0.75):
    intr = np.abs(((t - T0 + 0.5*P) % P) - 0.5*P) < 0.75*dur_d
    ff = f.copy(); 
    trend = flatten(t, np.where(intr, np.nan, f), method='biweight', window_length=window, return_trend=True)[1]
    # interpolate trend over masked transits
    ok = np.isfinite(trend); trend = np.interp(t, t[ok], trend[ok])
    fl = f/trend; m = np.isfinite(fl)
    t, fl = t[m], fl[m]
    oot = ~(np.abs(((t - T0 + 0.5*P) % P) - 0.5*P) < 0.75*dur_d)
    med = np.median(fl[oot]); sd = 1.4826*np.median(np.abs(fl[oot]-med))
    keep = (fl < med + 4*sd) & (fl > med - 8*sd) | ~oot
    keep &= (fl > med - 20*sd)
    return t[keep], fl[keep]
def phase(t, P, T0): return ((t - T0 + 0.5*P) % P) - 0.5*P
def boxdepth(t, f, P, T0, dur_d, frac=0.6):
    ph = phase(t, P, T0); i = np.abs(ph) < frac*dur_d/2; o = (np.abs(ph) > dur_d) & (np.abs(ph) < 4*dur_d)
    sd = np.std(f[o]); d = np.median(f[o]) - np.mean(f[i]); e = sd/np.sqrt(i.sum())
    return d, e, int(i.sum())
def tpf(epic, seg='c102'):
    import lightkurve as lk
    p = [x for x in _f(epic, '*lpd-targ.fits.gz') if seg in x or seg == '']
    return lk.read(p[-1], quality_bitmask='default')
def fresh(epic, P, T0, dur_d, seg='c102', thresh=3, windows=20):
    """Fresh aperture photometry from the TPF: threshold aperture, SAP, lightkurve SFFCorrector with transits masked."""
    import lightkurve as lk
    tp = tpf(epic, seg)
    ap = tp.create_threshold_mask(threshold=thresh, reference_pixel='center')
    lc = tp.to_lightcurve(aperture_mask=ap).remove_nans()
    lc = lc[np.isfinite(lc.centroid_col.value) & np.isfinite(lc.centroid_row.value)]
    t = lc.time.value
    intr = np.abs(phase(t, P, T0)) < 0.75*dur_d
    corr = lk.SFFCorrector(lc).correct(windows=windows, bins=15, cadence_mask=~intr, restore_trend=True)
    f = corr.flux.value; m = np.isfinite(f)
    return t[m], f[m]/np.nanmedian(f[m]), ap
