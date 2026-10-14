"""Copied from ../../220179502/scripts/pix_null.py; adds an optional second-planet mask (POT P2 T02 T14_2) so the other planet's transits are excluded from in-transit and reference cadences. Null test: centroid offset of each no-transit fake-epoch difference image (same estimator) vs the real one.
usage: pix2.py TPF P T0 T14_H RA DEC OUTJSON [P2 T02 T14_2_H]"""
import sys, json, numpy as np, warnings; warnings.filterwarnings('ignore')
import lightkurve as lk, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
tp_path, P, T0, T14, RA, DEC, OUT = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])/24, float(sys.argv[5]), float(sys.argv[6]), sys.argv[7]
tp = lk.read(tp_path, quality_bitmask='hard')
P2, T02, T142 = (float(sys.argv[8]), float(sys.argv[9]), float(sys.argv[10])/24) if len(sys.argv) > 10 else (1e9, 0.0, 0.0)
t = tp.time.value; F = tp.flux.value; pc1 = tp.pos_corr1; pc2 = tp.pos_corr2
ok = np.isfinite(t) & np.isfinite(pc1) & np.isfinite(pc2) & np.all(np.isfinite(F.reshape(len(t), -1)) | True, axis=1)
ok &= np.abs(((t - T02 + 0.5*P2) % P2) - 0.5*P2) > T142
t, F, pc1, pc2 = t[ok], F[ok], pc1[ok], pc2[ok]
# remove slow trend per pixel? local matching within 1.5 d keeps it small
def diffimg(T0x):
    ph = ((t - T0x + 0.5*P) % P) - 0.5*P
    phr = ((t - T0 + 0.5*P) % P) - 0.5*P
    intr = np.where(np.abs(ph) < 0.35*T14)[0]
    oot_ok = (np.abs(phr) > 0.75*T14)
    D = []; O = []
    for i in intr:
        cand = np.where(oot_ok & (np.abs(t - t[i]) < 1.5) & (np.abs(ph) > 0.75*T14))[0]
        if len(cand) < 8: continue
        d = np.hypot(pc1[cand]-pc1[i], pc2[cand]-pc2[i]); sel = cand[np.argsort(d)[:8]]
        o = np.nanmean(F[sel], axis=0); D.append(o - F[i]); O.append(o)
    if len(D) < 3: return None, None, 0
    return np.nanmean(D, axis=0), np.nanmean(O, axis=0), len(D)
diff, oot, n = diffimg(T0)
rng = np.random.default_rng(1); fakes = []
for k in range(200):
    off = rng.uniform(0.15, 0.85)*P
    if abs(((off + 0.5*P) % P) - 0.5*P) < 3*T14: continue
    d, _, nn = diffimg(T0 + off)
    if d is not None and nn >= 0.6*n: fakes.append(d)
    if len(fakes) >= 60: break
fakes = np.array(fakes); sig = np.nanstd(fakes, axis=0)
snr = diff/sig; depth_map = diff/oot
# aperture (threshold) and target pixel from WCS
ap = tp.create_threshold_mask(threshold=3, reference_pixel='center')
from astropy.wcs import WCS
w = tp.wcs; cx, cy = w.all_world2pix([[RA, DEC]], 0)[0]   # column, row in cutout coordinates
yy, xx = np.mgrid[:diff.shape[0], :diff.shape[1]]
def cen(img, m):
    im = np.where(m, img, 0); im = np.clip(im, 0, None); s = im.sum(); return (xx*im).sum()/s, (yy*im).sum()/s
m3 = (np.hypot(xx-cx, yy-cy) < 3.5)
dcx, dcy = cen(diff, m3); ocx, ocy = cen(oot, m3)
# bootstrap centroid error from fakes added to diff
bc = []
for fk in fakes: bc.append(cen(diff + (fk - np.nanmean(fakes, axis=0)), m3))
real = np.hypot(dcx-ocx, dcy-ocy)*3.98
bc = np.array(bc); e_dc = np.nanstd(bc, axis=0)
fo = np.array([np.hypot(*(np.array(cen(fk, m3)) - [ocx, ocy]))*3.98 for fk in fakes])
fs = np.array([np.nansum(fk[ap])/np.nansum(oot[ap])*1e6 for fk in fakes])
r = dict(real_offset_arcsec=float(real), real_depth_ap_ppm=float(np.nansum(diff[ap])/np.nansum(oot[ap])*1e6), n_fake=int(len(fakes)), fake_offset_arcsec_pct=dict(zip(['p16','p50','p84','p95'], np.percentile(fo, [16,50,84,95]).tolist())),
         frac_fake_offset_ge_real=float(np.mean(fo >= real)), fake_depth_ap_std_ppm=float(np.std(fs)), n_in_cadences=int(n), target_pix_colrow=[float(cx), float(cy)], oot_centroid=[float(ocx), float(ocy)], diff_centroid=[float(dcx), float(dcy)], offset_pix=[float(dcx-ocx), float(dcy-ocy)], diff_centroid_err_pix=e_dc.tolist(), offset_sigma=float(np.hypot((dcx-ocx)/e_dc[0], (dcy-ocy)/e_dc[1])), offset_from_target_pix=[float(dcx-cx), float(dcy-cy)], depth_ap_snr=float(np.nansum(diff[ap])/np.nansum(oot[ap])*1e6/np.std(fs)), npix_ap=int(ap.sum()), max_pixel_snr=float(np.nanmax(snr)))
np.savez(OUT.replace('.json', '.npz'), diff=diff, oot=oot, sig=sig, snr=snr, ap=ap, cx=cx, cy=cy)
json.dump(r, open(OUT, 'w'), indent=1); print(r)
