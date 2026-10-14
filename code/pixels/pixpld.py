"""Per-pixel transit amplitude map by linear regression (PLD-style, Deming+2015 / EVEREST idea).
Each pixel: flux = transit_template*A + B-spline(time, 1.0 d knots) + PCA(pixel fractions, 15) + squares of top 5 PCs.
usage: pixpld.py TPF P T0 T14_H RpRs aRs b RA DEC OUT"""
import sys, json, numpy as np, warnings; warnings.filterwarnings('ignore')
import lightkurve as lk, batman, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.interpolate import BSpline
tpp, P, T0, T14, k, aR, b, RA, DEC, OUT = sys.argv[1], *map(float, sys.argv[2:10]), sys.argv[10]
tp = lk.read(tpp, quality_bitmask='default')
t = tp.time.value; F = tp.flux.value.astype(float)
good = np.isfinite(t) & np.all(np.isfinite(F.reshape(len(t), -1)[:, np.isfinite(np.nanmedian(F, 0)).ravel()]), axis=1)
t, F = t[good], F[good]
med = np.nanmedian(F, 0); pix = np.isfinite(med) & (med > 0)
big = tp.create_threshold_mask(threshold=1, reference_pixel='center') & pix
Fb = F[:, big]; frac = Fb/Fb.sum(1, keepdims=True); frac = (frac - frac.mean(0))/frac.std(0)
U, S, Vt = np.linalg.svd(frac, full_matrices=False); pcs = U[:, :15]
reg = np.column_stack([pcs, pcs[:, :5]**2])
# time spline
kn = np.arange(t.min()-3, t.max()+3, 1.0); deg = 3
Bs = np.column_stack([BSpline(kn, np.eye(len(kn)-deg-1)[i], deg, extrapolate=False)(t) for i in range(len(kn)-deg-1)])
Bs = np.nan_to_num(Bs); Bs = Bs[:, Bs.sum(0) > 0.5]
pm = batman.TransitParams(); pm.t0=T0; pm.per=P; pm.rp=k; pm.a=aR; pm.inc=np.degrees(np.arccos(b/aR)); pm.ecc=0; pm.w=90; pm.limb_dark='quadratic'; pm.u=[0.3,0.3]
tm = batman.TransitModel(pm, t, supersample_factor=15, exp_time=1765.5/86400).light_curve(pm) - 1
tmpl = tm/abs(tm.min())   # unit-depth template (negative)
X = np.column_stack([tmpl, Bs, reg])
A = np.full(med.shape, np.nan); E = np.full(med.shape, np.nan)
for (r, c) in zip(*np.where(pix)):
    y = F[:, r, c]
    beta, *_ = np.linalg.lstsq(X, y, rcond=None); res = y - X@beta
    # robust: clip 5 sigma and refit
    s = 1.4826*np.median(np.abs(res)); m = np.abs(res) < 5*s
    beta, *_ = np.linalg.lstsq(X[m], y[m], rcond=None); res = y[m] - X[m]@beta; s2 = np.var(res)
    cov = s2*np.linalg.pinv(X[m].T@X[m]); A[r, c] = -beta[0]; E[r, c] = np.sqrt(cov[0, 0])
# A: transit flux deficit per pixel (e-/s). Compare to PSF image
from astropy.wcs import WCS
cx, cy = tp.wcs.all_world2pix([[RA, DEC]], 0)[0]
yy, xx = np.mgrid[:med.shape[0], :med.shape[1]]
sn = A/E
# fit A = s*med (target PSF scaled); chi2 of that vs free
mm = pix & (np.hypot(xx-cx, yy-cy) < 4)
w = 1/E[mm]**2; sfit = np.sum(w*A[mm]*med[mm])/np.sum(w*med[mm]**2); e_s = 1/np.sqrt(np.sum(w*med[mm]**2))
chi2_psf = float(np.sum(w*(A[mm]-sfit*med[mm])**2)); chi2_null = float(np.sum(w*A[mm]**2)); dof = int(mm.sum()-1)
# centroid of A via weighted fit of shifted PSF: grid search of integer-free shift using interpolation of med
from scipy.ndimage import shift as ndshift
best = (1e99, 0, 0)
for dx in np.arange(-2, 2.01, 0.1):
    for dy in np.arange(-2, 2.01, 0.1):
        mt = ndshift(np.nan_to_num(med), (dy, dx), order=1)
        ss = np.sum(w*A[mm]*mt[mm])/np.sum(w*mt[mm]**2); ch = np.sum(w*(A[mm]-ss*mt[mm])**2)
        if ch < best[0]: best = (ch, dx, dy)
chis = []
grid = []
for dx in np.arange(-2, 2.01, 0.1):
    for dy in np.arange(-2, 2.01, 0.1):
        mt = ndshift(np.nan_to_num(med), (dy, dx), order=1)
        ss = np.sum(w*A[mm]*mt[mm])/np.sum(w*mt[mm]**2); grid.append((dx, dy, np.sum(w*(A[mm]-ss*mt[mm])**2)))
grid = np.array(grid); in1 = grid[grid[:, 2] < best[0] + 2.30]   # 1-sigma 2D region
ap = tp.create_threshold_mask(threshold=3, reference_pixel='center') & pix
res = dict(tpf=tpp, n_cad=int(len(t)), depth_ap_ppm=float(np.nansum(A[ap])/np.nansum(med[ap])*1e6),
           depth_ap_err_ppm=float(np.sqrt(np.nansum(E[ap]**2))/np.nansum(med[ap])*1e6),
           psf_scale_depth_ppm=float(sfit*1e6), psf_scale_err_ppm=float(e_s*1e6), chi2_psf=chi2_psf, chi2_null=chi2_null, dof=dof,
           best_shift_pix=[float(best[1]), float(best[2])], shift_1sig_range_x=[float(in1[:, 0].min()), float(in1[:, 0].max())],
           shift_1sig_range_y=[float(in1[:, 1].min()), float(in1[:, 1].max())], best_shift_arcsec=float(np.hypot(best[1], best[2])*3.98),
           max_pixel_sn=float(np.nanmax(sn)), argmax_sn_colrow=[int(np.unravel_index(np.nanargmax(sn), sn.shape)[1]), int(np.unravel_index(np.nanargmax(sn), sn.shape)[0])],
           target_colrow=[float(cx), float(cy)], brightest_colrow=[int(np.unravel_index(np.nanargmax(med), med.shape)[1]), int(np.unravel_index(np.nanargmax(med), med.shape)[0])])
json.dump(res, open(OUT + '.json', 'w'), indent=1); np.savez(OUT + '.npz', A=A, E=E, med=med, ap=ap)
fig, ax = plt.subplots(1, 4, figsize=(17, 4))
for a, im, ti in zip(ax, [med, A, sn, np.where(sn > 1, A/med*1e6, np.nan)], ['median flux', 'transit deficit A (e/s)', 'A / sigma', 'per-pixel depth ppm (A/sig>1)']):
    c = a.imshow(im, origin='lower'); plt.colorbar(c, ax=a, fraction=0.046); a.plot(cx, cy, 'r*', ms=12)
    a.plot(cx + best[1], cy + best[2], 'wx', ms=10); a.contour(ap, levels=[0.5], colors='w', linewidths=0.8); a.set_title(ti)
plt.suptitle(f"{OUT.split('/')[-1]}: depth(ap) {res['depth_ap_ppm']:.0f}+-{res['depth_ap_err_ppm']:.0f} ppm; PSF shift {res['best_shift_arcsec']:.1f}\"")
plt.tight_layout(); plt.savefig(OUT + '.png', dpi=90)
for kk, v in res.items(): print(kk, v)
