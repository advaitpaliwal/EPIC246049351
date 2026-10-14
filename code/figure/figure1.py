# Three-panel figure for EPIC 246049351 (folded / unfolded campaign / four single transits), 300 dpi, single-column width.
# Same layout as the EPIC 220179502 figure: the raw light curve gets its own strip with its own scale (it ran off the
# axes), transit ticks sit below the detrended data, single-transit panels are scaled to their data and labelled above the axes.
# Data: lc/everest.npz (raw EVEREST t_raw,f_raw; wotan biweight 0.75 d, transits masked -> t,f). Model: median posterior, fit/result_prior_everest.json.
# usage: python -I figure1.py K2DEEP_DIR OUT_PNG
import sys, json, numpy as np, batman, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
D, OUT = sys.argv[1], sys.argv[2]; EPIC, CAMP = "246049351", 12
R = json.load(open(f'{D}/fit/result_prior_everest.json')); z = np.load(f'{D}/lc/everest.npz')
P, T0, T14 = R['P_d'][0], R['T0_BKJD'][0], R['T14_h'][0]/24
q1, q2 = R['q1'][0], R['q2'][0]
p = batman.TransitParams(); p.t0 = T0; p.per = P; p.rp = R['RpRs'][0]; p.a = R['a_Rs'][0]; p.inc = R['inc_deg'][0]
p.ecc = 0.; p.w = 90.; p.limb_dark = 'quadratic'; p.u = [2*np.sqrt(q1)*q2, np.sqrt(q1)*(1-2*q2)]
def model(t): return batman.TransitModel(p, t, supersample_factor=15, exp_time=1765.5/86400).light_curve(p)
tr_, fr = z['t_raw'], z['f_raw']; t, fd = z['t'], z['f']
trend = np.interp(t, tr_, fr)/fd
oot = np.abs(((t - T0 + 0.5*P) % P) - 0.5*P) > 0.75*T14
sig = 1.4826*np.median(np.abs(fd[oot] - np.median(fd[oot])))
def per_transit():
    rows = []
    for n in np.unique(np.round((t - T0)/P)):
        tc = T0 + n*P; w = np.abs(t - tc) < 0.3; inn = w & (np.abs(t - tc) < T14/2)
        if inn.sum() < 3 or (w & ~inn).sum() < 10: continue
        tm = model(t[w]) - 1; y = fd[w] - 1; y = y - np.median(y[np.abs(t[w] - tc) > T14])
        A = np.sum(tm**2); dep = np.sum(tm*y)/A*(-tm.min()); err = sig/np.sqrt(A)*(-tm.min())
        rows.append((n, tc, dep*1e6, err*1e6, dep/err, inn.sum()))
    return np.array(rows)
pt = per_transit(); sel = pt[np.isin(pt[:, 0], [-7, -3, 4, 9])]  # three typical transits and the shallow epoch -3 (thruster firing 1.5 h before mid-transit)
print('sigma per cadence (MAD, OOT) ppm %.0f' % (sig*1e6)); print('full-coverage transits', len(pt))
for r in pt: print('epoch %d tc %.4f depth %.0f err %.0f S/N %.2f n_in %d' % tuple(r))
print('zoomed epochs', sel[:, 0].astype(int))
plt.rcParams.update({'font.size': 7, 'axes.linewidth': 0.6, 'xtick.direction': 'in', 'ytick.direction': 'in',
                     'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.top': True, 'ytick.right': True})
dep = R['depth_ppm'][0]
def binned(x, y, wd):
    e = np.arange(x.min(), x.max() + wd, wd); i = np.digitize(x, e); u = np.unique(i)
    return np.array([x[i == j].mean() for j in u]), np.array([y[i == j].mean() for j in u]), np.array([sig/np.sqrt((i == j).sum()) for j in u])
def lab(ax, s): ax.text(-0.2, 1.0, s, transform=ax.transAxes, fontsize=8, fontweight='bold', va='bottom')
fig = plt.figure(figsize=(3.45, 6.6))
gs = fig.add_gridspec(7, 1, height_ratios=[1.25, 0.42, 0.62, 0.78, 0.36, 1.6, 0.08], hspace=0.0)
# (a) folded
ax = fig.add_subplot(gs[0]); ph = ((t - T0 + 0.5*P) % P - 0.5*P)*24; y = (fd - 1)*1e6; w = np.abs(ph) < 6
hb, yb, eb = binned(ph[w], y[w], 0.5)
ax.plot(ph[w], y[w], '.', c='0.78', ms=1.5, rasterized=True, label='29.4-min cadences')
ax.errorbar(hb, yb, eb, fmt='o', c='k', ms=2, lw=0.6, capsize=0, label='30-min bins')
tt = np.linspace(-0.25, 0.25, 800); p.t0 = 0.; ax.plot(tt*24, (model(tt) - 1)*1e6, c='C1', lw=1.1, label='median model'); p.t0 = T0
ax.set_xlim(-6, 6); ax.set_ylim(-5600, 3400); ax.set_ylabel('relative flux [ppm]'); ax.set_xlabel('hours from mid-transit', labelpad=1)
ax.legend(fontsize=5.6, loc='lower left', frameon=False, handletextpad=0.3, borderaxespad=0.3)
ax.set_title('EPIC %s, K2 C%d EVEREST, P = %.5f d' % (EPIC, CAMP, P), fontsize=6.5, pad=3); lab(ax, 'a')
# (b) unfolded campaign: raw + trend (top strip), detrended (bottom strip)
ax1 = fig.add_subplot(gs[2]); ax2 = fig.add_subplot(gs[3], sharex=ax1)
ax1.plot(tr_, (fr - 1)*100, '.', c='0.75', ms=0.6, rasterized=True); ax1.plot(t, (trend - 1)*100, c='C1', lw=0.6)
rl = np.percentile((fr - 1)*100, [0.2, 99.8]); ax1.set_ylabel('raw [%]'); ax1.set_ylim(rl[0] - 0.1*np.ptp(rl), rl[1] + 0.35*np.ptp(rl)); ax1.set_yticks([-1, 0, 1]); plt.setp(ax1.get_xticklabels(), visible=False)
ax1.text(0.02, 0.92, 'EVEREST and 0.75-d biweight trend', transform=ax1.transAxes, va='top', fontsize=5.6)
ax2.plot(t, (fd - 1)*1e6, '.', c='0.78', ms=0.6, rasterized=True)
xa, ya, _ = binned(t, fd - 1, 2/24); ax2.plot(xa, ya*1e6, '.', c='k', ms=1.2)
tc = T0 + np.arange(np.ceil((t.min() - T0)/P), np.floor((t.max() - T0)/P) + 1)*P
y0, y1 = -6600, -5300
for c in tc: ax2.plot([c, c], [y0, y1], c='C3', lw=0.6)
for c in sel[:, 1]: ax2.plot([c, c], [y0, y1], c='C3', lw=1.8)
ax2.set_ylim(-7000, 7000); ax2.set_yticks([-4000, 0, 4000]); ax2.set_xlim(tr_.min() - 0.5, tr_.max() + 0.5); ax2.set_ylabel('flux [ppm]')
ax2.text(0.02, 0.95, 'detrended; black: 2-h bins; ticks: transits (thick: panel c)', transform=ax2.transAxes, va='top', fontsize=5.6)
ax2.set_xlabel('time [BKJD = BJD$_{\\rm TDB}$ $-$ 2454833]', labelpad=1); lab(ax1, 'b')
ax1.set_title('Campaign %d, 29.4-min cadence' % CAMP, fontsize=6.5, pad=3)
# (c) four single transits
sub = gs[5].subgridspec(2, 2, hspace=0.3, wspace=0.07)
axs = [fig.add_subplot(sub[i//2, i % 2]) for i in range(4)]
win = []
for row in sel: w = np.abs(t - row[1]) < 7/24; win.append((fd[w] - 1)*1e6)
lo = min(np.min(v) for v in win) - 1.1*sig*1e6; hi = max(np.max(v) for v in win) + 1.1*sig*1e6
for i, (ax, row) in enumerate(zip(axs, sel)):
    c0 = row[1]; w = np.abs(t - c0) < 0.3; h = (t[w] - c0)*24; yy = (fd[w] - 1)*1e6
    ax.errorbar(h, yy, sig*1e6, fmt='o', c='k', ms=1.8, lw=0.5, capsize=0)
    tt = np.linspace(-0.3, 0.3, 600); p.t0 = 0.; ax.plot(tt*24, (model(tt) - 1)*1e6, c='C1', lw=0.9); p.t0 = T0
    ax.set_xlim(-7, 7); ax.set_ylim(lo, hi); ax.axhline(0, c='0.5', lw=0.4, ls=':')
    ax.set_title('BKJD %.2f, S/N %.1f%s' % (c0, row[4], ' (shallow)' if row[0] == -3 else ''), fontsize=5.8, pad=2)
    if i % 2: ax.set_yticklabels([])
    else: ax.set_ylabel('flux [ppm]')
    if i < 2: ax.set_xticklabels([])
lab(axs[0], 'c')
fig.subplots_adjust(left=0.19, right=0.97, top=0.965, bottom=0.05)
x0 = axs[2].get_position().x0; x1 = axs[3].get_position().x1
fig.text((x0 + x1)/2, axs[-1].get_position().y0 - 0.032, 'hours from mid-transit', ha='center', va='top', fontsize=7)
fig.savefig(OUT, dpi=300); print('wrote', OUT, 'size in', fig.get_size_inches())
