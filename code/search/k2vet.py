"""Gate vetting for K2 TLS detections. Usage: k2vet.py EPIC CAMP PERIOD T0 DURATION [Rstar Mstar Kp]
Gates in order, stop at first failure:
 a) present in EVEREST and K2SFF with consistent depths, and in both campaign halves
 b) odd/even within 3 sigma, no secondary, planet shape, stellar density consistent
 c) no Gaia DR3 neighbour within 12" bright enough to produce the depth
"""
import sys, os, json, urllib.request, urllib.parse, shutil, warnings
import numpy as np
from astropy.io import fits
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(__file__))
from k2tls_hunt import load_everest, clean, fetch, LC

SRC = os.path.expanduser("k2-hunt/_sources")
G = 6.674e-8

def sff_url(epic, c):
    e = str(epic); cs = "%02d" % int(c)
    cs = {"10": "102", "11": "112"}.get(cs, cs)  # C10/C11: second (longer) segment
    return f"https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HLSP/k2sff/c{cs}/{e[:4]}00000/{e[4:]}/hlsp_k2sff_k2_lightcurve_{e}-c{cs}_kepler_v1_llc.fits"

def load_sff(epic, c):
    p = fetch(sff_url(epic, c), f"{LC}/sff_{epic}_c{c}.fits")
    d = fits.getdata(p, 1)  # BESTAPER
    t, f = d["T"], d["FCOR"]
    m = np.isfinite(t) & np.isfinite(f) & (d["MOVING"] == 0)
    return t[m], f[m] / np.nanmedian(f[m])

def phase(t, P, t0):
    return ((t - t0 + 0.5 * P) % P) - 0.5 * P

def depth(t, f, P, t0, dur, sel=None):
    ph = phase(t, P, t0)
    de = max(dur, 0.06)
    if sel is None: sel = np.ones(len(t), bool)
    inn = sel & (np.abs(ph) < 0.4 * de)
    out = sel & (np.abs(ph) > 1.0 * de) & (np.abs(ph) < 4.0 * de)
    allout = np.abs(ph) > 1.0 * de
    if inn.sum() < 3 or out.sum() < 10: return np.nan, np.nan
    sd = 1.4826 * np.median(np.abs(f[allout] - np.median(f[allout])))
    return float(np.median(f[out]) - np.mean(f[inn])), float(sd * np.sqrt(1.0 / inn.sum() + 1.0 / out.sum()))

def cons(d1, e1, d2, e2):
    return abs(d1 - d2) / np.hypot(e1, e2)

def trapfit(t, f, P, t0, dur):
    from scipy.optimize import curve_fit
    ph = phase(t, P, t0); m = np.abs(ph) < 2.5 * dur
    x, y = ph[m], f[m]
    def trap(x, dep, T14, frac, x0):
        T12 = frac * T14 / 2
        a = np.abs(x - x0); out = np.ones_like(a)
        full = a <= T14 / 2 - T12
        ing = (a > T14 / 2 - T12) & (a < T14 / 2)
        out[full] = 1 - dep
        out[ing] = 1 - dep * (T14 / 2 - a[ing]) / max(T12, 1e-6)
        return out
    best = None
    for fr in (0.1, 0.4, 0.8):
        try:
            p, _ = curve_fit(trap, x, y, p0=[max(1 - np.min(y), 1e-4) / 2, dur, fr, 0],
                             bounds=([0, 0.3 * dur, 0.01, -dur / 2], [0.5, 3 * dur, 1.0, dur / 2]))
            chi = np.sum((y - trap(x, *p)) ** 2)
            if best is None or chi < best[1]: best = (p, chi)
        except Exception:
            pass
    return best[0] if best else None

def gaia_neighbours(ra, dec, epic):
    q = (f"SELECT source_id,ra,dec,phot_g_mean_mag,bp_rp,parallax,pmra,pmdec,"
         f"DISTANCE(POINT('ICRS',ra,dec),POINT('ICRS',{ra},{dec}))*3600 AS sep "
         f"FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},{30/3600}))")
    url = "https://gea.esac.esa.int/tap-server/tap/sync?" + urllib.parse.urlencode(
        dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=q))
    out = f"{SRC}/gaia_{epic}.csv"
    if not os.path.exists(out):
        with urllib.request.urlopen(url, timeout=120) as r, open(out, "wb") as fo:
            shutil.copyfileobj(r, fo)
        open(out + ".url", "w").write(url + "\n")
    import csv
    return list(csv.DictReader(open(out)))

def vet(epic, c, P, t0, dur, R=None, M=None, ra=None, dec=None, gmag=None, gaia_id=None):
    r = dict(epic=epic, camp=c, P=P, t0=t0, dur=dur)
    te, fe = clean(*load_everest(epic, c))
    de, ee = depth(te, fe, P, t0, dur)
    r.update(depth_ev=de, err_ev=ee, snr_ev=de / ee)
    # gate a
    try:
        ts, fs = clean(*load_sff(epic, c))
        ds, es = depth(ts, fs, P, t0, dur)
    except Exception as ex:
        r.update(gate="a", fail=f"K2SFF unavailable: {ex}"[:120]); return r
    r.update(depth_sff=ds, err_sff=es)
    if not (ds / es > 3):
        r.update(gate="a", fail=f"not in K2SFF ({ds/es:.1f} sigma)"); return r
    if cons(de, ee, ds, es) > 3 and abs(ds / de - 1) > 0.3:
        r.update(gate="a", fail=f"EVEREST/K2SFF depth mismatch {cons(de,ee,ds,es):.1f} sigma"); return r
    mid = 0.5 * (te.min() + te.max())
    d1, e1 = depth(te, fe, P, t0, dur, te < mid); d2, e2 = depth(te, fe, P, t0, dur, te >= mid)
    r.update(depth_h1=d1, depth_h2=d2, err_h1=e1, err_h2=e2)
    if not (d1 / e1 > 2.5 and d2 / e2 > 2.5) or (cons(d1, e1, d2, e2) > 3 and abs(d1 / d2 - 1) > 0.3):
        r.update(gate="a", fail=f"halves h1 {d1/e1:.1f}s h2 {d2/e2:.1f}s diff {cons(d1,e1,d2,e2):.1f}s"); return r
    # gate a (cont.): epoch-to-epoch consistency, errors from empirical single-window scatter off transit
    w = max(dur, 0.06)
    def ep_depth(tc):
        dt = te - tc; i = np.abs(dt) < 0.4 * w; o = (np.abs(dt) > w) & (np.abs(dt) < 4 * w)
        return (np.median(fe[o]) - fe[i].mean()) if (i.sum() >= 2 and o.sum() > 5) else np.nan
    rng = np.random.default_rng(1)
    ph_all = phase(te, P, t0)
    rnd = [ep_depth(x) for x in rng.uniform(te.min(), te.max(), 300) if abs(phase(np.array([x]), P, t0)[0]) > 2 * w]
    sig1 = 1.4826 * np.nanmedian(np.abs(np.array(rnd) - np.nanmedian(rnd)))
    eps = np.arange(np.ceil((te.min() - t0) / P), np.floor((te.max() - t0) / P) + 1)
    dd = np.array([ep_depth(t0 + k * P) for k in eps]); dd = dd[np.isfinite(dd)]
    if len(dd) < 3:
        r.update(gate="a", fail=f"fewer than 3 measurable transit epochs ({len(dd)})"); return r
    if True:
        chi = np.sum((dd - np.mean(dd)) ** 2) / sig1 ** 2 / (len(dd) - 1)
        nlow = int(np.sum(dd < np.mean(dd) - 3 * sig1))
        r.update(epoch_depths=[round(x * 1e6) for x in dd], epoch_sig=sig1 * 1e6, epoch_redchi=chi)
        nweak = int(np.sum((dd < 0.25 * np.mean(dd)) & ((np.mean(dd) - dd) / sig1 > 3)))
        r.update(epoch_nweak=nweak)
        if chi > 6 or nweak > 0:
            r.update(gate="a", fail=f"epoch depths inconsistent (red. chi2 {chi:.1f}, {nlow} epochs >3 sigma shallow)"); return r
    # gate b
    n = np.round((te - t0) / P)
    do, eo = depth(te, fe, P, t0, dur, n % 2 == 1); dv, ev = depth(te, fe, P, t0, dur, n % 2 == 0)
    r.update(depth_odd=do, depth_even=dv, oe_sigma=cons(do, eo, dv, ev))
    if cons(do, eo, dv, ev) > 3:  # strict 3 sigma
        r.update(gate="b", fail=f"odd/even {cons(do,eo,dv,ev):.1f} sigma"); return r
    sec = []
    for ph0 in np.linspace(0.3, 0.7, 41):
        dsec, esec = depth(te, fe, P, t0 + ph0 * P, dur)
        if np.isfinite(dsec): sec.append((dsec / esec, ph0, dsec))
    smax = max(sec) if sec else (0, 0, 0)
    r.update(sec_sigma=smax[0], sec_phase=smax[1], sec_depth=smax[2])
    if smax[0] > 4:
        r.update(gate="b", fail=f"secondary {smax[0]:.1f} sigma at phase {smax[1]:.2f}"); return r
    p = trapfit(te, fe, P, t0, dur)
    if p is None:
        r.update(gate="b", fail="shape fit failed"); return r
    dep, T14, frac, x0 = p
    r.update(fit_depth=dep, T14=T14, ingress_frac=frac)
    k = np.sqrt(dep)
    if frac > 0.85 and k > 0.2:
        r.update(gate="b", fail=f"V-shaped (T12/T14 frac {frac:.2f}) with k={k:.2f}"); return r
    if T14 * 24 < 0.6:
        r.update(gate="b", fail=f"duration {T14*24:.2f} h shorter than ~1 long cadence (not resolved)"); return r
    if frac > 0.85 and T14 * 24 < 1.0:
        r.update(gate="b", fail=f"V-shaped short event (frac {frac:.2f}, {T14*24:.2f} h)"); return r
    if dep > 0.3:
        r.update(gate="b", fail=f"depth {dep:.2f} too deep for planet"); return r
    # density from duration assuming b=0, circular: upper bound
    aR = P / (np.pi * T14) * (1 + k)
    rho_max = 3 * np.pi * aR ** 3 / (G * (P * 86400) ** 2)
    if R and M:
        rho_star = 1.41 * M / R ** 3
        r.update(rho_transit_max=rho_max, rho_star=rho_star)
        if rho_max > 30 * rho_star:
            r.update(gate="b", fail=f"duration too short: rho_max {rho_max:.0f} vs star {rho_star:.1f} g/cc (needs b~1 grazing)"); return r
        if rho_max < 0.25 * rho_star:
            r.update(gate="b", fail=f"duration too long: rho_max {rho_max:.2f} vs star {rho_star:.2f} g/cc"); return r
    # gate c
    if ra is not None:
        nb = gaia_neighbours(ra, dec, epic)
        tgt = None
        for s in nb:
            if gaia_id and str(s["source_id"]) == str(gaia_id): tgt = s
        if tgt is None and nb:
            tgt = min(nb, key=lambda s: float(s["sep"]))
        gt = float(tgt["phot_g_mean_mag"])
        bad = []
        for s in nb:
            if s is tgt or not s["phot_g_mean_mag"]: continue
            sep = np.hypot((float(s["ra"]) - float(tgt["ra"])) * np.cos(np.radians(dec)), float(s["dec"]) - float(tgt["dec"])) * 3600
            if sep > 12: continue
            fr = 10 ** (-0.4 * (float(s["phot_g_mean_mag"]) - gt))
            maxdep = fr / (1 + fr)
            if maxdep >= dep:
                bad.append(dict(id=s["source_id"], sep=round(sep, 2), G=float(s["phot_g_mean_mag"]), maxdepth=round(maxdep, 4)))
        r.update(neighbours_12=bad, target_gaia=tgt["source_id"])
        if bad:
            r.update(gate="c", fail=f"{len(bad)} Gaia neighbour(s) within 12\" can produce depth"); return r
    r.update(gate="pass")
    return r

if __name__ == "__main__":
    a = sys.argv[1:]
    print(json.dumps(vet(int(a[0]), int(a[1]), float(a[2]), float(a[3]), float(a[4]))))
