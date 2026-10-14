"""TLS search of EVEREST K2 light curves. Usage: search.py targets.csv out.jsonl worker_idx nworkers"""
import sys, os, json, warnings, urllib.request
import numpy as np
from astropy.io import fits
warnings.filterwarnings("ignore")
from transitleastsquares import transitleastsquares
from wotan import flatten

LC = os.path.expanduser("k2-hunt/lc")

def ev_url(epic, c):
    e = str(epic)
    cs = "%02d" % int(c) if int(c) < 100 else str(c)
    cs = {"10": "102", "11": "112"}.get(cs, cs)  # C10/C11: second (longer) segment
    return f"https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HLSP/everest/v2/c{cs}/{e[:4]}00000/{e[4:]}/hlsp_everest_k2_llc_{e}-c{cs}_kepler_v2.0_lc.fits"

def fetch(url, path):
    if not os.path.exists(path):
        import subprocess
        tmp = path + ".part%d" % os.getpid()
        subprocess.run(["curl", "-sf", "--max-time", "120", "-o", tmp, url], check=True)
        os.rename(tmp, path)
    return path

def load_everest(epic, c):
    p = fetch(ev_url(epic, c), f"{LC}/ev_{epic}_c{c}.fits")
    d = fits.getdata(p, 1)
    t, f, q = d["TIME"], d["FLUX"], d["QUALITY"]
    m = np.isfinite(t) & np.isfinite(f) & ((q & ~((1 << 23) | (1 << 24))) == 0)  # keep EVEREST-outlier-flagged cadences (they include transits)
    # EVEREST also flags bad cadences in BADBIT-like mask; quality==0 keeps clean ones
    return t[m], f[m] / np.nanmedian(f[m])

def clean(t, f, window=0.4):
    fl = flatten(t, f, method="biweight", window_length=window)
    m = np.isfinite(fl)
    t, fl = t[m], fl[m]
    med, sd = np.median(fl), 1.4826 * np.median(np.abs(fl - np.median(fl)))
    m = (fl < med + 3 * sd) & (fl > med - 20 * sd)  # clip upward outliers and extreme low ones
    return t[m], fl[m]

def run(row):
    epic, c = int(row["EPIC"]), int(row["Camp"])
    t, f = load_everest(epic, c)
    if len(t) < 500:
        return dict(epic=epic, camp=c, status="short")
    t, f = clean(t, f)
    R = float(row["Rstar"]) if row["Rstar"] else 0.4
    M = float(row["Mstar"]) if row["Mstar"] else 0.4
    out = dict(epic=epic, camp=c, n=len(t), cdpp=float(np.std(f) * 1e6))
    res_all = []
    tt, ff = t.copy(), f.copy()
    for k in range(2):  # up to two signals per star
        model = transitleastsquares(tt, ff)
        r = model.power(period_min=0.5, period_max=min(30.0, (tt.max() - tt.min()) / 2.0),
                        R_star=R, M_star=M, R_star_min=0.08, R_star_max=1.0,
                        M_star_min=0.08, M_star_max=1.0, oversampling_factor=2,
                        duration_grid_step=1.15, use_threads=1, show_progress_bar=False)
        if not np.isfinite(r.SDE):
            break
        s = dict(period=float(r.period), t0=float(r.T0), depth=float(1 - r.depth),
                 duration=float(r.duration), sde=float(r.SDE), snr=float(r.snr),
                 ntr=int(r.distinct_transit_count), odd_even_mismatch=float(r.odd_even_mismatch),
                 rp_rs=float(r.rp_rs))
        res_all.append(s)
        if r.SDE < 9:
            break
        intr = np.abs(((tt - r.T0 + 0.5 * r.period) % r.period) - 0.5 * r.period) < r.duration
        tt, ff = tt[~intr], ff[~intr]
    out["signals"] = res_all
    return out

if __name__ == "__main__":
    import csv
    tf, of, wi, nw = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    rows = list(csv.DictReader(open(tf)))[wi::nw]
    done = set()
    import glob
    for g in glob.glob(os.path.join(os.path.dirname(of), "w*.jsonl")):
        for l in open(g):
            j = json.loads(l); done.add((j["epic"], j["camp"]))
    rows = [r for r in rows if (int(r["EPIC"]), int(r["Camp"])) not in done]
    from concurrent.futures import ThreadPoolExecutor
    pool = ThreadPoolExecutor(4)
    def pf(r):
        try: fetch(ev_url(r["EPIC"], r["Camp"]), f"{LC}/ev_{int(r['EPIC'])}_c{int(r['Camp'])}.fits")
        except Exception: pass
    futs = {}
    with open(of, "a") as fo:
        for i, row in enumerate(rows):
            for r2 in rows[i:i + 6]:
                k2 = (r2["EPIC"], r2["Camp"])
                if k2 not in futs: futs[k2] = pool.submit(pf, r2)
            futs[(row["EPIC"], row["Camp"])].result()
            key = (int(row["EPIC"]), int(row["Camp"]))
            if key in done:
                continue
            try:
                o = run(row)
            except Exception as e:
                o = dict(epic=key[0], camp=key[1], status="err", err=str(e)[:200])
            fo.write(json.dumps(o) + "\n"); fo.flush()
            # keep files only for strong signals
            sig = o.get("signals") or []
            if not any(s["sde"] >= 9 and s["snr"] >= 7 for s in sig):
                p = f"{LC}/ev_{key[0]}_c{key[1]}.fits"
                if os.path.exists(p):
                    os.remove(p)
