# Adapted from ../../248780202/fpp_joint/scripts/aggregate.py (single planet). Per variant: FPP = 1 - sum P(TP, PTP, DTP); NFPP = sum P(NTP, NEB, NEBx2P),
# recomputed from each probs csv and asserted equal to TRICERATOPS' own numbers.
import json, hashlib, numpy as np, pandas as pd
from pathlib import Path
D = Path('k2-deep/246049351/fpp_joint')
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
out = {'input_sha256': {f: sha(D / f) for f in ['out/target_246049351.pkl', 'molusc_out/molusc_gr_for_trc_N1e6.csv', 'molusc_out/molusc_grs_for_trc_N1e6.csv', 'cc832_arcsec_dmag.csv']}}
for v, nr in [('gr', 10), ('grs', 10), ('none', 1), ('cc', 1)]:
    rows, sc = [], []
    for k in range(nr):
        f = D / f'out/res_{v}_{k}.json'
        if not f.exists(): continue
        r = json.loads(f.read_text()); p = pd.read_csv(D / f'out/probs_{v}_{k}.csv')
        assert np.isclose(p.prob.sum(), 1, atol=1e-9) and (p.prob >= 0).all()
        fpp = 1 - p.loc[p.scenario.isin(['TP', 'PTP', 'DTP']), 'prob'].sum(); nfpp = p.loc[p.scenario.isin(['NTP', 'NEB', 'NEBx2P']), 'prob'].sum()
        assert np.isclose(fpp, r['FPP'], atol=1e-9) and np.isclose(nfpp, r['NFPP'], atol=1e-9)
        rows.append({x: r[x] for x in r if x != 'aperture'}); sc.append(p.groupby('scenario').prob.sum().rename(k))
    if not rows: continue
    d = pd.DataFrame(rows); s = pd.concat(sc, axis=1).fillna(0)
    st = pd.DataFrame({'mean': s.mean(axis=1), 'sd': s.std(axis=1, ddof=1), 'min': s.min(axis=1), 'max': s.max(axis=1)})
    st.to_csv(D / f'scenario_summary_{v}.csv'); d.to_csv(D / f'run_results_{v}.csv', index=False)
    out[v] = {'n_runs': len(d), 'seeds': d.seed.tolist(), 'err_ppm': float(d.err_ppm.iloc[0]), 'nbins': int(d.nbins.iloc[0]), 'depth_ppm': float(d.depth_ppm.iloc[0]), 'n_stars': int(d.n_stars.iloc[0]),
              'FPP': {'mean': d.FPP.mean(), 'sd': d.FPP.std(ddof=1) if len(d) > 1 else None, 'min': d.FPP.min(), 'max': d.FPP.max(), 'runs_le_0.035': int((d.FPP <= 0.035).sum()), 'all': d.FPP.tolist()},
              'NFPP': {'mean': d.NFPP.mean(), 'max': d.NFPP.max()}, 'scenarios': {i: {c: float(st.loc[i, c]) for c in st.columns} for i in st.index if st.loc[i, 'mean'] > 1e-6}}
(D / 'summary.json').write_text(json.dumps(out, indent=1, default=float) + '\n'); print(json.dumps(out, indent=1, default=float))
