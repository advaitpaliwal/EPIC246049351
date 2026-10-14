# EPIC 246049351.01: a 1.9 Earth-radius planet candidate on a 4.3-day orbit around an M dwarf from K2 Campaign 12

Advait Paliwal (independent researcher, advaitpaliwal@gmail.com), 2026-10-13.

A dip in brightness of about 0.26% repeats every 4.33 days in K2 Campaign 12 (2016–2017) photometry of EPIC 246049351, an M dwarf 175 pc away. Seventeen transits fall in the campaign data. If it is a planet, it is about 1.9 times the size of Earth and receives about 15 times the sunlight Earth does.

| | |
|---|---|
| Host | EPIC 246049351 = TIC 49678675 = Gaia DR3 2438563103592977408 = 2MASS J23255947-0816135, M dwarf (TIC Teff 3544 ± 157 K, R 0.377 ± 0.012 R☉, M 0.366 M☉), Gaia DR3 parallax 5.727 mas, G = 15.96, Kp = 16.04, RUWE 1.04 |
| Period | 4.333377 (+0.000262/−0.000271) d |
| Mid-transit | 2457774.4261 ± 0.0014 BJD_TDB |
| Depth / duration | 2609 (+189/−176) ppm / 1.58 (+0.05/−0.09) h (EVEREST light curve) |
| Radius | 1.94 ± 0.11 R⊕ (for R★ = 0.377 R☉) |
| Equilibrium temperature / insolation | 544 ± 24 K / 14.6 (+2.8/−2.4) S⊕ |
| Data | K2 C12 (2016 December – 2017 March), 30-minute cadence; 17 transits |
| High-resolution imaging | Gemini-South Zorro speckle, 562 and 832 nm, 2021-09-23 (PI Howell; https://exofop.ipac.caltech.edu/tess/target.php?id=49678675): no companion; Δm 4.3 at 0.1″ and 5.3 at 0.5″ (832 nm) |
| False-positive probability (TRICERATOPS, no neighbours removed, 10 runs each) | 0.0098 ± 0.0019 with MOLUSC Gaia limits (range 0.0069–0.0124); 0.0071 ± 0.0013 with the 832 nm speckle curve added (range 0.0050–0.0088); 0.050 with neither; NFPP 0 |

Status: **planet candidate**. Not a TOI, CTOI or K2 candidate in the NASA Exoplanet Archive or ExoFOP, and not in SIMBAD, as of 2026-10-13. Odd and even transits agree (0.8σ), there is no secondary eclipse (−118 ± 131 ppm at phase 0.5), the seventeen transit depths are consistent (χ² = 15.1 for 16 degrees of freedom), and the fitted depths from EVEREST (2609 ppm), K2SFF (2671 ppm) and new aperture photometry (3061 ppm) agree within 0.7σ; no K2SC product exists. One transit (BKJD 2928.43) is shallow (1186 ± 585 ppm in EVEREST, absent in K2SFF) and follows a thruster firing 1.5 h before mid-transit; leaving it out changes the box depth from 2315 to 2405 ppm. The difference-image centroid is 1.5″ (1.3σ) from the out-of-transit centroid, which 63% of no-transit difference images exceed. TESS (Sectors 42, 70, 92) cannot test the transit at T = 14.85, and the single Campaign 19 transit is inconclusive.

## Contents

- `note/figure1_EPIC246049351.png` — (a) C12 EVEREST photometry folded on the fitted ephemeris with 30-minute bins and the median model; (b) the whole campaign before and after detrending; (c) four single transits, including the shallow one
- `data/lc/` — the detrended light curve with phase and model (`everest_detrended_c12.csv`), the undetrended EVEREST light curve, and box depths per reduction
- `data/fit/` — posterior summary (`posterior_summary_everest.csv`, 40 walkers × 40,000 steps) and fit outputs for EVEREST (with and without the density prior), K2SFF and new aperture photometry; the Campaign 19 check
- `data/vetting/`, `data/injection/` — odd/even, secondary-eclipse, phase-scan and per-transit checks; per-epoch depths, quality flags and thruster firings; injection-recovery and null searches
- `data/pixels/` — difference-image centroid and its null test, depth versus aperture, per-pixel fit, with the same tests on K2-325 b for reference
- `data/fpp/` — false-positive probability runs (per-run results, scenario tables, MOLUSC summaries, TRICERATOPS star list, the 832 nm contrast curve as used)
- `data/tess/`, `data/catalogue/`, `data/search/` — TESS coverage and sensitivity, catalogue queries with their URLs (including the speckle sensitivity files and the 2026-10-13 recheck), and the search detection
- `code/` — search, light-curve, fit, vetting, pixel, FPP, TESS, catalogue and figure scripts

Data: K2 (NASA, via MAST), EVEREST and K2SFF high-level science products, TESS (TESScut), Gemini-South Zorro speckle data (via ExoFOP), ESA Gaia DR3, TIC v8.2. Software: lightkurve, batman, emcee, transitleastsquares, wotan, TRICERATOPS, MOLUSC, tess-point, astropy, matplotlib. Licences: code MIT; text, figures and data CC BY 4.0.

Scripts are kept as run. File locations inside them are given relative to a single working directory (for example `k2-deep/246049351/`); set that directory up, or edit the paths, to rerun them.
