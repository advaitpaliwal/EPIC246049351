"""Per-transit depths in EVEREST and K2SFF for one detection. Usage: EPIC CAMP P T0 DUR"""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from k2tls_hunt import load_everest, clean
from k2vet import load_sff
e, c, P, t0, dur = int(sys.argv[1]), int(sys.argv[2]), *map(float, sys.argv[3:6])
for lab, (t, f) in [("EV", clean(*load_everest(e, c))), ("SFF", clean(*load_sff(e, c)))]:
    n = np.round((t - t0) / P); sd = 1.4826 * np.median(np.abs(f - np.median(f)))
    for k in np.unique(n):
        dt = t - (t0 + k * P); i = np.abs(dt) < 0.4 * max(dur, 0.06); o = (np.abs(dt) > max(dur, .06)) & (np.abs(dt) < 4 * max(dur, .06))
        if i.sum() and o.sum() > 5:
            print(lab, int(k), round(t0 + k * P, 3), i.sum(), round((np.median(f[o]) - f[i].mean()) * 1e6), "+/-", round(sd / np.sqrt(i.sum()) * 1e6))
