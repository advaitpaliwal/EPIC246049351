# Search every eumemic results_public file for the target's identifiers. usage: python -I eumemic_check.py DIR FILELIST
import sys, os, gzip
D, L = sys.argv[1], sys.argv[2]
ids = ['49678675', '246049351', '2438563103592977408']
tot = 0
for f in sorted(os.listdir(D)):
    p = os.path.join(D, f); op = gzip.open if f.endswith('.gz') else open
    with op(p, 'rt', errors='replace') as h: txt = h.read()
    hits = {i: txt.count(i) for i in ids}; tot += sum(hits.values()); print(f, os.path.getsize(p), len(txt.splitlines()), hits)
names = open(L).read().splitlines()
print('files', len(os.listdir(D)), 'total hits', tot, 'paths containing an ID', [n for n in names if any(i in n for i in ids)])
print('candidate folders', len({n.split('/')[1] for n in names if n.startswith('candidates/') and n.count('/') >= 1}))
