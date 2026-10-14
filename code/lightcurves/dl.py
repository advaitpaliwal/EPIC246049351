import sys, os, warnings; warnings.filterwarnings('ignore')
import lightkurve as lk
epic, camp = sys.argv[1], int(sys.argv[2])
D = os.path.expanduser(f"k2-deep/{epic}/lc")
sr = lk.search_lightcurve(f"EPIC {epic}", mission="K2")
print(sr)
open(f"{D}/search_lc.txt","w").write(str(sr)+"\nquery: lightkurve.search_lightcurve('EPIC %s', mission='K2') -> MAST\n"%epic)
st = lk.search_targetpixelfile(f"EPIC {epic}", mission="K2")
print(st); open(f"{D}/search_tpf.txt","w").write(str(st))
for i,r in enumerate(st):
    f = r.download(download_dir=D); print('tpf', f.path if hasattr(f,'path') else f)
for i,r in enumerate(sr):
    a = str(r.author[0]).lower()
    if a in ('k2sc','k2sff','everest','k2','k2varcat','polar'):
        try:
            f = r.download(download_dir=D); print('lc', a, getattr(f,'meta',{}).get('FILENAME'))
        except Exception as ex: print('fail', a, ex)
