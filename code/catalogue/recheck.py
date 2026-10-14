# Fresh catalogue recheck for EPIC 246049351 / TIC 49678675. usage: python recheck.py OUTDIR
import os, sys, json, urllib.parse, subprocess, math, datetime
D = sys.argv[1]
ra, dec = 351.497778491632, -8.27044515731495; q = urllib.parse.quote
def get(name, url):
    p = os.path.join(D, name); open(p + '.url', 'w').write(url + '\nfetched ' + datetime.datetime.now(datetime.UTC).isoformat() + '\n')
    r = subprocess.run(['curl', '-sL', '--max-time', '300', '-o', p, url]); print(name, r.returncode, os.path.getsize(p) if os.path.exists(p) else -1, flush=True)
c = abs(math.cos(math.radians(dec)))
box = f"ra between {ra-2/60/c} and {ra+2/60/c} and dec between {dec-2/60} and {dec+2/60}"
nea = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?format=csv&query="
get('nea_ps_box.csv', nea + q(f"select pl_name,hostname,ra,dec,pl_orbper,disc_facility,default_flag from ps where {box}"))
get('nea_k2pandc_box.csv', nea + q(f"select pl_name,epic_hostname,epic_candname,ra,dec,pl_orbper,disposition,disp_refname from k2pandc where {box}"))
get('nea_k2pandc_epic.csv', nea + q("select pl_name,epic_hostname,epic_candname,pl_orbper,disposition,disp_refname from k2pandc where epic_hostname like '%246049351%'"))
get('nea_toi_box.csv', nea + q(f"select toi,tid,ra,dec,pl_orbper,tfopwg_disp from toi where {box}"))
get('nea_k2names.csv', nea + q("select * from k2names where epic_id like '%246049351%'"))
sq = f"SELECT main_id,ra,dec,otype,DISTANCE(POINT('ICRS',ra,dec),POINT('ICRS',{ra},{dec}))*3600 AS sep FROM basic WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},{2/60}))=1"
get('simbad_2arcmin.csv', "https://simbad.cds.unistra.fr/simbad/sim-tap/sync?request=doQuery&lang=adql&format=csv&query=" + q(sq))
for nm, i in [('epic', 'EPIC 246049351'), ('tic', 'TIC 49678675'), ('gaia', 'Gaia DR3 2438563103592977408'), ('2mass', '2MASS J23255947-0816135')]:
    get(f'simbad_id_{nm}.txt', "https://simbad.cds.unistra.fr/simbad/sim-id?output.format=ASCII&Ident=" + q(i))
get('simbad_coo_10as.txt', f"https://simbad.cds.unistra.fr/simbad/sim-coo?output.format=ASCII&Coord={ra}%20{dec}&Radius=10&Radius.unit=arcsec")
get('exofop_toi.csv', 'https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv')
get('exofop_ctoi.csv', 'https://exofop.ipac.caltech.edu/tess/download_ctoi.php?sort=ctoi&output=csv')
get('exofop_tic49678675.json', 'https://exofop.ipac.caltech.edu/tess/target.php?id=49678675&json')
get('eumemic_tree.json', 'https://api.github.com/repos/eumemic/exoplanets/git/trees/HEAD?recursive=1')
get('eumemic_commits.json', 'https://api.github.com/repos/eumemic/exoplanets/commits?per_page=5')
sha = json.load(open(os.path.join(D, 'eumemic_commits.json')))[0]['sha']
tree = json.load(open(os.path.join(D, 'eumemic_tree.json')))['tree']
files = [t['path'] for t in tree if t['type'] == 'blob' and t['path'].startswith('results_public/')]
E = os.path.join(D, f'eumemic_{sha[:7]}'); os.makedirs(E, exist_ok=True)
open(E + '.url', 'w').write(f'https://github.com/eumemic/exoplanets/tree/{sha}/results_public fetched {datetime.datetime.now(datetime.UTC).isoformat()}\n')
open(os.path.join(D, 'eumemic_files.txt'), 'w').write('\n'.join(t['path'] for t in tree) + '\n')
for f in files:
    subprocess.run(['curl', '-sL', '--max-time', '300', '-o', os.path.join(E, os.path.basename(f)), f'https://raw.githubusercontent.com/eumemic/exoplanets/{sha}/{f}'])
print('eumemic', sha, len(files), 'files')
