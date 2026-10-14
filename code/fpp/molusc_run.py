# Run MOLUSC CLI with its Gaia TAP queries routed to VizieR I/355/gaiadr3 (same DR3 astrometry as gaiadr3/gaiaedr3.gaia_source).
import sys, os, re, numpy as np, astropy.units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
import astroquery.gaia
class _Job:
    def __init__(s,t): s.t=t
    def get_results(s): return s.t
def _launch(q, **kw):
    ra,dec=map(float,re.search(r"POINT\('ICRS', ([-\d.]+), ([-\d.]+)\)",q).groups())
    c=SkyCoord(ra*u.deg,dec*u.deg)
    t=Vizier(columns=['**','_r'],row_limit=-1,timeout=300).query_region(c,radius=30*u.arcsec,catalog='I/355/gaiadr3')[0]
    t.sort('_r'); t=t[:10]
    o=Table(); o['dist']=np.asarray(t['_r'],float)/3600; o['ra']=t['RA_ICRS']; o['dec']=t['DE_ICRS']
    o['parallax']=t['Plx']; o['parallax_error']=t['e_Plx']; o['phot_g_mean_mag']=t['Gmag']; o['bp_rp']=t['BP-RP']
    o['astrometric_n_good_obs_al']=t['NgAL']; o['astrometric_chi2_al']=t['chi2AL']; o['ruwe']=t['RUWE']; o['astrometric_params_solved']=t['Solved']
    print('VizieR Gaia: nearest', t['Source'][0], 'RUWE', t['RUWE'][0], flush=True)
    return _Job(o)
astroquery.gaia.Gaia.launch_job=_launch; astroquery.gaia.Gaia.launch_job_async=_launch
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)),'molusc','code'))
sys.path.insert(0,'.')
import BinaryStarGUI as B
B.Gaia.launch_job=_launch
B.Application.star_jitter = -1  # upstream CLI bug: attribute never initialised (no stellar-jitter test requested)
np.random.seed(int(os.environ.get('SEED','1')))
B.Application(sys.argv).start()
