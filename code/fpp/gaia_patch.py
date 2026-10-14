# Route TRICERATOPS' Gaia DR3 field-star query through VizieR (I/355/gaiadr3); ESA Gaia TAP was hanging.
import re, numpy as np, astropy.units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
import astroquery.gaia
class _Job:
    def __init__(s,t): s.t=t
    def get_results(s): return s.t
def _launch(adql, verbose=False, **kw):
    m=re.search(r"CIRCLE\('ICRS',\s*([-\d.e]+),\s*([-\d.e]+),\s*([-\d.e]+)\)\)\s*AND phot_g_mean_mag < ([\d.]+)",adql)
    ra,dec,r,gl=map(float,m.groups())
    v=Vizier(columns=['Gmag','BPmag','RPmag','Plx','e_Plx'],column_filters={'Gmag':f'<{gl}'},row_limit=-1,timeout=600)
    t=v.query_region(SkyCoord(ra*u.deg,dec*u.deg),radius=r*u.deg,catalog='I/355/gaiadr3')[0]
    print('VizieR Gaia field stars:',len(t),flush=True)
    o=Table(); o['phot_g_mean_mag']=t['Gmag']; o['phot_bp_mean_mag']=t['BPmag']; o['phot_rp_mean_mag']=t['RPmag']
    o['parallax']=t['Plx']; o['parallax_over_error']=np.asarray(t['Plx'],float)/np.asarray(t['e_Plx'],float)
    return _Job(o)
astroquery.gaia.Gaia.launch_job=_launch; astroquery.gaia.Gaia.launch_job_async=_launch
