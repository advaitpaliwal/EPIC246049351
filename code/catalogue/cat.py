"""Catalogue queries for one EPIC. usage: cat.py EPIC RA DEC GAIA"""
import sys, os, urllib.parse, subprocess, json
epic, ra, dec, gaia = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
D = os.path.expanduser(f"k2-deep/{epic}/catalogue"); os.makedirs(D, exist_ok=True)
def get(name, url, data=None):
    p = os.path.join(D, name)
    open(p + ".url", "w").write(url + ("\nPOST " + data if data else "") + "\n")
    cmd = ["curl", "-sL", "--max-time", "120", "-o", p, url]
    if data: cmd[1:1] = ["--data", data]
    r = subprocess.run(cmd); print(name, r.returncode, os.path.getsize(p) if os.path.exists(p) else -1, flush=True)
q = urllib.parse.quote
nea = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?format=csv&query="
box = f"ra between {ra-2/60/abs(__import__('math').cos(dec*3.14159/180))} and {ra+2/60/abs(__import__('math').cos(dec*3.14159/180))} and dec between {dec-2/60} and {dec+2/60}"
get("nea_ps_box.csv", nea + q(f"select pl_name,hostname,ra,dec,pl_orbper,disc_facility,default_flag from ps where {box}"))
get("nea_k2pandc_box.csv", nea + q(f"select pl_name,epic_hostname,epic_candname,ra,dec,pl_orbper,disposition,disp_refname from k2pandc where {box}"))
get("nea_toi_box.csv", nea + q(f"select toi,tid,ra,dec,pl_orbper,tfopwg_disp from toi where {box}"))
get("nea_k2names.csv", nea + q(f"select * from k2names where epic_id like '%{epic}%'"))
sq = f"SELECT main_id,ra,dec,otype,DISTANCE(POINT('ICRS',ra,dec),POINT('ICRS',{ra},{dec}))*3600 AS sep FROM basic WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},{2/60}))=1"
get("simbad_2arcmin.csv", "https://simbad.cds.unistra.fr/simbad/sim-tap/sync?request=doQuery&lang=adql&format=csv&query=" + q(sq))
get("simbad_ids.txt", "https://simbad.cds.unistra.fr/simbad/sim-id?output.format=ASCII&Ident=" + q(f"EPIC {epic}"))
gq = f"SELECT source_id,ra,dec,phot_g_mean_mag,bp_rp,parallax,ruwe,non_single_star,ipd_frac_multi_peak,ipd_gof_harmonic_amplitude,astrometric_excess_noise,astrometric_excess_noise_sig,phot_variable_flag,radial_velocity,radial_velocity_error,rv_nb_transits,rv_chisq_pvalue,rv_renormalised_gof,DISTANCE(POINT('ICRS',ra,dec),POINT('ICRS',{ra},{dec}))*3600 AS sep FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{dec},{60/3600}))"
get("gaia_dr3_60as.csv", "https://gea.esac.esa.int/tap-server/tap/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY=" + q(gq))
get("gaia_nss.csv", "https://gea.esac.esa.int/tap-server/tap/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY=" + q(f"SELECT * FROM gaiadr3.nss_two_body_orbit WHERE source_id={gaia}"))
get("gaia_vari.csv", "https://gea.esac.esa.int/tap-server/tap/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&QUERY=" + q(f"SELECT * FROM gaiadr3.vari_summary WHERE source_id={gaia}"))
get("mast_tic.json", "https://mast.stsci.edu/api/v0/invoke", "request=" + q(json.dumps({"service":"Mast.Catalogs.Tic.Cone","params":{"ra":ra,"dec":dec,"radius":2/60},"format":"json","pagesize":200})))
get("exofop_k2.html", f"https://exofop.ipac.caltech.edu/k2/edit_target.php?id={epic}")
get("vizier_all_5as.tsv", f"https://vizier.cds.unistra.fr/viz-bin/votable?-c={ra}%20{dec:+f}&-c.rs=5&-out.max=20&-out.add=_r")
