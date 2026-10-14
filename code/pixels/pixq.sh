#!/bin/sh
# waits for the vet chain, then: reference download + ephemeris, difference-image null test (pix2.py) on target and K2-325 b, depth vs aperture, PLD per-pixel map.
cd k2-deep/246049351
while [ ! -f fit/vet_done ]; do sleep 5; done
TPF=$(ls lc/mastDownload/K2/ktwo246049351-c12_lc/*lpd-targ.fits.gz)
python -I -u scripts/ref_k2325.py > pix/ref/ref.log 2>&1
python -I -u scripts/pix2.py $TPF 4.3333771136882735 2941.4260584848926 1.5751261977658904 351.4979894572662 -8.270597422518295 pix/pix_target.json > pix/pix_target.log 2>&1
RT=$(ls pix/ref/mastDownload/K2/*/*lpd-targ.fits.gz); RP=$(python -c "import json;r=json.load(open('pix/ref/ref_ephem.json'));print(r['P'],r['T0_BKJD'],r['T14_h'])")
RC=$(python -c "import json;d=json.load(open('catalogue/nea_c12_confirmed.csv.url')) " 2>/dev/null; echo "353.9169507 -7.6867213")
python -I -u scripts/pix2.py $RT $RP $RC pix/pix_ref_K2-325b.json > pix/pix_ref.log 2>&1
python -I -u ../scripts/apdep.py $TPF 4.3333771136882735 2941.4260584848926 1.5751261977658904 pix/apdep.json > pix/apdep.log 2>&1
python -I -u ../scripts/pixpld.py $TPF 4.3333771136882735 2941.4260584848926 1.5751261977658904 0.04707 21.196 0.277 351.4979894572662 -8.270597422518295 pix/pld_target > pix/pld.log 2>&1
python -I -u scripts/tess3.py pix2x2 > tess/tess_pix2x2.log 2>&1
echo D > pix/done
