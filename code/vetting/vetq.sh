#!/bin/sh
cd k2-deep/246049351
for s in everest k2sff fresh; do python -I -u ../scripts/vet.py 246049351 $s 4.3333771136882735 2941.4260584848926 1.5751261977658904 > fit/vet_$s.log 2>&1; done
python -I -u scripts/epoch_k2sff.py > fit/epoch_k2sff.log 2>&1
python -I -u scripts/thruster3.py > fit/thruster_alias.log 2>&1
for s in everest k2sff fresh; do python -I -u scripts/pertr_null3.py $s > fit/pertr_null_$s.log 2>&1; done
echo D > fit/vet_done
