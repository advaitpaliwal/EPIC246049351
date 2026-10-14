#!/bin/sh
# init target (VizieR Gaia shim), then 10 runs gr (MOLUSC RUWE+Gaia), 10 runs grs (+speckle in MOLUSC and TRICERATOPS), 1 none, 1 cc reference. One worker.
D=k2-deep/246049351/fpp_joint
PY=python
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
GAIAPATCH=1 $PY -I -u $D/scripts/run_joint.py init 0 > $D/out/init.log 2>&1
for v in "gr 0" "gr 1" "gr 2" "gr 3" "gr 4" "gr 5" "gr 6" "gr 7" "gr 8" "gr 9" "none 0" "grs 0" "grs 1" "grs 2" "grs 3" "grs 4" "grs 5" "grs 6" "grs 7" "grs 8" "grs 9" "cc 0"; do
  set -- $v; $PY -I -u $D/scripts/run_joint.py $1 $2 > $D/out/run_$1_$2.log 2>&1
done
echo ALL_DONE > $D/out/done
