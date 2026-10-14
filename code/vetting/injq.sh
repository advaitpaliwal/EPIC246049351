#!/bin/sh
cd k2-deep/246049351
export OPENBLAS_NUM_THREADS=1
python -I -u scripts/inj3.py 246049351 everest 0.04707 21.196 4.33338 1 0 inj/real.json real > inj/real.log 2>&1
python -I -u scripts/inj3.py 246049351 everest 0.04707 21.196 4.33338 50 11 inj/inj.json inj > inj/inj.log 2>&1
python -I -u scripts/inj3.py 246049351 everest 0.04707 21.196 4.33338 10 12 inj/null.json null > inj/null.log 2>&1
echo D > inj/done
