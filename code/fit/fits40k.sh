#!/bin/sh
# Final chains: 40 walkers x 40000 steps (first 1/3 discarded, thin 5) for EVEREST prior and free; K2SFF and fresh prior kept at 40000 too.
cd k2-deep/246049351
python -I -u scripts/fit3.py 246049351 $1 4.33316777 2906.7617 1.69 6.80341 0.259179 0.377432 0.0119343 3544 157 40000 $2 > fit/run_$2_$1.log 2>&1
