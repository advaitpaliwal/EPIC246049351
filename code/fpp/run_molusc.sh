#!/bin/zsh
# MOLUSC (https://github.com/woodml/MOLUSC @ b26efbf) 2 x 100k companions per variant, M = 0.365799 Msun (TIC 49678675 mass), Gaia DR3 position (epoch 2016).
# Variant gr: --ruwe --gaia (joint analysis). Variant grs: + Gemini-S/Zorro 832 nm speckle curve (ExoFOP, 2021-09-23) as --ao with --filter I.
# Copied from ../../248780202/fpp_joint/scripts/run_molusc.sh.
F=k2-deep/246049351/fpp_joint
PYM=python
cd $F
for s in 21 22; do
  SEED=$s $PYM -I -W ignore molusc_run.py -v --ruwe --gaia -a $F/molusc_out/gr_s$s 23h25m59.5175s -- -08d16m14.151s 100000 0.365799 > molusc_out/gr_s$s.log 2>&1
  SEED=$s $PYM -I -W ignore molusc_run.py -v --ruwe --gaia --ao $F/cc832_molusc_mas.txt --filter I -a $F/molusc_out/grs_s$s 23h25m59.5175s -- -08d16m14.151s 100000 0.365799 > molusc_out/grs_s$s.log 2>&1
done
echo MOLUSC_DONE > molusc_out/done
