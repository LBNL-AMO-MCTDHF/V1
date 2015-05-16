#!/bin/bash


 export OMP_NUM_THREADS=1
 mpirun  -n 256 -npernode 9 --bind-to-core --bysocket --cpus-per-proc 1 --report-bindings  chmctdhf_sinc Inp=Input.Inp.Neutral Debug=42 |tee Outs/Out.neutral



