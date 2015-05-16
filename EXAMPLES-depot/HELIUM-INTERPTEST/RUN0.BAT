#!/bin/bash


export OMP_NUM_THREADS=3
mpirun  -n 5 -npernode 5 --bind-to-core --bysocket --cpus-per-proc 3 --report-bindings chmctdhf_sinc Inp=Input.Inp.States0 |tee Outs/Out.states0



