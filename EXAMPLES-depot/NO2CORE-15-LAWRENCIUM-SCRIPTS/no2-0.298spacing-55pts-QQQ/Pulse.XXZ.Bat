#!/bin/bash

threads=1
procs=121
NPN=18

export OMP_NUM_THREADS=$threads

mpirun -n $procs -npernode $NPN --bind-to-core --bysocket --cpus-per-proc $threads --report-bindings chmctdhf_sinc Inp=Input.Inp.Pulse.xxz |tee Outs/Out.pulse.xxz

