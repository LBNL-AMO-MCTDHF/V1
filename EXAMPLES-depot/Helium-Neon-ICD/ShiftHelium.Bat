#!/bin/bash

threads=8
procs=8
NPN=2

export OMP_NUM_THREADS=$threads
mpirun  -n $procs -npernode $NPN --bind-to-core --bysocket --cpus-per-proc $threads --report-bindings chmctdhf_sinc $1 Inp=Input.Inp.ShiftHelium |tee Outs/Out.shifthelium



