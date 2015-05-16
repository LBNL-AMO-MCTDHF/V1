#!/bin/bash

threads=1
procs=64
NPN=16

export OMP_NUM_THREADS=$threads
mpirun  -n $procs -npernode $NPN --bind-to-core --bysocket --cpus-per-proc $threads --report-bindings chmctdhf_sinc $1 Inp=Input.Inp.Butadiene |tee Outs/Out.butadiene



