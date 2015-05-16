#!/bin/bash

threads=8
procs=16

export OMP_NUM_THREADS=$threads
mpirun  -n $procs -npernode 2 --bind-to-core --bysocket --cpus-per-proc $threads --report-bindings chmctdhf_sinc.debug $1 Inp=Input.Inp.States1 |tee Outs/Out.states1



