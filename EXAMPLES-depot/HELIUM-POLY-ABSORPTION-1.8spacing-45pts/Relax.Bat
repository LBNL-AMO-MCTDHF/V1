#!/bin/bash

threads=1
procs=25
NPN=13

export OMP_NUM_THREADS=$threads

mpirun -n $procs -npernode $NPN --bind-to-core --bysocket --cpus-per-proc $threads --report-bindings chmctdhf_sinc Inp=Input.Inp.Relax |tee Outs/Out.relax

