#!/bin/bash

threads=1
procs=49
NPN=17

export OMP_NUM_THREADS=$threads

mpirun -n $procs -npernode $NPN --bind-to-core --bysocket --cpus-per-proc $threads --report-bindings chmctdhf_sinc Inp=Input.Inp.Pulse.x.df |tee Outs/Out.pulse.x.df

