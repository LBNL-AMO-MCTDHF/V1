#!/bin/bash

source ./RUNS.H

export OMP_NUM_THREADS=$threads
mpirun  -n $procs -npernode $NPN --bind-to-core --bysocket --cpus-per-proc $threads --report-bindings chmctdhf_diatom.debug $OPTION Inp=Input.Inp.States0$1 |tee Outs/Out.states0$1



