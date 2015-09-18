#!/bin/bash

threads=1
procs=9

export OMP_NUM_THREADS=$threads

mpirun -n $procs chmctdhf_sinc.mpi Inp=Input.Inp.Pulse |tee Outs/Out.pulse

