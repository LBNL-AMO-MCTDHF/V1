
threads=8
procs=16

export OMP_NUM_THREADS=$threads
mpirun  -n $procs -npernode 2 --bind-to-core --bysocket --cpus-per-proc $threads --report-bindings chmctdhf_sinc $1 Inp=Input.Inp.States2 |tee Outs/Out.states2



