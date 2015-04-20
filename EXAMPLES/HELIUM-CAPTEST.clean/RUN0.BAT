
threads=2
procs=21

export OMP_NUM_THREADS=$threads
mpirun  -n $procs -npernode 8 --bind-to-core --bysocket --cpus-per-proc $threads --report-bindings chmctdhf_sinc Inp=Input.Inp.States0 |tee Outs/Out.states0



