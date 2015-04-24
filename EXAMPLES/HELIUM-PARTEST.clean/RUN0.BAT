
threads=5
procs=7

export OMP_NUM_THREADS=$threads
mpirun  -n $procs -npernode 4 --bind-to-core --bysocket --cpus-per-proc $threads --report-bindings chmctdhf_sinc.debug Inp=Input.Inp.States0 |tee Outs/Out.states0



