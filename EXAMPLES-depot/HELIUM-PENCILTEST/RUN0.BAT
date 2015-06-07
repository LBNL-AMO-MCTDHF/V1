
threads=5
procs=9

export OMP_NUM_THREADS=$threads
mpirun  -n $procs -npernode 4 --bind-to-core --bysocket --cpus-per-proc $threads --report-bindings chmctdhf_sinc.debug $1 Inp=Input.Inp.States0 MPIOUT |tee Outs/Out.states0



