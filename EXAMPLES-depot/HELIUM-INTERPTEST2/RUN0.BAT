

export OMP_NUM_THREADS=2
mpirun  -n 31 -npernode 8 --bind-to-core --bysocket --cpus-per-proc 2 --report-bindings chmctdhf_sinc Inp=Input.Inp.States0 |tee Outs/Out.states0



