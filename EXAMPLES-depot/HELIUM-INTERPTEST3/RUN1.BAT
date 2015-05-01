


export OMP_NUM_THREADS=9
mpirun  -n 7 -npernode 2 --bind-to-core --bysocket --cpus-per-proc 7 --report-bindings chmctdhf_sinc Inp=Input.Inp.States |tee Outs/Out.states


