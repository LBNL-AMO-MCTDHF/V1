
export OMP_NUM_THREADS=5
mpirun  -n 15 -npernode 2 --bind-to-core --bysocket --cpus-per-proc 5 --report-bindings chmctdhf_sinc Inp=Input.Inp.States |tee Outs/Out.states

