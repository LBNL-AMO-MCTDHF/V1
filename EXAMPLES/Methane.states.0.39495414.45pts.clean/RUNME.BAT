
export OMP_NUM_THREADS=1
mpirun  -n 15 -npernode 15 --bind-to-core --bysocket --cpus-per-proc 1 --report-bindings chmctdhf_sinc Inp=Input.Inp.States |tee Outs/Out.states

