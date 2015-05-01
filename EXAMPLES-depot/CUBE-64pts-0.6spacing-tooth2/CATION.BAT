

export OMP_NUM_THREADS=8
mpirun  -n 8 -npernode 2 --bind-to-core --bysocket --cpus-per-proc 8 --report-bindings chmctdhf_sinc Inp=Input.Inp.Cation |tee Outs/Out.cation


