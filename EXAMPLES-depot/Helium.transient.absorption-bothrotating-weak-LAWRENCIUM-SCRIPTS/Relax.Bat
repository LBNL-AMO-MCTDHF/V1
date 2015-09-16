
export OMP_NUM_THREADS=18

mpirun -n 1 -npernode 18 --bind-to-core --cpus-per-proc 18 --report-bindings chmctdhf_atom Inp=Input.Inp.Relax |tee Outs/Out.relax

