
threads=1
procs=8

export OMP_NUM_THREADS=$threads
mpirun -n $procs chmctdhf_sinc.mpi Inp=Input.Inp.Relax |tee Outs/Out.relax

