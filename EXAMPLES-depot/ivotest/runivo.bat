
procs=16
threads=1
NPN=16

 export OMP_NUM_THREADS=$threads
 mpirun  -n $procs -npernode $NPN --bind-to-core --bysocket --cpus-per-proc $threads --report-bindings  chmctdhf_sinc.debug Inp=Input.Inp.ivo |tee Outs/Out.ivo






