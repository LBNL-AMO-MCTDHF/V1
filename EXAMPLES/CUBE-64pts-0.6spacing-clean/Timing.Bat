

proc=16
ppn=4
threads=4

ext=${proc}proc.${ppn}ppn
export OMP_NUM_THREADS=$threads
rm -r timing.relax0 timing.relax0.$ext
mpirun -n $proc -npernode $ppn chmctdhf_sinc Inp=Input.Inp.Relax0 MPIOUT Debug=42 Step=0.05 |tee Outs/Out.relax0.timing.$ext
mv timing.relax0 timing.relax0.$ext

proc=4
ppn=1
threads=16

ext=${proc}proc.${ppn}ppn
export OMP_NUM_THREADS=$threads
rm -r timing.relax0 timing.relax0.$ext
mpirun -n $proc -npernode $ppn chmctdhf_sinc Inp=Input.Inp.Relax0 MPIOUT Debug=42 Step=0.05 |tee Outs/Out.relax0.timing.$ext
mv timing.relax0 timing.relax0.$ext


exit

proc=64
ppn=16
threads=1

ext=${proc}proc.${ppn}ppn
export OMP_NUM_THREADS=$threads
rm -r timing.relax0 timing.relax0.$ext
mpirun -n $proc -npernode $ppn chmctdhf_sinc Inp=Input.Inp.Relax0 MPIOUT Debug=42 Step=0.05 |tee Outs/Out.relax0.timing.$ext
mv timing.relax0 timing.relax0.$ext

exit

