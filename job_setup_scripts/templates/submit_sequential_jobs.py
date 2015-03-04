from glob import glob
from subprocess import call
from os import getcwd,chdir

scriptlist=glob("*/mpijob.seq.slurm")
currentdir=getcwd()
for scriptname in scriptlist:
    dirname=scriptname.split("/")[0]+"/"
    #print("dirname\t"+dirname)
    chdir(dirname)
    print("submitting job in "+dirname)
    retval=call("sbatch mpijob.seq.slurm",shell=True)
    if(retval!=0):
        print("problem in "+dirname+"!\t"+str(retval))
    chdir(currentdir)
    #print("returned to "+getcwd())
