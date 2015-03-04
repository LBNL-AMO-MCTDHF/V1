import glob
import sys
from shutil import copyfile
from math import ceil, floor


###################
def templatereplace(oldfilename,newfilename,dictlist):
    oldfile=open(oldfilename)
    newfilestring=oldfile.read()
    
    for dict in dictlist:
        newfilestring=dictionaryreplace(newfilestring,dict)
        
    
    newfile=open(newfilename,'w+')
    newfile.write(newfilestring)
    #os.chmod(newfilename,stat. S_IXUSR)
    oldfile.close()
    newfile.close()

def maketaskfile(taskfilename, ncalcs):
    newfile=open(taskfilename,'w+')
    for i in range(ncalcs):
        newfile.write('bash runjob.sh '+str(i)+"\n")
    newfile.close()
def walltimestr(walltime):
    hrs=int(floor(walltime))
    mins=int(floor(60*(walltime-hrs)))
    return str(hrs)+":"+str(mins).zfill(2)+":00"

def jobfiledictionary(jobname,nnodes,corespernode,walltime):
    dict={"#jobname#":jobname, "#nnodes#":str(nnodes), "#corespernode#":str(corespernode), "#walltime#":walltimestr(walltime)}
    return dict

#def jobfiledictionary(nodenum,walltime,nmin,nmax):
#    
#    dict={"#nodeid#":str(nodenum), "#walltime#":walltimestr(walltime), "#nmin#":str(int(nmin)), "#nmax#":str(int(nmax))}
#    return dict

def dictionaryreplace(text,dictionary):
    newtext=text
    for i,j in dictionary.items():#.iteritems():
        newtext=newtext.replace(i,j)
    return newtext

def makejobfile(resultdir):
    if(resultdir[-1]!="/"):
        resultdir=resultdir+"/"
    tmpjobname=jobname+"."+resultdir[:-1]
    templatefilename="templates/htjobscript.slurm.template"

    masterscriptfile=resultdir+"queuesubmission.sh"
    ms=open(masterscriptfile,'w')

    dirlist=glob.glob(resultdir+"[0-9]*/")
    ncalcs=len(dirlist)
    print("ncalcs\t"+str(ncalcs))
    print("nnodes, nprocs\t"+str(nnodes)+"\t"+str(nnodes*corespernode))
    
    calcspernode=(ncalcs/nnodes)
    print("calcs per node\t"+str(calcspernode))
    print("calcs per processor\t"+str(calcspernode*1./corespernode))
    print("hours per run\t"+str(hrsperrun))
    print("safety factor\t"+str(safetyfactor))
    walltime=ncalcs*hrsperrun*safetyfactor/(1.*nnodes*corespernode)
    print("estimated walltime\t"+str(walltime))
    
    jobdict=jobfiledictionary(tmpjobname,nnodes,corespernode,walltime)
    newfilename=resultdir+"htjobscript.slurm"
    taskfilename=resultdir+"taskfile.txt"
    print("new filename\t"+newfilename)
    templatereplace(templatefilename,newfilename,[jobdict])
    maketaskfile(taskfilename,ncalcs)
    copyfile("templates/runjob.sh",resultdir+"runjob.sh")



################################
#Main Program

jobname="MCTDHF"
nnodes=8
corespernode=16
hrsperrun=4.#15/60.#1.
safetyfactor=2.#overestimate time for calculation to avoid running out of time

print("argv\t"+str(sys.argv))
#resultdir=sys.argv[-1]#"results_tmp/"
for resultdir in sys.argv[1:]:
    copyfile("./createHTjobfiles.slurm.py", resultdir+"createHTjobfiles.slurm.py")
    makejobfile(resultdir)

#for nodenum in range(nnodes):
#    nmin=nodenum*calcspernode
#    nmax=min((nodenum+1)*calcspernode,ncalcs)
#    walltime=(nmax-nmin)*(hrsperrun/corespernode)*safetyfactor
#    nodedict=jobfiledictionary(nodenum,walltime,nmin,nmax)
#    newfilename=resultdir+"multiprocessingjob."+str(nodenum)+".pbs"
#    print("new filename\t"+newfilename)
#    ms.write("qsub multiprocessingjob."+str(nodenum)+".pbs\n")
#    templatereplace(templatefilename,newfilename,[nodedict])
#ms.close()
    
