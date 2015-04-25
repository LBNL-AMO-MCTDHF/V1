from numpy import *
import sys
import shutil
import subprocess
import os
from itertools import product as iterproduct
aut=24.2#attoseconds per atomic unit of time
Hrt=27.21

#this script is designed to write input files for Dan Haxton's MCTDHF program
#corresponding to cosine squared pulses.  An abstract pulse class is defined
#which is capable of calculating the various quantities needed to create the
#input file


################################################################################
#pulse class contains all information needed to characterize a pulse and create
#the relevant input files

class pulse:
    def __init__(self, I0=1.e-4, pulsestrength=None, w0=.0565, wenv=.014,
                 CEPtheta=0, jtheta=0, jphiz=0, propeuler0=0, propeuler1=-pi/2,
                 propeuler2=-pi/2, tcenter=0, tstart=None , pulsetype=2,
                 chirp=0, ramp=0, longstep=0, rightcircular=False,
                 leftcircular=False):
        self.I0=I0#intensity divided by 10^16 W/cm^2
        if(pulsestrength==None):
            self.pulsestrength=sqrt(I0/3.51)/w0#definition in MCTDHF code
            #(hard coded 3.51 not explained in MCTDHF code  )
        else:
            self.pulsestrength=pulsestrength
            self.I0=pow(w0*pulsestrength,2)*3.51
        self.w0=w0#central frequency in atomic units
        self.wenv=wenv#envelope frequency in atomic units
        self.CEPtheta=CEPtheta#phase at peak of envelope
        self.jtheta=jtheta# jones theta: vec{E} = E0(zhat cos(jtheta) exp(i jphiz)+
                                #xhat sin(jtheta))
        self.jphiz=jphiz#jones phiz: vec{E} = E0(zhat cos(jtheta) exp(i jphiz)+
                                # xhat sin(jtheta))
        #given euler angles, define vectors for propagation direction and 2
        #polarization directions s.t. vpol1 x vpol2 = vprop                        

        #default values for right- or left-circularly polarized pulses
        if(rightcircular):
            self.jtheta=pi/4
            self.jphiz=-pi/2
        if(leftcircular):
            self.jtheta=pi/4
            self.jphiz=pi/2
        self.propeuler0=propeuler0
        self.propeuler1=propeuler1
        self.propeuler2=propeuler2
        self.initpolvectors()
        self.tcenter=tcenter
        if(tstart!=None):
            self.settstart(tstart)
        self.pulsetype=pulsetype
        self.chirp=chirp
        self.ramp=ramp
        self.longstep=longstep
        self.initMCTDHF()

    def initpolvectors(self):
        vprop=[0,0,1]
        vpol1=[1,0,0]
        vpol2=[0,1,0]
        rotmat0=Rx(self.propeuler0)
        rotmat1=Rz(self.propeuler1)
        rotmat2=Rx(self.propeuler2)
        self.vprop=dot(rotmat2,dot(rotmat1,vprop))
        self.vpol1=dot(rotmat2,dot(rotmat1,vpol1))
        self.vpol2=dot(rotmat2,dot(rotmat1,vpol2))
#        print("self.vprop\t"+str(self.vprop))
#        print("self.vpol1\t"+str(self.vpol1))
#        print("self.vpol2\t"+str(self.vpol2))

    def initMCTDHF(self):
        #set up two MCTDHF pulses containing the same information
        E0=self.E0()
        E1=E0*cos(self.jtheta)
        I1=E1**2
        E2=E0*sin(self.jtheta)
        I2=E2**2
        pulsestrength1=self.pulsestrength*cos(self.jtheta)
        pulsestrength2=self.pulsestrength*sin(self.jtheta)
        
        #MCTDHF takes I1, I2 as inputs rather than E1, E2.  Negative signs for
        #E1, E2 #correspond to a phase shift of pi
        if(sign(E1)<0):
            signphase1=pi
        else:
            signphase1=0
        if(sign(E2)<0):
            signphase2=pi
        else:
            signphase2=0
        pulser1, pulsetheta1, pulsephi1=topolar(self.vpol1)
        pulser2, pulsetheta2, pulsephi2=topolar(self.vpol2)
        
        strtphase=self.startphase()
        strtphase1=mod(strtphase+self.jphiz+signphase1, 2*pi)
        strtphase2=mod(strtphase+signphase2, 2*pi)
#apologies for using I1, I2, etc to define MCTDHFpulse0, MCTDHFpulse1, but it's
#better to name the pulses by their positions in the pulse list than otherwise
        #print ("self.pulsetype\t"+str(self.pulsetype))
        self.MCTDHFpulse0=MCTDHFpulse(pulsetype=self.pulsetype, w0=self.w0,
                                      wenv=self.wenv, I0=I1,
                                      pulsestrength=pulsestrength1,
                                      pulsetheta=pulsetheta1,
                                      tstart=self.tstart(),
                                      phaseshift=strtphase1, chirp=self.chirp,
                                      ramp=self.ramp, longstep=self.longstep)
        self.MCTDHFpulse1=MCTDHFpulse(pulsetype=self.pulsetype, w0=self.w0,
                                      wenv=self.wenv, I0=I2,
                                      pulsestrength=pulsestrength2,
                                      pulsetheta=pulsetheta2,
                                      tstart=self.tstart(),
                                      phaseshift=strtphase2, chirp=self.chirp,
                                      ramp=self.ramp, longstep=self.longstep)
        self.MCTDHFpulselist=[self.MCTDHFpulse0, self.MCTDHFpulse1]
        


    def twidth(self):
        return pi/self.wenv

    def settstart(self, tstart):
        self.tcenter=tstart+self.twidth()/2
    
    def E0(self):
        #print("self.I0\t"+str(self.I0))
        return sqrt(self.I0)
    
    def startphase(self):
        #MCTDHF code has pulse of the form A(t) sin(wenv t)^2
        #sin(wosc*t+startphase) however, it's often more convenient to
        #think about a pulse of the form 
        #A(t) cos(wenv(t-tmid))^2 cos(w(t-tmid)+centerphase).

        #This subroutine calculates startphase corresponding to a chosen
        #centerphase -- ie, so that 
        #sin(wosc*(t-tstart)+startphase)= cos(wosc*(t-tcenter)+centerphase),
        #where tcenter=tstart+(pi/2)/wenv
        accumulatedphase=self.twidth()/2*self.w0
        starttheta=mod(self.CEPtheta-accumulatedphase,2*pi)
        return starttheta
 
    def tstart(self):
        return self.tcenter-self.twidth()/2

    def tostringarray(self):
        #return list containing variable names and strings corresponding to their values
        return array([['I0', str(self.I0)], ['w0', str(self.w0)],
                      ['wenv',str(self.wenv)], ['CEPtheta',str(self.CEPtheta)],
                      ['jtheta',str(self.jtheta)], ['jphiz',str(self.jphiz)],
                      ['propeuler0',str(self.propeuler0)],
                      ['propeuler1',str(self.propeuler1)],
                      ['propeuler2',str(self.propeuler2)],
                      ['vprop',str(self.vprop)],
                      ['vpol1',str(self.vpol1)], 
                      ['vpol2',str(self.vpol2)], 
                      ['tcenter',str(self.tcenter)], 
                      ['pulsetype',str(self.pulsetype)], 
                      ['chirp',str(self.chirp)], 
                      ['ramp',str(self.ramp)], 
                      ['longstep',str(self.longstep)], 
                     ])
#    def todict(self):
#        return {
#            'I0':str(self.I0),
#            'wenv':str(self.wenv),
#            'jtheta':str(self.jtheta),
#            'propeuler0':str(self.propeuler0),
#            'propeuler1':str(self.propeuler1),
#            'propeuler2':str(self.propeuler2),
#            'vprop':str(self.vprop),
#            'vpol1':str(self.vpol1),
#            'vpol2':str(self.vpol2),
#            'tcenter':str(self.tcenter),
#            'pulsetype':str(self.pulsetype),
#            'chirp':str(self.chirp),
#            'ramp':str(self.ramp),
#            'longstep':str(self.longstep),
#            'pulsestrength':str(self.pulsestrength)
#        }
    
    def todict(self):
        dictlist=[]
        for i in range(len(self.MCTDHFpulselist)):
            MCTDHFpulse=self.MCTDHFpulselist[i]
            dictlist.append(MCTDHFpulse.todict()) 
        retdict=dictjoin(dictlist)
        return retdict


    def tostring(self):
        tmparray=self.tostringarray()
        return "\t".join(tmparray[:,1])

    def stringkey(self, pulseindx=0):
        tmparray=self.tostringarray()
        joinstr="_"+str(pulseindx)+"\t"
        return joinstr.join(tmparray[:,0])+"_"+str(pulseindx)

################################################################################
#exponential pulse class gives a pulse of the form 
#cos(wenv t)**2 e^(1j w2 t)=cos(wenv t)**2 * (cos(w2 t)+1j*sin(w2 t))
class exponential_pulse:
    def __init__(self, **kwargs):
#cosine component
        self.pulse0=pulse(**kwargs)
#sine component
        CEPtheta=self.pulse0.CEPtheta
        pulsestrength=self.pulse0.pulsestrength

        self.pulse1=pulse(**kwargs)
        self.pulse1.CEPtheta-=pi/2
        self.pulse1.CEPtheta=mod(self.pulse1.CEPtheta, 2*pi)
        self.pulse1.pulsestrength*=1j
        self.pulse1.initMCTDHF()

#make list of component pulses, MCTDHF pulses
        self.pulselist=[self.pulse0, self.pulse1]
        self.MCTDHFpulselist=self.pulse0.MCTDHFpulselist+self.pulse1.MCTDHFpulselist

    def todict(self):
        dictlist=[]
        for i in range(len(self.MCTDHFpulselist)):
            MCTDHFpulse=self.MCTDHFpulselist[i]
            dictlist.append(MCTDHFpulse.todict()) 
        retdict=dictjoin(dictlist)
        return retdict



################################################################################
#MCTDHF pulse class contains all information needed to write an MCTDHF input file
class MCTDHFpulse:
    def __init__(self, pulsetype=2, w0=.0565, wenv=.014, I0=1.e-4, pulsestrength=0, pulsetheta=0,
                 tstart=0., phaseshift=0, chirp=0, ramp=0, longstep=0):
#this class contains the information that can go into an MCTDHF pulse and
#returns a dictionary which maps the names of the variables to strings
#containing their value
        self.pulsetype=pulsetype
        self.pulsestrength=pulsestrength
        self.w0=w0
        self.wenv=wenv
        self.I0=I0
        self.pulsetheta=pulsetheta
        self.tstart=tstart
        self.phaseshift=phaseshift
        self.chirp=chirp
        self.longstep=longstep
    
    def todict(self):
        return {'#pulsetypelist#':str(self.pulsetype), '#omega2list#':str(self.w0),
                '#pulsestartlist#':str(self.tstart),
                '#omegalist#':str(self.wenv), '#intensitylist#':str(self.I0),
                '#pulsethetalist#':str(self.pulsetheta),
                '#CEPlist#':str(self.phaseshift), '#chirplist#':str(self.chirp),
                '#longsteplist#':str(self.longstep),
                '#pulsestrengthlist#':fortran_cmplx_str(self.pulsestrength)}


################################################################################
#Helper functions

def fortran_cmplx_str(zinp):
    rval=real(zinp)
    imval=imag(zinp)
    return "("+str(rval)+", "+str(imval)+")"

def FWHM_to_wenv(FWHM_fs):
    return pi/(FWHM_fs*2/aut)

def wenv_to_tmid(wenv):
    return (pi/2)/wenv

def phasearray(nphasepoints):
    return arange(0.,2.,2./nphasepoints)

#functions to choose dt,npts for desired energy resolution
def tparams(dE, maxE):
    #dE is desired energy resolution, maxE is maximum energy we wish to resolve
    dEHrt=dE/Hrt
    maxEHrt=maxE/Hrt
    ErangeHrt=2*maxEHrt
    npts=int(floor(ErangeHrt/dEHrt))

    #critical sampling is 2 samples per period
    dt=(2*pi)/(ErangeHrt)
    Tmax=dt*npts
    return dt,Tmax,npts



def pulselistkey(pulselist):
    retstr=''
    for i in range(len(pulselist)):
        retstr+=pulselist[i].stringkey(pulseindx=i)+"\t"
    return retstr+"\n"
def printloopkey(resultdir, loopnamelist, looplist, filename='loopkey.txt'):
    tmpfile=open(resultdir+filename,'w')
    tmpfile.write("arrays looped over in this calculation (outer loops first, inner loops last)\n")
    for i in range(len(loopnamelist)):
#        tmpfile.write(loopnamelist[i]+"\t"+str(looplist[i])+"\n")
        tmpfile.write(loopnamelist[i]+"\t"+"\t".join(map(str,looplist[i]))+"\n")
    tmpfile.close


def printparameterkey_full(resultdir, pulsearray, filename='parameterkey_full.txt'):
    tmpfile=open(resultdir+filename,'w')
    tmpfile.write(pulselistkey(pulsearray[0]))
    for i in range(len(pulsearray)):
        tmpfile.write(pulselisttostring(pulsearray[i]))
    tmpfile.close


def dictjoin(dictlist, joinstr=", "):
#join the strings saved in different pulse dictionaries
    retdict={}
    for key in dictlist[0]:
        tmpstr=''
        for i in range(len(dictlist)):
            tmpstr+=str(dictlist[i][key])+", "
        retdict[key]=tmpstr[:-2]
    return retdict

def pulselistdictionary(pulselist):
#first, make a list of the MCTHDF pulse list dictionaries
    pulsedictlist=[]
    for i in range(len(pulselist)):
        for j in range(len(pulselist[i].MCTDHFpulselist)):
            pulsedictlist+=[pulselist[i].MCTDHFpulselist[j].todict()]
    dict0=dictjoin(pulsedictlist)
    return dict0

def makeinputfiles(pulsearray, resultdir):
    masterscriptfile="fullcalculation.sh"
    for calcindx in range(len(pulsearray)):
        pulselist=pulsearray[calcindx]
        dirstr=resultdir+str(calcindx)+"/"
        subprocess.call(["mkdir", dirstr])
        
        initialrelaxation(masterscriptfile,dirstr)
        multiplepulses(masterscriptfile,dirstr,1,pulselist)

def initialrelaxation(masterscriptfile,dirstr):
    bashtemplate="templates/Relax.Bat.template"
    inputtemplate="templates/Input.Inp.Relax.template"

    steproot=stepstr(0)
    
    bashfilename=steproot+".Bat"
    inpfilename="Input.Inp."+steproot
    outfilename="Out.States."+steproot

    #list of replacement rules
    dictlist=[inpoutdictionary(steproot,steproot),filedictionary(inpfilename,outfilename)]

    #use dictlist to convert template files to usable input & script files
    templatereplace(bashtemplate,dirstr+bashfilename,dictlist)
    templatereplace(inputtemplate,dirstr+inpfilename,dictlist)

    masterscript=open(dirstr+masterscriptfile,'w')
    masterscript.write("bash "+bashfilename+"\n")
    masterscript.close()



def multiplepulses(masterscriptfile, dirstr, stepindx, pulselist):
    dict0=pulselistdictionary(pulselist)

    
    #next, make a dictionary corresponding to the values which aren't contained
    #in the previous dictionary
    numpulses=sum(map(lambda pls: len(pls.MCTDHFpulselist), pulselist))#2*len(pulselist)
    if(velocitygauge):
        velstr="1"
    else:
        velstr="0"
    dict1={'#numpulses#':str(numpulses), '#timestep#':str(timestep),
           '#measurementtimestep#':str(measurementtimestep),
           '#tfinal#':str(tfinal), "#velflag#":velstr}
     
    inpstr=stepstr(stepindx-1)
    outstr=stepstr(stepindx)
    steproot=stepstr(stepindx)

    bashtemplate="templates/Multiplepulses.Bat.template"
    inputtemplate="templates/Input.Inp.Multiplepulses.template"

    bashfilename=steproot+".Bat"
    inpfilename="Input.Inp."+steproot
    outfilename="Out.States."+steproot

    dictlist=[dict0, dict1, inpoutdictionary(inpstr, outstr),
              filedictionary(inpfilename, outfilename)]
    #use dictlist to convert template files to usable input & script files
    templatereplace(bashtemplate,dirstr+bashfilename,dictlist)
    templatereplace(inputtemplate,dirstr+inpfilename,dictlist)

    masterscript=open(dirstr+masterscriptfile,'a')
    masterscript.write("bash "+bashfilename+"\n")
    masterscript.close()


def inpoutdictionary(inpstr,outstr):
    dict={'#inpstr#':inpstr,'#outstr#':outstr}
    return dict

def filedictionary(inpfilename,outfilename):
    dict={'#inpfilename#':inpfilename,'#outfilename#':outfilename,'#gridstr#':gridstring,'#atomstr#':atomstring}
    return dict




def replacementdict(pulselist):
    dictlist=[]
    for i in range(len(pulselist)):
        for j in range(len(pulselist[i].MCTDHFpulselist)):
            dictlist+=pulselist[i].MCTDHFpulselist[j].todict()

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
    
def stepstr(n):
    return "step_"+str(n)



def Rx(theta):
#rotation matrix about x axis
    return [[1,0,0],[0,cos(theta),-sin(theta)],[0,sin(theta),cos(theta)]]

def Rz(theta):
    #rotation matrix about z axis
    return [[cos(theta),-sin(theta),0],[sin(theta),cos(theta),0],[0,0,1]]

def Ry(theta):
#rotation matrix about y axis
    return [[cos(theta),0,sin(theta)],[0,1,0],[-sin(theta),0,cos(theta)]]

def topolar(cartvec):
#convert cartesian vector to polar coords
    vx=cartvec[0]
    vy=cartvec[1]
    vz=cartvec[2]
    r=sqrt(vx**2+vy**2+vz**2)
    theta=arccos(vz/r)
    phi=angle(vx+1j*vy)
    return r, theta, phi

def nm_to_eV(nm):
    hc=1239.841#planck's constant in eV*nm
    eV=hc/nm
    return eV

def eV_to_nm(eV):
    hc=1239.841#planck's constant in eV*nm
    nm=hc/eV
    return nm

def dictionaryreplace(text,dictionary):
    newtext=text
    for i,j in dictionary.items():#.iteritems():
        newtext=newtext.replace(i,j)
    return newtext

################################################################################
#Main program

##test initiation
#pls1=pulse()
#print("check to see that pls1 is initiated")
#print("pls1 vprop\t"+str(pls1.vprop))
#print("pls1 vpol1\t"+str(pls1.vpol1))
#print("pls1 vpol2\t"+str(pls1.vpol2))


evstr=sys.argv[-1]
if(len(sys.argv)<2):
    evstr="22"
evfloat=float(evstr)
resultdir="results_"+str(evstr)+"eV/"
#resultdir=sys.argv[-2]#"results_tmp/"
if(resultdir[-1]!="/"):
    resultdir=resultdir+"/"

#copy this program into resultdir
subprocess.call(['mkdir',resultdir])
shutil.copy('./corotating_counterrotating_streaking.py',resultdir)

#copy templates directory to resultdir
templatetargetstr=resultdir+"templates"
if(os.path.exists(templatetargetstr)):
    shutil.rmtree(templatetargetstr)
shutil.copytree("./templates",templatetargetstr)

#read grid parameters from file
tmpfile=open("templates/gridfile.txt")
gridstring="".join(tmpfile.readlines())
tmpfile.close()
#read atom parameters from file
tmpfile=open("templates/atomfile.txt")
atomstring="".join(tmpfile.readlines())
tmpfile.close()

################################################################################
#main program

velocitygauge=False#True#True: use velocity gauge.  False: use length gauge

timestep=0.05
wir=nm_to_eV(800)/Hrt
wxuv=evfloat/Hrt

Iir=3e-4
Ixuv=1.e-6
IRduration=11e3#pulse FWHM in attoseconds        
IRwenv=FWHM_to_wenv(IRduration)
IRtmid_atomic=wenv_to_tmid(IRwenv)#midpoint of IR pulse in atomic units
XUVduration=11e3
XUVwenv=FWHM_to_wenv(XUVduration)



dE=.25#.08#desired energy resolution in eV
maxE=40.#20.#desired maximum energy range in eV
dtdip,Tdipmax,ndtdip=tparams(dE,maxE)
dtdip=ceil(dtdip/timestep)*timestep#make interpulse delays an integral
                               #number of timesteps
measurementtimestep=dtdip#interval between calculating dipole
tfinal=Tdipmax#total interval to propagate when calculating dipole



#make array corresponding to delay between the centers of the two pulses
dtstart=0.
dtstop=0.
deltadt=100#
dtarray=arange(dtstart,dtstop+deltadt,deltadt)/aut
ndt=len(dtarray)

#number of phase points for XUV and IR pulses
nxuvphase=5
xuvphasearray=phasearray(nxuvphase)

ncophase=5
cophasearray=phasearray(ncophase)

ncounterphase=5
counterphasearray=phasearray(ncounterphase)

paramindx=0
pulsearray=[]
loopnamelist=['dtarray','xuvphasearray','cophasearray','counterphasearray']
looplist=[dtarray, xuvphasearray, cophasearray, counterphasearray]
printloopkey(resultdir,loopnamelist, looplist)
for i,j,k,l in iterproduct(range(ndt), range(nxuvphase), range(ncophase),
                             range(ncounterphase)):
    dt=dtarray[i]
    xuvphase=xuvphasearray[j]*pi
    cophase=cophasearray[k]*pi
    counterphase=counterphasearray[l]*pi
    #print("sumphase\t"+str(sumphasearray[l])+"\t"+str(sumphase))
    #IR pulse is split into two pulses: E(t)= 1/sqrt(2)*E0*(cos(wt+IRphase+sumphase)+cos(wt+IRphase-sumphase))
    pulselist=[exponential_pulse(I0=Iir/sqrt(2), w0=wir, wenv=IRwenv,
                                 CEPtheta=cophase,tstart=0),
               exponential_pulse(I0=Iir/sqrt(2), w0=-wir, wenv=IRwenv,
                                 CEPtheta=counterphase,tstart=0),
               pulse(I0=Ixuv, w0=wxuv, wenv=XUVwenv, CEPtheta=xuvphase,
                     tcenter=IRtmid_atomic+dt)]
#    pulselist=[pulse(I0=Iir/sqrt(2),  w0=wir, wenv=IRwenv,
#                     CEPtheta=IRphase+sumphase, tstart=0),
#               pulse(I0=Iir/sqrt(2),  w0=wir, wenv=IRwenv,
#                     CEPtheta=IRphase-sumphase, tstart=0),
#               pulse(I0=Ixuv, w0=wxuv, wenv=XUVwenv,
#                     CEPtheta=xuvphase, tcenter=IRtmid_atomic+dt)]
    pulsearray.append(pulselist) 

#make input files corresponding to this calculation
makeinputfiles(pulsearray, resultdir)
