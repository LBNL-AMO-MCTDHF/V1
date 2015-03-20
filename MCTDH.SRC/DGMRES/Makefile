
include ../Makefile.header

#dgmrcso = \
### 121212       dgmres.o   dorth.o   dpigmr.o   drlcal.o   dxlcal.o   isdgmr.o
#       dgmres.o   dorth.o   dpigmr.o   drlcal.o   dxlcal.o   isdgmr.o

dgmrcso = \
       dgmres.o   dorth.o   dpigmr.o   drlcal.o   dxlcal.o   isdgmr.o \
	dgmres-parallel.o dnrm2-parallel.o dorth-parallel.o dpigmr-parallel.o isdgmr-parallel.o




default: go

clean:
	rm *.o *.mod *~ *genmod*

go: $(dgmrcso)


%.o : %.f ../Makefile.header ; $(FC) $(F77FLAGS) -c $<

#.f.o: 
#	$(FC) $(F77FLAGS) -c $<



