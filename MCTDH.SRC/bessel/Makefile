
include ../Makefile.header

bessrcso = \
	./dgamln.o ./fdump.o ./i1mach.o ./j4save.o \
	./xercnt.o ./xerhlt.o ./xermsg.o \
	./xerprn.o ./xersve.o ./xgetua.o ./zabs.o \
	./zacai.o ./zairy.o ./zasyi.o \
	./zbesy.o ./zbesh.o ./zbesj.o ./zbinu.o \
	./zbknu.o ./zbuni.o ./zdiv.o \
	./zexp.o ./zkscl.o ./zlog.o ./zmlri.o \
	./zmlt.o ./zrati.o ./zs1s2.o ./d1mach.o \
	./zseri.o ./zshch.o ./zsqrt.o ./zuchk.o \
	./zunhj.o ./zuni1.o ./zuni2.o \
	./zunik.o ./zuoik.o ./zwrsk.o ./djairy.o \
	./dbesj.o ./dasyjy.o ./initds.o \
	./dlngam.o ./dgamma.o ./dgamln.o ./dgamlm.o \
	./dcsevl.o ./d9lgmc.o ./zacon.o \
	./zbunk.o ./zunk1.o ./zunk2.o ./drc3jj.o

default: go

clean:
	rm *.o *.mod *~  *genmod*

go: $(bessrcso)



.f.o: 
	$(FC) $(F77FLAGS) -c $<




