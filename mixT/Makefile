#include Makefile.$(OSTYPE)
include Makefile.linux


OBJ = mk_cal.o sigism.o readrmf.o calutil.o mekal.o util.o


all: librp.a mk_cal weightT

librp.a:
	(cd ReadPar; make)

mk_cal: $(OBJ) librp.a 
	$(FC) -o mk_cal $(OBJ) $(LIBFITS) -L./ -lrp

weightT:  weightT.o util2.o
	$(FC) -o weightT $^ -L./ -lrp


# If have troubles compiling mekal.o, setup the rule something like this:
#
# .F.o:
# 	$(FC) $(FFLAGS) $(CPPFLAGS) -F -o /tmp/fcc.gov.f $<
# 	$(FC) $(FFLAGS) -c -o $@  /tmp/fcc.gov.f


clean:
	rm ReadPar/*.o *.o *.a tcal.dat

