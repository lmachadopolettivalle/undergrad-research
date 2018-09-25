include ../Makefile.linux

CPPFLAGS = $(READLINE) -DEXPANDVARS=1 $(SHAREDFLAGS)

OBJ = recordtype.o remleadblanks.o parstring.o irafstring.o parse_p_str.o \
	chomp.o getarg.o lcase.o commandline.o logics.o parse_cont.o \
	parfile.o parse_group.o readparf77_c.o readparf77.o \
	get_pv_def.o expand_vars.o cl_r_par.o stringpar.o numbered_par.o

all: ../librp.a
#all: ../librp.a ../librp.so

../librp.a: ../librp.a($(OBJ))
	ranlib ../librp.a

../librp.so: $(OBJ)
	$(SHAREDLINK) -o ../librp.so $(READLINELIBS) $(OBJ)

distclean:
	rm -f ../librp.a ../librp.so *.o
