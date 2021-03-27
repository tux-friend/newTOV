include ./make.inc

SOURCES=eosmodule.F90 readtable.F90 nuc_eos.F90 bisection.F90 findtemp.F90 findrho.F90 linterp_many.F90
FSOURCES=linterp.f

CLEANSTUFF=rm -rf *.o *.mod *.a

OBJECTS=$(SOURCES:.F90=.o)
FOBJECTS=$(FSOURCES:.f=.o)

EXTRADEPS=

MODINC=$(HDF5INCS)

all: nuc_eos.a newTON Testing

newTON: nuc_eos.a newTON.F90
	$(F90) $(F90FLAGS) -o newTON newTON.F90 nuc_eos.a $(HDF5LIBS)

Testing: nuc_eos.a Testing.F90
	$(F90) $(F90FLAGS) -o Testing Testing.F90 nuc_eos.a $(HDF5LIBS)
       
nuc_eos.a: $(OBJECTS) $(FOBJECTS)
	ar r nuc_eos.a *.o

$(OBJECTS): %.o: %.F90 $(EXTRADEPS)
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) -c $< -o $@

$(FOBJECTS): %.o: %.f $(EXTRADEPS)
	$(F90) $(F90FLAGS) $(DEFS) $(MODINC) -c $< -o $@

clean: 
	$(CLEANSTUFF)
