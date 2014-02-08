NDEBUG :=  
NDEBUG := t
MPI    :=
OMP    :=

#COMP := Intel
COMP := gfortran

MKVERBOSE := t

FPARALLEL := ../../../..

include $(FPARALLEL)/mk/GMakedefs.mak

#F90FLAGS += -assume byterecl

f90sources += testeos.f90

include $(FPARALLEL)/boxlib/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/boxlib

include $(FPARALLEL)/extern/EOS/collapse/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/extern/EOS/collapse
FINCLUDE_LOCATIONS += $(FPARALLEL)/extern/EOS/collapse

include $(FPARALLEL)/extern/networks/collapse/GPackage.mak
VPATH_LOCATIONS += $(FPARALLEL)/extern/networks/collapse

#include $(FPARALLEL)/extern/VODE/GPackage.mak
#VPATH_LOCATIONS += $(FPARALLEL)/extern/VODE

testeos.$(suf).exe: $(objects)
	$(LINK.f90) -o testeos.$(suf).exe $(objects) $(libraries)

include $(FPARALLEL)/mk/GMakerules.mak

