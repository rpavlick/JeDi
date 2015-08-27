BLDDIR = $(abspath .)
include $(BLDDIR)/make.inc

PPOPTS=$(MPIOPT) $(CMP_OPT) $(DEBUGOPT) $(NETCDFOPT)

# *****************************
# * generate a version string *
# *****************************
dummy:=$(shell rm -f version.f90)

GIT_VERSION := $(shell git describe --dirty --long --always --all)
DATE:=$(shell date)
UNIX:=$(shell uname -n; echo -n ": "; uname -s; uname -r; uname -m; uname -p)

dummy:=$(shell echo creating version.f90 >&2)
dummy:=$(shell echo 'subroutine globe_show_version ()' >version.f90)
dummy:=$(shell echo 'write(*,*) "JeDi-DGVM: git revision $(GIT_VERSION)"' >>version.f90)
dummy:=$(shell echo 'write(*,*) "compiled on: $(DATE)"' >>version.f90)
dummy:=$(shell echo 'write(*,*) "compiled at host: $(UNIX)"' >>version.f90)
dummy:=$(shell echo 'write(*,*) "compiler used:"' >>version.f90)
dummy:=$(shell $(FC) $(VOPT) 2>&1 | awk '{print "write(*,*) \"" $$0 "\""}' >>version.f90)
dummy:=$(shell echo 'end subroutine globe_show_version' >>version.f90)

dummy:=$(shell $(FC) $(FOPTS) $(REAL8OPT) -c -o version.o version.f90)
VERSION=version.o

# ********************
# * model file lists *
# ********************

GLOBE_FILES=globe_mod_func.F90 globe_mod_flux.F90 globe_mod_stat.F90 \
						globe_mod.F90 globe.F90 globe_sub.F90 globe_fio.F90 \
						globe_surf.F90 globe_functions.F90  \
						globe_mpimod.F90
JAM_FILES=jam_mod.F90 jam.F90 jam_fio.F90 jam_sub.F90 jam_globe.F90 \
						jam_carbon.F90
JEDI_FILES=jedi_mod.F90 jedi_mod_dyn.F90 jedi.F90 jedi_sub.F90 \
						jedi_fio.F90 jedi_globe.F90 jedi_dyn.F90 \
						jedi_opti.F90 jedi_spec.F90 jedi_kristin.f90

# ***********************************************************
# * model configuration                                     *
# ***********************************************************
PROG_NAME=jedi.x

DIR_GLOBE=./globe/
SRC_GLOBE=$(GLOBE_FILES)
DIR_ATMOS=./jam/
SRC_ATMOS=$(JAM_FILES)
DIR_SOIL=./nosoil/
SRC_SOIL=nosoil_globe.f90
DIR_VEG=./jedi/
SRC_VEG=$(JEDI_FILES)

# setup the make variables
SRC=$(SRC_GLOBE) $(SRC_ATMOS) $(SRC_VEG) $(SRC_SOIL)
TMP=$(SRC:.f90=.o)
OBJ=$(TMP:.F90=.o) $(VERSION)
VPATH=$(DIR_GLOBE) $(DIR_ATMOS) $(DIR_VEG) $(DIR_SOIL) .

# list all needed object files and make the binary
jedi.x: $(OBJ)
	@echo linking $@
	@$(FC) $(FOPTS) $(FPPOPT) $(REAL8OPT) -o $@ $(OBJ) $(MPIOPT) $(LIBNETCDF) $(LIBMATH)

# make all listed object files
%.o : %.f90
	@echo compiling $<
	@$(FC) $(FOPTS) $(FPPOPT) $(INCNETCDF) $(REAL8OPT) -c -o $@ $<

%.o : %.F90
	@echo compiling $<
	@$(FC) $(FOPTS) $(FPPOPT) $(INCNETCDF) $(REAL8OPT) $(PPOPTS) -c -o $@ $<
