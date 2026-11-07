#==============================================================================
# jupiter7 viztools  Makefile
#
#   make                       # gfortran + OpenMP; bin2hdf with raw output only
#   make COMPILER=intel        # ifx (add FC=ifort for ifort)
#   make COMPILER=nec          # NEC SX-Aurora nfort
#   make COMPILER=fujitsu      # Fujitsu frtpx
#   make COMPILER=nvidia       # nvfortran
#   make OMP=no                # without OpenMP (serial)
#   make HDF5=yes HDF5_DIR=/path/to/hdf5   # enable HDF5 output
#   make test                  # regression tests on small synthetic data (gfortran, python3+numpy)
#==============================================================================
COMPILER ?= gnu
OMP      ?= yes
HDF5     ?= no
HDF5_DIR ?= /usr

ifeq ($(COMPILER),gnu)
  FC        = gfortran
  FFLAGS   ?= -O2
  OMPFLAG   = -fopenmp
  MODFLAG   = -J
endif
ifeq ($(COMPILER),intel)
  # -fp-model precise: the Intel default (fast) optimizes assuming there are no NaNs, so the
  # missing-value (NaN) test ieee_is_nan may stop working
  # (kept separate in FPFLAGS so that it is added even when FFLAGS is given)
  FC        = ifx
  FFLAGS   ?= -O2
  FPFLAGS   = -fp-model precise
  OMPFLAG   = -qopenmp
  MODFLAG   = -module
endif
ifeq ($(COMPILER),nec)
  FC        = nfort
  FFLAGS   ?= -O3
  OMPFLAG   = -fopenmp
  MODFLAG   = -module
endif
ifeq ($(COMPILER),fujitsu)
  FC        = frtpx
  FFLAGS   ?= -Kfast
  OMPFLAG   = -Kopenmp
  MODFLAG   = -M
endif
ifeq ($(COMPILER),nvidia)
  FC        = nvfortran
  FFLAGS   ?= -O2
  OMPFLAG   = -mp
  MODFLAG   = -module
endif
ifeq ($(OMP),no)
  OMPFLAG   =
endif

ifeq ($(HDF5),yes)
  # HDF5 include / library locations (change for your system)
  H5INC   ?= -I$(HDF5_DIR)/include
  H5LIB   ?= -L$(HDF5_DIR)/lib -lhdf5_fortran -lhdf5 -lz
  H5DEF    = -DUSE_HDF5
endif

SRC  = src
OBJD = build
BIN  = bin
ALLFLAGS = $(FFLAGS) $(FPFLAGS) $(OMPFLAG) $(MODFLAG) $(OBJD) -I$(OBJD)

PROGS = $(BIN)/combine3d $(BIN)/vorticity $(BIN)/slice_plot $(BIN)/bin2hdf $(BIN)/polar_map $(BIN)/mag3d

all: $(PROGS)

$(OBJD) $(BIN):
	mkdir -p $@

$(OBJD)/mhd_common.o: $(SRC)/mhd_common.f90 | $(OBJD)
	$(FC) $(ALLFLAGS) -c $< -o $@

$(OBJD)/mhdunits.o: $(SRC)/mhdunits.f90 | $(OBJD)
	$(FC) $(ALLFLAGS) -c $< -o $@

$(OBJD)/psplot.o: $(SRC)/psplot.f90 $(OBJD)/mhd_common.o
	$(FC) $(ALLFLAGS) -c $< -o $@

$(BIN)/combine3d: $(SRC)/combine3d.f90 $(OBJD)/mhd_common.o | $(BIN)
	$(FC) $(ALLFLAGS) $< $(OBJD)/mhd_common.o -o $@

$(BIN)/vorticity: $(SRC)/vorticity.f90 $(OBJD)/mhd_common.o | $(BIN)
	$(FC) $(ALLFLAGS) $< $(OBJD)/mhd_common.o -o $@

$(BIN)/slice_plot: $(SRC)/slice_plot.f90 $(OBJD)/psplot.o $(OBJD)/mhdunits.o $(OBJD)/mhd_common.o | $(BIN)
	$(FC) $(ALLFLAGS) $< $(OBJD)/psplot.o $(OBJD)/mhdunits.o $(OBJD)/mhd_common.o -o $@

$(OBJD)/mhdfield.o: $(SRC)/mhdfield.f90 $(OBJD)/mhd_common.o $(OBJD)/mhdunits.o
	$(FC) $(ALLFLAGS) -c $< -o $@

$(OBJD)/scene3d.o: $(SRC)/scene3d.f90 $(OBJD)/mhd_common.o
	$(FC) $(ALLFLAGS) -c $< -o $@

$(OBJD)/isosurf.o: $(SRC)/isosurf.f90 $(OBJD)/mhd_common.o
	$(FC) $(ALLFLAGS) -c $< -o $@

MAG3D_OBJ = $(OBJD)/mhdfield.o $(OBJD)/scene3d.o $(OBJD)/isosurf.o $(OBJD)/psplot.o $(OBJD)/mhdunits.o $(OBJD)/mhd_common.o
$(BIN)/mag3d: $(SRC)/mag3d.f90 $(MAG3D_OBJ) | $(BIN)
	$(FC) $(ALLFLAGS) $< $(MAG3D_OBJ) -o $@

$(BIN)/polar_map: $(SRC)/polar_map.f90 $(OBJD)/psplot.o $(OBJD)/mhdunits.o $(OBJD)/mhd_common.o | $(BIN)
	$(FC) $(ALLFLAGS) $< $(OBJD)/psplot.o $(OBJD)/mhdunits.o $(OBJD)/mhd_common.o -o $@

$(BIN)/bin2hdf: $(SRC)/bin2hdf.F90 $(OBJD)/mhdunits.o $(OBJD)/mhd_common.o | $(BIN)
	$(FC) $(ALLFLAGS) $(H5DEF) $(H5INC) $< $(OBJD)/mhdunits.o $(OBJD)/mhd_common.o $(H5LIB) -o $@

test: all
	cd test && ./run_tests.sh

clean:
	rm -rf $(OBJD) $(BIN)

.PHONY: all clean test
