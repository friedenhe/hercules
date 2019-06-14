# Makefile for HERCULES_V1.1.0

SRCS_f90d1 = \
module_boundary.f90 \
module_cnvdiff.f90 \
module_spectral.f90 \
module_variables.f90 \
module_navier_stokes.f90 \
module_poisson_solver.f90 \
module_parameters.f90 \
module_io.f90 \
main.f90 \
module_tools.f90 

OBJS_f90d1 = \
module_boundary.o \
module_cnvdiff.o \
module_spectral.o \
module_variables.o \
module_navier_stokes.o \
module_poisson_solver.o \
module_parameters.o \
module_io.o \
main.o \
module_tools.o 

SRC_DIR_f90d1 = 
OBJS_DIR = obj/
EXE_DIR = ./

EXE = hercules.exe
FC = mpif90
IDIR = -Iinclude
CFLAGS = -O3 -DDOUBLE_PREC -J$(OBJS_DIR) $(IDIR)
LFLAGS = 
LIBS = lib/lib2decomp_fft.a lib/libfftw3.a 

VPATH = $(SRC_DIR_f90d1):$(OBJS_DIR)
OBJS = $(addprefix $(OBJS_DIR), $(OBJS_f90d1))

all : $(EXE)

$(EXE) : $(OBJS_f90d1)
	@mkdir -p $(EXE_DIR)
	$(FC) -o $(EXE_DIR)$(EXE) $(OBJS) $(LFLAGS) $(LIBS)

$(OBJS_f90d1):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90d1)$(@:.o=.f90) -o $(OBJS_DIR)$@

clean :
	rm -f $(OBJS_DIR)*.*
	rm -f $(EXE_DIR)$(EXE)

# Dependencies of files
module_boundary.o: \
    module_boundary.f90 \
    module_parameters.o
module_cnvdiff.o: \
    module_cnvdiff.f90 \
    module_parameters.o
module_spectral.o: \
    module_spectral.f90 \
    module_parameters.o \
    module_tools.o \
    module_variables.o
module_variables.o: \
    module_variables.f90 \
    module_parameters.o
module_navier_stokes.o: \
    module_navier_stokes.f90 \
    module_boundary.o \
    module_cnvdiff.o \
    module_io.o \
    module_parameters.o \
    module_poisson_solver.o \
    module_spectral.o \
    module_tools.o \
    module_variables.o
module_poisson_solver.o: \
    module_poisson_solver.f90 \
    module_parameters.o \
    module_tools.o \
    module_variables.o
module_parameters.o: \
    module_parameters.f90
module_io.o: \
    module_io.f90 \
    module_parameters.o \
    module_tools.o
main.o: \
    main.f90 \
    module_navier_stokes.o \
    module_parameters.o
module_tools.o: \
    module_tools.f90 \
    module_parameters.o \
    module_variables.o

