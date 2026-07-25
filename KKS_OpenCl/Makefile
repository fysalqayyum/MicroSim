FDIR=functions
SDIR=solverloop
TDIR=tdbs

# Site configuration. All values may be overridden on the make command line.
# Example:
#   make microsim_kks_opencl CC=/opt/mpi/bin/mpicc \
#        GSL_ROOT=/opt/gsl CUDA_ROOT=/usr/local/cuda
GSL_ROOT ?= $(HOME)/gsl_2_7
CUDA_ROOT ?= /usr/local/cuda
HDF5_INCLUDE ?= /usr/include/hdf5/openmpi

CC=mpicc
HC=h5pcc
CFLAGS=-I. -I$(GSL_ROOT)/include -I$(CUDA_ROOT)/include -I$(HDF5_INCLUDE)

LDIR=-L$(CUDA_ROOT)/lib64 -L$(GSL_ROOT)/lib

_FUNCTION_DEPS = global_vars.h CL_global_vars.h CL_initialize_variables.h      \
                 CL_device_kernel_build.h CL_buffer_allocation.h CL_create_kernel_args.h  \
                 CL_Initialize_domain.h CL_kernel_init_temperature.h functions.h matrix.h  \
                 utility_functions.h functionH.h functionF_01.h functionF_02.h functionF_03.h  \
                 CL_functionF_03_HostUpdate.h functionF_04.h CL_functionF_04_HostUpdate.h  \
                 initialize_variables.h filling.h reading_input_parameters.h read_boundary_conditions.h  \
                 free_variables.h fill_domain.h CL_DeviceToHost.h Temperature_gradient.h shift.h CL_Shift.h CL_functionF_03_HostUpdate.h 
                 
DEPS = $(patsubst %,$(FDIR)/%,$(_FUNCTION_DEPS))

_SOLVERLOOP_DEPS = serialinfo_xy.h  file_writer.h initialize_functions_solverloop.h mpi_xy.h \
			CL_Solve_phi_com_Function_F_2.h CL_Solve_phi_com_Function_F_3.h CL_Solve_phi_com_Function_F_4.h \
			CL_Solve_phi_com.h CL_Update_Temperature.h defines.h defines1.h Function_Thermo_calls.h FunctionDefines_CL.h \
			FunctionF_2_c_mu_NR_GEData2.h functionH_CL_5th.h functionW_02_CL.h GibbsEnergyData2.h matrix_nsolCL.h FunctionF_4_SplineCPU.h CL_Kim_Kernel.cl ThermoCL_2.c ThermoCL_2.h 

DEPS += $(patsubst %,$(SDIR)/%,$(_SOLVERLOOP_DEPS))

_tdb_DEPS = Thermo.c Thermo.h

DEPS += $(patsubst %,$(TDIR)/%,$(_tdb_DEPS))

LIBS=-lm -lgsl -lgslcblas -lOpenCL
WARNFLAGS ?= -Wall -Wextra

all : microsim_kks_opencl reconstruct

microsim_kks_opencl : microsim_kks_opencl.o
	$(CC) -g -o microsim_kks_opencl microsim_kks_opencl.o $(CFLAGS) $(IDIR) $(LDIR) $(LIBS) $(WARNFLAGS)

microsim_kks_opencl.o : $(DEPS)


_HEADERS_RECONSTRUCT = global_vars.h CL_global_vars.h functions.h matrix.h utility_functions.h reading_input_parameters.h


DEPS_RECONSTRUCT = $(patsubst %,$(FDIR)/%,$(_HEADERS_RECONSTRUCT))



# echo "Making reconstruct for collating different processor files into consolidated .vtk files. Valid for non .h5 files"

reconstruct : reconstruct.o
	$(HC) -o reconstruct reconstruct.o $(CFLAGS) $(LDIR) $(LIBS)

reconstruct.o : $(DEPS_RECONSTRUCT)

.PHONY : clean

clean :
	rm -f microsim_kks_opencl.o microsim_kks_opencl reconstruct *.out *.o  *.dat *.gz *.vtk *.txt *.mv *.log *.err *.epsc *.eps *.mat  *.png *.op *.jpg *.jpeg *.csl *.bd   *.bd *.out
