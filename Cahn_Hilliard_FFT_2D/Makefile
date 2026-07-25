FDIR=functions
SDIR=solverloop
CC=gcc

# GSL (gsl_rng) and FFTW3 are required but not always on the default
# include/lib search path. On macOS neither ships with the system, so
# default to the Homebrew prefixes; override GSL_PREFIX/FFTW_PREFIX on
# the command line if your libraries live elsewhere (e.g. a conda env).
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
GSL_PREFIX  ?= $(shell brew --prefix gsl 2>/dev/null)
FFTW_PREFIX ?= $(shell brew --prefix fftw 2>/dev/null)
CFLAGS = -I. -I$(GSL_PREFIX)/include -I$(FFTW_PREFIX)/include
RPATH_LIBS = -Wl,-rpath,$(GSL_PREFIX)/lib -Wl,-rpath,$(FFTW_PREFIX)/lib
else
GSL_PREFIX  ?=
FFTW_PREFIX ?=
CFLAGS = -I.
RPATH_LIBS =
endif

_FUNCTION_DEPS = global_vars.h functions.h matrix.h utility_functions.h  \
                 filling.h reading_input_parameters.h free_variables.h fill_domain.h

DEPS = $(patsubst %,$(FDIR)/%,$(_FUNCTION_DEPS))

_SOLVERLOOP_DEPS = serialinfo_xy.h GibbsEnergyData.h functions_fftw.h  file_writer.h

DEPS += $(patsubst %,$(SDIR)/%,$(_SOLVERLOOP_DEPS))


LIBS = -L$(GSL_PREFIX)/lib -L$(FFTW_PREFIX)/lib -lfftw3 -lgsl -lm

microsim_ch_fft : microsim_ch_fft.o
	$(CC) -o microsim_ch_fft microsim_ch_fft.o $(CFLAGS) $(LIBS) $(RPATH_LIBS)

microsim_ch_fft.o : $(DEPS)       

.PHONY : clean

clean :

	rm -r microsim_ch_fft.o microsim_ch_fft   


