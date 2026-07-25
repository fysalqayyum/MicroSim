SHELL = /bin/sh
# Override on the command line if a stale nvcc precedes the intended
# toolkit on PATH, e.g.: make NVCC=/usr/local/cuda/bin/nvcc
NVCC ?= nvcc
COMPOPS = -Xptxas -O3
#COMPOPS = -Xptxas -g
# Default targets P100/V100 (sm_60/sm_70) per this solver's README.
# Override for other GPUs, e.g.: make ARCH=-arch=sm_80   (A100)
ARCH ?= -arch=sm_70
LIBS =  -lcufft -lcuda -lcudart
#INCS =  -I/usr/local/cuda-11.2/samples/common/inc
HEADERS = functions/global_vars.h

SRCDIR = .
BINDIR = DATA
SRCFILES = $(SRCDIR)/KKS_cuFFT.cu
all: build

build: build_dir $(HEADERS)
	$(NVCC) $(COMPOPS) $(ARCH) $(SRCFILES) -o microsim_kks_cufft $(LIBS)

build_dir:
	mkdir -p $(BINDIR)

clear:
	-\rm $(BINDIR)/*.*

clean:
	-\rm *microsim_kks_cufft
	-\rm -r $(BINDIR)
#	@\rm -f job.*.*
## End of the makefile
