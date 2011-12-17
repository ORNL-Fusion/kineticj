NAME := bin/kineticj

# Defaults for dlg-hp.ornl.gov
GCCDIR := /home/dg6/code/gcc/gcc-${GNUVER}/bin
ALGLIBDIR := /home/dg6/code/alglib/cpp/src
NETCDFDIR := /home/dg6/code/netcdf/gnu_${GNUVER}# must be an --enable-cxx-4 dist
CUDADIR := /home/dg6/code/cuda/4.1/cuda
CUDALIBDIR = ${CUDADIR}/lib64
CUDA_ARCH := sm_13
CUDA_SDK_DIR := /home/dg6/code/cuda/NVIDIA_GPU_Computing_SDK
LIBCONFIGDIR := /home/dg6/code/libconfig
GOOGLE_PERF_DIR := /home/dg6/code/google-perftools
PAPI_DIR :=/home/dg6/code/papi/gnu_${GNUVER}

# Catch for greendl.* (my laptop)
ifeq ($(findstring greendl,$(HOSTNAME_OSX)),greendl)
GCCDIR := /opt/local/bin
ALGLIBDIR := /home/dg6/code/alglib/cpp/src
NETCDFDIR := /opt/local
CUDADIR := /usr/local/cuda
CUDALIBDIR := ${CUDADIR}/lib
CUDA_ARCH := sm_11
CUDA_SDK_DIR := /Developer/GPU\ Computing
endif

# Catch for xbmc.* (my tv)
MACHINE:=$(shell hostname)
ifeq ($(MACHINE),xbmc-desktop)
GCCDIR := /usr/bin
ALGLIBDIR := 
NETCDFDIR := /home/xbmc/code/netcdf# must be an --enable-cxx-4 dist
CUDADIR := 
CUDALIBDIR = ${CUDADIR}/lib64
CUDA_ARCH := sm_13
CUDA_SDK_DIR := 
LIBCONFIGDIR := 
endif

CUDA_SDK_INC := $(CUDA_SDK_DIR)/C/common/inc

CC := $(GCCDIR)/gcc
CPP := $(GCCDIR)/g++
NVCC := $(CUDADIR)/bin/nvcc

MODULES := src include

INCLUDEFLAGS := -I${PAPI_DIR}/include -I$(NETCDFDIR)/include -I$(GOOGLE_PERF_DIR)/include #-I$(ALGLIBDIR) -I$(CUDA_SDK_DIR)  -I$(CUDA_SDK_INC) -I$(LIBCONFIGDIR)/include
CFLAGS := 
CPPFLAGS := -g -pg -O3
NVCCFLAGS := --compiler-bindir $(GCCDIR) -arch $(CUDA_ARCH) --ptxas-options=-v #-g -G 
LFLAGS := -L${PAPI_DIR}/lib -L$(NETCDFDIR)/lib -L/home/dg6/code/google-perftools/lib #-L$(CUDALIBDIR) -L$(LIBCONFIGDIR)/lib
LIBS := -lnetcdf_c++4 -lpapi #-lnetcdf #-lprofiler #$(ALGLIBDIR)/*.o -lcuda -lcudart -lconfig++

USECUDA:=0
DEBUG:=0
PITCH_SCATTERING:=0
ENERGY_SCATTERING:=0
GOOSE:=0
SAVE_ORBITS:=0

LINK := $(CPP) ${CPPFLAGS} 

# You shouldn't have to go below here
#
# DLG: 	Added the -x c to force c file type so that 
# 		the .cu files will work too :)

DIRNAME = `dirname $1`
MAKEDEPS = $(CC) -MM -MG $2 -x c $3 | sed -e "s@^\(.*\)\.o:@.dep/$1/\1.d obj/$1/\1.o:@"

.PHONY : all

all : $(NAME)

# look for include files in each of the modules
INCLUDEFLAGS += $(patsubst %, -I%, $(MODULES))

CFLAGS += $(INCLUDEFLAGS)
CPPFLAGS += $(INCLUDEFLAGS) -DDEBUGLEVEL=$(DEBUG) \
			-DPITCH_SCATTERING=$(PITCH_SCATTERING) \
			-DENERGY_SCATTERING=$(ENERGY_SCATTERING) \
			-DGOOSE=$(GOOSE) \
			-D__SAVE_ORBITS__=$(SAVE_ORBITS)
NVCCFLAGS += $(INCLUDEFLAGS) -DDEBUGLEVEL=$(DEBUG) \
			 -DPITCH_SCATTERING=$(PITCH_SCATTERING) \
			 -DENERGY_SCATTERING=$(ENERGY_SCATTERING) \
			 -DGOOSE=$(GOOSE) \
			 -D__SAVE_ORBITS__=$(SAVE_ORBITS)

# determine the object files
SRCTYPES := c cpp 
ifeq ($(USECUDA),1)
SRCTYPES += cu
CPPFLAGS += -DUSECUDA  #-D__PROFILING__
NVCCFLAGS += #-D__PROFILING__
endif
OBJ := $(foreach srctype, $(SRCTYPES), $(patsubst %.$(srctype), obj/%.o, $(wildcard $(patsubst %, %/*.$(srctype), $(MODULES)))))

# link the program
$(NAME) : $(OBJ)
	$(LINK) $(LFLAGS) -o $@ $(OBJ) $(LIBS)

# calculate include dependencies
.dep/%.d : %.cpp
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.d$$||'`
	$(call MAKEDEPS,$(call DIRNAME, $<), $(INCLUDEFLAGS), $<) > $@

obj/%.o : %.cpp
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.o$$||'`
	$(CPP) $(CPPFLAGS) -c -o $@ $<

.dep/%.d : %.c
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.d$$||'`
	$(call MAKEDEPS,$(call DIRNAME, $<), $(CFLAGS), $<) > $@

obj/%.o : %.c
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.o$$||'`
	$(CC) $(CFLAGS) -c -o $@ $<

.dep/%.d : %.cu
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.d$$||'`
	$(call MAKEDEPS,$(call DIRNAME, $<), $(INCLUDEFLAGS), $<) > $@

obj/%.o : %.cu
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.o$$||'`
	$(NVCC) $(NVCCFLAGS) -c -o $@ $<


# include the C include dependencies
DEP := $(patsubst obj/%.o, .dep/%.d, $(OBJ))

ifneq ($(MAKECMDGOALS),clean)
include $(DEP)
endif

clean:
	-@rm $(NAME) $(OBJ) $(DEP) .dep/src/*

allclean: 
	-@rm $(NAME) $(OBJ) $(DEP) .dep/src/* webFace.wt src_webFace/*.o

webFace.wt: src_webFace/webFaceApp.o
	g++ -o $@ $< -L ~/code/wt/lib -lwt -lwthttp -lboost_signals-mt -lboost_filesystem-mt -lboost_system-mt \
			-L$(LIBCONFIGDIR)/lib -lconfig++
	@echo 'Run webApp using ...'
	@echo 'WT_TMP_DIR=/home/dg6/code/sMC/tmp ./webFace.wt --docroot ./ --http-address 0.0.0.0 --http-port 8080 -c ./wt_config.xml'

src_webFace/webFaceApp.o: src_webFace/webFaceApp.cpp
	g++ -c $< -I ~/code/wt/include -o $@ -I$(LIBCONFIGDIR)/include

