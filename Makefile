.SUFFIXES:
.SUFFIXES: .c .cpp .cu

NAME := bin/kineticj

# Defaults to dlg-hp.ornl.gov

GCCDIR :=  
ALGLIBDIR := ${HOME}/code/alglib/cpp/src
NETCDFDIR := ${HOME}/code/netcdf/gnu_4.7.0# must be an --enable-cxx-4 dist
CUDADIR := ${HOME}/code/cuda/4.1/cuda
CUDALIBDIR = ${CUDADIR}/lib64
CUDA_ARCH := sm_13
CUDA_SDK_DIR := ${HOME}/cuda/NVIDIA_GPU_Computing_SDK
LIBCONFIGDIR := ${HOME}/code/libconfig
GOOGLE_PERF_DIR := ${HOME}/code/google-perftools
PAPI_DIR := ${HOME}/code/papi/gnu_${GNUVER}

# Catch for greendl.* (my laptop)
ifeq ($(findstring greendl,$(HOSTNAME_OSX)),greendl)
GCCDIR := /opt/local/bin
ALGLIBDIR := ${HOME}/code/alglib/cpp/src
NETCDFDIR := ${HOME}/code/netcdf/gnu_4.3.6
LIBCONFIGDIR := ${HOME}/code/libconfig/gnu_4.3.6
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

CC := gcc
CPP := g++

NVCC := $(CUDADIR)/bin/nvcc

MODULES := src include

INCLUDEFLAGS := -I$(LIBCONFIGDIR)/include  \
		-I$(NETCDFDIR)/include #-I$(GOOGLE_PERF_DIR)/include -I${PAPI_DIR}/include
OPENMPFLAGS := #-fopenmp
DEBUGFLAGS := #-g -pg
OPTFLAGS := -O3
CFLAGS := 
CXXFLAGS := ${OPENMPFLAGS} ${DEBUGFLAGS} ${OPTGLAGS} 
CPPFLAGS :=
#NVCCFLAGS := --compiler-bindir $(GCCDIR) -arch $(CUDA_ARCH) --ptxas-options=-v #-g -G 
NVCCFLAGS := -arch $(CUDA_ARCH) --ptxas-options=-v #-g -G 
LFLAGS := -L$(NETCDFDIR)/lib -L$(LIBCONFIGDIR)/lib #-L${PAPI_DIR}/lib -L${HOME}/code/google-perftools/lib #-L$(CUDALIBDIR) 
LIBS := -lnetcdf_c++4 -lconfig++ #-lpapi -lnetcdf #-lprofiler #$(ALGLIBDIR)/*.o -lcuda -lcudart

CPPFLAGS += -DDEBUGLEVEL=0
CPPFLAGS += -DUSEPAPI=0
CPPFLAGS += -D__SAVE_ORBITS__=1
CPPFLAGS += -DLOWMEM=0
CPPFLAGS += -D_PARTICLE_BOUNDARY=3 # 1 = particle absorbing walls, 2 = periodic, 3 = reflective

LINK := $(CPP) ${CXXFLAGS} 

# You shouldn't have to go below here
#
# DLG: 	Added the -x c to force c file type so that 
# 		the .cu files will work too :)

DIRNAME = `dirname $1`
#MAKEDEPS = $(GCCDIR)/gcc -MM -MG $2 -x c $3 | sed -e "s@^\(.*\)\.o:@.dep/$1/\1.d obj/$1/\1.o:@"
MAKEDEPS = ${CC} -MM -MG $2 -x c $3 | sed -e "s@^\(.*\)\.o:@.dep/$1/\1.d obj/$1/\1.o:@"

.PHONY : all

all : $(NAME)

# look for include files in each of the modules
INCLUDEFLAGS += $(patsubst %, -I%, $(MODULES))

CFLAGS += $(INCLUDEFLAGS)
CXXFLAGS += $(INCLUDEFLAGS) 
NVCCFLAGS += $(INCLUDEFLAGS) 

# determine the object files
SRCTYPES := c cpp 
ifeq ($(USECUDA),1)
SRCTYPES += cu
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
	$(CPP) $(CXXFLAGS) ${CPPFLAGS} -c -o $@ $<

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
	$(NVCC) $(NVCCFLAGS) ${CPPFLAGS} -c -o $@ $<


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

