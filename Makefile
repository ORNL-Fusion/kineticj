.SUFFIXES:
.SUFFIXES: .c .cpp .cu

CPPFLAGS :=

# Select only one.

USE_SERIAL:= 0
USE_CUDA  := 0
USE_OPENMP:= 1

NAME := bin/kineticj

LIBS :=  
INCLUDEFLAGS :=  

GOOGLE_PERF_DIR := ${HOME}/code/google-perftools
PAPI_DIR := ${HOME}/code/papi/gnu_${GNUVER}

CC := gcc
CPP := g++
GCC := gcc
VENDOR := GCC_

ThisMachine := $(shell uname -n)

ifneq (,$(findstring titan,$(ThisMachine)))
ThisMachine := titan
endif

ifneq (,$(findstring edison,$(ThisMachine)))
ThisMachine := edison 
endif

include machine-makefiles/Makefile.$(ThisMachine)
include Makefile.flags

ifeq ($(USE_SERIAL),1)
THRUST_POLICY:=THRUST_DEVICE_SYSTEM_CPP
NVCC := ${CPP} -O3 -DTHRUST_DEVICE_SYSTEM=${THRUST_POLICY} -I ${CUDA_INCDIR}
NVCCFLAGS:= -x c++ 
CPPFLAGS+= -D__THRUST
endif

ifeq ($(USE_CUDA),1)
THRUST_POLICY:=THRUST_DEVICE_SYSTEM_CUDA
NVCC := nvcc -O3 -DTHRUST_DEVICE_SYSTEM=${THRUST_POLICY} 
NVCCFLAGS:= -dc --expt-relaxed-constexpr -Xcompiler 
endif


ifeq ($(USE_OPENMP),1)
THRUST_POLICY:=THRUST_DEVICE_SYSTEM_OMP
NVCC := ${CPP} -O3  -DTHRUST_DEVICE_SYSTEM=${THRUST_POLICY} -I ${CUDA_INCDIR}
NVCCFLAGS := -x c++ -fopenmp
LFLAGS += -fopenmp -lgomp
CPPFLAGS+= -D__THRUST
endif


MODULES := src include

OPENMPFLAGS := $($(VENDOR)OPENMPFLAGS)
DEBUGFLAGS := $($(VENDOR)DEBUGFLAGS)
OPTFLAGS := $($(VENDOR)OPTFLAGS)

CFLAGS := 

CXXFLAGS := ${OPENMPFLAGS} ${DEBUGFLAGS} ${OPTFLAGS} 

CPPFLAGS += -DDEBUGLEVEL=0
CPPFLAGS += -DDEBUG_LINES=0
CPPFLAGS += -DDEBUG_INTERP=0
CPPFLAGS += -DUSEPAPI=0
CPPFLAGS += -DLOWMEM_USEPAPI=0
CPPFLAGS += -D_PARTICLE_BOUNDARY=1 # 1 = particle absorbing walls, 2 = periodic, 3 = reflective
CPPFLAGS += -DCOMPLEX_WRF=0
CPPFLAGS += -DDEBUG_GC=0
CPPFLAGS += -DDEBUG_EVAL_VGC=0
CPPFLAGS += -DDEBUG_EVAL_APAR=0
CPPFLAGS += -DCLOCK=1
CPPFLAGS += -DPRINT_INFO=1
CPPFLAGS += -DGC_ORBITS=0
CPPFLAGS += -DDEBUG_MAXWELLIAN=0
CPPFLAGS += -DDEBUG_FORCE_TERM=0
CPPFLAGS += -DDEBUG_MOVE=0
CPPFLAGS += -DDEBUG_ROTATION=0
CPPFLAGS += -std=c++11
CPPFLAGS += -DDO_CPU_ITERATOR_APPROACH=0
CPPFLAGS += -DDO_CPU_APPROACH=0
CPPFLAGS += -DLOWMEM_ORBIT_WRITE=0
CPPFLAGS += -DF1_WRITE=0
CPPFLAGS += -DDEBUG_INTVECARRAY=0
CPPFLAGS += -DDEBUG_READ_E_FIELD=0
CPPFLAGS += -DCYLINDRICAL_INPUT_FIELDS=0


# You shouldn't have to go below here
#
# DLG: 	Added the -x c to force c file type so that 
# 		the .cu files will work too :)

DIRNAME = `dirname $1`
MAKEDEPS = ${GCC} -MM -MG $2 -x c $3 | sed -e "s@^\(.*\)\.o:@.dep/$1/\1.d obj/$1/\1.o:@"

.PHONY : all

all : $(NAME)

# look for include files in each of the modules
INCLUDEFLAGS += $(patsubst %, -I%, $(MODULES))

CFLAGS += $(INCLUDEFLAGS)
CXXFLAGS += $(INCLUDEFLAGS) 
NVCCFLAGS += $(INCLUDEFLAGS) 

# determine the object files
SRCTYPES := c cpp cu 
LINK := $(CPP) $(CXXFLAGS) $(LFLAGS) 

LINK := $(NVCC) $(LFLAGS)  

OBJ := $(foreach srctype, $(SRCTYPES), $(patsubst %.$(srctype), obj/%.o, $(wildcard $(patsubst %, %/*.$(srctype), $(MODULES)))))

# link the program
$(NAME) : $(OBJ)
	$(LINK) -o $@ $(OBJ) $(LIBS)

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

