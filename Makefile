.SUFFIXES:
.SUFFIXES: .c .cpp .cu

USECUDA:=1

NAME := bin/kineticj

LIBS :=  
INCLUDEFLAGS :=  

GOOGLE_PERF_DIR := ${HOME}/code/google-perftools
PAPI_DIR := ${HOME}/code/papi/gnu_${GNUVER}

GCCDIR :=  
PGIDIR := 

VENDOR := GCC_
CC := gcc
CPP := g++
#NVCC := nvcc -g -G
NVCC := nvcc -O3 


#VENDOR := PGI_
#CC := pgcc
#CPP := pgcpp

ThisMachine := $(shell uname -n)

ifneq (,$(findstring titan,$(ThisMachine)))
ThisMachine := titan
endif

include machine-makefiles/Makefile.$(ThisMachine)
include Makefile.flags

MODULES := src include

OPENMPFLAGS := $($(VENDOR)OPENMPFLAGS)
DEBUGFLAGS := $($(VENDOR)DEBUGFLAGS)
OPTFLAGS := $($(VENDOR)OPTFLAGS)

CFLAGS := 

CXXFLAGS := ${OPENMPFLAGS} ${DEBUGFLAGS} ${OPTFLAGS} 

CPPFLAGS :=
CPPFLAGS += -DDEBUGLEVEL=0
CPPFLAGS += -DDEBUG_LINES=0
CPPFLAGS += -DDEBUG_INTERP=0
CPPFLAGS += -DUSEPAPI=0
CPPFLAGS += -DLOWMEM_USEPAPI=0
CPPFLAGS += -D_PARTICLE_BOUNDARY=2 # 1 = particle absorbing walls, 2 = periodic, 3 = reflective
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
CPPFLAGS += -DLOWMEM_ORBIT_WRITE=0
CPPFLAGS += -DDEBUG_ROTATION=0
CPPFLAGS += -DF1_WRITE=0
CPPFLAGS += -std=c++11
CPPFLAGS += -DDO_CPU_ITERATOR_APPROACH=0
CPPFLAGS += -DDO_CPU_APPROACH=0
CPPFLAGS += -DDEBUG_INTVECARRAY=0
CPPFLAGS += -DDEBUG_READ_E_FIELD=0
CPPFLAGS += -DCYLINDRICAL_INPUT_FIELDS=0

NVCCFLAGS := -dc --expt-relaxed-constexpr -Xcompiler -fopenmp

# You shouldn't have to go below here
#
# DLG: 	Added the -x c to force c file type so that 
# 		the .cu files will work too :)

DIRNAME = `dirname $1`
MAKEDEPS = gcc -MM -MG $2 -x c $3 | sed -e "s@^\(.*\)\.o:@.dep/$1/\1.d obj/$1/\1.o:@"
#MAKEDEPS = ${CC} -MM -MG $2 -x c $3 | sed -e "s@^\(.*\)\.o:@.dep/$1/\1.d obj/$1/\1.o:@"

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

ifeq ($(USECUDA),1)
LINK := $(NVCC) $(LFLAGS) -lcuda -Xcompiler -fopenmp 
else
NVCCFLAGS += --x c++
endif

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

