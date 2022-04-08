ENABLE_OMP_OFFLOAD=0
CXX = clang++
ifeq ($(ENABLE_OMP_OFFLOAD),1)
#CXX = CC
CXX = clang++
#CXX = nvc++
#CXX = xlc++
#SM=80
endif

OPTFLAGS = -O3 -DPRINT_EXTRA_NEDGES -std=c++11 #-DCHECK_RESULTS #-DPRINT_RESULTS
CXXFLAGS = -g -I. $(OPTFLAGS)
ifeq ($(ENABLE_OMP_OFFLOAD),1)
#CXXFLAGS_THREADS += -mp=gpu -Minfo=mp -DUSE_OMP_OFFLOAD
#CXXFLAGS_THREADS += -fopenmp-targets=nvptx64 -Xopenmp-target=nvptx64 -march=sm_${SM} -DUSE_OMP_OFFLOAD
#CXXFLAGS_THREADS += -qsmp=omp -qoffload -DUSE_OMP_OFFLOAD
CXXFLAGS_THREADS +=  -DGRAPH_FT_LOAD=4 -DUSE_OMP_OFFLOAD -fopenmp -I${ROCM_PATH}/include -fopenmp-assume-no-thread-state -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn-amd-amdhsa -march=gfx90a
#CXXFLAGS_THREADS += -mp=gpu -Minfo=mp -DUSE_OMP_ACCELERATOR
#CXXFLAGS_THREADS += -fopenmp-targets=nvptx64 -Xopenmp-target=nvptx64 -march=sm_${SM} -DUSE_OMP_ACCELERATOR
#PE_MPICH_GTL_DIR_amd_gfx90a="-L${CRAY_MPICH_ROOTDIR}/gtl/lib"
#PE_MPICH_GTL_LIBS_amd_gfx90a="-lmpi_gtl_hsa"
LDFLAGS = -L${ROCM_PATH}/lib -lamdhip64
else
CXXFLAGS_THREADS = -fopenmp -DGRAPH_FT_LOAD=4 #-I/usr/lib/gcc/x86_64-redhat-linux/4.8.5/include/
endif
OBJ = main.o
SRC = main.cpp
TARGET = matching_omp

OBJS = $(OBJ) 
TARGETS = $(TARGET)

all: $(TARGETS)

$(TARGET):  $(OBJ)
	$(LDAPP) $(CXX) $(CXXFLAGS_THREADS) -o $@ $+ $(LDFLAGS) $(CXXFLAGS) 

$(OBJ): $(SRC)
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(CXXFLAGS_THREADS) -c $< -o $@

.PHONY: clean 

clean:
	rm -rf *~ *.dSYM nc.vg.* $(OBJS) $(TARGETS)
