EXEC   = cholla

OPTIMIZE =  -O2  

DIR = ./src
CFILES = $(wildcard $(DIR)/*.c)
CPPFILES = $(wildcard $(DIR)/*.cpp)
CUDAFILES = $(wildcard $(DIR)/*.cu)

OBJS = $(subst .c,.o,$(CFILES)) $(subst .cpp,.o,$(CPPFILES)) $(subst .cu,.o,$(CUDAFILES))


#To use GPUs, CUDA must be turned on here
#Optional error checking can also be enabled
CUDA = -DCUDA -DCUDA_ERROR_CHECK


#To use MPI, MPI_FLAGS must be set to -DMPI_CHOLLA
#otherwise gcc/g++ will be used for serial compilation
MPI_FLAGS =  -DMPI_CHOLLA

ifdef MPI_FLAGS
  CC	= mpicc
  CXX = mpicxx

  #MPI_FLAGS += -DSLAB
  MPI_FLAGS += -DBLOCK

else
  CC  = cc
  CXX = c++
endif


#define the NVIDIA CUDA compiler
NVCC	= nvcc

.SUFFIXES : .c .cpp .cu .o

#PRECISION = -DPRECISION=1
PRECISION = -DPRECISION=2

#OUTPUT = -DBINARY
OUTPUT = -DHDF5

#RECONSTRUCTION = -DPCM
#RECONSTRUCTION = -DPLMP
#RECONSTRUCTION = -DPLMC
#RECONSTRUCTION = -DPPMP
RECONSTRUCTION = -DPPMC

SOLVER = -DEXACT
#SOLVER = -DROE
#SOLVER = -DHLLC

#INTEGRATOR = -DCTU 
INTEGRATOR = -DVL

#COOLING = -DCOOLING_CPU
#COOLING = -DCOOLING_GPU

INCL   = -I./ -I/usr/local/cuda/include -I/usr/local/include
NVLIBS = -L/usr/local/cuda/lib -lcuda -lcudart
LIBS   = -L/usr/local/lib -lm -lhdf5 #-lgsl


FLAGS = $(CUDA) $(PRECISION) $(OUTPUT) $(RECONSTRUCTION) $(SOLVER) $(INTEGRATOR) $(COOLING)
CFLAGS 	  = $(OPTIMIZE) $(FLAGS) $(MPI_FLAGS)
CXXFLAGS  = $(OPTIMIZE) $(FLAGS) $(MPI_FLAGS)
NVCCFLAGS = $(FLAGS) -fmad=false -m64 -Xcompiler -arch -Xcompiler x86_64 -gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30 -Wno-deprecated-gpu-targets
ALL_LDFLAGS	 += -Xlinker -F/Library/Frameworks -Xlinker -framework -Xlinker CUDA


%.o:	%.c
		$(CC) $(CFLAGS)  $(INCL)  -c $< -o $@ 

%.o:	%.cpp
		$(CXX) $(CXXFLAGS)  $(INCL)  -c $< -o $@ 

%.o:	%.cu
		$(NVCC) $(NVCCFLAGS)  $(INCL)  -c $< -o $@ 

$(EXEC): $(OBJS) 
	 	 $(CXX) $(OBJS) $(LIBS) $(NVLIBS) -o $(EXEC) $(INCL) $(ALL_LDFLAGS)

#$(OBJS): $(INCL) 

.PHONY : clean

clean:
	 rm -f $(OBJS) $(EXEC)

