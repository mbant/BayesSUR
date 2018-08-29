CC=g++
CFLAGS= -c -Wall -std=c++11 -fopenmp
LDFLAGS= -larmadillo -lpthread -llapack -lopenblas -fopenmp
#LDFLAGS= -larmadillo -lpthread -llapack -lnvblas -fopenmp

SOURCES_BVS=global.cpp utils.cpp distr.cpp junction_tree.cpp HESS_Chain.cpp SSUR_Chain.cpp ESS_Sampler.h drive.cpp 
#ESS_Atom.h interface only
OBJECTS_BVS=$(SOURCES_BVS:.cpp=.o)

all:$(SOURCES_BVS) BVS

BVS: OPTIM_FLAGS := -O3
BVS: $(OBJECTS_BVS)
	@echo [Linking and producing executable]:
	$(CC) $(OBJECTS_BVS) -o BVS_Reg $(LDFLAGS)
	
BVS_DEBUG: $(OBJECTS_BVS)
	@echo [Linking and producing executable]:
	$(CC) $(OBJECTS_BVS) -o BVS_DEBUG_Reg $(LDFLAGS) -ggdb3

%.o: %.cpp
	@echo [Compiling]: $<
	$(CC) $(CFLAGS) $(OPTIM_FLAGS) -o $@ -c $<	
	

clean: 
	@echo [Cleaning: ]
	rm *.o 
	rm *_Reg 
