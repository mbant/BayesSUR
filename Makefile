SOURCE_DIR=R2SSUR/src
VPATH=$(SOURCE_DIR)

CC=g++
CFLAGS= -c -Wall -Wno-reorder -std=c++11 -fopenmp -I$(SOURCE_DIR)/ -DCCODE

OPENLDFLAGS= -larmadillo -lpthread -lopenblas -fopenmp
NVLDFLAGS= -larmadillo -lpthread -lnvblas -fopenmp

SOURCES_BVS=$(SOURCE_DIR)/global.cpp $(SOURCE_DIR)/utils.cpp $(SOURCE_DIR)/distr.cpp $(SOURCE_DIR)/junction_tree.cpp $(SOURCE_DIR)/HESS_Chain.cpp $(SOURCE_DIR)/dSUR_Chain.cpp $(SOURCE_DIR)/SSUR_Chain.cpp $(SOURCE_DIR)/drive.cpp $(SOURCE_DIR)/main.cpp 
#ESS_Atom.h interface only
OBJECTS_BVS=$(SOURCES_BVS:.cpp=.o)

all:$(SOURCES_BVS) BVS_NONVIDIA

BVS_NVIDIA: OPTIM_FLAGS := -O3
BVS_NVIDIA: $(OBJECTS_BVS)
	@echo [Linking and producing executable]:
	if [ -e /usr/lib/libnvblas.so -o -e /usr/lib64/libnvblas.so ] ; then $(CC) $(OBJECTS_BVS) -o BVS_Reg $(NVLDFLAGS) ; else $(CC) $(OBJECTS_BVS) -o BVS_Reg $(OPENLDFLAGS) ; fi
# I hate that this is searching only in the standard path, the correct way is with automake ./config, but that will come later...

BVS_NONVIDIA: OPTIM_FLAGS := -O3
BVS_NONVIDIA: $(OBJECTS_BVS)
	@echo [Linking and producing executable]:
	$(CC) $(OBJECTS_BVS) -o BVS_Reg $(OPENLDFLAGS)

BVS_DEBUG: $(OBJECTS_BVS)
	@echo [Linking and producing executable]:
	$(CC) $(OBJECTS_BVS) -o BVS_DEBUG_Reg $(LDFLAGS) -ggdb3

%.o: %.cpp
	@echo [Compiling]: $<
	$(CC) $(CFLAGS) $(OPTIM_FLAGS) -o $@ -c $<	

clean: 
	@echo [Cleaning: ]
	rm $(SOURCE_DIR)/*.o ; rm *_Reg ; rm simul* ; rm blocks* ; rm structureGraph* ; rm call.sh ; rm -R results ;

remake: 
	@echo [Cleaning compilation objets only: ]
	rm $(SOURCE_DIR)/*.o ; rm *_Reg ;
