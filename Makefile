SOURCE_DIR=R2SSUR/src
VPATH=$(SOURCE_DIR)

CC=g++
CFLAGS= -c -Wall -Wno-reorder -std=c++11 -I$(SOURCE_DIR)/ -DCCODE -fopenmp

OPENLDFLAGS= -larmadillo -lpthread -lopenblas -fopenmp
NVLDFLAGS= -larmadillo -lpthread -lnvblas -fopenmp

SOURCES_BVS=$(SOURCE_DIR)/global.cpp $(SOURCE_DIR)/utils.cpp $(SOURCE_DIR)/distr.cpp $(SOURCE_DIR)/junction_tree.cpp $(SOURCE_DIR)/HESS_Chain.cpp $(SOURCE_DIR)/SUR_Chain.cpp $(SOURCE_DIR)/drive.cpp main.cpp 
#ESS_Atom.h and Parameters_type.h are interface only
OBJECTS_BVS=$(SOURCES_BVS:.cpp=.o)

SOURCES_XML=$(SOURCE_DIR)/pugixml.cpp
OBJECTS_XML=$(SOURCES_XML:.cpp=.o)

all:$(SOURCES_XML) $(SOURCES_BVS) BVS_NONVIDIA

BVS_NVIDIA: OPTIM_FLAGS := -O3
BVS_NVIDIA: $(OBJECTS_BVS)
	@echo [Linking and producing executable]:
	if [ -e /usr/lib/libnvblas.so -o -e /usr/lib64/libnvblas.so ] ; then $(CC) $(OBJECTS_BVS) -o BVS_Reg $(NVLDFLAGS) ; else $(CC) $(OBJECTS_BVS) -o BVS_Reg $(OPENLDFLAGS) ; fi
# I hate that this is searching only in the standard path, the correct way is with automake ./config, but that will come later...
# if you compile against nvblas run with LD_PRELOAD=/usr/lib64/libnvblas.so or similar path to libnvblas.so

BVS_NONVIDIA: OPTIM_FLAGS := -O3
BVS_NONVIDIA: $(OBJECTS_XML) $(OBJECTS_BVS)
	@echo [Linking and producing executable]:
	$(CC) $(OBJECTS_XML) $(OBJECTS_BVS) -o BVS_Reg $(OPENLDFLAGS)

BVS_DEBUG: OPTIM_FLAGS := -O0 #maybe O1 ? -pg ? linker as well?
BVS_DEBUG: $(OBJECTS_XML) $(OBJECTS_BVS)
	@echo [Linking and producing executable]:
	$(CC) $(OBJECTS_XML) $(OBJECTS_BVS) -o BVS_DEBUG_Reg $(OPENLDFLAGS) -ggdb3 -g -lprofiler 

%.o: %.cpp
	@echo [Compiling]: $<
	$(CC) $(CFLAGS) $(OPTIM_FLAGS) -o $@ -c $<	

clean: 
	@echo [Cleaning: ]
	rm *.o; rm $(SOURCE_DIR)/*.o ; rm *_Reg ; rm data* ; rm blocks* ; rm structureGraph* ; rm call.sh ; rm -R results ;

remake: 
	@echo [Cleaning compilation objets only: ]
	rm *.o; rm $(SOURCE_DIR)/*.o ; rm *_Reg ;
