CC=g++
CFLAGS= -c -Wall -std=c++11 -fopenmp
LDFLAGS= -larmadillo -lpthread -llapack -lopenblas -fopenmp
#LDFLAGS= -larmadillo -lpthread -llapack -lnvblas -fopenmp

SOURCES_HESS=global.cpp utils.cpp distr.cpp HESS.cpp drive_HESS.cpp
OBJECTS_HESS=$(SOURCES_HESS:.cpp=.o)

SOURCES_SSUR=global.cpp utils.cpp distr.cpp junction_tree.cpp SSUR_Chain.cpp drive_SSUR.cpp
OBJECTS_SSUR=$(SOURCES_SSUR:.cpp=.o)

SOURCES_SSUR_DEBUG=global.cpp utils.cpp distr.cpp junction_tree.cpp SSUR_Chain.cpp drive_SSUR.cpp
OBJECTS_SSUR_DEBUG=$(SOURCES_SSUR_DEBUG:.cpp=.o)

all: $(SOURCES_HESS) $(SOURCES_SSUR) HESS SSUR

$(OBJECTS_HESS): EXTRA_FLAGS := -O3

$(OBJECTS_SSUR): EXTRA_FLAGS := -O3

$(OBJECTS_SSUR_DEBUG): EXTRA_FLAGS :=

HESS: $(OBJECTS_HESS)
	@echo [Linking and producing executable]:
	$(CC) $(OBJECTS_HESS) -o HESS_Reg $(LDFLAGS)

SSUR: $(OBJECTS_SSUR)
	@echo [Linking and producing executable]:
	$(CC) $(OBJECTS_SSUR) -o SSUR_Reg $(LDFLAGS)
	
SSUR_DEBUG: $(OBJECTS_SSUR_DEBUG)
	@echo [Linking and producing executable]:
	$(CC) $(OBJECTS_SSUR_DEBUG) -o SSUR_Reg $(LDFLAGS) -ggdb3

%.o: %.cpp
	@echo [Compiling]: $<
	$(CC) $(CFLAGS) $(EXTRA_FLAGS) -o $@ -c $<	
	
clean:
	rm *.o
	rm *_Reg
