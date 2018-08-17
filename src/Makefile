CC=g++
CFLAGS= -c -Wall -std=c++11 -fopenmp
LDFLAGS= -larmadillo -lpthread -llapack -lopenblas -fopenmp
#LDFLAGS= -larmadillo -lpthread -llapack -lnvblas -fopenmp

SOURCES_HESS=global.cpp utils.cpp distr.cpp HESS.cpp drive_HESS.cpp
OBJECTS_HESS=$(SOURCES_HESS:.cpp=.o)

SOURCES_SSUR=global.cpp utils.cpp distr.cpp junction_tree.cpp SSUR_Chain.cpp ESS_Sampler.cpp drive_SSUR.cpp
OBJECTS_SSUR=$(SOURCES_SSUR:.cpp=.o)

all: $(SOURCES_HESS) $(SOURCES_SSUR) HESS SSUR

HESS: OPTIM_FLAGS := -O3
HESS: $(OBJECTS_HESS)
	@echo [Linking and producing executable]:
	$(CC) $(OBJECTS_HESS) -o HESS_Reg $(LDFLAGS)

SSUR: OPTIM_FLAGS := -O3
SSUR: $(OBJECTS_SSUR)
	@echo [Linking and producing executable]:
	$(CC) $(OBJECTS_SSUR) -o SSUR_Reg $(LDFLAGS)

%.o: %.cpp
	@echo [Compiling]: $<
	$(CC) $(CFLAGS) $(OPTIM_FLAGS) -o $@ -c $<

clean:
	rm *.o
	rm *_Reg
