CC=gcc-7
CXX=g++

OBJDIR=./obj
SRCDIR=./src
INCDIR =  -I~/usr/include -I ./include -I/usr/local/include -I$(SRCDIR)
LIBDIR  =  -L~/usr/lib  -L/usr/local/lib 

DEBUG=
OPTIM=-O3
#BLASOARMA=
BLASOARMA=-DARMA_DONT_USE_WRAPPER -llapack -larmadillo -lopenblas -fopenmp
CFLAGS= $(DEBUG) $(OPTIM) -std=c++11 -fpermissive -Wall -w $(BLASOARMA)

LIBS= $(BLASOARMA) -lm -larmadillo -w 

_DEPS = Solver.hpp DataProcess.hpp Config.hpp Projection.hpp Gibbs.hpp Alloy.hpp Connectivity.hpp ChronoP.hpp KMax.hpp GaussianCompare.cpp
DEPS = $(patsubst %,$(SRCDIR)/%,$(_DEPS))

_OBJ = Solver.o DataProcess.o Config.o Projection.o Gibbs.o Alloy.o Connectivity.o ChronoP.o KMax.o GaussianCompare.o
OBJ = $(patsubst %,$(OBJDIR)/%,$(_OBJ))

NAME=solver

all: main

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CFLAGS)

main: $(OBJ) $(OBJDIR)/main.o
	$(CXX) -o $(NAME) $(OBJ) $(OBJDIR)/main.o $(CFLAGS) $(LIBS)

run:
	./solver

.PHONY: clean

clean:
	rm -f $(OBJDIR)/*.o *~ $(NAME)
