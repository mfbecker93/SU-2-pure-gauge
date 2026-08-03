# TODO: for HPC replace with mpicxx
CXX = clang++

GRID = $(HOME)/PhD/Grid/grid_install_SU2

TARGET = build/pure_gauge

SRC = src/main.cc

CXXFLAGS = -std=c++17 \
           $(shell $(GRID)/bin/grid-config --cxxflags) \
           -I/opt/homebrew/opt/hdf5/include \
           -I/opt/homebrew/opt/openssl@3/include\
		   -I/opt/homebrew/opt/fftw/include 


LDFLAGS = -L$(GRID)/lib \
		  -L$(HOME)/PhD/Grid/c-lime/lib \
          -L/opt/homebrew/opt/hdf5/lib \
          -L/opt/homebrew/opt/openssl@3/lib\
		  -L/opt/homebrew/opt/fftw/lib

LIBS = $(shell $(GRID)/bin/grid-config --libs)

all:
	mkdir -p build
	$(CXX) $(SRC) $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $(TARGET)

clean:
	rm -f $(TARGET)