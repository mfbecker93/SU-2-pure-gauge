CXX = mpicxx

GRID = $(HOME)/local/grid-su2

TARGET = build/pure_gauge

SRC = src/main.cc

CXXFLAGS = $(shell $(GRID)/bin/grid-config --cxxflags) \
           -I/usr/include/hdf5/serial

LDFLAGS = -L$(GRID)/lib \
          -L$(HOME)/local/lib \
          -L/usr/lib/x86_64-linux-gnu/hdf5/serial \
          -L/usr/lib/llvm-18/lib

LIBS = $(shell $(GRID)/bin/grid-config --libs) \
       -lunwind-x86_64 \
       -lunwind


all:
	$(CXX) $(SRC) $(CXXFLAGS) $(LDFLAGS) $(LIBS) -o $(TARGET)


clean:
	rm -f $(TARGET)
