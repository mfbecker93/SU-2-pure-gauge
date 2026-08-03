# README
## 1. General Information
This repository contains the code for a SU(2) pure gauge theory on the lattice realised with Grid v0.7.0. Grid itself is not included in this repository. For Grid related information refer to section 2. From the Grid log files (not inlcuded in this repo) observables have been extracted und analysed. The details of the analysis can be found in section 3 (TODO).

## 2. Grid

### Grid Instalation
Grid was installed following the instructions on: https://github.com/telos-collaboration/Grid 

Grid commit: 2279f4c2

### Grid configuration

Grid was built with:
- Nc = 2
- SIMD = AVX2
- Threading enabled
- MPI3 communication
- FFTW support enabled
- HDF5 support enabled
- LIME support enabled

To achieve this, use following configure command after adjusting the paths to your system's paths and adapting -march=native to your architecture since this flag was used for local optimisation.

../configure \
    --prefix=/home/max/local/grid-su2 \
    --enable-Nc=2 \
    --enable-simd=AVX2 \
    --enable-comms=mpi-auto \
    --with-lime=/home/max/local \
    --with-hdf5=/usr \
    CC=clang \
    CXX=clang++ \
    'CPPFLAGS=-I/home/max/local/include -I/usr/include/hdf5/serial' \
    'LDFLAGS=-L/home/max/local/lib -L/usr/lib/x86_64-linux-gnu/hdf5/serial -L/usr/lib/llvm-18/lib' \
    'CXXFLAGS=-O3 -march=native -fopenmp' \
    LIBS=-lomp


This configuration leads to the following specifications:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Summary of configuration for Grid v0.7.0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
----- GIT VERSION -------------------------------------\
commit: 2279f4c2\
branch: develop\
date  : 2026-06-15\
----- PLATFORM ----------------------------------------\
architecture (build)        : x86_64\
os (build)                  : linux-gnu\
architecture (target)       : x86_64\
os (target)                 : linux-gnu\
compiler vendor             : clang\
compiler version            : 18.1.3\
----- BUILD OPTIONS -----------------------------------\
Nc                          : 2\
SIMD                        : AVX2\
Threading                   : yes\
Acceleration                : none\
Unified virtual memory      : yes\
Communications type         : mpi3\
Shared memory allocator     : no\
Shared memory mmap path     : /var/lib/hugetlbfs/global/pagesize-2MB/\
Default precision           :\
Software FP16 conversion    : yes\
RNG choice                  : sitmo\
GMP                         : no\
LAPACK                      : no\
FFTW                        : yes\
LIME (ILDG support)         : yes\
HDF5                        : yes\
build DOXYGEN documentation : no\
Sp2n                        : no

### Dependencies

Required software:
- Grid v0.7.0
- C++ compiler (clang 18.1.3 used)
- MPI (OpenMPI, MPI3)
- FFTW3
- HDF5
- LIME (ILDG support)
- OpenMP
- libunwind

## 3. Analysis

TODO