
# SMEAGOL library version 1.2

This library is based on unaltered SMEAGOL/SIESTA v1.2 source code.

The file `src/NEGF/Develop/ScissorOperator/scissor.F90` was replaced with
its dummy version, as the original one needs a SIESTA-specific library
(fdf) to parse input files.

The directory `src/SIESTA_DUMMY` contains source files which offer
functionality originally provided by SIESTA.

There are patch files in `patches/` directory which can be (partially)
applied to the original SMEAGOL source code in any order.

To compile the library

  0. (optional step) patch the source code

     patch -p1 < patches/file.patch

  1. provide the name of a Fortran compiler and build flags via arch.make file.
     A typical arch.make file for Classic Intel Fortran compiler is the following

```
# Fortran compiler
FC                = mpif90
# Utility to create a static library
AR                = ar -r
# Fortran compiler's flags
FCFLAGS           = -DMPI -qopenmp -xHost -O2 -g
# Flags to specify the layout used by source files
FCFLAGS_FIXEDFORM = -fixed
FCFLAGS_FREEFORM  = -free
```

     Although it is optional, produce debugging information flag (-g)
     allows CP2K to print accurate backtrace of the entire stack if
     it happens to crash within a SMEAGOL library's routine.
     
  2. build the library using GNU make. Parallel build is supported.
```
     make -j 8
```
Then Fortran module files and static library can be found in `obj/` and `lib/`
directories respectively. The command
```
     make clean
```
simply removes these `obj/` and `lib/` directories.

