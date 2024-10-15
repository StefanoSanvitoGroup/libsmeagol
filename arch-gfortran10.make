# Fortran compiler
FC                = mpif90

# Utility to create a static library
AR                = ar -r

# Fortran compiler's flags.
# The source code is not Fortran 95/2003/2008 standard compliant due to DOUBLE COMPLEX variables.
# Therefore GNU extensions need to be enabled.
FCFLAGS           = -DMPI -fopenmp -march=native -O3 -g -std=gnu -fallow-argument-mismatch \
                    -fexternal-blas -fblas-matmul-limit=0 -fno-omit-frame-pointer -funroll-loops

# Flags to specify the layout used by source files.
#   ifort determines a source file's layout based on the filename's suffix.
#   gfortran expects Fortran 95 source files to have free-form layout, which is not always the case.
FCFLAGS_FIXEDFORM = -ffixed-form
FCFLAGS_FREEFORM  = -ffree-form -ffree-line-length-none
