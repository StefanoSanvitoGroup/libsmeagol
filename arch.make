# Fortran compiler
FC                = mpif90

# Utility to create a static library
AR                = ar -r

# Fortran compiler's flags
FCFLAGS           = -DMPI -qopenmp -xHost -O2 -g -fno-omit-frame-pointer

# Flags to specify the layout used by source files
FCFLAGS_FIXEDFORM = -fixed
FCFLAGS_FREEFORM  = -free

