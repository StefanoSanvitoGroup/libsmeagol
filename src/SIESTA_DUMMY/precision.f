      MODULE precision
C
C  One can rely on REAL32 and REAL64 constants from ISO_FORTRAN_ENV module
C  introduced in Fortran 2008 standard.
C
      IMPLICIT NONE
      PRIVATE

      integer, parameter, public :: sp = selected_real_kind(6,30)
      integer, parameter, public :: dp = selected_real_kind(14,100)

      END MODULE precision
