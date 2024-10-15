      SUBROUTINE f77flush
C
C  This version of the subroutine requires a Fortran 2003 compiler
C
      USE ISO_FORTRAN_ENV, ONLY : output_unit
      IMPLICIT NONE

      FLUSH(output_unit)
      END SUBROUTINE f77flush
