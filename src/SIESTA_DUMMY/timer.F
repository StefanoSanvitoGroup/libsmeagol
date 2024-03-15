      SUBROUTINE TIMER (PROG,IOPT)
      implicit none
      character(len=*) :: prog
      integer :: iopt

C  The aim of this subroutine is to measure and print both CPU time
C  and wall-time spent in a particular routine.
C
C  ARGUMENTS:
C    PROG: the name of the routine to benchmark
C    IOPT: operation code
C       0 => reset all measured times
C       1 => start counting time spent in the routine
C       2 => stop counting time spent in the routine
C       3 => print time spent in the routine
C
C  SMEAGOL uses this routine to measure performance of its routines
C  (iopt = 1 and 2) but newer prints it (iopt = 3). SMEAGOL expects
C  SIESTA to report the measured times at the end of its run.
C
C  Other DFT programs usually have their own benchmarking routines.
C  Calling a SIESTA-specific routine from their code base does
C  not seem to be a good idea.
C
C  SMEAGOL also has its own benchmarking facility which can be enabled
C  by setting emtimings @ negfmod = .TRUE. global variable.

      END SUBROUTINE TIMER
