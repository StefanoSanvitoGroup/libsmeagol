c
c Minified version of the original SIESTA's io.f module.
c
c SMEAGOL only relies on io_assign() and io_close() routines
c from this module. All other routines were removed.
c
c
c
c This module implements an interface to the FORTRAN logical unit
c system. Based on code by Richard Maine.
c
c Alberto Garcia, December 30, 1996
c Rewritten as a single subroutine
c with multiple entry points, March 7, 1998
c
c This scheme is actually the closest in spirit to f90 modules, but
c in perfectly legal f77.
c
c---------------------------------------------------------------
c
      subroutine io
c
c     Logical unit management. Units 0 to min_lun-1 are "reserved",
c     since most of the "typical" files (output, etc) use them.
c
c     Logical units min_lun to min_max are managed by this module.

      implicit none
c
c----------------------------------------------------------------
c     Module variables
c
      integer min_lun, max_lun, nunits
      parameter (min_lun=10, max_lun=99, nunits=max_lun-min_lun+1)
      logical lun_is_free(min_lun:max_lun)

      save lun_is_free
c-----------------------------------------------------------------
c
c     Internal and dummy variables
c
      integer i, lun, stat
      logical used
c
c-----------------------------------------------------------------
c     Initialization section
c
c
      data lun_is_free /nunits*.true./
c
c------------------------------------------------------------------
c
c     Logical unit management
c
      entry io_assign(lun)
c
c     Looks for a free unit and assigns it to lun
c
      do lun= min_lun, max_lun
         if (lun_is_free(lun)) then
            inquire(unit=lun, opened=used, iostat=stat)
            if (stat .ne. 0) used = .true.
            lun_is_free(lun) = .false.
            if (.not. used) return
         endif
      enddo
      CALL stop_not_implemented("[STOP] No luns available in io_assign")
      return
c
c===
c
      entry io_close(lun)
c
c     Use this routine instead of a simple close!!
c
      close(lun)
      if (lun .ge. min_lun .and. lun .le. max_lun)
     $                     lun_is_free(lun) = .true.
      return

      end subroutine io
