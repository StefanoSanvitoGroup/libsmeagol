      module atmfuncs
      implicit none
      private

      public :: symfio,labelfis,cnfigfio,lofio,mofio,zetafio

      contains

      FUNCTION symfio(is, io)
      integer, intent(in) :: is, io
      character(len=20) :: symfio

      symfio = ""
      CALL stop_not_implemented("[STOP] symfio() is not implemented")

      END FUNCTION symfio

      FUNCTION labelfis(is)
      integer, intent(in) :: is
      character(len=20) ::  labelfis

      labelfis = ""
      CALL stop_not_implemented("[STOP] labelfis() is not implemented")

      end function labelfis

      FUNCTION cnfigfio(is, io)
      integer, intent(in) :: is, io
      integer :: cnfigfio

      cnfigfio = 0
      CALL stop_not_implemented("[STOP] cnfigfio() is not implemented")

      END FUNCTION cnfigfio

      FUNCTION lofio(is, io)
      integer, intent(in) :: is, io
      integer :: lofio

      lofio = 0
      CALL stop_not_implemented("[STOP] lofio() is not implemented")

      END FUNCTION lofio

      FUNCTION mofio(is, io)
      integer, intent(in) :: is, io
      integer :: mofio

      mofio = 0
      CALL stop_not_implemented("[STOP] mofio() is not implemented")

      END FUNCTION mofio

      FUNCTION zetafio(is, io)
      integer, intent(in) :: is, io
      integer :: zetafio

      zetafio = 0
      CALL stop_not_implemented("[STOP] zetafio() is not implemented")

      END FUNCTION zetafio

      end module atmfuncs
