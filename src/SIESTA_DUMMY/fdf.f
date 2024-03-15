      module fdf
      implicit none
      private

      public :: fdf_block

      contains

      function fdf_block(label,unit) result(is_block)
      character(len=*), intent(in) :: label
      integer :: unit
      logical :: is_block

      is_block = .false.
      CALL stop_not_implemented("[STOP] fdf_block() is not implemented")

      end function fdf_block

      end module fdf
