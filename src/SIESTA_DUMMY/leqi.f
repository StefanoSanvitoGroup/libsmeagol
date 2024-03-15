      FUNCTION LEQI(str1, str2)
      IMPLICIT NONE
      character(len=*) :: str1, str2
      logical :: leqi
C
C  Case-insensitive lexical equal-to comparison
C
      character(len=LEN_TRIM(str1)) :: str1_uc
      character(len=LEN_TRIM(str2)) :: str2_uc

      leqi = .FALSE.

      IF (LEN(str1_uc) == LEN(str2_uc)) THEN
         str1_uc = str1
         CALL toupper_str(str1_uc)

         str2_uc = str2
         CALL toupper_str(str2_uc)

         leqi = (str1_uc == str2_uc)
      END IF

      CONTAINS

      ELEMENTAL SUBROUTINE toupper_str(str)
      character(len=*), intent(inout)                    :: str

      integer, parameter :: iachar_a = IACHAR('a')
      integer, parameter :: iachar_z = IACHAR('z')
      integer, parameter :: iachar_a_minus_A = IACHAR('a')-IACHAR('A')

      integer :: i, iascii

      DO i = 1, LEN(str)
         iascii = IACHAR(str(i:i))
         IF ((iascii >= iachar_a) .AND. (iascii <= iachar_z)) THEN
            str(i:i) = ACHAR(iascii - iachar_a_minus_A)
         END IF
      END DO

      END SUBROUTINE toupper_str

      END FUNCTION LEQI
