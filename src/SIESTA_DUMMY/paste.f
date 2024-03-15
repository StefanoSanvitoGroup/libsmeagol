      FUNCTION PASTE(str1, str2)
      IMPLICIT NONE
      character(len=*) :: str1, str2
      character(len=*) :: paste

C Concatenates two strings removing all trailing blank characters
      paste = TRIM(str1)//TRIM(str2)

      END FUNCTION PASTE

      FUNCTION PASTEB(str1, str2)
      IMPLICIT NONE
      character(len=*) :: str1, str2
      character(len=*) :: pasteb

C Concatenates two strings leaving one blank character between them
      pasteb = TRIM(str1)//' '//TRIM(str2)

      END FUNCTION PASTEB
