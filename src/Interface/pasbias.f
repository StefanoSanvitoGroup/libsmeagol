C Subroutine pasbias
C Written by V. M. Garcia-Suarez. October 2003.

      character*(*) function pasbias( str1, str2 )

      character*(*) str1, str2
      m = len(str1)
      do 10 l = 1, m
        if (str1(l:l) .ne. ' ') then
          n = l
          goto 20
        endif
   10 continue
   20 pasbias = str1(n:m-1)//str2
      end

