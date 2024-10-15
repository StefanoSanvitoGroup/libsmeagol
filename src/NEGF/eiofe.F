! 
! Copyright (c) Authors:
! Ivan Rungger and Stefano Sanvito
! Trinity College Dublin, Ireland
! October 2008 
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! THE SUBROUTINES
!                   REDUCEENE,
!                   GET_ENEIOFEALL,
!                   GET_NENEUNIQUE,
!                   RANKDOUBLE,
!                   EIOFE  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!

      subroutine reduceene(eimin,
     &            ene1d,eiene1d,nene,
     &            ene1dr,eiene1dr,nener)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      IMPLICIT NONE

      INTEGER nene,nener,i
      double precision ene1d(nene),eimin
      double complex eiene1d(nene)
      double precision, pointer :: ene1dr(:)
      double complex, pointer :: eiene1dr(:)

      nener=0
      do i=1,nene
        if(DABS(DIMAG(eiene1d(i))).lt.eimin)nener=nener+1
      enddo
      allocate(eiene1dr(nener),ene1dr(nener))
      nener=0
      do i=1,nene
        if(DABS(DIMAG(eiene1d(i))).lt.eimin)then
          nener=nener+1
          eiene1dr(nener)=eiene1d(i)
          ene1dr(nener)=ene1d(i)
        endif
      enddo

      end subroutine reduceene

      subroutine get_eneiofeall(e,nrank,neitot,eis,nei,nene)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      integer neitot,nene,eis(nene),nei(nene),i,nrank(neitot),neib,ie
      double precision e(neitot),eb

      eb=e(nrank(1))
      neib=1
      ie=1
      eis(ie)=1
      do i=2,neitot
        if(eb.ne.e(nrank(i)))then
          nei(ie)=i-eis(ie)
          ie=ie+1
          eis(ie)=i
          eb=e(nrank(i))
        endif
      enddo
      nei(nene)=neitot-eis(nene)+1
      end

      subroutine get_neneunique(e,nrank,neitot,nene)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      integer nene,neitot,i,nrank(neitot)
      double precision e(neitot),eb

      nene=1
      eb=e(nrank(1))
      do i=2,neitot
        if(eb.ne.e(nrank(i)))then
          nene=nene+1
          eb=e(nrank(i))
        endif
      enddo


      end

      subroutine rankdouble(e,nrank,n)
      
! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      implicit none
      integer n,nrank(n)
      double precision e(n)
     
      call mrgrnk2(e, nrank,n)
      
      end

      subroutine eiofe(e,ei,neitot,eizeron,ezero,nzero)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      implicit none

      integer i,j,neitot,nene,i2,i2min,i2max,imin,imax,
     .    nzero,izero,nzerostep
      integer, allocatable :: nrank(:),nrankei(:),eies(:),neie(:),
     .    neigh(:)
      double precision dx,ene,ene2,de,e(neitot)
      double precision, allocatable :: ebuf(:),eibufr(:)
      double complex ei(neitot),zi
      double complex, allocatable :: eibuf(:),
     .    eibuf2(:)
      double precision, pointer :: ezero(:,:)
      double complex, pointer :: eizeron(:,:)
      double precision, allocatable :: ezerobuf(:,:)
      double complex, allocatable :: eizerobuf(:,:),eizero(:,:)

      PARAMETER(zi=(0.0d0,1.0d0))


      allocate(ebuf(neitot), nrank(neitot))
      call rankdouble(e,nrank,neitot)

      call get_neneunique(e,nrank,neitot,nene)
      write(*,*)"#nene=",nene
      allocate(neie(nene),eies(nene))
      call get_eneiofeall(e,nrank,neitot,eies,neie,nene)

!!      do i=1,nene
!!!!        write(*,*)i,"neie(i)=",eies(i),neie(i)
!!        do j=0,neie(i)-1
!!          write(*,'(2i,3f)')i,neie(i),e(nrank(eies(i)+j)),
!!     .        DREAL(ei(nrank(eies(i)+j))),
!!     .        DIMAG(ei(nrank(eies(i)+j)))
!!        enddo
!!      enddo


      do i=1,nene
!!        write(*,*)i,"neie(i)=",eies(i),neie(i)
        allocate(eibuf(neie(i)),eibufr(neie(i)),nrankei(neie(i)))
        do j=1,neie(i)
          eibufr(j)=ei(nrank(eies(i)+j-1))
        enddo
        call rankdouble(eibufr,nrankei,neie(i))
        do j=1,neie(i)
          eibuf(j)=ei(nrank(eies(i)+nrankei(j)-1))
        enddo
        do j=1,neie(i)
          ei(nrank(eies(i)+j-1))=eibuf(j)
        enddo
!        do j=1,neie(i)
!          write(*,'(2i,3f)')i,neie(i),e(nrank(eies(i)+j-1)),
!     .        DREAL(ei(nrank(eies(i)+j-1))),
!     .        DIMAG(ei(nrank(eies(i)+j-1)))
!        enddo
!!!        do j=1,neie(i)
!!!          write(*,'(4f)')e(nrank(eies(i)+j-1)),e(nrank(eies(i)+j-1)),
!!!     .        DREAL(ei(nrank(eies(i)+j-1))),
!!!     .        DIMAG(ei(nrank(eies(i)+j-1)))
!!!        enddo
        deallocate(eibuf,eibufr,nrankei)
      enddo
!      stop

      nzerostep=2
      nzero=nzerostep
      izero=1
      allocate(eizero(nzero,2),ezero(nzero,2))
      do i=1,nene-1
        ene=e(nrank(eies(i)))
        ene2=e(nrank(eies(i+1)))
        de=ene2-ene
        allocate(eibuf(neie(i)),eibuf2(neie(i+1)),neigh(neie(i)))
        neigh=0

        do j=1,neie(i)
          eibuf(j)=ei(nrank(eies(i)+j-1))
        enddo
        do j=1,neie(i+1)
          eibuf2(j)=ei(nrank(eies(i+1)+j-1))
        enddo

        imin=1
        imax=neie(i)
        i2min=1
        i2max=neie(i+1)

        loopei: do j=1,neie(i)
          dx=1D100

          do i2=i2min,i2max
            if(ABS(DREAL(eibuf(j))-DREAL(eibuf2(i2))).gt.dx)then
              i2min=i2-1
              exit
            endif
            neigh(j)=i2
            dx=ABS(DREAL(eibuf(j))-DREAL(eibuf2(i2)))
            if(i2.eq.i2max) exit loopei
          enddo

          if(j.lt.neie(i))then
            if((ABS(DREAL(eibuf(j+1)-eibuf2(neigh(j)))).lt.dx)
     .          .or.(ABS(DREAL(eibuf(j)-eibuf2(neigh(j)))).gt.5 * de))
     .          then
!              write(*,'("de",i,3f)')j,DREAL(eibuf(j)),DREAL(eibuf(j+1)),
!     .       DREAL(eibuf2(neigh(j)))
              neigh(j)=-1
            
!             stop
            endif
          endif

        enddo loopei

        imin=1
        imax=neie(i)
        i2min=1
        i2max=neie(i+1)

        loopei2: do j=1,neie(i)
          dx=1D100
          if(neigh(j).gt.0)cycle
          do i2=i2min,i2max
            if(ABS(DREAL(eibuf(j)-eibuf2(i2))).gt.dx)then
              i2min=i2-1
              exit
            endif
            if(ABS(DREAL(eibuf(j)-eibuf2(i2))).lt.2D0 * de)then
              neigh(j)=i2
            endif
            dx=ABS(DREAL(eibuf(j))-DREAL(eibuf2(i2)))
            if(i2.eq.i2max) exit loopei2
          enddo

        enddo loopei2


!!        do j=1,neie(i)
!!!          if(neigh(j).ne.0)
!!          if(neigh(j).gt.0)then
!!            write(*,'(3(i," "),3(f," "))')
!!     .        i,j,neigh(j),ene,DREAL(eibuf(j)),
!!     .        DREAL(eibuf2(neigh(j)))
!!          else
!!            write(*,'(3(i," "),f)')i,j,neigh(j),DREAL(eibuf(j))
!!          endif
!!        enddo

        do j=1,neie(i)
          if(neigh(j).gt.0)then
            if(((ene.le.DREAL(eibuf(j))).and.
     .          (ene2.ge.DREAL(eibuf2(neigh(j))))).or.
     .          ((ene.ge.DREAL(eibuf(j))).and.
     .          (ene2.le.DREAL(eibuf2(neigh(j))))))then
              if(izero.gt.nzero)then
!                write(*,*)"reallocating",nzero
                allocate(eizerobuf(nzero,2),ezerobuf(nzero,2))
                eizerobuf=eizero
                ezerobuf=ezero
                deallocate(eizero,ezero)
                allocate(eizero(nzero+nzerostep,2),
     .              ezero(nzero+nzerostep,2))
                eizero(1:nzero,:)=eizerobuf
                ezero(1:nzero,:)=ezerobuf
                nzero=nzero+nzerostep
                deallocate(eizerobuf,ezerobuf)
              endif
              ezero(izero,1)=ene 
              ezero(izero,2)=ene2
              eizero(izero,1)=eibuf(j)
              eizero(izero,2)=eibuf2(neigh(j))
              izero=izero+1

!              write(*,'("zero",4(f," "))')ene,ene2,DREAL(eibuf(j)),
!     .          DREAL(eibuf2(neigh(j)))

            endif
          endif

        enddo
        nzero=izero-1
        
        deallocate(neigh,eibuf,eibuf2)
      enddo
      allocate(eizeron(nzero,3))
      eizeron=0D0
      do i=1,nzero
        eizeron(i,1)=eizero(i,1)
        eizeron(i,2)=eizero(i,2)
        eizeron(i,3)=(ezero(i,1) * DREAL(eizero(i,2)) - ezero(i,2) *
     &      DREAL(eizero(i,1)))/
     &      (ezero(i,1)-ezero(i,2)-
     &      DREAL(eizero(i,1))+DREAL(eizero(i,2)))
        eizeron(i,3)=eizeron(i,3)+ zi *(
     &      ((DIMAG(eizero(i,2))-DIMAG(eizero(i,1)))/
     &      (ezero(i,2)-ezero(i,1))) *
     &      (DREAL(eizeron(i,3))-ezero(i,1))+
     &      DIMAG(eizero(i,1)))
      enddo
      deallocate(eizero)

      end subroutine eiofe

