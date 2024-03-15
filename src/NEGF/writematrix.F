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
!                   WRITEPHIASCII,
!                   WRITEPHIBIN,
!                   WRITEMAT2,
!                   WRITEMAT7C,
!                   WRITEMAT7B,
!                   WRITEMATCOMPLEX,
!                   WRITEMATREAL,
!                   WRITEMAT7,
!                   WRITEMAT6,
!                   WRITEMAT8,
!                   WRITEMAT9,
!                   WRITEMAT3,
!                   WRITEMAT4,
!                   WRITEMAT5,
!                   WRITEINDEX  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!

      subroutine writephiascii(N1,psi,psit,ene,ibs)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************


      use atmfuncs, only: symfio,labelfis,cnfigfio
      use negfmod

      implicit none

      integer iuo,N1,ibs, iuo1

      character sym_label(em_nuo)*20, atm_label(em_nuo)*20
      DOUBLE COMPLEX :: psi(N1),psit(N1),ene
      logical noncollinear



      if (N1 .eq. 2*em_nuo) then
          noncollinear =.true.
      else
          noncollinear =.false.
      endif
!      write(*,*)"nuon1=",em_nuo,n1
            write(12347,*)
            write(12347,'(a22,2x,i6)') 'Wavefunction = ', ibs
            write(12347,'(a22,2x,f10.6)') 'Energy (eV) = ', 
     .                                  DREAL(ene)*13.60580d0
        write(12347,'(a72)')   ' ---------------------------------------
     .--------------------------------'
        write(12347,'(a72)')  '  Atom  Species Orb-global  Orb-in-atom
     . Orb-type      Re(psi)  Im(psi)'


      do iuo1 = 1, N1
        if (iuo1 .gt. em_nuo) then
            iuo =iuo1 -em_nuo
        else
            iuo =iuo1
        endif
        sym_label(iuo)=symfio(em_isa(em_iaorb(iuo)),em_iphorb(iuo))
        atm_label(iuo)=labelfis(em_isa(em_iaorb(iuo)))
        write(12347,'(i6,5x,a3,1x,i10,8x,i3,7x,i1,a7,1x,4(f10.6))')
     .       em_iaorb(iuo),atm_label(iuo),iuo,
     .       em_iphorb(iuo),
     .       cnfigfio(em_isa(em_iaorb(iuo)),em_iphorb(iuo)),
     .       sym_label(iuo),DBLE(psi(iuo1)),DIMAG(psi(iuo1)),
     .       DBLE(psit(iuo1)),DIMAG(psit(iuo1))
      enddo

!            do j = 1,ntot
!              write(6,'(i6,5x,a3,1x,i10,8x,i3,7x,i1,a7,1x,2(f10.6))') 
!     .                 iaorb(j),labelfis(isa(iaorb(j))),j,
!     .                 iphorb(j), cnfigfio(isa(iaorb(j)),iphorb(j)),
!     .                 symfio(isa(iaorb(j)),iphorb(j)),
!     .                 aux(1,j), aux(2,j)
!            enddo

        write(12347,'(a72)')   ' ---------------------------------------
     .--------------------------------'

      end subroutine writephiascii



      subroutine writephibin(N1,psi,ene,ibs)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      use atmfuncs, only: symfio,labelfis,cnfigfio
      use negfmod

      implicit none

      integer iuo,N1,ibs, iuo1

      character sym_label(em_nuo)*20, atm_label(em_nuo)*20
      DOUBLE COMPLEX :: psi(N1),ene
      logical noncollinear


      if (N1 .eq. 2*em_nuo) then
          noncollinear =.true.
      else
          noncollinear =.false.
      endif

            write(12348) ibs
            write(12348) DREAL(ene)*13.60580d0

      do iuo1 = 1, N1
        if (iuo1 .gt. em_nuo) then
            iuo =iuo1 -em_nuo
        else
            iuo =iuo1
        endif
        sym_label(iuo)=symfio(em_isa(em_iaorb(iuo)),em_iphorb(iuo))
        atm_label(iuo)=labelfis(em_isa(em_iaorb(iuo)))
!        write(12348)
!     .       em_iaorb(iuo),atm_label(iuo),iuo,
!     .       em_iphorb(iuo),
!     .       cnfigfio(em_isa(em_iaorb(iuo)),em_iphorb(iuo)),
!     .       sym_label(iuo),DREAL(psi(iuo)),DIMAG(psi(iuo))
        write(12348)
     .       em_iaorb(iuo),
     .       labelfis(em_isa(em_iaorb(iuo))),iuo,
     .       em_iphorb(iuo),
     .       cnfigfio(em_isa(em_iaorb(iuo)),em_iphorb(iuo)),
     .       symfio(em_isa(em_iaorb(iuo)),em_iphorb(iuo)),
     .       DBLE(psi(iuo1)),DIMAG(psi(iuo1))
!        write(12348)
!     .       7,labelfis(em_isa(em_iaorb(iuo))),2,
!     .       3,
!     .       4,
!     .       symfio(em_isa(em_iaorb(iuo)),em_iphorb(iuo)),
!     .       DREAL(psi(iuo)),DIMAG(psi(iuo))
      enddo

!            do j = 1,ntot
!              write(6,'(i6,5x,a3,1x,i10,8x,i3,7x,i1,a7,1x,2(f10.6))') 
!     .                 iaorb(j),labelfis(isa(iaorb(j))),j,
!     .                 iphorb(j), cnfigfio(isa(iaorb(j)),iphorb(j)),
!     .                 symfio(isa(iaorb(j)),iphorb(j)),
!     .                 aux(1,j), aux(2,j)
!            enddo

      end subroutine writephibin



      subroutine writemat2(mat,n,m,nam)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      IMPLICIT NONE
      INTEGER n,m,i,jm
      DOUBLE COMPLEX, DIMENSION(n,m) :: mat
      CHARACTER(LEN=*) :: nam

      i=n
      write(*,*)"start ",nam
        do jm=1,m
      do i=1,n
          write(*,*)"[",i,",",jm,"]=",mat(i,jm)
        enddo
      enddo
      write(*,*)"end ",nam


      END subroutine writemat2

      subroutine writemat7c(ene,mat,n,m,nam)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      IMPLICIT NONE
      INTEGER n,m,i,jm
      DOUBLE PRECISION ene
      DOUBLE COMPLEX, DIMENSION(n,m) :: mat
      CHARACTER(LEN=*) :: nam

      do i=1,n
        do jm=1,m
          write(12347,'(a,"[",i9,",",i9,"]=",3e15.6)')nam,i,jm,ene,
     .              DREAL(mat(i,jm)),DIMAG(mat(i,jm))
        enddo
      enddo

      END subroutine writemat7c


      subroutine writemat7b(ene,mat,n,m,nam)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      IMPLICIT NONE
      INTEGER n,m,i,jm
      DOUBLE PRECISION ene
      DOUBLE COMPLEX, DIMENSION(n,m) :: mat
      CHARACTER(LEN=*) :: nam

      return
      do i=1,n
        do jm=1,m
          write(12347,'(a,"[",i9,",",i9,"]=",3e15.6)')nam,i,jm,ene,
     .              DREAL(mat(i,jm)),DIMAG(mat(i,jm))
        enddo
      enddo

      END subroutine writemat7b

      subroutine writematcomplex(mat,n,m,nam)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      IMPLICIT NONE
      INTEGER n,m,i,jm
      DOUBLE COMPLEX, DIMENSION(n,m) :: mat
      CHARACTER(LEN=*) :: nam

      do i=1,n
        do jm=1,m
          write(*,'(a,"[",i9,",",i9,"]=",3e15.6)')nam,i,jm,
     .              dreal(mat(i,jm)),dimag(mat(i,jm))
        enddo
      enddo

      END subroutine writematcomplex



      subroutine writematreal(mat,n,m,nam)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      IMPLICIT NONE
      INTEGER n,m,i,jm
      DOUBLE PRECISION, DIMENSION(n,m) :: mat
      CHARACTER(LEN=*) :: nam

      do i=1,n
        do jm=1,m
          write(12347,'(a,"[",i9,",",i9,"]=",3e15.6)')nam,i,jm,
     .              mat(i,jm)
        enddo
      enddo

      END subroutine writematreal



      subroutine writemat7(ene,mat,n,m,nam)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      IMPLICIT NONE
      INTEGER n,m,i,jm
      DOUBLE PRECISION ene
      DOUBLE COMPLEX, DIMENSION(n,m) :: mat
      CHARACTER(LEN=*) :: nam

      do i=1,n
        do jm=1,m
          write(12347,'(a,"[",i9,",",i9,"]=",3e15.6)')nam,i,jm,ene,
     .              DREAL(mat(i,jm)),DIMAG(mat(i,jm))
        enddo
      enddo

      END subroutine writemat7


      subroutine writemat6(ene,mat,n,m,nam)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      IMPLICIT NONE
      INTEGER n,m,i,jm
      DOUBLE PRECISION ene
      DOUBLE COMPLEX, DIMENSION(n,m) :: mat
      CHARACTER(LEN=*) :: nam

      do jm=1,m
        do i=1,n
          if(ABS(mat(i,jm)).gt.1D0)
     .    write(12346,'(a,"[",i9,",",i9,"]=",3e15.6)')nam,i,jm,ene,
     .              DREAL(mat(i,jm)),DIMAG(mat(i,jm))
        enddo
      enddo

      END subroutine writemat6

      subroutine writemat8(ene,mat,n,m,tol,nam)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      IMPLICIT NONE
      INTEGER n,m,i,j
      DOUBLE PRECISION ene,tol
      DOUBLE COMPLEX, DIMENSION(n,m) :: mat
      CHARACTER(LEN=*) :: nam

      if(tol .lt. 0.0D0) return

      do i=1,n
        do j=1,m
          if(ABS(mat(i,j)).ge.tol)
     .    write(12346,'(a,"[",i9,",",i9,"]=",3e15.6)')nam,i,j,ene,
     .              DREAL(mat(i,j)),DIMAG(mat(i,j))
        enddo
      enddo

      END subroutine writemat8

      subroutine writemat9(ik,ene,mat,n,m,tol,nam)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      IMPLICIT NONE
      INTEGER n,m,i,jm,ik
      DOUBLE PRECISION ene,tol
      DOUBLE COMPLEX, DIMENSION(n,m) :: mat
      CHARACTER(LEN=*) :: nam
      CHARACTER(LEN=100) :: nam2
     

      write(12346,*)"start ",nam," tolerance=",tol
        do jm=1,m
      do i=1,n
          if(i==jm)then
            nam2=TRIM(nam // "_D")
          else
            nam2=TRIM(nam)
          endif
          if((ABS(mat(i,jm)).ge.tol).and.(tol.ge.0D0))
!          if(.true.)
     .        write(12346,'(a,"[",i9,",",i9,"]=",i7,3e15.6)')
     .        TRIM(nam2),i,jm,ik,ene,DREAL(mat(i,jm)),DIMAG(mat(i,jm))
        enddo
      enddo
      write(12346,*)"end",nam

      END subroutine writemat9






      subroutine writemat3(ene,mat,n,m,nam)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      IMPLICIT NONE
      INTEGER n,m,i,jm
      DOUBLE PRECISION ene
      DOUBLE COMPLEX, DIMENSION(n,m) :: mat
      CHARACTER(LEN=*) :: nam

      do jm=1,m
        do i=1,n
          if(ABS(mat(i,jm)).gt.1D0)
     .    write(*,'(a,"[",i9,",",i9,"]=",3e15.6)')nam,i,jm,ene,
     .              DREAL(mat(i,jm)),DIMAG(mat(i,jm))
        enddo
      enddo

      END subroutine writemat3

      subroutine writemat4(ene,mat,n,m,nam)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      IMPLICIT NONE
      INTEGER n,m,i,jm
      DOUBLE PRECISION ene
      DOUBLE COMPLEX, DIMENSION(n,m) :: mat
      CHARACTER(LEN=*) :: nam
      integer   el,sr


      el=80
      sr=145

      do jm=1,m
        do i=1,n
          if(((jm.le.el).and.(i.ge.sr)).or.((i.le.el).and.(jm.ge.sr)))
     .       then
          if(ABS(mat(i,jm)).gt.0.0)
     .    write(*,'(a,"[",i9,",",i9,"]=",3e15.6)')nam,i,jm,ene,
     .              DREAL(mat(i,jm)),DIMAG(mat(i,jm))
          endif
        enddo
      enddo

      END subroutine writemat4


      subroutine writemat5(ene,mat,n,m,nam)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      IMPLICIT NONE
      INTEGER n,m,i,jm
      DOUBLE PRECISION ene
      DOUBLE COMPLEX, DIMENSION(n,m) :: mat
      CHARACTER(LEN=*) :: nam

      do jm=1,m
        do i=1,n
          if(jm.ne.i.or..true.)then
            if(ABS(mat(i,jm)).gt.1D-10)
     .        write(*,'(a,"[",i9,",",i9,"]=",3e15.6)')nam,i,jm,ene,
     .                DREAL(mat(i,jm)),DIMAG(mat(i,jm))
          elseif(.true.)then
            write(*,'("diag",a,"[",i9,",",i9,"]=",3e15.6)')nam,i,jm,ene,
     .                DREAL(mat(i,jm)),DIMAG(mat(i,jm))
          endif
        enddo
      enddo

      END subroutine writemat5



      subroutine writeindex(nau,nuo,iaorb,iphorb,isa)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      use atmfuncs, only: symfio,labelfis

      integer nau,nuo,iuo,isa(nau),iaorb(nuo),iphorb(nuo)

      character sym_label(nuo)*7, atm_label(nuo)*20


      do iuo = 1, nuo
             sym_label(iuo)=symfio(isa(iaorb(iuo)), iphorb(iuo))
             atm_label(iuo)=labelfis(isa(iaorb(iuo)))
             write(*,*)"label:",iuo,atm_label(iuo),
     .          iaorb(iuo),sym_label(iuo)
      enddo


      end subroutine writeindex


