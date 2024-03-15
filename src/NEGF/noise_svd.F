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
!                   ADDNOISE,
!                   ADDNOISER,
!                   SVDMK,
!                   SVDMTEST,
!                   SVDM,
!                   SVDM2,
!                   SVDM3  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!


      subroutine addnoise(mat,n,delta)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      IMPLICIT NONE

      INTEGER n,i,j
      double complex mat(n,n)
      double complex noise
      double PRECISION delta
      REAL  rnum,inum
      logical, SAVE :: initseed=.true.

      if(initseed)then
!        write(12347,*)"initseed"
        call random_seed
        initseed=.false.
      endif

      do i=1,n
        do j=1,n
          call random_number(rnum)
          call random_number(inum)
!          write(*,*)"rnum=",rnum
          noise=(rnum-0.5D0+ (0D0,1D0) * (inum-0.5D0)) * delta
!          if(mat(i,j).ne.0D0)then
            mat(i,j)=mat(i,j)+noise
!          endif
        enddo
      enddo


      end subroutine addnoise
     
      subroutine addnoiser(mat,n,delta)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      use negfmod, only: outinfo
      IMPLICIT NONE

      INTEGER n,i,j
      double complex mat(n,n)
      double complex noise
      double PRECISION delta
      REAL  rnum,inum
      logical, SAVE :: initseed=.true.

      if(initseed)then
          if (outinfo) write(12347,*)"initseed"
        call random_seed
        initseed=.false.
      endif

      do i=1,n
        do j=1,n
          call random_number(rnum)
          call random_number(inum)
!          write(*,*)"rnum=",rnum
          noise=(rnum-0.5D0+ (0D0,1D0) * (inum-0.5D0)) * delta
          mat(i,j)=mat(i,j) * ( 1D0 + noise )
        enddo
      enddo


      end subroutine addnoiser
     



      subroutine svdmk(
     .    matin,invmat,n,svdtol,noisemax,fmatinv,faddnoise,kappa,s1,sn)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      IMPLICIT NONE

      integer n,m,nwork,INFO,i,j,i1
      double COMPLEX matin(n,n),invmat(n,n),u(n,n),
     .    vt(n,n)
      double PRECISION svdtol,sv(n),svmax,condn,noisemax,
     .    kappa,sn,s1
      double COMPLEX, ALLOCATABLE :: work(:)
      DOUBLE PRECISION, ALLOCATABLE :: rwork(:)
      DOUBLE COMPLEX,PARAMETER :: alpha=(1.D0,0D0)
      logical faddnoise,fmatinv


      nwork=8
      allocate(work(nwork * n))
      allocate(rwork(5 * n))


!== begin old SVD

      CALL ZGESVD( 'A', 'A', n, n, matin, n,sv, u, n, vt,n,
     $                   work, nwork * n, rwork, INFO )
      if(info.ne.0)write(*,*)"infosvd=",info
      s1=sv(1)
      sn=sv(n)
      kappa=sv(1)/sv(n)
      return

      svmax=sv(1)
      m=n
      do i=1,n
        condn=svmax/sv(i)
!        write(*,*)"sigmak1=",DREAL(ene),i,sv(i),sv(i)/svmax,condn
        if(sv(i)/svmax.lt.svdtol)then
          m=i-1
          exit
        endif
      enddo

!      write(12347,*)"m,svmax=",m,svmax,svdtol
!      write(*,*)"m,svmax=",m,svmax,svdtol


      do i=1,n
        if(sv(i)/svmax.lt.svdtol)then
          sv(i)= svdtol * svmax 
        endif
!        sv(i)= sv(i) + svdtol * svmax
!        if(sv(i).lt.1D-15)then
!          sv(i)= sv(i) /1D0
!        endif
      enddo


      matin=0D0
      do i=1,n
        do j=1,n
          do i1=1,n
            matin(i,j)=matin(i,j)+u(i,i1) * sv(i1)  * vt(i1,j)
          enddo
        enddo
      enddo

!      k1o=0D0
!      do i=1,n
!        do j=1,n
!          do i1=1,n
!            k1o(i,j)=k1o(i,j)+
!     .          DCONJG(vt(i1,i)) *  DCONJG(u(j,i1)) / sv(i1) 
!          enddo
!        enddo
!      enddo


!      invmat=0D0
!      do i=1,n
!        invmat(i,i)=1D0/sv(i)
!      enddo
!      k1ob=matmul(invmat,k1o)

!      call geigenvalues2(k1ob,zv,n,k1o,info)


!      do i=1,n
!        call alphaofz(alphak,lz,zv(i))
!        write(*,*)"kapprox=",DREAL(ene),lz,alphak, svdtol
!      enddo


      if(fmatinv)then
        invmat=0D0
        do i=1,n
          do j=1,n
            do i1=1,n
              invmat(i,j)=invmat(i,j)+
     .            DCONJG(vt(i1,i)) *  DCONJG(u(j,i1)) / sv(i1) 
            enddo
          enddo
        enddo
      endif

      if(faddnoise)then
        call addnoise(matin,n,noisemax)
        call addnoise(invmat,n,noisemax)
      endif
!== end old SVD

      end subroutine svdmk




      subroutine svdmtest(ene,
     .    matin,invmat,n,svdtol,noisemax,fmatinv,faddnoise)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      use negfmod, only: outinfo
      IMPLICIT NONE

      integer n,m,nwork,INFO,i,j,i1
      double COMPLEX matin(n,n),invmat(n,n),u(n,n),
     .    vt(n,n),ene
      double PRECISION svdtol,sv(n),svmax,condn,noisemax
      double COMPLEX, ALLOCATABLE :: work(:)
      DOUBLE PRECISION, ALLOCATABLE :: rwork(:)
      DOUBLE COMPLEX,PARAMETER :: alpha=(1.D0,0D0)
      logical faddnoise,fmatinv


      nwork=8
      allocate(work(nwork * n))
      allocate(rwork(5 * n))


!== begin old SVD

      CALL ZGESVD( 'A', 'A', n, n, matin, n,sv, u, n, vt,n,
     $                   work, nwork * n, rwork, INFO )
      if(info.ne.0)write(*,*)"infosvd=",info

      svmax=sv(1)
      m=n
      do i=1,n
        condn=svmax/sv(i)
!        write(*,*)"sigmak1=",DREAL(ene),i,sv(i),sv(i)/svmax,condn
        if(sv(i)/svmax.lt.svdtol)then
          m=i-1
          exit
        endif
      enddo

      if (outinfo) 
     .  write(12347,*)"m,svmaxt=",DREAL(ene),svmax,sv(n),sv(n-1),
     .    sv(n-2),sv(n-3), svdtol,m
!      write(*,*)"m,svmax=",m,svmax,svdtol


      do i=1,n
        if(sv(i)/svmax.lt.svdtol)then
          sv(i)= svdtol * svmax 
        endif
!        sv(i)= sv(i) + svdtol * svmax
!        if(sv(i).lt.1D-15)then
!          sv(i)= sv(i) /1D0
!        endif
      enddo


      matin=0D0
      do i=1,n
        do j=1,n
          do i1=1,n
            matin(i,j)=matin(i,j)+u(i,i1) * sv(i1)  * vt(i1,j)
          enddo
        enddo
      enddo

!      k1o=0D0
!      do i=1,n
!        do j=1,n
!          do i1=1,n
!            k1o(i,j)=k1o(i,j)+
!     .          DCONJG(vt(i1,i)) *  DCONJG(u(j,i1)) / sv(i1) 
!          enddo
!        enddo
!      enddo


!      invmat=0D0
!      do i=1,n
!        invmat(i,i)=1D0/sv(i)
!      enddo
!      k1ob=matmul(invmat,k1o)

!      call geigenvalues2(k1ob,zv,n,k1o,info)


!      do i=1,n
!        call alphaofz(alphak,lz,zv(i))
!        write(*,*)"kapprox=",DREAL(ene),lz,alphak, svdtol
!      enddo


      if(fmatinv)then
        invmat=0D0
        do i=1,n
          do j=1,n
            do i1=1,n
              invmat(i,j)=invmat(i,j)+
     .            DCONJG(vt(i1,i)) *  DCONJG(u(j,i1)) / sv(i1) 
            enddo
          enddo
        enddo
      endif

      if(faddnoise)then
        call addnoise(matin,n,noisemax)
        call addnoise(invmat,n,noisemax)
      endif
!== end old SVD

      end subroutine svdmtest



      subroutine svdm(
     .    matin,invmat,n,svdtol,noisemax,fmatinv,faddnoise)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      IMPLICIT NONE

      integer n,m,nwork,INFO,i,j,i1
      double COMPLEX matin(n,n),invmat(n,n),u(n,n),
     .    vt(n,n)
      double PRECISION svdtol,sv(n),svmax,condn,noisemax
      double COMPLEX, ALLOCATABLE :: work(:)
      DOUBLE PRECISION, ALLOCATABLE :: rwork(:)
      DOUBLE COMPLEX,PARAMETER :: alpha=(1.D0,0D0)
      logical faddnoise,fmatinv


      nwork=8
      allocate(work(nwork * n))
      allocate(rwork(5 * n))


!== begin old SVD

      CALL ZGESVD( 'A', 'A', n, n, matin, n,sv, u, n, vt,n,
     $                   work, nwork * n, rwork, INFO )
      if(info.ne.0)write(*,*)"infosvd=",info

      svmax=sv(1)
      m=n
      do i=1,n
        condn=svmax/sv(i)
!        write(*,*)"sigmak1=",DREAL(ene),i,sv(i),sv(i)/svmax,condn
        if(sv(i)/svmax.lt.svdtol)then
          m=i-1
          exit
        endif
      enddo

!      write(12347,*)"m,svmax=",m,svmax,svdtol
!      write(*,*)"m,svmax=",m,svmax,svdtol


      do i=1,n
        if(sv(i)/svmax.lt.svdtol)then
          sv(i)= svdtol * svmax 
        endif
!        sv(i)= sv(i) + svdtol * svmax
!        if(sv(i).lt.1D-15)then
!          sv(i)= sv(i) /1D0
!        endif
      enddo


      matin=0D0
      do i=1,n
        do j=1,n
          do i1=1,n
            matin(i,j)=matin(i,j)+u(i,i1) * sv(i1)  * vt(i1,j)
          enddo
        enddo
      enddo

!      k1o=0D0
!      do i=1,n
!        do j=1,n
!          do i1=1,n
!            k1o(i,j)=k1o(i,j)+
!     .          DCONJG(vt(i1,i)) *  DCONJG(u(j,i1)) / sv(i1) 
!          enddo
!        enddo
!      enddo


!      invmat=0D0
!      do i=1,n
!        invmat(i,i)=1D0/sv(i)
!      enddo
!      k1ob=matmul(invmat,k1o)

!      call geigenvalues2(k1ob,zv,n,k1o,info)


!      do i=1,n
!        call alphaofz(alphak,lz,zv(i))
!        write(*,*)"kapprox=",DREAL(ene),lz,alphak, svdtol
!      enddo


      if(fmatinv)then
        invmat=0D0
        do i=1,n
          do j=1,n
            do i1=1,n
              invmat(i,j)=invmat(i,j)+
     .            DCONJG(vt(i1,i)) *  DCONJG(u(j,i1)) / sv(i1) 
            enddo
          enddo
        enddo
      endif

      if(faddnoise)then
        call addnoise(matin,n,noisemax)
        call addnoise(invmat,n,noisemax)
      endif
!== end old SVD

      end subroutine svdm



      subroutine svdm2(mat,matinv,n,svdtol)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      IMPLICIT NONE

      integer n,nwork,INFO,i,j
      double COMPLEX mat(n,n),matinv(n,n),k1o(n,n),u(n,n),
     .    vt(n,n)
      double PRECISION svdtol,sv(n),InvTolNorm,svdtolint
      double COMPLEX, ALLOCATABLE :: work(:)
      DOUBLE PRECISION, ALLOCATABLE :: rwork(:)
      DOUBLE COMPLEX,PARAMETER :: alpha=(1.D0,0D0)


      nwork=8
      allocate(work(nwork * n))
      allocate(rwork(5 * n))


!== begin old SVD

      k1o=mat
      CALL ZGESVD( 'A', 'A', n, n, k1o, n,sv, u, n, vt,n,
     $                   work, nwork * n, rwork, INFO )
      write(*,*)"infosvdm2=",INFO

       SVDTOLint=svdtol
       InvTolNorm=1D+3

       Do I=1,n
        Do J=1,n
          If (sv(J) .GT. SVDTolint) then
            if (sv(J).GT. InvTolNorm) then
              write(*,*)"gtinvtolnorm"
              k1o(I,J)=DCONJG(vt(J,I))/InvTolNorm
              matinv(I,J)=u(I,J)*InvTolNorm
            else
              k1o(I,J)=DCONJG(vt(J,I))/sv(J)
              matinv(I,J)=u(I,J)*sv(J)
            endif
          else    
            k1o(I,J)=DCONJG(vt(J,I))/(SVDTolint)
            matinv(I,J)=u(I,J)*(SVDTolint)
          EndIf
        EndDo
       EndDo
       
       CALL ZGEMM('N','N',n,n,n,(1.D0,0.D0),matinv,n,vt,n,
     &  (0.D0,0.D0),mat,n)   
       CALL ZGEMM('N','C',n,n,n,(1.D0,0.D0),k1o,n,u,n,
     &  (0.D0,0.D0),matinv,n)
       






!-------      sigma=0D0
!-------      do i=1,n
!-------        if(sv(i).lt.svdtol)then
!-------          sigma(i,i)= svdtol * svmax 
!-------        else
!-------          sigma(i,i)= sv(i) 
!-------        endif
!-------      enddo
!-------
!-------      call ZGEMM('N', 'N', n, n, n,
!-------     .            alpha, u, n, sigma, n, 0, k1ob, n)
!-------      call ZGEMM('N', 'N', n, n, n,
!-------     .            alpha, k1ob, n, vt, n, 0, k1o, n)
!-------      mat=k1o
!-------
!-------      sigma=0D0
!-------      do i=1,n
!-------        if(sv(i)/svmax.lt.svdtol)then
!-------          sigma(i,i)=1D0/( svdtol * svmax )
!-------        else
!-------!          sigma(i,i)=0D0
!-------          sigma(i,i)=1D0/( sv(i) )
!-------        endif
!-------      enddo
!-------
!-------      call ZGEMM('N', 'N', n, n, n,
!-------     .        alpha, DCONJG(TRANSPOSE(vt)), n, sigma, n, 0, k1ob, n)
!-------      call ZGEMM('N', 'N', n, n, n,
!-------     .        alpha, k1ob, n, DCONJG(TRANSPOSE(u)), n, 0, k1o, n)
!-------      matinv=k1o


!== end old SVD

      end subroutine svdm2



      subroutine svdm3(mat,matinv,n,svdtol)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      use negfmod, only: outinfo
      IMPLICIT NONE

      integer n,m,nwork,INFO,i
      double COMPLEX mat(n,n),matinv(n,n),k1ob(n,n),u(n,n),
     .    vt(n,n),sigma(n,n)
      double PRECISION svdtol,sv(n),svmax,condn
      double COMPLEX, ALLOCATABLE :: work(:)
      DOUBLE PRECISION, ALLOCATABLE :: rwork(:)
      DOUBLE COMPLEX,PARAMETER :: alpha=(1.D0,0D0)


      nwork=8
      allocate(work(nwork * n))
      allocate(rwork(5 * n))


!== begin old SVD

      CALL ZGESVD( 'A', 'A', n, n, mat, n,sv, u, n, vt,n,
     $                   work, nwork * n, rwork, INFO )
      if(info.ne.0)write(*,*)"infosvd=",info

      svmax=sv(1)
      m=n
      do i=1,n
        condn=svmax/sv(i)
!        write(*,*)"sigmak1=",DREAL(ene),sv(i),sv(i)/svmax,condn
        if(sv(i)/svmax.lt.svdtol)then
          m=i-1
          exit
        endif
      enddo

      if (outinfo) write(12347,*)"m,svmax=",m,svmax,svdtol
!      call f77flush

      sigma=0D0
      do i=1,n
        if(sv(i)/svmax.lt.svdtol)then
          sigma(i,i)= svdtol * svmax 
        else
          sigma(i,i)= sv(i) 
        endif
      enddo

      call ZGEMM('N', 'N', n, n, n,
     .            alpha, u, n, sigma, n, 0, k1ob, n)
      call ZGEMM('N', 'N', n, n, n,
     .            alpha, k1ob, n, vt, n, 0, mat, n)


      sigma=0D0
      do i=1,n
        if(sv(i)/svmax.lt.svdtol)then
          sigma(i,i)= 1D0/(svdtol * svmax )
        else
          sigma(i,i)= 1D0/sv(i) 
        endif
      enddo

      call ZGEMM('N', 'N', n, n, n,
     .    alpha, DCONJG(TRANSPOSE(vt)), n, sigma, n, 0, k1ob, n)
      call ZGEMM('N', 'N', n, n, n,
     .    alpha, k1ob, n, DCONJG(TRANSPOSE(u)), n, 0, matinv, n)

!== end old SVD

      end subroutine svdm3




