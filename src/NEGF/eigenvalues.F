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
!                   GEIGENVALUESLR,
!                   GEIGENVALUES3,
!                   GEIGENVALUES2,
!                   GEIGENVALUES1,
!                   GEIGENVALUES,
!                   GEIGENVALUESHERMITIAN  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!

      subroutine geigenvalueslr(mat,mat2,walpha,wbeta,n,vrg,vlg,info)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      use negfmod, only: outinfo
      IMPLICIT NONE

      INTEGER n,LWORK,INFO
      double complex mat(n,n),mat2(n,n),walpha(n),wbeta(n)
      double complex matbuf(n,n),mat2buf(n,n)
      double complex vlg(n,n),vrg(n,n),work(4 * n)
      double PRECISION rwork(8 * n)

      LWORK=4 * n
      
      matbuf=mat
      mat2buf=mat2
      call ZGGEV( 'V' , 'V', n, matbuf ,n, mat2buf, n,
     &    walpha, wbeta, vlg, n, vrg, n, WORK, LWORK ,
     &    RWORK, INFO )
      if(info.ne.0 .and. outinfo) write(12347,*)"info_c2=",info


      end subroutine geigenvalueslr
     





      subroutine geigenvalues3(mat,zv,n,vrg,vlg,info)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      IMPLICIT NONE

      INTEGER n,LWORK,INFO,ilo,ihi
      double complex mat(n,n),matbuf(n,n)
      double PRECISION scaleev(n),ABNRM,RCONDE(n),RCONDV(n)
      double complex vlg(n,n),vrg(n,n),work(4 * n),zv(n)
      double PRECISION rwork(2 * n)

      LWORK=4 * n
      
      matbuf=mat
      call ZGEEVX( 'P', 'V', 'V','N', N, matbuf, N, zv,
     &    vlg, N, vrg, N,ilo,ihi,scaleev,abnrm,
     &    rconde,rcondv, WORK, LWORK , RWORK, INFO )


      end subroutine geigenvalues3
     

      subroutine geigenvalues2(mat,zv,n,vrg,info)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      IMPLICIT NONE

      INTEGER n,LWORK,INFO,ilo,ihi
      double complex mat(n,n),matbuf(n,n)
      double PRECISION scaleev(n),ABNRM,RCONDE(n),RCONDV(n)
      double complex vlg(n,n),vrg(n,n),work(4 * n),zv(n)
      double PRECISION rwork(2 * n)

      LWORK=4 * n
      
      matbuf=mat
      call ZGEEVX( 'P', 'N', 'V','N', N, matbuf, N, zv,
     &    vlg, N, vrg, N,ilo,ihi,scaleev,abnrm,
     &    rconde,rcondv, WORK, LWORK , RWORK, INFO )


      end subroutine geigenvalues2



      subroutine geigenvalues1(mat,zv,n)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      use negfmod, only: outinfo
      IMPLICIT NONE

      INTEGER n,LWORK,INFO,ilo,ihi
      double complex mat(n,n),matbuf(n,n)
      double PRECISION scaleev(n),ABNRM,RCONDE(n),RCONDV(n)
      double complex vlg(n,n),vrg(n,n),work(4 * n),zv(n)
      double PRECISION rwork(2 * n)

      LWORK=4 * n
      

      matbuf=mat
      call ZGEEVX( 'P', 'N', 'N','N', N, matbuf, N, zv,
     &    vlg, N, vrg, N,ilo,ihi,scaleev,abnrm,
     &    rconde,rcondv, WORK, LWORK , RWORK, INFO )
      if(info.ne.0 .and. outinfo) write(12347,*)"info_ce1=",info


      end subroutine geigenvalues1
     



      subroutine geigenvalues(mat,mat2,walpha,wbeta,n,vrg,info)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      use negfmod, only: outinfo
      IMPLICIT NONE

      INTEGER n,LWORK,INFO
      double complex mat(n,n),mat2(n,n),walpha(n),wbeta(n)
      double complex matbuf(n,n),mat2buf(n,n)
      double complex vlg(n,n),vrg(n,n),work(4 * n)
      double PRECISION rwork(8 * n)

      LWORK=4 * n
      
      matbuf=mat
      mat2buf=mat2
      call ZGGEV( 'N' , 'V', n, matbuf ,n, mat2buf, n,
     &    walpha, wbeta, vlg, n, vrg, n, WORK, LWORK ,
     &    RWORK, INFO )
      if(info.ne.0 .and. outinfo) write(12347,*)"info_c2b=",info


      end subroutine geigenvalues
     

      subroutine geigenvalueshermitian(mat,eig,n,vr,info)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

      IMPLICIT NONE

      integer, intent(in) :: n
      double complex, intent(in) :: mat(n,n)
      double precision, intent(out) :: eig(n)
      double complex, intent(out) :: vr(n,n)
      integer, intent(out) :: info
     
      integer lwork
      double complex work(4 * n)
      double precision rwork(3 * n-2)

      LWORK=4 * n
      
      vr=mat
      call ZHEEV('V', 'L' , n, vr, n, eig, WORK, LWORK, RWORK,
     $                  INFO )

      end subroutine geigenvalueshermitian



