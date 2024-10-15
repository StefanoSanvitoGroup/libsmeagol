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
!                   UPDATERHONEQ,
!                   UPDATERHONEQ_NC,
!                   UPDATERHO_NC,
!                   UPDATERHOSPARSE,
!                   UPDATERHODENSE_NC,
!                   UPDATERHODENSE,
!                   UPDATERHO1  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
subroutine UpdateRhoNEQ(nnz,n1,nl,nr,nlead,q,j,b,gf1,gf2,const, set_rho_boundary)

   implicit none

   integer,intent(in):: nnz,n1,nl,nr,q(n1+1),j(nnz),nlead
   double complex, intent(in):: const
   double complex, intent(in):: gf1(nlead,n1),gf2(nlead,n1)
   logical, intent(in):: set_rho_boundary
   double complex, intent(inout):: b(nnz)

   integer ii,ind,jj,i1
   double complex gfadd

!$omp parallel do default(shared) private(ii,ind,jj,gfadd) schedule(dynamic)
   do ii=1,n1
     do ind= q(ii),q(ii+1)-1
       jj=j(ind)
       if ((((II .GT. NL) .AND. (II .LE. N1-NR)) .OR.((JJ .GT. NL) .AND. (JJ .LE. N1-NR))).or.set_rho_boundary) THEN

         gfadd=0D0
         do i1=1,nlead
           gfadd=gfadd+gf2(i1,ii)*gf1(i1,jj)
         enddo

         b(ind)=b(ind)+ const * gfadd
       
       endif
     enddo
   enddo
!$omp end parallel do

end subroutine UpdateRhoNEQ

subroutine UpdateRhoNEQ_nc(nnz,n12,nl,nr,nlead,q,j,b1,b2,b3,b4,gf1,gf2,const, set_rho_boundary)

   implicit none

   integer nnz,n12,nl,nr,q(n12/2+1),j(nnz),nlead
   double complex const
   double complex, intent(inout) :: b1(nnz),b2(nnz),b3(nnz),b4(nnz)
   double complex gf1(nlead,n12),gf2(nlead,n12)
   logical, intent(in):: set_rho_boundary

   integer ii,ind,jj,i1,n1
   double complex gfadd

   n1=n12/2

   do ii=1,n1
     do ind= q(ii),q(ii+1)-1


!       if(s(ind)==0.0D0)cycle

       jj=j(ind)

       if ((((II .GT. NL) .AND. (II .LE. N1-NR)) .OR.((JJ .GT. NL) .AND. (JJ .LE. N1-NR))).or.set_rho_boundary) THEN

         gfadd=0D0
         do i1=1,nlead
           gfadd=gfadd+gf2(i1,ii)*gf1(i1,jj)
         enddo
         b1(ind)=b1(ind)+ const * gfadd


         gfadd=0D0
         do i1=1,nlead
           gfadd=gfadd+gf2(i1,n1+ii)*gf1(i1,n1+jj)
         enddo
         b2(ind)=b2(ind)+ const * gfadd


         gfadd=0D0
         do i1=1,nlead
           gfadd=gfadd+gf2(i1,ii)*gf1(i1,n1+jj)
         enddo
         b3(ind)=b3(ind)+ const * gfadd


         gfadd=0D0
         do i1=1,nlead
           gfadd=gfadd+gf2(i1,n1+ii)*gf1(i1,jj)
         enddo
         b4(ind)=b4(ind)+ const * gfadd
       
       endif
     enddo
   enddo

end subroutine UpdateRhoNEQ_nc



  SUBROUTINE updaterho_nc(rhogeneral,ematgeneral,emforces,ispin,nspin,gf,nl,nr,gfmattype,weightc,cl,cr,weightrho,ene, set_rho_boundary)

    use mTypes

    implicit none
    logical, intent(in):: set_rho_boundary
    integer, intent(in) :: ispin,nspin
    type(matrixTypeGeneral), intent(inout) :: rhogeneral(nspin),ematgeneral(nspin)
    logical, intent(in) :: emforces
    type(matrixTypeGeneral),intent(in) :: gf
    integer nl,nr,gfmattype,ii,jj,n1,ind2,ind
    double complex weightc,cl,cr,ene
    double precision weightrho

    if(gfmattype.eq.0)then
      if(nspin<=2)then
        call updaterhodense(rhogeneral(ispin),ematgeneral(ispin),emforces,gf,nl,nr,gfmattype,weightc,cl,cr,weightrho,ene,set_rho_boundary)
      else
        call updaterhodense_nc(rhogeneral,ematgeneral,emforces,nspin,gf,rhogeneral(1)%irows,nl/2,nr/2,gfmattype,weightc,cl,cr,weightrho,ene,set_rho_boundary)
      endif
    elseif(gfmattype.eq.2)then
      call updaterhosparse(rhogeneral(ispin),ematgeneral(ispin),emforces,gf,nl,nr,gfmattype,weightc,cl,cr,weightrho,ene,set_rho_boundary)
    endif


  end SUBROUTINE updaterho_nc


  SUBROUTINE updaterhosparse(rhogenerals,ematgenerals,emforces,gf,nl,nr,gfmattype,weightc,cl,cr,weightrho,ene, set_rho_boundary)

    use mMatrixUtil
    use mTypes

    implicit none
    logical, intent(in):: set_rho_boundary
    type(matrixTypeGeneral) :: rhogenerals,ematgenerals,gf,gfdagger
    logical, intent(in) :: emforces
    integer nl,nr,gfmattype,ii,jj,n1,ind2,ind
    double complex weightc,cl,cr,ene
!    double complex mat1(rhogenerals%irows,rhogenerals%irows),mat2(rhogenerals%irows,rhogenerals%irows)
    double complex gfij,drhoij,gfji
    double precision weightrho
    type(ioType) :: io
    DOUBLE COMPLEX, PARAMETER :: zi=(0.D0,1.D0)
    DOUBLE PRECISION, PARAMETER :: PI=3.141592654D0
    double complex w(rhogenerals%iCols)
    integer idxrow(rhogenerals%iCols),nj

    io%isDebug=.false.

    n1=rhogenerals%irows

    call AllocateMatrixGeneral(n1,n1,gf%matSparse%nnz,gfmattype,gfdagger,"updaterhosparse", io)
    call mathermitianCRS(gf%MatSparse,gfdagger%MatSparse)

!    mat1=0D0
!    do ii=1,n1
!      do ind=gf%matSparse%q(ii),gf%matSparse%q(ii+1)-1
!        jj=gf%matSparse%j(ind)
!        mat1(ii,jj)=gf%matSparse%b(ind)
!      enddo
!    enddo
!
!    mat2=0D0
!    do ii=1,n1
!      do ind=gfdagger%matSparse%q(ii),gfdagger%matSparse%q(ii+1)-1
!        jj=gfdagger%matSparse%j(ind)
!        mat2(ii,jj)=gfdagger%matSparse%b(ind)
!      enddo
!    enddo
!
!    write(*,*)"gd-gd=",maxval(abs(DCONJG(TRANSPOSE(mat1))-mat2))


    w=0D0
    DO II=1,N1

      nj=0
      do ind=gf%matSparse%q(ii),gf%matSparse%q(ii+1)-1
        jj=gf%matSparse%j(ind)
!        w(jj)=gf%matSparse%b(ind)
        w(jj)=(-zi/(2.0D0*PI))*weightc*((1D0-weightrho)*cl +   weightrho * cr)*gf%matSparse%b(ind)

        nj=nj+1
        idxrow(nj)=jj
      enddo

      do ind=rhogenerals%matSparse%q(ii),rhogenerals%matSparse%q(ii+1)-1
        jj=rhogenerals%matSparse%j(ind)
        if ((((II .GT. NL) .AND. (II .LE. N1-NR)) .OR.((JJ .GT. NL) .AND. (JJ .LE. N1-NR))).or.set_rho_boundary) THEN
          rhogenerals%matSparse%b(ind)=rhogenerals%matSparse%b(ind)+w(jj)
          if(emforces)then
            ematgenerals%matSparse%b(ind)=ematgenerals%matSparse%b(ind)+w(jj)*ene
          endif
        endif
      enddo

      do jj=1,nj
        w(idxrow(jj))=0D0
      enddo

      nj=0
      do ind=gfdagger%matSparse%q(ii),gfdagger%matSparse%q(ii+1)-1
        jj=gfdagger%matSparse%j(ind)
        w(jj)=(-zi/(2.0D0*PI))*DCONJG(weightc)*((1D0-weightrho)*(DCONJG(cl)) + weightrho *(DCONJG(cr)))*(-gfdagger%matSparse%b(ind))

        nj=nj+1
        idxrow(nj)=jj
      enddo

      do ind=rhogenerals%matSparse%q(ii),rhogenerals%matSparse%q(ii+1)-1
        jj=rhogenerals%matSparse%j(ind)
        if ((((II .GT. NL) .AND. (II .LE. N1-NR)) .OR.((JJ .GT. NL) .AND. (JJ .LE. N1-NR))).or.set_rho_boundary) THEN
          rhogenerals%matSparse%b(ind)=rhogenerals%matSparse%b(ind)+w(jj)
          if(emforces)then
            ematgenerals%matSparse%b(ind)=ematgenerals%matSparse%b(ind)+w(jj)*DCONJG(ene)
          endif
        endif
      enddo

      do jj=1,nj
        w(idxrow(jj))=0D0
      enddo

    ENDDO
 
    call DestroyMatrixGeneral(gfdagger,"updaterhosparse",io)

  end SUBROUTINE updaterhosparse

  SUBROUTINE updaterhodense_nc(rhogeneral,ematgeneral,emforces,nspin,gf,n1,nl,nr,gfmattype,weightc,cl,cr,weightrho,ene, set_rho_boundary)

    use mTypes
    use mMatrixUtil, only:WriteMatrixSparse

    implicit none
    logical, intent(in):: set_rho_boundary
    integer, intent(in) :: nspin,n1,nl,nr
    type(matrixTypeGeneral) :: rhogeneral(nspin),ematgeneral(nspin),gf
    logical, intent(in) :: emforces
    integer gfmattype,ii,jj,ind2,ind
    double complex weightc,cl,cr,ene
    double complex gfij,drhoij,gfji,c1,c2
    double precision weightrho
    DOUBLE COMPLEX, PARAMETER :: zi=(0.D0,1.D0)
    DOUBLE PRECISION, PARAMETER :: PI=3.141592654D0
    double complex, allocatable :: mat(:,:)


    c1=(-zi/(2.0D0*PI))*weightc*((1D0-weightrho)*(cl) + weightrho *(cr))
    c2=-DCONJG(c1)

    do ii=1,n1
      do ind=rhogeneral(1)%matSparse%q(ii),rhogeneral(1)%matSparse%q(ii+1)-1
        jj =rhogeneral(1)%matSparse%j(ind)
        if ((((II .GT. NL) .AND. (II .LE. N1-NR)) .OR.((JJ .GT. NL) .AND. (JJ .LE. N1-NR))).or.set_rho_boundary) THEN

         gfij=gf%matdense%a(II,JJ)
         gfji=gf%matdense%a(JJ,II)
         rhogeneral(1)%matSparse%b(ind)=rhogeneral(1)%matSparse%b(ind)+c1*gfij-c2*DCONJG(gfji)
         if(emforces)then
           ematgeneral(1)%matSparse%b(ind)=ematgeneral(1)%matSparse%b(ind)+c1*ene*gfij-c2*DCONJG(ene)*DCONJG(gfji)
         endif

         gfij=gf%matdense%a(n1+II,n1+JJ)
         gfji=gf%matdense%a(n1+JJ,n1+II)
         rhogeneral(2)%matSparse%b(ind)=rhogeneral(2)%matSparse%b(ind)+c1*gfij-c2*DCONJG(gfji)
         if(emforces)then
           ematgeneral(2)%matSparse%b(ind)=ematgeneral(2)%matSparse%b(ind)+c1*ene*gfij-c2*DCONJG(ene)*DCONJG(gfji)
         endif


         gfij=gf%matdense%a(II,n1+JJ)
         gfji=gf%matdense%a(n1+JJ,II)
         rhogeneral(3)%matSparse%b(ind)=rhogeneral(3)%matSparse%b(ind)+c1*gfij-c2*DCONJG(gfji)
         if(emforces)then
           ematgeneral(3)%matSparse%b(ind)=ematgeneral(3)%matSparse%b(ind)+c1*ene*gfij-c2*DCONJG(ene)*DCONJG(gfji)
         endif


         gfij=gf%matdense%a(n1+II,JJ)
         gfji=gf%matdense%a(JJ,n1+II)
         rhogeneral(4)%matSparse%b(ind)=rhogeneral(4)%matSparse%b(ind)+c1*gfij-c2*DCONJG(gfji)
         if(emforces)then
           ematgeneral(4)%matSparse%b(ind)=ematgeneral(4)%matSparse%b(ind)+c1*ene*gfij-c2*DCONJG(ene)*DCONJG(gfji)
         endif

        ENDIF
      ENDDO
    ENDDO

!  call WriteMatrixSparse(rhogeneral(1),"rho1")
!  call WriteMatrixSparse(rhogeneral(2),"rho2")
!  call WriteMatrixSparse(rhogeneral(3),"rho3")
!  call WriteMatrixSparse(rhogeneral(4),"rho4")

  end SUBROUTINE updaterhodense_nc


  SUBROUTINE updaterhodense(rhogenerals,ematgeneral,emforces,gf,nl,nr,gfmattype,weightc,cl,cr,weightrho,ene, set_rho_boundary)

    use mTypes

    implicit none
    logical, intent(in):: set_rho_boundary
    logical, intent(in) :: emforces
    type(matrixTypeGeneral) :: rhogenerals,ematgeneral,gf
    integer nl,nr,gfmattype,ii,jj,n1,ind2,ind, nnz
    double complex weightc,cl,cr,ene
    double complex gfij,drhoij,denematij,gfji,c1,c2
    double precision weightrho
    DOUBLE COMPLEX, PARAMETER :: zi=(0.D0,1.D0)
    DOUBLE PRECISION, PARAMETER :: PI=3.141592654D0

    n1=rhogenerals%irows
    nnz=rhogenerals%matSparse%nnz
!    write(*,*)"weights=",weightc,cl,cr
!    write(*,*)"weightcl=",cl
!    write(*,*)"weightcr=",cr
    c1=(-zi/(2.0D0*PI))*weightc*((1D0-weightrho)*(cl) + weightrho *(cr))
    c2=-DCONJG(c1)

    do ii=1,n1
      do ind=rhogenerals%matSparse%q(ii),rhogenerals%matSparse%q(ii+1)-1
        JJ=rhogenerals%matSparse%j(ind)
        if ((((II .GT. NL) .AND. (II .LE. N1-NR)) .OR.((JJ .GT. NL) .AND. (JJ .LE. N1-NR))).or.set_rho_boundary) THEN

          gfij=gf%matdense%a(II,JJ)
          gfji=gf%matdense%a(JJ,II)
          rhogenerals%matSparse%b(ind)=rhogenerals%matSparse%b(ind)+c1*gfij-c2*DCONJG(gfji)

          if(emforces)then
            ematgeneral%matSparse%b(ind)=ematgeneral%matSparse%b(ind)+c1*ene * gfij -c2*DCONJG(ene) * DCONJG(gfji)
          endif
        ENDIF
      enddo
    enddo
 

  end SUBROUTINE updaterhodense



