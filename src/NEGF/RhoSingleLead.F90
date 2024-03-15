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
!                   RHOSINGLELEAD,
!                   RHOSINGLELEADGAMMA  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!

  subroutine RhoSingleLead(ei,GF_iter1l,n1,nl,nr,nlead,sigmal,rhogeneralp,ematgeneralp,nspin,weightc,ispin)

    use mMPI_NEGF
    use negfmod
    use mTypes
    use mMatrixUtil
    use mONInterface

    implicit none

      
  integer, intent(in) :: n1,nlead,nl,nr,nspin,ispin
  type(matrixTypeGeneral), intent(inout) :: rhogeneralp(nspin)
  type(matrixTypeGeneral), intent(inout) :: ematgeneralp(nspin)
  DOUBLE COMPLEX,intent(in) :: sigmal(nlead,nlead),weightc
  DOUBLE COMPLEX, intent(in) :: GF_iter1l(n1,nlead),ei

  DOUBLE COMPLEX, ALLOCATABLE :: gammal(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: GF_iter2l(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: gf1(:,:),gf2(:,:)
  double complex, allocatable :: b(:)
  double complex, allocatable :: bnc(:,:)
  DOUBLE PRECISION, PARAMETER :: PI=3.141592654D0
  DOUBLE COMPLEX, PARAMETER :: zi=(0.D0,1.D0)
  integer*4:: sc_0,sc_1,sc_r,sc_m
  integer is

  if(emtimings)then
    CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
  endif



  ALLOCATE(gammal(nlead,nlead))
  gammal=zi*(sigmal-DCONJG(TRANSPOSE(sigmal)))

  ALLOCATE(GF_iter2l(N1,nlead))
  CALL ZGEMM('N','N',N1,nlead,nlead,(1.D0,0.D0),GF_iter1l,N1, gammal,nlead,(0.D0,0.D0),GF_iter2l,N1)
  deallocate(gammal)
            
  ALLOCATE(gf2(nlead,N1))
  gf2=transpose(GF_iter2l)
  deallocate(GF_iter2l)
  ALLOCATE(gf1(nlead,N1))
  gf1=DCONJG(transpose(GF_iter1l))

  if(emtimings)then
    CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
    write(12347,'(A,f12.6)') 'TbeforeUpdateRho',(sc_1-sc_0)*1.0d0/sc_r
    CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
  endif

  if(nspin<=2)then
    allocate(b(rhogeneralp(1)%matSparse%nnz))
    b=0.0D0

    call UpdateRhoNEQ(rhogeneralp(1)%matSparse%nnz,n1,nl,nr,nlead,rhogeneralp(1)%matSparse%q,rhogeneralp(1)%matSparse%j,b,gf1,gf2, weightc,.true.)

    rhogeneralp(ispin)%matSparse%b(:)=rhogeneralp(ispin)%matSparse%b(:)+b(:)
    if(emforces) ematgeneralp(ispin)%matSparse%b(:)= ematgeneralp(ispin)%matSparse%b(:)+Ei * b(:)
    deallocate(b)
  else
    allocate(bnc(rhogeneralp(1)%matSparse%nnz,4))
    bnc=0.0D0

    call UpdateRhoNEQ_nc(rhogeneralp(1)%matSparse%nnz,n1, nl/2,nr/2,nlead, rhogeneralp(1)%matSparse%q,rhogeneralp(1)%matSparse%j,bnc(:,1),bnc(:,2), bnc(:,3),bnc(:,4), gf1,gf2,weightc ,.true.)

    do is=1,nspin
      rhogeneralp(is)%matSparse%b= rhogeneralp(is)%matSparse%b+bnc(:,is)
    enddo

    if(emforces) then
      do is=1,nspin
        ematgeneralp(is)%matSparse%b= ematgeneralp(is)%matSparse%b+ei* bnc(:,is)
      enddo
    endif

    deallocate(bnc)
  endif
  deallocate(gf1,gf2)

  if(emtimings)then
    CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
    write(12347,'(A,f12.6)') 'TUpdateRho',(sc_1-sc_0)*1.0d0/sc_r
  endif

  end subroutine RhoSingleLead


  subroutine RhoSingleLeadGamma(GF_iter1l,n1,nl,nr,nlead,gammal,rhogeneralp,nspin,weightc,ispin)

    use mMPI_NEGF
    use negfmod
    use mTypes
    use mMatrixUtil
    use mONInterface

    implicit none

      
  integer, intent(in) :: n1,nlead,nl,nr,nspin,ispin
  type(matrixTypeGeneral), intent(inout) :: rhogeneralp(nspin)
  DOUBLE COMPLEX,intent(in) :: gammal(nlead,nlead),weightc
  DOUBLE COMPLEX, intent(in) :: GF_iter1l(n1,nlead)

  DOUBLE COMPLEX, ALLOCATABLE :: GF_iter2l(:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: gf1(:,:),gf2(:,:)
  double complex, allocatable :: b(:)
  DOUBLE PRECISION, PARAMETER :: PI=3.141592654D0
  DOUBLE COMPLEX, PARAMETER :: zi=(0.D0,1.D0)


  ALLOCATE(GF_iter2l(N1,nlead))
  CALL ZGEMM('N','N',N1,nlead,nlead,(1.D0,0.D0),GF_iter1l,N1, gammal,nlead,(0.D0,0.D0),GF_iter2l,N1)
            
  ALLOCATE(gf2(nlead,N1))
  gf2=transpose(GF_iter2l)
  deallocate(GF_iter2l)
  ALLOCATE(gf1(nlead,N1))
  gf1=DCONJG(transpose(GF_iter1l))

  if(nspin<=2)then
    allocate(b(rhogeneralp(1)%matSparse%nnz))
    b=0.0D0

    call UpdateRhoNEQ(rhogeneralp(1)%matSparse%nnz,n1,nl,nr,nlead,rhogeneralp(1)%matSparse%q,rhogeneralp(1)%matSparse%j,b,gf1,gf2, weightc * (1.D0,0D0)/(2.D0*PI),.true.)

    rhogeneralp(ispin)%matSparse%b(:)=rhogeneralp(ispin)%matSparse%b(:)+b(:)
!    if(emforces) ematgeneralp(ispin)%matSparse%b(:)= ematgeneralp(ispin)%matSparse%b(:)+Ei * b(:)
    deallocate(b)
  else
    call UpdateRhoNEQ_nc(rhogeneralp(1)%matSparse%nnz,n1, nl/2,nr/2,nlead, rhogeneralp(1)%matSparse%q,rhogeneralp(1)%matSparse%j,rhogeneralp(1)%matSparse%b,rhogeneralp(2)%matSparse%b, rhogeneralp(3)%matSparse%b,rhogeneralp(4)%matSparse%b, gf1,gf2,weightc * (1.D0,0D0)/(2.D0*PI),.true.)
  endif
  deallocate(gf1,gf2)
 
  end subroutine RhoSingleLeadGamma

