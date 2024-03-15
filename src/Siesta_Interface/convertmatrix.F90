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
! THE SUBROUTINE
!                   CONVERTMATRIXSIESTATOSMEAGOLGENERAL  
! IN THIS FILE IS LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
subroutine convertmatrixsiestatosmeagolgeneral(H,S,maxnh,numh,listhptr,listh,n1loc,xij,NspinRealInputMatrix,NspinComplexMatrix,N1,nmat,indxuo,no,kpoint,hgeneral,sgeneral,rhogeneral,ematgeneral,emforces,maxnelerow,nProcs,iProc,mattype)

  use mMatrixUtil
  use mTypes
  use mMPI_NEGF

  implicit none

  integer, intent(in)::nProcs,iProc,mattype
  integer maxnh,NspinRealInputMatrix,NspinComplexMatrix,n1,n1loc,nmat
  double precision, intent(in) :: H(maxnh,NspinRealInputMatrix),S(maxnh),xij(3,maxnh)
  type(matrixTypeGeneral),intent(inout) :: rhogeneral(NspinComplexMatrix),ematgeneral(NspinComplexMatrix)
  type(matrixTypeGeneral),intent(inout) :: hgeneral(NspinComplexMatrix),sgeneral
  logical, intent(in) :: emforces

  integer io,j,no,ind,jo,iuo,juo,indxuo(no),i,numh(n1loc),listhptr(n1loc),listh(maxnh),iio,ii,ispin
  integer ind2,is,maxnelerow,ind3,ind4
  double precision kxij,kpoint(3)
  double complex, allocatable :: hrow(:,:),srow(:)
  

  double complex sckxij
  DOUBLE COMPLEX, PARAMETER :: zi=(0.D0,1.D0)
  integer  mynode,nnodes,bnode,MPIerror,neletotal
  integer, allocatable :: nstartrow(:),nelerow(:),listj(:),nelenz(:)
  type(ioType) :: iout
  integer n1l,i1,i2,nele,nnz,istart,iend,nrows,neletarget,jproc,mpi_group

  if(mattype==2)then
    call convertmatrixsiestatosmeagolparallel(H,S,maxnh,numh,listhptr,listh,n1loc,xij,NspinRealInputMatrix,NspinComplexMatrix,N1,nmat,indxuo,no,kpoint,hgeneral,sgeneral,rhogeneral,ematgeneral,emforces,maxnelerow)
  elseif(mattype==3)then
    call convertmatrixsiestatosmeagolparallel2(H,S,maxnh,numh,listhptr,listh,n1loc,xij,NspinRealInputMatrix,NspinComplexMatrix,N1,indxuo,no,kpoint,hgeneral,sgeneral,rhogeneral,ematgeneral,emforces,maxnelerow,nProcs,iProc)
  endif

end subroutine convertmatrixsiestatosmeagolgeneral


