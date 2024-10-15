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
!                   MERGE_RHOK  
! IN THIS FILE IS LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!

subroutine merge_rhok(NspinRealInputMatrix,rhogeneral)
  use mMatrixUtil
  use mTypes
  use mMPI_NEGF

  implicit none

  integer, intent(in)::NspinRealInputMatrix
  type(matrixTypeGeneral),intent(inout) :: rhogeneral(NspinRealInputMatrix)
  double complex, allocatable :: b(:)
  double complex b2

  integer i,ispin,MPIerror


#ifdef MPI
  do ispin=1,NspinRealInputMatrix


!    do i=1,rhogeneral(1)%matSparseP%matSparse%nnz
!      rhogeneral(ispin)%matSparseP%matSparse%b(i)=mynode_negf
!    enddo


!    write(12347,*)"ispino=",ispin
    if(mynode_groupk==0)then
      allocate(b(rhogeneral(1)%matSparseP%matSparse%nnz))
      b=rhogeneral(ispin)%matSparseP%matSparse%b
      call MPI_REDUCE(b,rhogeneral(ispin)%matSparseP%matSparse%b,rhogeneral(1)%matSparseP%matSparse%nnz,DAT_dcomplex,MPI_SUM,0,groupk_comm,MPIerror)
      deallocate(b)
    else
!      allocate(b(1))
      call MPI_REDUCE(rhogeneral(ispin)%matSparseP%matSparse%b,b2,rhogeneral(1)%matSparseP%matSparse%nnz,DAT_dcomplex,MPI_SUM,0,groupk_comm,MPIerror)
!      deallocate(b)
    endif
!    do i=1,rhogeneral(1)%matSparseP%matSparse%nnz
!      write(12347,*)"rhogeneral(ispin)%matSparseP%matSparse%b(i)=",rhogeneral(ispin)%matSparseP%matSparse%b(i),mynode_groupk,mynode_negf
!    enddo



  enddo
#endif

 
end subroutine merge_rhok




