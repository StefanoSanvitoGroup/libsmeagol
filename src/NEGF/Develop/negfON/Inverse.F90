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
!                   ALLOCATESPACE,
!                   DEALLOCATESPACE,
!                   GAUSSIANFORWARDELIMINATION,
!                   GAUSSIANFORWARDELIMINATIONSPARSE,
!                   GAUSSIANBACKWARDELIMINATION,
!                   GAUSSIANBACKWARDELIMINATIONSPARSESTEP,
!                   GAUSSIANBACKWARDELIMINATIONSPARSE,
!                   INVERSEDIAGONALBLOCKS,
!                   INVERSEDIAGONALBLOCKSSPARSE,
!                   INVERSEDIAGONALOFFDIAGONALBLOCKSSPARSE,
!                   INVERSEOFFDIAGONALBLOCKS,
!                   INVERSEOFFDIAGONALBLOCKSSPARSE,
!                   INVERSEFIRSTCOLUMNBLOCKS,
!                   INVERSE1NBLOCKSSPARSE,
!                   INVERSEFIRSTCOLUMNBLOCKSSPARSE,
!                   INVERSELASTCOLUMNBLOCKS,
!                   INVERSELASTCOLUMNBLOCKSSPARSE,
!                   CHECKIFIRSTCOLUMNBLOCKS,
!                   CHECKILASTCOLUMNBLOCKS  
! AND
! THE MODULE
!                   MINVERSE  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
!> \brief functions and routines needed to compute the sparse inverse
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 1st of April 2009
!> \remarks
!> \todo
!> \bug
module mInverse
  use mConstants
  use mTypes
  use mMatrixUtil
  implicit none
  private
!
  public :: InverseDiagonalBlocks
  public :: InverseOffDiagonalBlocks
  public :: InverseFirstColumnBlocks
  public :: InverseLastColumnBlocks
  public :: AllocateSpace
  public :: DeallocateSpace
  public :: CheckIFirstColumnBlocks
  public :: CheckILastColumnBlocks
  public :: InverseDiagonalBlocksSparse
  public :: InverseDiagonalOffdiagonalBlocksSparse
  public :: InverseOffDiagonalBlocksSparse
  public :: InverseLastColumnBlocksSparse
  public :: InverseFirstColumnBlocksSparse
  public :: Inverse1NBlocksSparse
!
 contains

!> \brief allocate
!> \details job selects what space should be allocated according to the calculations that will be requested later \n
!> job is the decimal for the following binary shorthand "FLD" D-Diagonal (always present), L- last column blocks, F First column blocks \n
!> 1: 001 only diagonal \n
!> 5: 101  first column and diagonal
!> 3: 011 last columns and diagonal
!> 7: 111 first column, last column and diagonal
  subroutine AllocateSpace(iBlocks,sigmaL,sigmaR,g0,g1,gm1,h0,h1,hm1,job,io,gi1,gin,ml,mm)
    character (len=*), parameter :: sMyName = "AllocateSpace"
    type(matrixType), intent(inout) :: h0(:),hm1(:),h1(:)
    type(matrixType), intent(inout), allocatable :: sigmaL(:),sigmaR(:),g0(:),g1(:),gm1(:)
    type(matrixType), intent(inout), optional, allocatable :: gi1(:),gin(:),ml(:),mm(:)
    type(ioType), intent(inout) :: io
    integer, intent(in) :: iBlocks,job

    integer :: i

    call AllocateArray(iBlocks,sigmaL,sMyName,io)
    call AllocateArray(iBlocks,sigmaR,sMyName,io)
    call AllocateArray(iBlocks,g0,sMyName,io)

    if ((job==7) .or. (job==5)) then
      call AllocateArray(iBlocks,gi1,sMyName,io)
      call AllocateArray(iBlocks-1,mm,sMyName,io)
    endif
    if ((job==7) .or. (job==3)) then
      call AllocateArray(iBlocks,gin,sMyName,io)
      call AllocateArray(iBlocks-1,ml,sMyName,io)
    endif
    if ((job==7) .or. (job==4)) then
      call AllocateArray(iBlocks-1,g1,sMyName,io)
      call AllocateArray(iBlocks-1,gm1,sMyName,io)
    endif

    do i=1,iBlocks
      call AllocateMatrix(h0(i)%iRows,h0(i)%iCols,h0(i)%iHorz,h0(i)%iVert,sigmaL(i),sMyName,io)
      sigmaL(i)%a = kczero
      call AllocateMatrix(h0(i)%iRows,h0(i)%iCols,h0(i)%iHorz,h0(i)%iVert,sigmaR(i),sMyName,io)
      sigmaR(i)%a = kczero
      call AllocateMatrix(h0(i)%iRows,h0(i)%iCols,h0(i)%iHorz,h0(i)%iVert,g0(i),sMyName,io)
      g0(i)%a = kczero
      if ((job==7) .or. (job==5)) then
        call AllocateMatrix(h0(i)%iRows,h0(1)%iCols,h0(1)%iHorz,h0(i)%iVert,gi1(i),sMyName,io)
        g0(i)%a = kczero
        if (i<=iBlocks-1) then
          call AllocateMatrix(h0(i)%iRows,h0(i+1)%iCols,h0(i)%iHorz,h0(1)%iVert,mm(i),sMyName,io)
          mm(i)%a = kczero
        endif
      endif

      if ((job==7) .or. (job==3)) then
        call AllocateMatrix(h0(i)%iRows,h0(iBlocks)%iCols,h0(iBlocks)%iHorz,h0(i)%iVert,gin(i),sMyName,io)
        g0(i)%a = kczero
        if (i>=2) then
          call AllocateMatrix(h0(i)%iRows,h0(i-1)%iCols,h0(i)%iHorz,h0(1)%iVert,ml(i-1),sMyName,io)
          ml(i-1)%a = kczero
        endif
      endif
      if ((job==7) .or. (job==4)) then
        if(i<iBlocks)then
          call AllocateMatrix(h1(i)%iRows,h1(i)%iCols,h1(i)%iHorz,h1(i)%iVert,g1(i),sMyName,io)
          call AllocateMatrix(hm1(i)%iRows,hm1(i)%iCols,hm1(i)%iHorz,hm1(i)%iVert,gm1(i),sMyName,io)
        endif
      endif
    enddo
  end Subroutine AllocateSpace

!> \brief deallocate
!> \details job selects what space should be deallocated according to the calculations that where requested before \n
!> job is the decimal for the following binary shorthand "FLD" D-Diagonal (always present), L- last column blocks, F First column blocks \n
!> 1: 001 only diagonal \n
!> 5: 101  first column and diagonal
!> 3: 011 last columns and diagonal
!> 7: 111 first column, last column and diagonal
  subroutine DeallocateSpace(iBlocks,sigmaL,sigmaR,g0,g1,gm1,h0,h1,hm1,job,io,gi1,gin,ml,mm)
    character (len=*), parameter :: sMyName = "DeallocateSpace"
    type(matrixType), intent(inout),allocatable :: h0(:),hm1(:),h1(:)
    type(matrixType), intent(inout),allocatable :: sigmaL(:),sigmaR(:),g0(:),g1(:),gm1(:)
    type(matrixType), intent(inout),allocatable, optional :: gi1(:),gin(:),ml(:),mm(:)
    type(ioType), intent(inout) :: io
    integer, intent(in) :: iBlocks,job

    integer :: i

    do i=1,iBlocks-1
      call DestroyMatrix(h0(i),sMyName,io)
      call DestroyMatrix(h1(i),sMyName,io)
      call DestroyMatrix(hm1(i),sMyName,io)
      call DestroyMatrix(g0(i),sMyName,io)
      call DestroyMatrix(sigmaL(i),sMyName,io)
      call DestroyMatrix(sigmaR(i),sMyName,io)
      if ((job==7) .or. (job==5)) then
        call DestroyMatrix(gi1(i),sMyName,io)
        call DestroyMatrix(mm(i),sMyName,io)
      endif
      if ((job==7) .or. (job==3)) then
        call DestroyMatrix(gin(i),sMyName,io)
        call DestroyMatrix(ml(i),sMyName,io)
      endif
      if ((job==7) .or. (job==4)) then
        call DestroyMatrix(g1(i),sMyName,io)
        call DestroyMatrix(gm1(i),sMyName,io)
      endif
    enddo
    call DestroyMatrix(h0(iBlocks),sMyName,io)
    call DestroyMatrix(g0(iBlocks),sMyName,io)
    call DestroyMatrix(sigmaL(iBlocks),sMyName,io)
    call DestroyMatrix(sigmaR(iBlocks),sMyName,io)
    call DestroyArray(g0,sMyName,io)
    call DestroyArray(sigmaL,sMyName,io)
    call DestroyArray(sigmaR,sMyName,io)
    call DestroyArray(h0,sMyName,io)
    call DestroyArray(h1,sMyName,io)
    call DestroyArray(hm1,sMyName,io)
    if ((job==7) .or. (job==5)) then
      call DestroyMatrix(gi1(iBlocks),sMyName,io)
      call DestroyArray(gi1,sMyName,io)
      call DestroyArray(ml,sMyName,io)
    endif
    if ((job==7) .or. (job==3)) then
      call DestroyMatrix(gin(iBlocks),sMyName,io)
      call DestroyArray(mm,sMyName,io)
      call DestroyArray(gin,sMyName,io)
    endif
    if ((job==7) .or. (job==4)) then
      call DestroyArray(g1,sMyName,io)
      call DestroyArray(gm1,sMyName,io)
    endif
  end subroutine DeallocateSpace


!> \brief forward gaussian elimination
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 14th of April, 2009
!> \param h0(:),hm1(:),h1(:) matrixType the blocks of the tridiagonal matrix
!> \param sigmaL(:) matrixType block matrices needed for the forward elimination
!> \param io ioType, structure that has the I/O info
!> \param iBlocks integer number of diagonal blocks
  subroutine GaussianForwardElimination (h0,hm1,h1,sigmaL,iBlocks,io)
    character (len=*), parameter :: sMyName = "GaussianForwardElimination"
    type(matrixType), intent(inout) :: h0(:),hm1(:),h1(:),sigmaL(:)
    type(ioType), intent(inout) :: io
    integer, intent(in) :: iBlocks

    integer :: i,m,n
    integer :: t1,t2,ta1,ta2
    type(matrixType) :: m1,m2


!    call system_clock(t1)
!    call sleep(1)
!    call system_clock(t2)
!    write(*,*)"tsystem=",t1,t2
!    write(*,*)"tsystem2-tsystem1 (1 second)=",(1D0 * (t2-t1))/1d4

    call system_clock(ta1)
    do i=1,iBlocks-1

      call system_clock(t1)
!      write(*,*)"gfe=",i

      n=h0(i)%iRows
      m=h0(i+1)%iRows
!      write(*,*)"mn=",i,n,m
      call AllocateMatrix(n,n,m1,sMyName,io)
      call AllocateMatrix(m,n,m2,sMyName,io)
!       m1%a=h0(i)%a-sigmaL(i)%a
      call MatrixAdd(m1,kcone,h0(i),-kcone,sigmaL(i),io)       
      call system_clock(t2)
!      write(*,*)"time_setup_gfe=",(1D0 * (t2-t1))/1d4

      call system_clock(t1)
      call Inversematrix(m1,io)
      call system_clock(t2)
!      write(*,*)"time_inverse_gfe=",(1D0 * (t2-t1))/1d4

      call system_clock(t1)
      call ProductCeAxB(m2,kcone,hm1(i),m1,io)
      call system_clock(t2)
!      write(*,*)"time_mmult1_gfe=",(1D0 * (t2-t1))/1d4

      call system_clock(t1)
      call ProductCeAxB(sigmaL(i+1),kcone,m2,h1(i),io)
      call system_clock(t2)
!      write(*,*)"time_mmult2_gfe=",(1D0 * (t2-t1))/1d4

      call system_clock(t1)
      call DestroyMatrix(m1,sMyName,io)
      call DestroyMatrix(m2,sMyName,io)
      call system_clock(t2)
!      write(*,*)"time_destroym_gfe=",(1D0 * (t2-t1))/1d4

    enddo
    call system_clock(ta2)
!    write(*,*)"time_gfe=",(1D0 * (ta2-ta1))/1d4
  end subroutine GaussianForwardElimination

!> \brief forward gaussian elimination sprse
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 4th of March, 2010
!> \param h0(:),hm1(:),h1(:) matrixSparseType the blocks of the tridiagonal matrix
!> \param sigmaL(:) matrixType block matrices needed for the forward elimination
!> \param io ioType, structure that has the I/O info
!> \param iBlocks integer number of diagonal blocks
  subroutine GaussianForwardEliminationSparse (h0,hm1,h1,sigmaL,iBlocks,io)
    character (len=*), parameter :: sMyName = "GaussianForwardEliminationSparse"
    type(matrixSparseType), intent(inout) :: h0(:),hm1(:),h1(:)
    type(matrixType), intent(inout) :: sigmaL(:)
    type(ioType), intent(inout) :: io
    integer, intent(in) :: iBlocks

    integer :: i,m,n
    integer :: t1,t2,ta1,ta2
    type(matrixType) :: m1,m2


!    call system_clock(t1)
!    call sleep(1)
!    call system_clock(t2)
!    write(*,*)"tsystem=",t1,t2
!    write(*,*)"tsystem2-tsystem1 (1 second)=",(1D0 * (t2-t1))/1d4

    call system_clock(ta1)
    do i=1,iBlocks-1

!!!!      if(.false..and.maxval(abs(sigmaL(i)%a))< 0.0000000000001_kdp)then
!!!!!         some action to do if one of the self-energies is zero, i.e. there
!is no coupling between layers i and i-1; this is needed in case of numerical
!problems
!one possibility is to add a small negative imaginary part to sigma in such
!cases
!!!!        cycle
!!!!      endif

      call system_clock(t1)
!      write(*,*)"gfe=",i

      n=h0(i)%iRows
      m=h0(i+1)%iRows
!      write(*,*)"mn=",i,n,m
      call AllocateMatrix(n,n,m1,sMyName,io)
      call AllocateMatrix(m,n,m2,sMyName,io)
!       m1%a=h0(i)%a-sigmaL(i)%a
     call MatrixAddSparse(m1,-kcone,sigmaL(i),kcone,h0(i),io)
      call system_clock(t2)
!      write(*,*)"time_setup_gfe=",(1D0 * (t2-t1))/1d4

      call system_clock(t1)
      call Inversematrix(m1,io)
      call system_clock(t2)
!      write(*,*)"time_inverse_gfe=",(1D0 * (t2-t1))/1d4

      call system_clock(t1)
      call ProductCeAxB(m2,kcone,hm1(i),m1,io)
      call system_clock(t2)
!      write(*,*)"time_mmult1_gfe=",(1D0 * (t2-t1))/1d4

      call system_clock(t1)
      call ProductCeAxB(sigmaL(i+1),kcone,m2,h1(i),io)
      call system_clock(t2)
!      write(*,*)"time_mmult2_gfe=",(1D0 * (t2-t1))/1d4

      call system_clock(t1)
      call DestroyMatrix(m1,sMyName,io)
      call DestroyMatrix(m2,sMyName,io)
      call system_clock(t2)
!      write(*,*)"time_destroym_gfe=",(1D0 * (t2-t1))/1d4

    enddo
    call system_clock(ta2)
!    write(*,*)"time_gfe=",(1D0 * (ta2-ta1))/1d4
  end subroutine GaussianForwardEliminationSparse


!> \brief backward gaussian elimination
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 14th of April, 2009
!> \param h0(:),hm1(:),h1(:) matrixType the blocks of the tridiagonal matrix
!> \param sigmaR(:) matrixType block matrices needed for the backward elimination
!> \param io ioType, structure that has the I/O info
!> \param iBlocks integer number of diagonal blocks
  subroutine GaussianBackwardElimination(h0,hm1,h1,sigmaR,iBlocks,io)
    character (len=*), parameter :: sMyName = "GaussianBackwardElimination"
    type(matrixType), intent(inout) :: h0(:),hm1(:),h1(:),sigmaR(:)
    type(ioType), intent(inout) :: io
    integer, intent(in) :: iBlocks

    integer :: i,m,n
    type(matrixType) :: m1,m2

    do i=iBlocks,2,-1
!      write(*,*)"gbe=",i
      n=h0(i)%iRows
      m=h0(i-1)%iRows
      call AllocateMatrix(n,n,m1,sMyName,io)
      call AllocateMatrix(m,n,m2,sMyName,io)

!       m1%a=h0(i)%a-sigmaR(i)%a
      call MatrixAdd(m1,kcone,h0(i),-kcone,sigmaR(i),io)
      call Inversematrix(m1,io)
      call ProductCeAxB(m2,kcone,h1(i-1),m1,io)
      call ProductCeAxB(sigmaR(i-1),kcone,m2,hm1(i-1),io)

      call DestroyMatrix(m1,sMyName,io)
      call DestroyMatrix(m2,sMyName,io)
    enddo

  end subroutine GaussianBackwardElimination

!> \brief backward gaussian elimination sparse
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 14th of April, 2009
!> \param h0(:),hm1(:),h1(:) matrixSparseType the blocks of the tridiagonal matrix
!> \param sigmaR(:) matrixType block matrices needed for the backward elimination
!> \param io ioType, structure that has the I/O info
!> \param iBlocks integer number of diagonal blocks
  subroutine GaussianBackwardEliminationSparseStep(h0,hm1,h1,sigmaRin,sigmaRout,io)
    character (len=*), parameter :: sMyName = "GaussianBackwardEliminationSparseStep"
    type(matrixSparseType), intent(inout) :: h0,hm1,h1
    type(matrixType), intent(inout) :: sigmaRin
    type(matrixType), intent(inout) :: sigmaRout
    type(ioType), intent(inout) :: io

    integer :: m,n
    type(matrixType) :: m1,m2

    n=h0%iRows
    m=h1%iRows
    call AllocateMatrix(n,n,m1,sMyName,io)
    call AllocateMatrix(m,n,m2,sMyName,io)

    call MatrixAddSparse(m1,-kcone,sigmaRin,kcone,h0,io)
    call Inversematrix(m1,io)
    call ProductCeAxB(m2,kcone,h1,m1,io)
    call ProductCeAxB(sigmaRout,kcone,m2,hm1,io)

    call DestroyMatrix(m1,sMyName,io)
    call DestroyMatrix(m2,sMyName,io)

  end subroutine GaussianBackwardEliminationSparseStep


!> \brief backward gaussian elimination sparse
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 14th of April, 2009
!> \param h0(:),hm1(:),h1(:) matrixSparseType the blocks of the tridiagonal matrix
!> \param sigmaR(:) matrixType block matrices needed for the backward elimination
!> \param io ioType, structure that has the I/O info
!> \param iBlocks integer number of diagonal blocks
  subroutine GaussianBackwardEliminationSparse(h0,hm1,h1,sigmaR,iBlocks,io)
    character (len=*), parameter :: sMyName = "GaussianBackwardEliminationSparse"
    type(matrixSparseType), intent(inout) :: h0(:),hm1(:),h1(:)
    type(matrixType), intent(inout) :: sigmaR(:)
    type(ioType), intent(inout) :: io
    integer, intent(in) :: iBlocks

    integer :: i,m,n
    type(matrixType) :: m1,m2

    do i=iBlocks,2,-1
!      write(*,*)"gbe=",i
      n=h0(i)%iRows
      m=h0(i-1)%iRows
      call AllocateMatrix(n,n,m1,sMyName,io)
      call AllocateMatrix(m,n,m2,sMyName,io)

!       m1%a=h0(i)%a-sigmaR(i)%a
      call MatrixAddSparse(m1,-kcone,sigmaR(i),kcone,h0(i),io)
      call Inversematrix(m1,io)
      call ProductCeAxB(m2,kcone,h1(i-1),m1,io)
      call ProductCeAxB(sigmaR(i-1),kcone,m2,hm1(i-1),io)

      call DestroyMatrix(m1,sMyName,io)
      call DestroyMatrix(m2,sMyName,io)
    enddo

  end subroutine GaussianBackwardEliminationSparse

!> \brief computes the inverse blocks for the diagonal part
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 14th of April, 2009
!> \param h0,h1,hm1 matrixType contains blocks to be inversed
!> \param g0 matrixType contains the inverse blocks
!> \param sigmaR, sigmaL matrixType block matrices generated by the backward and forward gaussian elimination
!> \param io ioType, structure that has the I/O info
!> \param iBlocks integer number of diagonal blocks
  subroutine InverseDiagonalBlocks (h0,h1,hm1,g0,sigmaL,sigmaR,iBlocks,io)
    character (len=*), parameter :: sMyName = "InverseDiagonalBlocks"
    type(matrixType), intent(inout) :: h0(:),g0(:), sigmaL(:),sigmaR(:),h1(:),hm1(:)
    integer, intent(inout) :: iBlocks
    type(ioType), intent(inout) :: io

    integer :: i,t1,t2
  ! forward gaussian
!    write(*,*)"gbf="
    call system_clock(t1)
    call GaussianForwardElimination(h0,hm1,h1,sigmaL,iBlocks,io)
    call system_clock(t2)
!    write(*,*)"time_gaussianforwardelimination=",(1D0 * (t2-t1))/1d4
  ! backward gaussian
! Note: we might remove GaussianBackwardElimination by calculating sigmar on the fly,
! Note: however we need sigmal unless we have invertable h1 and hm1
! Note: one option to save memory would be to split teh sytem into sub-systems, and in this
!       way only store the self-energies of the sub-systems; e.g we might split teh sysetm into
!       two peaces, in which way we would save about half the memory
! Note: maybe we can split the system recursively into two sub-partg, in which way we would
!       reduce the memory usage drastically for very long systems
    call system_clock(t1)
    call GaussianBackwardElimination(h0,hm1,h1,sigmaR,iBlocks,io)
    call system_clock(t2)
!    write(*,*)"time_gaussianbackwardelimination=",(1D0 * (t2-t1))/1d4
    call system_clock(t1)
    do i=1,iBlocks
!      write(*,*)"iinverse=",i
!       g0(i)%a=h0(i)%a-sigmaL(i)%a-sigmaR(i)%a
      call MatrixAdd(g0(i),kcone,h0(i),-kcone,sigmaL(i),-kcone,sigmaR(i),io)
      call Inversematrix(g0(i),io)
    enddo
    call system_clock(t2)
!    write(*,*)"time_diag_inversions=",(1D0 * (t2-t1))/1d4
  end subroutine InverseDiagonalBlocks


!> \brief computes the inverse blocks for the diagonal part sparse format
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 4th of March, 2010
!> \param h0 matrixSparseType contains diagonal blocks to be inverted
!> \param h1 matrixSparseType contains the blocks above the diagonal to be inverted
!> \param hm1 matrixSparseType contains the blocks under the diagonal to be inverted
!> \param g0 matrixType contains the inverse blocks
!> \param sigmaR, sigmaL matrixType block matrices generated by the backward and forward gaussian elimination
!> \param io ioType, structure that has the I/O info
!> \param iBlocks integer number of diagonal blocks
  subroutine InverseDiagonalBlocksSparse (h0,h1,hm1,g0,sigmaL,sigmaR,iBlocks,io)
    character (len=*), parameter :: sMyName = "InverseDiagonalBlocksSparse"
    type(matrixSparseType), intent(inout) :: h0(:),h1(:),hm1(:)
    type(matrixType), intent(inout) :: g0(:), sigmaL(:),sigmaR(:)
    integer, intent(inout) :: iBlocks
    type(ioType), intent(inout) :: io

    integer :: i,t1,t2
  ! forward gaussian
!    write(*,*)"gbf="
    call system_clock(t1)
    call GaussianForwardEliminationSparse(h0,hm1,h1,sigmaL,iBlocks,io)
    call system_clock(t2)
!    write(*,*)"time_gaussianforwardelimination=",(1D0 * (t2-t1))/1d4
  ! backward gaussian
! Note: we might remove GaussianBackwardElimination by calculating sigmar on the fly,
! Note: however we need sigmal unless we have invertable h1 and hm1
! Note: one option to save memory would be to split teh sytem into sub-systems, and in this
!       way only store the self-energies of the sub-systems; e.g we might split teh sysetm into
!       two peaces, in which way we would save about half the memory
! Note: maybe we can split the system recursively into two sub-partg, in which way we would
!       reduce the memory usage drastically for very long systems
    call system_clock(t1)
    call GaussianBackwardEliminationSparse(h0,hm1,h1,sigmaR,iBlocks,io)
    call system_clock(t2)
!    write(*,*)"time_gaussianbackwardelimination=",(1D0 * (t2-t1))/1d4
    call system_clock(t1)
    do i=1,iBlocks
!      write(*,*)"inverse=",i
!       g0(i)%a=h0(i)%a-sigmaL(i)%a-sigmaR(i)%a
      call MatrixAddSparse(g0(i),kcone,h0(i),-kcone,sigmaL(i),-kcone,sigmaR(i),io)
      call Inversematrix(g0(i),io)
    enddo
    call system_clock(t2)
!    write(*,*)"time_diag_inversions=",(1D0 * (t2-t1))/1d4
  end subroutine InverseDiagonalBlocksSparse


!> \brief computes the inverse blocks for the diagonal part sparse format
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 4th of March, 2010
!> \param h0 matrixSparseType contains diagonal blocks to be inverted
!> \param h1 matrixSparseType contains the blocks above the diagonal to be inverted
!> \param hm1 matrixSparseType contains the blocks under the diagonal to be inverted
!> \param g0 matrixType contains the inverse blocks
!> \param sigmaR, sigmaL matrixType block matrices generated by the backward and forward gaussian elimination
!> \param io ioType, structure that has the I/O info
!> \param iBlocks integer number of diagonal blocks
  subroutine InverseDiagonalOffdiagonalBlocksSparse(h0,h1,hm1,g0,sigmaL,iBlocks,gfsparse,opindex,io)
    character (len=*), parameter :: sMyName = "InverseDiagonalOffdiagonalBlocksSparse"
    type(matrixSparseType), intent(inout) :: h0(:),h1(:),hm1(:)
    type(matrixType), intent(inout) :: g0(:), sigmaL(:)
    integer, intent(inout) :: iBlocks
    type(matrixTypeGeneral), intent(inout) :: gfsparse
    integer, intent(in) :: opindex
    type(ioType), intent(inout) :: io

    type(matrixType) :: m1,m2
    integer :: m,n
    type(matrixType) :: g1,gm1
    type(matrixType), allocatable :: sigmaR(:)

    integer :: i,t1,t2
  ! forward gaussian
!    write(*,*)"gbf="
    call system_clock(t1)
    call GaussianForwardEliminationSparse(h0,hm1,h1,sigmaL,iBlocks,io)
    call system_clock(t2)
!    write(*,*)"time_gaussianforwardelimination=",(1D0 * (t2-t1))/1d4
  ! backward gaussian
! Note: we might remove GaussianBackwardElimination by calculating sigmar on the fly,
! Note: however we need sigmal unless we have invertable h1 and hm1
! Note: one option to save memory would be to split teh sytem into sub-systems, and in this
!       way only store the self-energies of the sub-systems; e.g we might split teh sysetm into
!       two peaces, in which way we would save about half the memory
! Note: maybe we can split the system recursively into two sub-partg, in which way we would
!       reduce the memory usage drastically for very long systems


    call AllocateArray(iBlocks,sigmaR,sMyName,io)

    call system_clock(t1)
    do i=iBlocks,1,-1

      call AllocateMatrix(h0(i)%iRows,h0(i)%iCols,h0(i)%iHorz,h0(i)%iVert,sigmaR(i),sMyName,io)
      sigmaR(i)%a = kczero
      if(i<iBlocks)then
        call GaussianBackwardEliminationSparseStep(h0(i+1),hm1(i),h1(i),sigmaR(i+1),sigmaR(i),io)
!!        if(maxval(abs(sigmaR(i+1)%a)) >= 0.0000000000001_kdp)then
!!          call GaussianBackwardEliminationSparseStep(h0(i+1),hm1(i),h1(i),sigmaR(i+1),sigmaR(i),io)
!!        else
!!          sigmaR(i)%a=kczero ! some action (not set to 0, since that is
!incorrect!) to deal with cases where there is no coupling betweel layers i and i+1
!!        endif
        call DestroyMatrix(sigmaR(i+1),sMyName,io)
      endif

      call AllocateMatrix(h0(i)%iRows,h0(i)%iCols,h0(i)%iHorz,h0(i)%iVert,g0(i),sMyName,io)
      g0(i)%a = kczero
      call MatrixAddSparse(g0(i),kcone,h0(i),-kcone,sigmaL(i),-kcone,sigmaR(i),io)


      call Inversematrix(g0(i),io)
      call CopySparseBlocksSingle(gfsparse,g0(i))

      if(i<iBlocks.and.(opindex==1.or.opindex==3.or.opindex==5))then
!        write(*,*)"ioi=",i,iBlocks

        n=h0(i)%iRows
        m=h0(i+1)%iRows
!        write(*,*)"mnid=",i,n,m
        call AllocateMatrix(n,n,m1,sMyName,io)
        call AllocateMatrix(n,m,m2,sMyName,io)

!         m1%a=h0(i)%a-sigmaL(i)%a
        call MatrixAddSparse(m1,-kcone,sigmaL(i),kcone,h0(i),io)
        call Inversematrix(m1,io)

        call AllocateMatrix(h1(i)%iRows,h1(i)%iCols,h1(i)%iHorz,h1(i)%iVert,g1,sMyName,io)
        call ProductCeAxB(m2,kcone,h1(i),g0(i+1),io)
        call ProductCeAxB(g1,-kcone,m1,m2,io)

        call CopySparseBlocksSingle(gfsparse,g1)
        call DestroyMatrix(g1,sMyName,io)

        call DestroyMatrix(m2,sMyName,io)
        call AllocateMatrix(m,n,m2,sMyName,io)

        call AllocateMatrix(hm1(i)%iRows,hm1(i)%iCols,hm1(i)%iHorz,hm1(i)%iVert,gm1,sMyName,io)
        call ProductCeAxB(m2,kcone,g0(i+1),hm1(i),io)
        call ProductCeAxB(gm1,-kcone,m2,m1,io)

        call CopySparseBlocksSingle(gfsparse,gm1)
        call DestroyMatrix(gm1,sMyName,io)


        call DestroyMatrix(m1,sMyName,io)
        call DestroyMatrix(m2,sMyName,io)
      endif

      if(i<iBlocks-1)call DestroyMatrix(g0(i+1),sMyName,io)

    enddo
!
    call DestroyMatrix(sigmaR(1),sMyName,io)
    call DestroyArray(sigmaR,sMyName,io)
!


    call system_clock(t2)
!    write(*,*)"time_diag_inversions=",(1D0 * (t2-t1))/1d4
  end subroutine InverseDiagonalOffdiagonalBlocksSparse



  subroutine InverseOffDiagonalBlocks (h0,h1,hm1,g0,g1,gm1,sigmaL,iBlocks,io)
    character (len=*), parameter :: sMyName = "InverseOffDiagonalBlocks"
    type(matrixType), intent(inout) :: h0(:),g0(:), sigmaL(:),h1(:),hm1(:),g1(:),gm1(:)
    integer, intent(inout) :: iBlocks
    type(ioType), intent(inout) :: io
    type(matrixType) :: m1,m2
    integer :: m,n
    integer :: i,t1,t2

    do i=1,iBlocks-1
!      write(*,*)"ioi=",i,iBlocks

      n=h0(i)%iRows
      m=h0(i+1)%iRows
!      write(*,*)"mnid=",i,n,m
      call AllocateMatrix(n,n,m1,sMyName,io)
      call AllocateMatrix(n,m,m2,sMyName,io)

!       m1%a=h0(i)%a-sigmaL(i)%a
      call MatrixAdd(m1,kcone,h0(i),-kcone,sigmaL(i),io)
      call Inversematrix(m1,io)

      call ProductCeAxB(m2,kcone,h1(i),g0(i+1),io)
      call ProductCeAxB(g1(i),-kcone,m1,m2,io)
!       g1(i)%a=-g1(i)%a

      call DestroyMatrix(m2,sMyName,io)
      call AllocateMatrix(m,n,m2,sMyName,io)

      call ProductCeAxB(m2,kcone,g0(i+1),hm1(i),io)
      call ProductCeAxB(gm1(i),-kcone,m2,m1,io)
!       gm1(i)%a=-gm1(i)%a

      call DestroyMatrix(m1,sMyName,io)
      call DestroyMatrix(m2,sMyName,io)
    enddo
  end subroutine InverseOffDiagonalBlocks

subroutine InverseOffDiagonalBlocksSparse (h0,h1,hm1,g0,g1,gm1,sigmaL,iBlocks,io)
    character (len=*), parameter :: sMyName = "InverseOffDiagonalBlocksSparse"
    type(matrixSparseType), intent(inout) :: h0(:),h1(:),hm1(:)
    type(matrixType), intent(inout) :: g0(:), sigmaL(:),g1(:),gm1(:)
    integer, intent(inout) :: iBlocks
    type(ioType), intent(inout) :: io
    type(matrixType) :: m1,m2
    integer :: m,n
    integer :: i,t1,t2

    do i=1,iBlocks-1
!      write(*,*)"ioi=",i,iBlocks

      n=h0(i)%iRows
      m=h0(i+1)%iRows
!      write(*,*)"mnid=",i,n,m
      call AllocateMatrix(n,n,m1,sMyName,io)
      call AllocateMatrix(n,m,m2,sMyName,io)

!       m1%a=h0(i)%a-sigmaL(i)%a
      call MatrixAddSparse(m1,-kcone,sigmaL(i),kcone,h0(i),io)
      call Inversematrix(m1,io)

      call ProductCeAxB(m2,kcone,h1(i),g0(i+1),io)
      call ProductCeAxB(g1(i),-kcone,m1,m2,io)
!       g1(i)%a=-g1(i)%a

      call DestroyMatrix(m2,sMyName,io)
      call AllocateMatrix(m,n,m2,sMyName,io)

      call ProductCeAxB(m2,kcone,g0(i+1),hm1(i),io)
      call ProductCeAxB(gm1(i),-kcone,m2,m1,io)
!       gm1(i)%a=-gm1(i)%a

      call DestroyMatrix(m1,sMyName,io)
      call DestroyMatrix(m2,sMyName,io)
    enddo
  end subroutine InverseOffDiagonalBlocksSparse



  subroutine InverseFirstColumnBlocks(h0,h1,hm1,gi1,a,ml,iBlocks,io)
    character(len=*), parameter :: sMyName = "InverseFirstColumnBlocks"
    type(matrixType), intent(inout) :: h0(:),h1(:),hm1(:),gi1(:),ml(:)
    type(matrixType), intent(inout) :: a
    type(ioType),intent(inout)  ::  io
    integer, intent(inout) :: iBlocks

    integer :: i
    type(matrixType) :: m1

    call AllocateMatrix(h0(iBlocks)%iRows,h0(iBlocks)%iRows,m1,sMyName,io)
!     m1%a=-h0(iBlocks)%a
    call CopyMatrix(m1,-kcone,h0(iBlocks))
    call Inversematrix(m1,io)
    call ProductCeAxB(ml(iBlocks-1),kcone,m1,hm1(iBlocks-1),io)
    call DestroyMatrix(m1,sMyName,io)
    do i=iBlocks-1,2,-1
      call AllocateMatrix(h0(i)%iRows,h0(i)%iRows,m1,sMyName,io)
      call ProductCeAxB(m1,kcone,h1(i),ml(i),io)
!       m1%a=-h0(i)%a-m1%a
      call MatrixAdd(m1,-kcone,h0(i),-kcone,m1,io)
      call Inversematrix(m1,io)
      call ProductCeAxB(ml(i-1),kcone,m1,hm1(i-1),io)
      call DestroyMatrix(m1,sMyName,io)
    enddo
!     gi1(1)%a=a%a
    call CopyMatrix(gi1(1),kcone,a)
    do i=2,iBlocks
      call ProductCeAxB(gi1(i),kcone,ml(i-1),gi1(i-1),io)
    enddo

  end subroutine InverseFirstColumnBlocks

  subroutine Inverse1NBlocksSparse(h0,h1,hm1,g1n,a,ml,iBlocks,io)
    character(len=*), parameter :: sMyName = "Inverse1NBlocksSparse"
    type(matrixSparseType), intent(inout) :: h0(:),h1(:),hm1(:)
    type(matrixType), intent(inout) :: g1n,ml(:)
    type(matrixType), intent(inout) :: a
    type(ioType),intent(inout)  ::  io
    integer, intent(inout) :: iBlocks

    integer :: i
    type(matrixType) :: m1,m2

    call AllocateMatrix(h0(iBlocks)%iRows,h0(iBlocks)%iRows,m1,sMyName,io)
!     m1%a=-h0(iBlocks)%a
    call CopySparse2Dense(m1,-kcone,h0(iBlocks))
    call Inversematrix(m1,io)
    call ProductCeAxB(ml(iBlocks-1),kcone,m1,hm1(iBlocks-1),io)
    call DestroyMatrix(m1,sMyName,io)
    do i=iBlocks-1,2,-1
      call AllocateMatrix(h0(i)%iRows,h0(i)%iRows,m1,sMyName,io)
      call ProductCeAxB(m1,kcone,h1(i),ml(i),io)
!       m1%a=-h0(i)%a-m1%a
      call MatrixAddSparse(m1,kcone,h0(i),io)
      m1%a=-m1%a
      call Inversematrix(m1,io)
      call ProductCeAxB(ml(i-1),kcone,m1,hm1(i-1),io)
      call DestroyMatrix(m1,sMyName,io)
    enddo

    call AllocateMatrix(h0(1)%iRows,h0(1)%iCols,h0(1)%iHorz,h0(1)%iVert,m1,sMyName,io)
    call CopyMatrix(m1,kcone,a)
    do i=2,iBlocks
      call AllocateMatrix(h0(i)%iRows,h0(1)%iCols,h0(1)%iHorz,h0(i)%iVert,m2,sMyName,io)
      m2%a = kczero
      call ProductCeAxB(m2,kcone,ml(i-1),m1,io)
      call DestroyMatrix(m1,sMyName,io)
      call AllocateMatrix(h0(i)%iRows,h0(1)%iCols,h0(1)%iHorz,h0(i)%iVert,m1,sMyName,io)
      call CopyMatrix(m1,kcone,m2)
      call DestroyMatrix(m2,sMyName,io)
    enddo
    call CopyMatrix(g1n,kcone,m1)
    call DestroyMatrix(m1,sMyName,io)

  end subroutine Inverse1NBlocksSparse


  subroutine InverseFirstColumnBlocksSparse(h0,h1,hm1,gi1,a,ml,iBlocks,io)
    character(len=*), parameter :: sMyName = "InverseFirstColumnBlocksSparse"
    type(matrixSparseType), intent(inout) :: h0(:),h1(:),hm1(:)
    type(matrixType), intent(inout) :: gi1(:),ml(:)
    type(matrixType), intent(inout) :: a
    type(ioType),intent(inout)  ::  io
    integer, intent(inout) :: iBlocks

    integer :: i
    type(matrixType) :: m1

    call AllocateMatrix(h0(iBlocks)%iRows,h0(iBlocks)%iRows,m1,sMyName,io)
!     m1%a=-h0(iBlocks)%a
    call CopySparse2Dense(m1,-kcone,h0(iBlocks))
    call Inversematrix(m1,io)
    call ProductCeAxB(ml(iBlocks-1),kcone,m1,hm1(iBlocks-1),io)
    call DestroyMatrix(m1,sMyName,io)
    do i=iBlocks-1,2,-1
      call AllocateMatrix(h0(i)%iRows,h0(i)%iRows,m1,sMyName,io)
      call ProductCeAxB(m1,kcone,h1(i),ml(i),io)
!       m1%a=-h0(i)%a-m1%a
      call MatrixAddSparse(m1,kcone,h0(i),io)
      m1%a=-m1%a
      call Inversematrix(m1,io)
      call ProductCeAxB(ml(i-1),kcone,m1,hm1(i-1),io)
      call DestroyMatrix(m1,sMyName,io)
    enddo
!     gi1(1)%a=a%a
    call CopyMatrix(gi1(1),kcone,a)
    do i=2,iBlocks
      call ProductCeAxB(gi1(i),kcone,ml(i-1),gi1(i-1),io)
    enddo

  end subroutine InverseFirstColumnBlocksSparse

  subroutine InverseLastColumnBlocks(h0,h1,hm1,gin,a,mm,iBlocks,io)
    character(len=*), parameter :: sMyName = "InverseLastColumnBlocks"
    type(matrixType), intent(inout) :: h0(:),h1(:),hm1(:),gin(:),mm(:)
    type(matrixType), intent(inout) :: a
    type(ioType),intent(inout)  ::  io
    integer, intent(inout) :: iBlocks

    integer :: i
    type(matrixType) :: m1

    call AllocateMatrix(h0(1)%iRows,h0(1)%iRows,m1,sMyName,io)
!     m1%a=-h0(1)%a
    call CopyMatrix(m1,-kcone,h0(1))
    call Inversematrix(m1,io)
    call ProductCeAxB(mm(1),kcone,m1,h1(1),io)
    call DestroyMatrix(m1,sMyName,io)
    do i=2,iBlocks-1
      call AllocateMatrix(h0(i)%iRows,h0(i)%iRows,m1,sMyName,io)
      call ProductCeAxB(m1,kcone,hm1(i-1),mm(i-1),io)
!       m1%a=-h0(i)%a-m1%a
      call MatrixAdd(m1,-kcone,h0(i),-kcone,m1,io)
      call Inversematrix(m1,io)
      call ProductCeAxB(mm(i),kcone,m1,h1(i),io)
      call DestroyMatrix(m1,sMyName,io)
    enddo
!     gin(iBlocks)%a=a%a
    call CopyMatrix(gin(iBlocks),kcone,a)
    do i=iBlocks-1,1,-1
      call ProductCeAxB(gin(i),kcone,mm(i),gin(i+1),io)
    enddo
  end subroutine InverseLastColumnBlocks

 subroutine InverseLastColumnBlocksSparse(h0,h1,hm1,gin,a,mm,iBlocks,io)
    character(len=*), parameter :: sMyName = "InverseLastColumnBlocksSparse"
    type(matrixSparseType), intent(inout) :: h0(:),h1(:),hm1(:)
    type(matrixType), intent(inout) :: gin(:),mm(:)
    type(matrixType), intent(inout) :: a
    type(ioType),intent(inout)  ::  io
    integer, intent(inout) :: iBlocks

    integer :: i
    type(matrixType) :: m1

    call AllocateMatrix(h0(1)%iRows,h0(1)%iRows,m1,sMyName,io)
!     m1%a=-h0(1)%a
    call CopySparse2Dense(m1,-kcone,h0(1))
    call Inversematrix(m1,io)
    call ProductCeAxB(mm(1),kcone,m1,h1(1),io)
    call DestroyMatrix(m1,sMyName,io)
    do i=2,iBlocks-1
      call AllocateMatrix(h0(i)%iRows,h0(i)%iRows,m1,sMyName,io)
      call ProductCeAxB(m1,kcone,hm1(i-1),mm(i-1),io)
!       m1%a=-h0(i)%a-m1%a
      call MatrixAddSparse(m1,kcone,h0(i),io)
      call CopyMatrix(m1,-kcone,m1)
!       m1%a=-m1%a
      call Inversematrix(m1,io)
      call ProductCeAxB(mm(i),kcone,m1,h1(i),io)
      call DestroyMatrix(m1,sMyName,io)
    enddo
!     gin(iBlocks)%a=a%a
    call CopyMatrix(gin(iBlocks),kcone,a)
    do i=iBlocks-1,1,-1
      call ProductCeAxB(gin(i),kcone,mm(i),gin(i+1),io)
    enddo
  end subroutine InverseLastColumnBlocksSparse



  subroutine CheckIFirstColumnBlocks(iBlocks,h0,h1,hm1,gi1,io)
    character(len=*), parameter :: sMyName="CheckIFirstColumnBlocks"
    integer, intent(inout) :: iBlocks
    type(ioType),intent(inout)  ::  io
    type(matrixType), intent(inout) :: gi1(:),h0(:),h1(:),hm1(:)

    type(matrixType) :: iden1,iden2,iden3
    integer :: j
    call AllocateMatrix(gi1(1)%iRows,gi1(1)%iCols,gi1(1)%iHorz,gi1(1)%iVert,iden1,sMyName,io)
    call AllocateMatrix(gi1(1)%iRows,gi1(1)%iCols,gi1(1)%iHorz,gi1(1)%iVert,iden2,sMyName,io)
    call ProductCeAxB(iden1,kcone,h0(1),gi1(1),io)
    call ProductCeAxB(iden2,kcone,h1(1),gi1(2),io)
    iden1%a=iden1%a+iden2%a
    write(io%iout,'(a,g12.6)')"Norm one for the first column blocks (1) ", NormL1(iden1,io)
    call DestroyMatrix(iden1,sMyName,io)
    call DestroyMatrix(iden2,sMyName,io)
    do j=2,iBlocks-2
      call AllocateMatrix(gi1(j)%iRows,gi1(j)%iCols,gi1(j)%iHorz,gi1(j)%iVert,iden1,sMyName,io)
      call AllocateMatrix(gi1(j)%iRows,gi1(j)%iCols,gi1(j)%iHorz,gi1(j)%iVert,iden2,sMyName,io)
      call AllocateMatrix(gi1(j)%iRows,gi1(j)%iCols,gi1(j)%iHorz,gi1(j)%iVert,iden3,sMyName,io)
      call ProductCeAxB(iden1,kcone,hm1(j-1),gi1(j-1),io)
      call ProductCeAxB(iden2,kcone,h0(j),gi1(j),io)
      call ProductCeAxB(iden3,kcone,h1(j),gi1(j+1),io)

      iden1%a=iden1%a+iden2%a+iden3%a
      write(io%iout,'(a,i0,a,g12.6)')"Norm one for the first column blocks (",j,") ",&
          NormL1(iden1,io)
      call DestroyMatrix(iden1,sMyName,io)
      call DestroyMatrix(iden2,sMyName,io)
      call DestroyMatrix(iden3,sMyName,io)
    enddo
      j=iBlocks
      call AllocateMatrix(gi1(j)%iRows,gi1(j)%iCols,gi1(j)%iHorz,gi1(j)%iVert,iden1,sMyName,io)
      call AllocateMatrix(gi1(j)%iRows,gi1(j)%iCols,gi1(j)%iHorz,gi1(j)%iVert,iden2,sMyName,io)
      call ProductCeAxB(iden1,kcone,hm1(j-1),gi1(j-1),io)
      call ProductCeAxB(iden2,kcone,h0(j),gi1(j),io)
      iden1%a=iden1%a+iden2%a
      write(io%iout,'(a,i0,a,g12.6)')"Norm one for the first column blocks (",j,") ",&
          NormL1(iden1,io)
      call DestroyMatrix(iden1,sMyName,io)
      call DestroyMatrix(iden2,sMyName,io)
  end subroutine CheckIFirstColumnBlocks

  subroutine CheckILastColumnBlocks(iBlocks,h0,h1,hm1,gin,io)
    character(len=*), parameter :: sMyName="CheckILastColumnBlocks"
    integer, intent(inout) :: iBlocks
    type(ioType),intent(inout)  ::  io
    type(matrixType), intent(inout) :: gin(:),h0(:),h1(:),hm1(:)

    type(matrixType) :: iden1,iden2,iden3
    integer :: j

    call AllocateMatrix(gin(1)%iRows,gin(1)%iCols,gin(1)%iHorz,gin(1)%iVert,iden1,sMyName,io)
    call AllocateMatrix(gin(1)%iRows,gin(1)%iCols,gin(1)%iHorz,gin(1)%iVert,iden2,sMyName,io)
    call ProductCeAxB(iden1,kcone,h0(1),gin(1),io)
    call ProductCeAxB(iden2,kcone,h1(1),gin(2),io)
    iden1%a=iden1%a+iden2%a
    write(io%iout,'(a,g12.6)')"Norm one for the last column blocks (1) ", NormL1(iden1,io)
    call DestroyMatrix(iden1,sMyName,io)
    call DestroyMatrix(iden2,sMyName,io)
    do j=2,iBlocks-1
      call AllocateMatrix(gin(j)%iRows,gin(j)%iCols,gin(j)%iHorz,gin(j)%iVert,iden1,sMyName,io)
      call AllocateMatrix(gin(j)%iRows,gin(j)%iCols,gin(j)%iHorz,gin(j)%iVert,iden2,sMyName,io)
      call AllocateMatrix(gin(j)%iRows,gin(j)%iCols,gin(j)%iHorz,gin(j)%iVert,iden3,sMyName,io)
      call ProductCeAxB(iden1,kcone,hm1(j-1),gin(j-1),io)
      call ProductCeAxB(iden2,kcone,h0(j),gin(j),io)
      call ProductCeAxB(iden3,kcone,h1(j),gin(j+1),io)
      iden1%a=iden1%a+iden2%a+iden3%a
      write(io%iout,'(a,i0,a,g12.6)')"Norm one for the last column blocks (",j,") ",&
          NormL1(iden1,io)
      call DestroyMatrix(iden1,sMyName,io)
      call DestroyMatrix(iden2,sMyName,io)
      call DestroyMatrix(iden3,sMyName,io)
    enddo
      j=iBlocks
      call AllocateMatrix(gin(j)%iRows,gin(j)%iCols,gin(j)%iHorz,gin(j)%iVert,iden1,sMyName,io)
      call AllocateMatrix(gin(j)%iRows,gin(j)%iCols,gin(j)%iHorz,gin(j)%iVert,iden2,sMyName,io)
      call ProductCeAxB(iden1,kcone,hm1(j-1),gin(j-1),io)
      call ProductCeAxB(iden2,kcone,h0(j),gin(j),io)
      iden1%a=iden1%a+iden2%a
      write(io%iout,'(a,i0,a,g12.6)')"Norm one for the last column blocks (",j,") ",&
          NormL1(iden1,io)
      call DestroyMatrix(iden1,sMyName,io)
      call DestroyMatrix(iden2,sMyName,io)
  end subroutine CheckILastColumnBlocks

end module mInverse
