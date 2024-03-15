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
!                   COPYDENSEBLOCKS,
!                   COPYDENSEBLOCKSSHIFT,
!                   COPYSPARSEBLOCKSSINGLE,
!                   COPYSPARSEBLOCKS,
!                   SETVALUESDENSEMATRIX,
!                   MEMORYACCOUNTANCYADD,
!                   MEMORYACCOUNTANCYSUB,
!                   ALLOCATEMATRIXA,
!                   ALLOCATEMATRIXB,
!                   ALLOCATEMATRIXCCS,
!                   PRINTMATRIXCRS3VECTORSDOUBLE,
!                   PRINTMATRIXCRS,
!                   PRINTMATRIXCRSP,
!                   MAXDIFFERENCEMATRIXCRS_CRS,
!                   MAXDIFFERENCEMATRIXCRS_CRSP,
!                   ALLOCATEMATRIXCRSP,
!                   ALLOCATEMATRIXCRS2,
!                   ALLOCATEMATRIXCRS,
!                   FILLBLOCKV1,
!                   FILLBLOCKV2,
!                   FILLBLOCKSPARSE,
!                   FILLBLOCKSPARSE2SPARSE,
!                   INVERSEMATRIXV1,
!                   DESTROYMATRIXA,
!                   IDENTITYA,
!                   PRODUCTCEAXBV1,
!                   ALLOCATEMATRIXARRAYV1,
!                   ALLOCATEMATRIXARRAYV2,
!                   DESTROYMATRIXARRAYV1,
!                   DESTROYMATRIXARRAYV2,
!                   DESTROYMATRIXGENERAL,
!                   ALLOCATEARRAYV1,
!                   ALLOCATEARRAYV3,
!                   ALLOCATEARRAYV2,
!                   COLLECTMATRIXGENERAL,
!                   PRINTMATRIXGENERAL,
!                   ALLOCATEMATRIXGENERALPARALLEL,
!                   ALLOCATEMATRIXGENERALSERIAL,
!                   DESTROYMATRIXCRS,
!                   DESTROYMATRIXCCS,
!                   DESTROYMATRIXSPARSEP,
!                   DESTROYMATRIXSPARSE,
!                   DESTROYARRAYV1,
!                   DESTROYARRAYV3,
!                   DESTROYARRAYV2,
!                   INCREASESPARSE,
!                   DECREASESPARSE,
!                   SCATTER,
!                   MATRIXADDSPARSEV1,
!                   MATRIXADDSPARSEV2,
!                   MATRIXADDSPARSEV3,
!                   MATRIXADDSPARSEV4,
!                   PRODUCTCEAXBV2,
!                   PRODUCTCEAXBV3,
!                   PRODUCTCEAXBV4,
!                   ROW2COLS,
!                   COL2ROWS,
!                   COPYSPARSE2DENSE,
!                   COPYSPARSECRS2DENSE,
!                   MATRIXADDV2,
!                   MATRIXADDV1,
!                   COPYMATRIX,
!                   MATDENSETOMATSPARSE,
!                   WRITEMATRIXSPARSE,
!                   DUPLICATEMATRIX  
! AND
! THE MODULES
!                   MMATRIXUTIL  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
!> \brief routines for matrix manipulation
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 1st April 2009
!> \remarks
!> \todo
!> \bug
module mMatrixUtil
  use mConstants
  use mTypes
  use mUseful
  use negfmod,only: outinfo
!$ use omp_lib
  implicit none
  private
!
  public :: MaxDifferenceMatrixCRS_CRSP
  public :: MaxDifferenceMatrixCRS_CRS
  public :: PrintMatrixGeneral
  public :: PrintMatrixCCS
  public :: PrintMatrixCRS
  public :: PrintMatrixCRS3VectorsDouble
  public :: CollectMatrixGeneral
  public :: BlockMatToFullMat
  public :: RandomSparseMatrix
  public :: AllocateMatrix
  public :: AllocateMatrixGeneral
  public :: FillBlock
  public :: InverseMatrix
  public :: DestroyMatrix
  public :: DestroyMatrixGeneral
  public :: DestroyMatrixSparse
  public :: MakeIdentity
  public :: ProductCeAxB
  public :: AllocateArray
  public :: DestroyArray
  public :: NormL1
  public :: SetValuesDenseMatrix
  public :: CopyDenseBlocks
  public :: CopyDenseBlocksShift
  public :: CopySparseBlocks
  public :: CopySparseBlocksSingle
  public :: AllocateMatrixCRS
  public :: AllocateMatrixCCS
  public :: DestroyMatrixCCS
  public :: DestroyMatrixCRS
  public :: MatrixAddSparse
  public :: MatrixAdd
  public :: CopyMatrix
  public :: CopySparse2Dense
  public :: CopySparseCRS2Dense
  public :: Row2ColS
  public :: Col2RowS
  public :: MatDenseToMatSparse
  public :: WriteMatrixSparse
  public :: DuplicateMatrix

  interface InverseMatrix
    module procedure InverseMatrixv1
  end interface

  interface FillBlock
    module procedure FillBlockv1, FillBlockv2,FillBlockSparse, FillBlockSparse2Sparse
  end interface

  interface DestroyArray
    module procedure DestroyMatrixArrayv1, DestroyMatrixArrayv2,&
      DestroyArrayv2, DestroyArrayv1, DestroyArrayv3
  end interface

  interface DestroyMatrix
    module procedure DestroyMatrixA
  end interface

  interface AllocateMatrix
    module procedure AllocateMatrixA, AllocateMatrixB
  end interface

  interface MakeIdentity
    module procedure IdentityA
  end interface

  interface AllocateMatrixGeneral
    module procedure AllocateMatrixGeneralSerial,AllocateMatrixGeneralParallel
  end interface


  interface AllocateArray
    module procedure AllocateMatrixArrayv1,AllocateMatrixArrayv2,&
     AllocateArrayv2, AllocateArrayv1, AllocateArrayv3
  end interface

  interface ProductCeAxB
    module procedure ProductCeAxBv1, ProductCeAxBv2, ProductCeAxBv3, ProductCeAxBv4
  end interface

  interface MatrixAddSparse
    module procedure MatrixAddSparsev1, MatrixAddSparsev2,MatrixAddSparsev3,MatrixAddSparsev4
  end interface

  interface MatrixAdd
    module procedure MatrixAddv1, MatrixAddv2
  end interface
  contains

!> \brief Copies the values of a matrix into a matrix type matrix
!> \details this subroutine takes a N by N matrix A as first argument,
!> the second argument is the dimension N, and the third is the output matrix B.
!> \author Ivan Rungger
!> \date 2009
!> \param A complex(kdp) , the matrix from which we copy the data
!> \param N integer , the dimension of the matrix A
!> \param h0 matrixType an array of matrice
!> \param iBlocks integer number of blocks for h0
  subroutine CopyDenseBlocks(A,N,h0,iBlocks)
    character (len=*), parameter :: sMyName = "CopyDenseBlocks"
    integer, intent(in) :: N,iBlocks
    complex(kdp) :: A(:,:)
    type(matrixType), intent(in) :: h0(:)
    integer :: ierror,i,i1,i2
    real(kdp) :: aux

    do i=1,iBlocks
      do i1=1,h0(i)%iRows
        do i2=1,h0(i)%iCols
          A(h0(i)%iVert+i1-1,h0(i)%iHorz+i2-1)=h0(i)%a(i1,i2)
        enddo
      enddo
    enddo

  end subroutine CopyDenseBlocks

!> \brief Copies the values of a matrix into a matrix type matrix shifted by i0v, i0h
!> \author Ivan Rungger
!> \date 2009
!> \param A complex(kdp) , the matrix from which we copy the data
!> \param N,M integer , the dimensions of the matrix A
!> \param h0 matrixType array, the matrices from which we shift
!> \param iBlocks integer, number of blocks for the h0
!> \param i0v, i0h integer the shifts for vertical and horizontal
  subroutine CopyDenseBlocksShift(A,N,M,h0,iBlocks,i0v,i0h)
    character (len=*), parameter :: sMyName = "CopyDenseBlocksShift"
    integer, intent(in) :: N,M,iBlocks,i0v,i0h
    complex(kdp) :: A(:,:)
    type(matrixType), intent(in) :: h0(:)
    integer :: ierror,i,i1,i2
    real(kdp) :: aux

    do i=1,iBlocks
      do i1=1,h0(i)%iRows
        do i2=1,h0(i)%iCols
!          write(12347,*)"pos=",h0(i)%iVert,h0(i)%iHorz,h0(i)%iVert-i0v+i1-1,h0(i)%iHorz-i0h+i2-1
          A(h0(i)%iVert-i0v+i1-1,h0(i)%iHorz-i0h+i2-1)=h0(i)%a(i1,i2)
        enddo
      enddo
    enddo

  end subroutine CopyDenseBlocksShift


!> \brief Copies the values of a matrix into a matrix type matrix
!> \details this subroutine takes a N by N matrix A as first argument,
!> the second argument is the dimension N, and the third is the output matrix B.
!> \author Ivan Rungger
!> \date 2009
!> \param A complex(kdp) , the matrix from which we copy the data
!> \param h0 matrixType array, the matrices from which we shift
!> \param iBlocks integer , number of blocks for the h0
  subroutine CopySparseBlocksSingle(A,h0)
    character (len=*), parameter :: sMyName = "CopySparseBlocksSingle"
    type(matrixTypeGeneral), intent(inout) :: A
    type(matrixType), intent(in) :: h0
    integer :: ierror,i,i1,i2
    real(kdp) :: aux
    integer :: dims(2),istartv,istarth,ii,jj,ind

    dims(1) = h0%iRows
    dims(2) = h0%iCols

    istartv= h0%iVert
    istarth= h0%iHorz
!!!!! an option should be added so both CCS and CRS can be copied.
    do ii=istartv,istartv+dims(1)-1
      do ind=A%matSparse%q(ii),A%matSparse%q(ii+1)-1
        if((A%matSparse%j(ind)-istarth+1 <= dims(2)).and.(A%matSparse%j(ind)-istarth+1 >= 1))then
          A%matSparse%b(ind) = h0%a(ii-istartv+1,A%matSparse%j(ind)-istarth+1)
        endif
      enddo
    enddo

  end subroutine CopySparseBlocksSingle


!> \brief Copies the values of a matrix into a matrix type matrix
!> \details this subroutine takes a N by N matrix A as first argument,
!> the second argument is the dimension N, and the third is the output matrix B.
!> \author Ivan Rungger
!> \date 2009
!> \param A complex(kdp) , the matrix from which we copy the data
!> \param h0 matrixType array, the matrices from which we shift
!> \param iBlocks integer , number of blocks for the h0
  subroutine CopySparseBlocks(A,h0,iBlocks)
    character (len=*), parameter :: sMyName = "CopySparseBlocks"
    type(matrixTypeGeneral), intent(inout) :: A
    integer, intent(in) :: iBlocks
    type(matrixType), intent(in) :: h0(:)
    integer :: ierror,i,i1,i2
    real(kdp) :: aux
    integer :: dims(2),istartv,istarth,ii,jj,ind

    do i=1,iBlocks
      dims(1) = h0(i)%iRows
      dims(2) = h0(i)%iCols

      istartv= h0(i)%iVert
      istarth= h0(i)%iHorz
!!!!!!! an option should be added so both CCS and CRS can be copied.
      do ii=istartv,istartv+dims(1)-1
        do ind=A%matSparse%q(ii),A%matSparse%q(ii+1)-1
          if((A%matSparse%j(ind)-istarth+1 <= dims(2)).and.(A%matSparse%j(ind)-istarth+1 >= 1))then
            A%matSparse%b(ind) = h0(i)%a(ii-istartv+1,A%matSparse%j(ind)-istarth+1)
          endif
        enddo
      enddo

    enddo

  end subroutine CopySparseBlocks


!> \brief Copies the values of a matrix into a matrix type matrix
!> \details this subroutine takes a N by N matrix A as first argument,
!> the second argument is the dimension N, and the third is the output matrix B.
!> \author Ivan Rungger
!> \date 2009
!> \param A(N,N) complex(kdp) , the matrix from which we copy the data
!> \param N integer , the dimension of the matrix A
!> \param mat matrixType, the matrix in which we read the data
  subroutine SetValuesDenseMatrix(A,N,mat)
    character (len=*), parameter :: sMyName = "SetValuesDenseMatrix"
    integer, intent(in) :: N
    complex(kdp), intent(in) :: A(:,:)
    type(matrixType), intent(inout) :: mat
    integer :: ierror,i,j
    real(kdp) :: aux

    do j=1,N
      do i=1,N
        mat%a(i,j)=A(i,j)
      enddo
    enddo

  end subroutine SetValuesDenseMatrix


!> \brief keeps the track of the momery allocated up to this point
!> \details It will add the lMem quantity and display a small report
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date April 2009
!> \param io type(ioType) contains the memory counter that gets updated
!> \param lMem integer(kdp2) the ammount of memory that would be added
!> \param str character contains the label for the program unit that allocated the memory
  subroutine MemoryAccountancyAdd(io,lMem,str)
    character(len=*), parameter :: sMyName="MemoryAccountancyAdd"
    type(ioType), intent(inout) :: io
    integer(kidp2), intent(in) :: lMem
    character(len=*), intent(in) ::  str

  character(len=kll) :: tmpStr, tmpStrb
    io%memCount=io%memCount + lMem
    if (lmem < 1024.0_kdp) then
      write(tmpStr,'(a,i0,a)')"Memory allocated this time: ",lmem,&
        "B by "//trim(str)//" up to this point: "
    elseif (lmem >= 1024.0_kdp) then
      write(tmpStr,'(a,f8.4,a)')"Memory allocated this time: ",lmem/1024.0_kdp,&
        "KiB by "//trim(str)//" up to this point: "
    elseif(lmem >= 1048576.0_kdp) then
      write(tmpStr,'(a,f8.4,a)')"Memory allocated this time: ",lmem/1048576.0_kdp,&
        "MiB by "//trim(str)//" up to this point: "
    elseif(lmem >= 1073741824.0_kdp) then
      write(tmpStr,'(a,f8.4,a)')"Memory allocated this time: ",lmem/1073741824.0_kdp,&
        "GiB by "//trim(str)//" up to this point: "
    elseif(lmem >= 1099511627776.0_kdp) then
      write(tmpStr,'(a,f8.4,a)')"Memory allocated this time: ",lmem/1099511627776.0_kdp,&
        "TiB by "//trim(str)//" up to this point: "
    endif
    if (io%memCount < 1024.0_kdp) then
      write(tmpStrb,'(i0,a)')int(io%memCount), "B"
    elseif (io%memCount >= 1024.0_kdp) then
      write(tmpStrb,'(f8.4,a)')io%memCount/1024.0_kdp, "KiB"
    elseif(io%memCount >= 1048576.0_kdp) then
      write(tmpStrb,'(f8.4,a)')io%memCount/1048576.0_kdp, "MiB"
    elseif(io%memCount >= 1073741824.0_kdp) then
      write(tmpStrb,'(f8.4,a)')io%memCount/1073741824.0_kdp, "GiB"
    elseif(io%memCount >= 1099511627776.0_kdp) then
      write(tmpStrb,'(f8.4,a)')io%memCount/1099511627776.0_kdp, "TiB"
    endif
    write(io%iout,'(a,a)')trim(tmpStr), trim(tmpStrb)
  end subroutine MemoryAccountancyAdd


!> \brief keeps the track of the momery allocated up to this point
!> \details It will substract the lMem quantity and display a small report
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date April 2009
!> \param io type(ioType) contains the memory counter that gets updated
!> \param lMem integer(kdp2) the ammount of memory that would be substracted
!> \param str character contains the label for the program unit that deallocated the memory
  subroutine MemoryAccountancySub(io,lMem,str)
    character(len=*), parameter :: sMyName="MemoryAccountancySub"
    type(ioType), intent(inout) :: io
    integer(kidp2), intent(in) :: lMem
    character(len=*), intent(in) ::  str

  character(len=kll) :: tmpStr, tmpStrb
    io%memCount=io%memCount - lMem
    if (lmem < 1024.0_kdp) then
      write(tmpStr,'(a,i0,a)')"Memory deallocated this time: ",lmem,&
        "B by "//trim(str)//" still allocated up to this point: "
    elseif (lmem >= 1024.0_kdp) then
      write(tmpStr,'(a,f8.4,a)')"Memory deallocated this time: ",lmem/1024.0_kdp,&
        "KiB by "//trim(str)//" still allocated up to this point: "
    elseif(lmem >= 1048576.0_kdp) then
      write(tmpStr,'(a,f8.4,a)')"Memory deallocated this time: ",lmem/1048576.0_kdp,&
        "MiB by "//trim(str)//" still allocated up to this point: "
    elseif(lmem >= 1073741824.0_kdp) then
      write(tmpStr,'(a,f8.4,a)')"Memory deallocated this time: ",lmem/1073741824.0_kdp,&
        "GiB by "//trim(str)//" still allocated up to this point: "
    elseif(lmem >= 1099511627776.0_kdp) then
      write(tmpStr,'(a,f8.4,a)')"Memory deallocated this time: ",lmem/1099511627776.0_kdp,&
        "TiB by "//trim(str)//" still allocated up to this point: "
    endif
    if (io%memCount < 1024.0_kdp) then
      write(tmpStrb,'(i0,a)')int(io%memCount), "B"
    elseif (io%memCount >= 1024.0_kdp) then
      write(tmpStrb,'(f8.4,a)')io%memCount/1024.0_kdp, "KiB"
    elseif(io%memCount >= 1048576.0_kdp) then
      write(tmpStrb,'(f8.4,a)')io%memCount/1048576.0_kdp, "MiB"
    elseif(io%memCount >= 1073741824.0_kdp) then
      write(tmpStrb,'(f8.4,a)')io%memCount/1073741824.0_kdp, "GiB"
    elseif(io%memCount >= 1099511627776.0_kdp) then
      write(tmpStrb,'(f8.4,a)')io%memCount/1099511627776.0_kdp, "TiB"
    endif
    write(io%iout,'(a,a)')trim(tmpStr), trim(tmpStrb)
  end subroutine MemoryAccountancySub



!> \brief allocates a matrix
!> \author Alin M Elena
!> \date 1st of April, 2009
!> \param rows interger no of rows
!> \param cols integer no of columns
!> \param mat matrixType the matrix
!> \param io ioType, structure that has the I/O info
!> \param substr character, calling routine
  subroutine AllocateMatrixA (rows,cols, mat,substr,io)
    character (len=*), parameter :: sMyName = "AllocateMatrixA"
    integer, intent(in) :: rows,cols
    type(matrixType), intent(inout):: mat
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr
    integer :: ierror
    integer(kidp2) :: localMem


    allocate(mat%a(rows,cols),stat=ierror)
#ifdef Memory
    call ErrorAllocate(ierror,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
#ifdef DebugNonStd
    localMem=sizeof(mat%a)
    call MemoryAccountancyAdd(io,localMem,trim(substr)//"/"//trim(sMyName))
#endif
    mat%iRows=rows
    mat%iCols=cols
  end subroutine AllocateMatrixA

!> \brief allocates a block matrix
!> \author Alin M Elena
!> \date 1st of April, 2009
!> \param rows interger no of rows
!> \param cols integer no of columns
!> \param horz, vert integer top left corner indices from the big matrix
!> \param blk blockType the matrix
!> \param io ioType, structure that has the I/O info
!> \param substr character, calling routine
  subroutine AllocateMatrixB(rows,cols,horz,vert,blk,substr,io)
    character (len=*), parameter :: sMyName = "AllocateMatrixB"
    integer, intent(in) :: rows,cols,horz,vert
    type(matrixType), intent(inout):: blk
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr
    integer :: ierror
    integer(kidp2) :: localMem

    allocate(blk%a(rows,cols),stat=ierror)
#ifdef Memory
    call ErrorAllocate(ierror,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
#ifdef DebugNonStd
    localMem=sizeof(blk%a)
    call MemoryAccountancyAdd(io,localMem,trim(substr)//"/"//trim(sMyName))
#endif
    blk%iRows=rows
    blk%iCols=cols
    blk%iHorz=horz
    blk%iVert=vert
  end subroutine AllocateMatrixB


!> \brief allocates a sparse matrix in CCS
!> \author Alin M Elena
!> \date 8th of January, 2010
!> \param rows integer no of rows
!> \param cols integer no of columns
!> \param mat matrixSparseType the matrix
!> \param nnz integer no of nonzero elements
!> \param io ioType, structure that has the I/O info
!> \param substr character, calling routine
  subroutine AllocateMatrixCCS(rows,cols, mat,nnz,substr,io)
    character (len=*), parameter :: sMyName = "AllocateMatrixCCS"
    integer, intent(in) :: rows,cols,nnz
    type(matrixSparseType), intent(inout):: mat
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr
    integer :: ierror
    integer(kidp2) :: localMem

    allocate(mat%a(nnz),stat=ierror)
#ifdef Memory
    call ErrorAllocate(ierror,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
    allocate(mat%i(nnz),stat=ierror)
#ifdef Memory
    call ErrorAllocate(ierror,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
    allocate(mat%p(cols+1),stat=ierror)
#ifdef Memory
    call ErrorAllocate(ierror,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
#ifdef DebugNonStd
    localMem=sizeof(mat%a)
    call MemoryAccountancyAdd(io,localMem,trim(substr)//"/"//trim(sMyName))
    localMem=sizeof(mat%i)
    call MemoryAccountancyAdd(io,localMem,trim(substr)//"/"//trim(sMyName))
    localMem=sizeof(mat%p)
    call MemoryAccountancyAdd(io,localMem,trim(substr)//"/"//trim(sMyName))
#endif
    mat%iRows=rows
    mat%iCols=cols
    mat%nnz=nnz
  end subroutine AllocateMatrixCCS

  subroutine PrintMatrixCRS3VectorsDouble(mat,a,n1,n2,substr,io)
  implicit none
    character (len=*), parameter :: sMyName = "PrintMatrixCRS"
    type(matrixSparseType), intent(in):: mat
    integer, intent(in) :: n1,n2
    real(kdp), intent(in):: a(n1,n2)
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr

    integer ii,jj,ind

    if (outinfo) then
    do ii=1,mat%iRows
      do ind=mat%q(ii),mat%q(ii+1)-1
        write(12347,*)substr,ii,mat%j(ind),DREAL(mat%b(ind)),DIMAG(mat%b(ind)),a(1,ind),a(2,ind),a(3,ind)
      enddo
    enddo
endif

  end subroutine PrintMatrixCRS3VectorsDouble


  subroutine PrintMatrixCRS(mat,substr,io)
  implicit none
    character (len=*), parameter :: sMyName = "PrintMatrixCRS"
    type(matrixSparseType), intent(in):: mat
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr

    integer ii,jj,ind

    if (outinfo) then
    do ii=1,mat%iRows
      do ind=mat%q(ii),mat%q(ii+1)-1
        write(12347,*)substr,ii, mat%j(ind),DREAL(mat%b(ind)),DIMAG(mat%b(ind))
      enddo
    enddo
endif

  end subroutine PrintMatrixCRS


  subroutine PrintMatrixCRSP(mat,substr,io)
  implicit none
    character (len=*), parameter :: sMyName = "PrintMatrixCRSP"
    type(matrixSparsePType), intent(in):: mat
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr

    integer ii,jj,ind

    if (outinfo) then
    do ii=1,mat%matSparse%iRows
      do ind=mat%matSparse%q(ii),mat%matSparse%q(ii+1)-1
        write(12347,*)substr,ii+mat%matSparse%iVert-1,mat%matSparse%j(ind)+mat%matSparse%iHorz-1,DREAL(mat%matSparse%b(ind)),DIMAG(mat%matSparse%b(ind))
      enddo
    enddo
endif

  end subroutine PrintMatrixCRSP

  subroutine MaxDifferenceMatrixCRS_CRS(matserial,matparallellocal,mynode_inverse,substr,io)
  implicit none
    character (len=*), parameter :: sMyName = "MaxDifferenceMatrixCRS_CRS"
    type(matrixTypeGeneral), intent(inout):: matserial,matparallellocal
    integer, intent(in) :: mynode_inverse
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr

    integer ii,jj,ind
    type(matrixSparseType) :: deltas

    if(mynode_inverse==0)then
!      call PrintMatrixGeneral(sgeneralpl,"sgeneralpl=",iout)

      call Row2Cols(matserial%matSparse,io)
      call Row2Cols(matparallellocal%matSparse,io)

      call MatrixAddSparse(deltas,(1.0D0,0.0D0),matserial%matSparse,(-1.0D0,0.0D0),matparallellocal%matSparse, io)
      if (outinfo)  write(12347,*)"deltamatrix=",substr,maxval(abs(deltas%a))
      call DestroyMatrixSparse(deltas,substr,io)

      call Col2Rows(matserial%matSparse,io)
      call Col2Rows(matparallellocal%matSparse,io)

    endif



  end subroutine MaxDifferenceMatrixCRS_CRS


  subroutine MaxDifferenceMatrixCRS_CRSP(matserial,matparallel,mynode_inverse,substr,io)
  implicit none
    character (len=*), parameter :: sMyName = "PrintMatrixCRSP"
    type(matrixTypeGeneral), intent(inout):: matserial,matparallel
    integer, intent(in) :: mynode_inverse
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr

    integer ii,jj,ind
    type(matrixTypeGeneral) :: matparallellocal
    type(matrixSparseType) :: deltas

    call CollectMatrixGeneral(matparallel,matparallellocal,substr,io)
    if(mynode_inverse==0)then
!      call PrintMatrixGeneral(sgeneralpl,"sgeneralpl=",iout)

      call Row2Cols(matserial%matSparse,io)
      call Row2Cols(matparallellocal%matSparse,io)

      call MatrixAddSparse(deltas,(1.0D0,0.0D0),matserial%matSparse,(-1.0D0,0.0D0),matparallellocal%matSparse, io)
      if (outinfo) write(12347,*)"deltamatrix=",substr,maxval(abs(deltas%a))
      call DestroyMatrixSparse(deltas,substr,io)
      call DestroyMatrixGeneral(matparallellocal,substr,io)

      call Col2Rows(matserial%matSparse,io)

    endif



  end subroutine MaxDifferenceMatrixCRS_CRSP



!> \brief allocates a sparse parallel matrix in CCS
!> \author Alin M Elena
!> \date 8th of January, 2010
!> \param rows integer no of rows
!> \param cols integer no of columns
!> \param mat matrixSparseType the matrix
!> \param nnz integer no of nonzero elements
!> \param io ioType, structure that has the I/O info
!> \param substr character, calling routine
  subroutine AllocateMatrixCRSP(mpi_group,mpi_comm,nProcs,iProc,rowsGlobal,colsGlobal,rows,cols,iHorz,iVert,nnz,mat,substr,io)
    character (len=*), parameter :: sMyName = "AllocateMatrixCRSP"
    integer, intent(in) :: mpi_group,mpi_comm
    integer, intent(in) :: nProcs,iProc
    integer, intent(in) :: rowsGlobal,colsGlobal
    integer, intent(in) :: rows,cols,nnz
    integer, intent(in) :: iHorz,iVert
    type(matrixSparsePType), intent(inout):: mat
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr
    integer :: ierror
    integer(kidp2) :: localMem


    mat%MPICommunicator=mpi_comm
    mat%MPIGroup=mpi_group
    mat%nProcs=nProcs
    mat%iProc=iProc
    mat%iRowsGlobal=rowsGlobal
    mat%iColsGlobal=colsGlobal

    call AllocateMatrixCRS2(rows,cols,iHorz,iVert,mat%matSparse,nnz,substr,io)


  end subroutine AllocateMatrixCRSP

!> \brief allocates a sparse matrix in CCS
!> \author Alin M Elena
!> \date 8th of January, 2010
!> \param rows integer no of rows
!> \param cols integer no of columns
!> \param mat matrixSparseType the matrix
!> \param nnz integer no of nonzero elements
!> \param io ioType, structure that has the I/O info
!> \param substr character, calling routine
  subroutine AllocateMatrixCRS2(rows,cols,iHorz,iVert,mat,nnz,substr,io)
    character (len=*), parameter :: sMyName = "AllocateMatrixCRS2"
    integer, intent(in) :: rows,cols,nnz,iHorz,iVert
    type(matrixSparseType), intent(inout):: mat
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr
    integer :: ierror
    integer(kidp2) :: localMem

    allocate(mat%b(nnz),stat=ierror)
#ifdef Memory
    call ErrorAllocate(ierror,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
    allocate(mat%j(nnz),stat=ierror)
#ifdef Memory
    call ErrorAllocate(ierror,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
    allocate(mat%q(rows+1),stat=ierror)
#ifdef Memory
    call ErrorAllocate(ierror,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
#ifdef DebugNonStd
    localMem=sizeof(mat%b)
    call MemoryAccountancyAdd(io,localMem,trim(substr)//"/"//trim(sMyName))
    localMem=sizeof(mat%j)
    call MemoryAccountancyAdd(io,localMem,trim(substr)//"/"//trim(sMyName))
    localMem=sizeof(mat%q)
    call MemoryAccountancyAdd(io,localMem,trim(substr)//"/"//trim(sMyName))
#endif
    mat%iRows=rows
    mat%iCols=cols
    mat%nnz=nnz
    mat%iHorz=iHorz
    mat%iVert=iVert
  end subroutine AllocateMatrixCRS2



!> \brief allocates a sparse matrix in CCS
!> \author Alin M Elena
!> \date 8th of January, 2010
!> \param rows integer no of rows
!> \param cols integer no of columns
!> \param mat matrixSparseType the matrix
!> \param nnz integer no of nonzero elements
!> \param io ioType, structure that has the I/O info
!> \param substr character, calling routine
  subroutine AllocateMatrixCRS(rows,cols, mat,nnz,substr,io)
    character (len=*), parameter :: sMyName = "AllocateMatrixCRS"
    integer, intent(in) :: rows,cols,nnz
    type(matrixSparseType), intent(inout):: mat
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr
    integer :: ierror
    integer(kidp2) :: localMem

    allocate(mat%b(nnz),stat=ierror)
#ifdef Memory
    call ErrorAllocate(ierror,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
    allocate(mat%j(nnz),stat=ierror)
#ifdef Memory
    call ErrorAllocate(ierror,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
    allocate(mat%q(rows+1),stat=ierror)
#ifdef Memory
    call ErrorAllocate(ierror,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
#ifdef DebugNonStd
    localMem=sizeof(mat%b)
    call MemoryAccountancyAdd(io,localMem,trim(substr)//"/"//trim(sMyName))
    localMem=sizeof(mat%j)
    call MemoryAccountancyAdd(io,localMem,trim(substr)//"/"//trim(sMyName))
    localMem=sizeof(mat%q)
    call MemoryAccountancyAdd(io,localMem,trim(substr)//"/"//trim(sMyName))
#endif
    mat%iRows=rows
    mat%iCols=cols
    mat%nnz=nnz
  end subroutine AllocateMatrixCRS



!> \brief copy the data into the block from the source.
!> \details The source location is the origin plus the size of the destinarion array
!> \author Michael Browne
!> \date 2009
!> \param matsource matrixType, source matrix
!> \param blkdest blockType, destination matrix
!> \remarks
  subroutine FillBlockv1(matsource,blkdest)
    character (len=*), parameter :: sMyName = "FillBlockv1"
    type(matrixType), intent(in) :: matsource
    type(matrixType), intent(inout) :: blkdest
    integer :: dims(2)

    ! get the shape of array we are copying into
!    write(*,*)"dims=",dims(1),dims(2)
    dims(1) = blkdest%iRows
!    write(*,*)"dims=",dims(1),dims(2)
    dims(2) = blkdest%iCols
!    write(*,*)"dims=",dims(1),dims(2)
    blkdest%a = matsource%a(blkdest%iVert : blkdest%iVert+dims(1)-1, blkdest%iHorz : blkdest%iHorz + dims(2)-1)

  end subroutine FillBlockv1

!> \brief copy the data into the block from the source.
!> \details The source location is the origin plus the size of the destinarion array
!> \author Alin M Elena
!> \date 2009
!> \param matsource complex(kdp), source matrix
!> \param blkdest blockType, destination matrix
!> \remarks
  subroutine FillBlockv2(matsource,blkdest)
    character (len=*), parameter :: sMyName = "FillBlockv2"
    complex(kdp), intent(in) :: matsource(:,:)
    type(matrixType), intent(inout) :: blkdest
    integer :: dims(2)

    ! get the shape of array we are copying into
    dims(1) = blkdest%iRows
    dims(2) = blkdest%iCols
    blkdest%a = matsource(blkdest%iVert : blkdest%iVert+dims(1)-1, blkdest%iHorz : blkdest%iHorz + dims(2)-1)

  end subroutine FillBlockv2


!> \brief copy the data into the block from the source.
!> \details The source location is the origin plus the size of the destinarion array
!> \author
!> \date 2009
!> \param matsource matrixType, source matrix
!> \param blkdest blockType, destination matrix
!> \remarks
  subroutine FillBlockSparse(matsource,blkdest)
    character (len=*), parameter :: sMyName = "FillBlockSparse"
    type(matrixTypeGeneral), intent(in) :: matsource
    type(matrixType), intent(inout) :: blkdest
    integer :: dims(2),istartv,istarth,ii,jj,ind,j

    dims(1) = blkdest%iRows
    dims(2) = blkdest%iCols

    istartv=blkdest%iVert
    istarth=blkdest%iHorz

! the same note as in routine CopyBlockSparse
    do ii=istartv,istartv+dims(1)-1
      do ind=matsource%matSparse%q(ii),matsource%matSparse%q(ii+1)-1
        if((matsource%matSparse%j(ind)-istarth+1 <= dims(2)).and.(matsource%matSparse%j(ind)-istarth+1 >= 1))then
          blkdest%a(ii-istartv+1,matsource%matSparse%j(ind)-istarth+1) = matsource%matSparse%b(ind)
        endif
      enddo
    enddo

    do ii=istartv,istartv+dims(1)-1
      do ind=matsource%matSparse%q(ii),matsource%matSparse%q(ii+1)-1
        if((matsource%matSparse%j(ind)-istarth+1 <= dims(2)).and.(matsource%matSparse%j(ind)-istarth+1 >= 1))then
          blkdest%a(ii-istartv+1,matsource%matSparse%j(ind)-istarth+1) = matsource%matSparse%b(ind)
        endif
      enddo
    enddo

  end subroutine FillBlockSparse

!> \brief copy the block of data into sparse CCS from a sparse CCS
!> \author Alin M Elena
!> \date 11th of February, 2010
!> \param A matrixSparseType, source matrix
!> \param B matrixSparseType, destination matrix
!> \param io ioType contains the i/o units and related
!> \remarks
  subroutine FillBlockSparse2Sparse(A,B,io)
    character (len=*), parameter :: sMyName = "FillBlockSparse"
    type(matrixSparseType), intent(in) :: A
    type(matrixSparseType), intent(inout) :: B
    type(ioType), intent(inout) :: io
    integer :: nnz,i,j
    nnz=0
    do j=B%iHorz,B%iHorz+B%iCols-1
       do i=A%p(j),A%p(j+1)-1
          if ((A%i(i)>=B%iVert).and.(A%i(i)<=B%iVert+B%iRows-1)) then
             nnz=nnz+1
          endif
       enddo
    enddo

    call AllocateMatrixCCS(B%iRows,B%iCols,B,nnz,sMyName,io)
    nnz=1
    do j=B%iHorz,B%iHorz+B%iCols-1
       B%p(j-B%iHorz+1)=nnz
       do i=A%p(j),A%p(j+1)-1
            if ((A%i(i)>=B%iVert).and.(A%i(i)<=B%iVert+B%iRows-1)) then
             B%a(nnz)=A%a(i)
             B%i(nnz)=A%i(i)-B%iVert+1
             nnz=nnz+1
          endif
       enddo
    enddo
    B%p(B%iCols+1)=nnz

  end subroutine FillBlockSparse2Sparse


!> \brief calculates the inverse \f$ A^{-1} \f$ of a matrix A
!> \details it uses lapack to compute the LU factorisation and then computes the inverse
!> \author Alin M Elena
!> \date 6th of April, 2009
!> \param mat matrixType, source matrix
!> \param io ioType, structure that has the I/O info
!> \remarks
  subroutine InverseMatrixv1 (mat,io)
  implicit none
    character (len=*), parameter :: sMyName = "InverseMatrixv1"
    type(matrixType), intent(inout) :: mat
    type(ioType), intent(inout) :: io

    integer :: lwork,ierror
    integer, allocatable :: ipiv(:)
    complex(kdp) :: dummy(1)
    complex(kdp), allocatable :: work(:)
    integer(kidp2) :: localMem
    integer :: t1,t2

!    write(12347,*)"inversematrixv1"
    call AllocateArray(mat%iRows,ipiv,sMyName,io)
    call system_clock(t1)
    call zgetrf(mat%iRows,mat%iRows,mat%a,mat%iRows,ipiv,ierror)
    call system_clock(t2)
!    write(*,*)"time_inverse_zgetrf=",(1D0 * (t2-t1))/1d4
    if (ierror /= 0 ) then
        if (outinfo) write(12347,'(a,a,i0)')trim(sMyName),": LU factorization failed with error: ", ierror
      stop
    endif
    if (io%isDebug) then
      write(io%iout,'(a,i0)')"info=",ierror
    endif
    lwork=-1
    call zgetri(mat%iRows,mat%a,mat%iRows,ipiv,dummy,lwork,ierror)
    lwork=dummy(1)
!    lwork=mat%iRows * 4
    call AllocateArray(lwork,work,sMyName,io)
    call system_clock(t1)
    call zgetri(mat%iRows,mat%a,mat%iRows,ipiv,work,lwork,ierror)
    call system_clock(t2)
!    write(*,*)"time_inverse_zgetri=",(1D0 * (t2-t1))/1d4
    if (ierror /= 0 ) then
      write(io%iout,'(a,a,i0)')trim(sMyName),": Inversion failed with error: ", ierror
      stop
    endif
    if (io%isDebug) then
      write(io%iout,'(a,i0)')"info=",ierror
    endif

    call DestroyArray(ipiv,sMyName,io)
    call DestroyArray(work,sMyName,io)
  end subroutine InverseMatrixv1

!> \brief deallocates the space used by a matrix of matrixType
!> \author Alin M Elena
!> \date 6th of April, 2009
!> \param mat matrixType, source matrix
!> \param io ioType, structure that has the I/O info
!> \param substr character, calling routine
  subroutine DestroyMatrixA(mat,substr,io)
    character (len=*), parameter :: sMyName = "DestroyMatrixA"
    type(matrixType), intent(inout) :: mat
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr

    integer :: info
    integer(kidp2) :: localMem

#ifdef DebugNonStd
    localMem=sizeof(mat%a)
    call MemoryAccountancySub(io,localMem,trim(substr)//"/"//trim(sMyName))
#endif
    deallocate(mat%a,stat=info)
#ifdef Memory
       call ErrorDeallocate(info,trim(substr)//"/"//trim(sMyName),io%iout)
#endif

  end subroutine DestroyMatrixA

!> \brief transforms a matrix to identity matrix
!> \author Alin M Elena
!> \date 6th of April 2009
!> \param mat matrixType, source matrix
!> \param io ioType, structure that has the I/O info
!> \remarks
  subroutine IdentityA(mat,io)
    character (len=*), parameter :: sMyName = "IdentityA"
    type(matrixType), intent(inout) :: mat
    type(ioType), intent(in) :: io

    integer :: i

    if (mat%iRows /= mat%iCols) then
      write(io%iout,'(a,a)')trim(sMyName), ": only square matrices can be made Identity matrix. Nothing will be done!!!"
    else
      mat%a=kczero
      do i=1,mat%iRows
        mat%a(i,i)=kcone
      enddo
    endif
  end subroutine IdentityA

!> \brief computes the product of two matrices \f$ C=\alpha AB \f$
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 7th of April, 2009
!> \param c matrixType result matrix
!> \param a,b matrixType matrices to be multiplied
!> \param alpha complex scaling factor 
!> \param io ioType, structure that has the I/O info
!> \remarks result is always of type matrixType
  subroutine ProductCeAxBv1 (c,alpha,a,b,io)
    character (len=*), parameter :: sMyName = "ProductCeAxBv1"
    type(matrixType), intent(inout) :: a,b,c
    type(ioType), intent(in) :: io
    complex(kdp), intent(in) :: alpha

    if ((a%iCols == b%iRows) .and. (c%iRows == a%iRows) .and. (c%iCols == b%iCols)) then
      call zgemm(kn, kn, a%iRows, b%iCols, a%iCols, alpha, a%a, a%iRows, b%a, b%iRows, kczero, c%a, c%iRows)
    else
      write(io%iout,'(a)')"The dimensions of matrices prevent multiplication. Please check them."
      stop
    endif
  end subroutine ProductCeAxBv1


 !> \brief allocates an array of matrices
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 7th od April 2009
!> \param array matrixType arrays of matrices to be allocated
!> \param iDim dimension of the array
!> \param io ioType, structure that has the I/O info
!> \param substr character, calling routine
  subroutine AllocateMatrixArrayv1(iDim,array,substr,io)
    character (len=*), parameter :: sMyName = "AllocateMatrixArrayv1"
    integer, intent(in) :: iDim
    type(matrixType), allocatable :: array(:)
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr

    integer :: ierror
    integer(kidp2) :: localMem
    allocate(array(iDim),stat=ierror)
#ifdef Memory
    call ErrorAllocate(ierror,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
#ifdef DebugNonStd
    localMem=sizeof(array)
    call MemoryAccountancyAdd(io,localMem,trim(substr)//"/"//trim(sMyName))
#endif
  end subroutine AllocateMatrixArrayv1

 !> \brief allocates an array of sparse matrices
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 2nd of February 2010
!> \param array matrixSparseType arrays of matrices to be allocated
!> \param iDim dimension of the array
!> \param io ioType, structure that has the I/O info
!> \param substr character, calling routine
  subroutine AllocateMatrixArrayv2(iDim,array,substr,io)
    character (len=*), parameter :: sMyName = "AllocateMatrixArrayv2"
    integer, intent(in) :: iDim
    type(matrixSparseType), allocatable :: array(:)
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr

    integer :: ierror
    integer(kidp2) :: localMem
    allocate(array(iDim),stat=ierror)
#ifdef Memory
    call ErrorAllocate(ierror,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
#ifdef DebugNonStd
    localMem=sizeof(array)
    call MemoryAccountancyAdd(io,localMem,trim(substr)//"/"//trim(sMyName))
#endif
  end subroutine AllocateMatrixArrayv2


!> \brief deallocates an array of matrices
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 7th of April, 2009
!> \param array matrixType arrays of matrices to be deallocated
!> \param io ioType, structure that has the I/O info
!> \param substr character, calling routine
  subroutine DestroyMatrixArrayv1(array,substr,io)
    character (len=*), parameter :: sMyName = "DestroyMatrixArrayv1"
    type(matrixType),allocatable,intent(inout) :: array(:)
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr
    integer :: ierror
    integer(kidp2) :: localMem

#ifdef DebugNonStd
    localMem=sizeof(array)
    call MemoryAccountancySub(io,localMem,trim(substr)//"/"//trim(sMyName))
#endif
    deallocate(array,stat=ierror)
#ifdef Memory
    call ErrorDeallocate(ierror,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
  end subroutine DestroyMatrixArrayv1

!> \brief deallocates an array of matrices
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 11th of February, 2010
!> \param array matrixType arrays of matrices to be deallocated
!> \param io ioType, structure that has the I/O info
!> \param substr character, calling routine
  subroutine DestroyMatrixArrayv2(array,substr,io)
    character (len=*), parameter :: sMyName = "DestroyMatrixArrayv2"
    type(matrixSparseType),allocatable,intent(inout) :: array(:)
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr
    integer :: ierror
    integer(kidp2) :: localMem

#ifdef DebugNonStd
    localMem=sizeof(array)
    call MemoryAccountancySub(io,localMem,trim(substr)//"/"//trim(sMyName))
#endif
    deallocate(array,stat=ierror)
#ifdef Memory
    call ErrorDeallocate(ierror,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
  end subroutine DestroyMatrixArrayv2

  subroutine DestroyMatrixGeneral(matg,substr,io)
    character (len=*), parameter :: sMyName = "DestroyMatrixGeneral"
    type(matrixTypeGeneral), intent(inout):: matg
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr

    if(matg%mattype == 0)then
      call DestroyMatrix(matg%matdense,substr,io)
    elseif(matg%mattype == 2)then
      call DestroyMatrixSparse(matg%matSparse,substr,io)
    elseif(matg%mattype == 3)then
      call DestroyMatrixSparseP(matg%matSparseP,substr,io)
    endif

    matg%mattype=-1
    matg%iRows=-1
    matg%iCols=-1
  end subroutine DestroyMatrixGeneral

!> \brief allocates an array of integers
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 14th of May 2009
!> \param array integer arrays of integers to be allocated
!> \param iDim dimension of the array
!> \param io ioType, structure that has the I/O info
!> \param substr character, calling routine
  subroutine AllocateArrayv1(iDim,array,substr,io)
    character (len=*), parameter :: sMyName = "AllocateArrayv1"
    integer, intent(in) :: iDim
    integer(kidp), allocatable :: array(:)
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr

    integer :: ierror
    integer(kidp2) :: localMem
    allocate(array(iDim),stat=ierror)
#ifdef Memory
    call ErrorAllocate(ierror,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
#ifdef DebugNonStd
    localMem=sizeof(array)
    call MemoryAccountancyAdd(io,localMem,trim(substr)//"/"//trim(sMyName))
#endif
  end subroutine AllocateArrayv1

!> \brief allocates an array of complex
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 14th of May 2009
!> \param array integer arrays of integers to be allocated
!> \param iDim dimension of the array
!> \param io ioType, structure that has the I/O info
!> \param substr character, calling routine
  subroutine AllocateArrayv3(iDim,array,substr,io)
    character (len=*), parameter :: sMyName = "AllocateArrayv3"
    integer, intent(in) :: iDim
    complex(kdp), allocatable :: array(:)
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr

    integer :: ierror
    integer(kidp2) :: localMem
    allocate(array(iDim),stat=ierror)
#ifdef Memory
    call ErrorAllocate(ierror,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
#ifdef DebugNonStd
    localMem=sizeof(array)
    call MemoryAccountancyAdd(io,localMem,trim(substr)//"/"//trim(sMyName))
#endif
  end subroutine AllocateArrayv3

!> \brief allocates an array of integers with idim rows and 2 columns
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 14th of May 2009
!> \param array integer arrays of integers to be allocated
!> \param iDim dimension of the array
!> \param io ioType, structure that has the I/O info
!> \param substr character, calling routine
  subroutine AllocateArrayv2(iDim,array,substr,io)
    character (len=*), parameter :: sMyName = "AllocateArrayv2"
    integer, intent(in) :: iDim
    integer(kidp), allocatable :: array(:,:)
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr

    integer :: ierror
    integer(kidp2) :: localMem
    allocate(array(iDim,1:2),stat=ierror)
#ifdef Memory
    call ErrorAllocate(ierror,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
#ifdef DebugNonStd
    localMem=sizeof(array)
    call MemoryAccountancyAdd(io,localMem,trim(substr)//"/"//trim(sMyName))
#endif
  end subroutine AllocateArrayv2

  subroutine CollectMatrixGeneral(matg,matg2,substr,io)

    use mMPI_NEGF

    character (len=*), parameter :: sMyName = "CollectMatrixGeneral"
    type(matrixTypeGeneral), intent(inout):: matg
    type(matrixTypeGeneral), intent(out):: matg2
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr

    integer mpi_comm,mynode,nnodes,jproc,nnz,MPIerror,i1,i2,i,ind,n1,ind2,nele
    integer, allocatable :: qGlobal(:),nRowsNode(:),iStartNode(:)
#ifdef MPI
      integer istatus(MPI_STATUS_SIZE)
#endif

#ifdef MPI

    mpi_comm=matg%matSparseP%MPICommunicator
    mynode=matg%matSparseP%iProc
    nnodes=matg%matSparseP%nProcs

    call MPI_REDUCE(matg%matSparseP%matSparse%nnz,nnz,1,MPI_integer,MPI_SUM,0,mpi_comm,MPIerror)
!    write(12347,*)"nnzglobal=",nnz,matg%matSparseP%matSparse%nnz

    allocate(nRowsNode(nnodes),iStartNode(nnodes))

    if(mynode==0)then
      nRowsNode(1)=matg%matSparseP%matSparse%iRows
      iStartNode(1)=matg%matSparseP%matSparse%iVert
      do i=1,nnodes-1
        call MPI_RECV(nRowsNode(i+1),1,MPI_integer, i, 1, mpi_comm, istatus, MPIerror)
        call MPI_RECV(iStartNode(i+1),1,MPI_integer, i, 1, mpi_comm, istatus, MPIerror)
        CALL MPI_BARRIER(mpi_comm, MPIerror)
      enddo

!      do i=1,nnodes
!        write(12347,*)"irows,istart=",i-1,iStartNode(i),nRowsNode(i)
!      enddo
    else
      do i=1,nnodes-1
        if(mynode.eq.i)then
          call MPI_SEND(matg%matSparseP%matSparse%iRows,1,MPI_integer,0,1,mpi_comm,MPIerror)
          call MPI_SEND(matg%matSparseP%matSparse%iVert,1,MPI_integer,0,1,mpi_comm,MPIerror)
        endif
        CALL MPI_BARRIER(mpi_comm, MPIerror)
      enddo
    endif

    
    n1=matg%iRows

    if(mynode==0)then

      call AllocateMatrixGeneral(0,0,0,0,n1,n1,n1,n1,1,1,nnz,2,matg2,"collectmatrix",io)

      do i1=iStartNode(1),iStartNode(1)+nRowsNode(1)
        matg2%matSparse%q(i1)=matg%matSparseP%matSparse%q(i1)
      enddo
      ind=matg2%matSparse%q(iStartNode(1))
      ind2=matg2%matSparse%q(iStartNode(1)+nRowsNode(1))-1
      matg2%matSparse%j(ind:ind2)=matg%matSparseP%matSparse%j(:)
      matg2%matSparse%b(ind:ind2)=matg%matSparseP%matSparse%b(:)
      do i=1,nnodes-1
        ind=matg2%matSparse%q(iStartNode(i+1))-1
        do i1=iStartNode(i+1),iStartNode(i+1)+nRowsNode(i+1)
          call MPI_RECV(matg2%matSparse%q(i1),1,MPI_integer, i, 1, mpi_comm, istatus, MPIerror)
          matg2%matSparse%q(i1)=matg2%matSparse%q(i1)+ind
        enddo

        ind=matg2%matSparse%q(iStartNode(i+1))
        ind2=matg2%matSparse%q(iStartNode(i+1)+nRowsNode(i+1))
        nele=ind2-ind
!        write(12347,*)"nnzrecv=",nele,i
        call MPI_RECV(matg2%matSparse%j(ind),nele,MPI_integer, i, 1, mpi_comm, istatus, MPIerror)
        call MPI_RECV(matg2%matSparse%b(ind),nele,DAT_dcomplex, i, 1, mpi_comm, istatus, MPIerror)

        CALL MPI_BARRIER(mpi_comm, MPIerror)
      enddo
!      do i=1,n1
!        write(12347,*)"qglobal=",i,matg2%matSparse%q(i)
!      enddo

    else
      do i=1,nnodes-1
        if(mynode.eq.i)then
          do i1=1,matg%matSparseP%matSparse%iRows+1
            call MPI_SEND(matg%matSparseP%matSparse%q(i1),1,MPI_integer,0,1,mpi_comm,MPIerror)
          enddo

!          write(12347,*)"nnzsend=",matg%matSparseP%matSparse%nnz
          call MPI_SEND(matg%matSparseP%matSparse%j(1),matg%matSparseP%matSparse%nnz,MPI_integer,0,1,mpi_comm,MPIerror)
          call MPI_SEND(matg%matSparseP%matSparse%b(1),matg%matSparseP%matSparse%nnz,DAT_dcomplex,0,1,mpi_comm,MPIerror)

        endif
        CALL MPI_BARRIER(mpi_comm, MPIerror)
      enddo
    endif

#endif


  end subroutine CollectMatrixGeneral



  subroutine PrintMatrixGeneral(matg,substr,io)
    character (len=*), parameter :: sMyName = "PrintMatrixGeneral"
    type(matrixTypeGeneral), intent(inout):: matg
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr

    if(matg%mattype == 2)then
      call PrintMatrixCRS(matg%matsparse,substr,io)
    elseif(matg%mattype == 3)then
      call PrintMatrixCRSP(matg%matsparsep,substr,io)
    endif

  end subroutine PrintMatrixGeneral


  subroutine AllocateMatrixGeneralParallel(mpi_group,mpi_comm,nProcs,iProc,rowsGlobal,colsGlobal,rows,cols,ihorz,ivert,nnz,mtype,matg,substr,io)
    character (len=*), parameter :: sMyName = "AllocateMatrixGeneralParallel"
    integer, intent(in) :: mpi_group,mpi_comm
    integer, intent(in) :: nProcs,iProc
    integer, intent(in) :: rowsGlobal,colsGlobal
    integer, intent(in) :: rows,cols,mtype
    integer, intent(in) :: ihorz,ivert
    integer, intent(in) :: nnz
    type(matrixTypeGeneral), intent(inout):: matg
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr
    integer :: ierror
    integer(kidp2) :: localMem

    matg%mattype=mtype
    matg%iRows=rowsGlobal
    matg%iCols=colsGlobal

    if(matg%mattype == 0)then
      call AllocateMatrix(rows,cols,matg%matdense,substr,io)
    elseif(matg%mattype == 2)then
      call AllocateMatrixCRS(rows,cols,matg%matsparse,nnz,substr,io)
    elseif(matg%mattype == 3)then
      call AllocateMatrixCRSP(mpi_group,mpi_comm,nProcs,iProc,rowsGlobal,colsGlobal,rows,cols,ihorz,ivert,nnz,matg%matsparsep,substr,io)
    endif

  end subroutine AllocateMatrixGeneralParallel


  subroutine AllocateMatrixGeneralSerial(rows,cols,nnz,mtype,matg,substr,io)
    character (len=*), parameter :: sMyName = "AllocateMatrixGeneralSerial"
    integer, intent(in) :: rows,cols,mtype
    integer, intent(in) :: nnz
    type(matrixTypeGeneral), intent(inout):: matg
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr
    integer :: ierror
    integer(kidp2) :: localMem

    matg%mattype=mtype
    matg%iRows=rows
    matg%iCols=cols

    if(matg%mattype == 0)then
      call AllocateMatrix(rows,cols,matg%matdense,substr,io)
    elseif(matg%mattype == 2)then
      call AllocateMatrixCRS(rows,cols,matg%matsparse,nnz,substr,io)
    endif

  end subroutine AllocateMatrixGeneralSerial

!> \brief deallocates the space used by a matrix of matrixSparseType CRS flavoue
!> \author Alin M Elena
!> \date 8th of January 2010
!> \param mat matrixSparseType, source matrix
!> \param io ioType, structure that has the I/O info
!> \param substr character, calling routine
  subroutine DestroyMatrixCRS(mat,substr,io)
    character (len=*), parameter :: sMyName = "DestroyMatrixCRS"
    type(matrixSparseType), intent(inout) :: mat
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr

    integer :: info
    integer(kidp2) :: localMem

#ifdef DebugNonStd
    localMem=sizeof(mat%b)
    call MemoryAccountancySub(io,localMem,trim(substr)//"/"//trim(sMyName))
    localMem=sizeof(mat%j)
    call MemoryAccountancySub(io,localMem,trim(substr)//"/"//trim(sMyName))
    localMem=sizeof(mat%q)
    call MemoryAccountancySub(io,localMem,trim(substr)//"/"//trim(sMyName))
#endif
    deallocate(mat%b, stat=info)
#ifdef Memory
     call ErrorDeallocate(info,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
     deallocate(mat%j, stat=info)
#ifdef Memory
      call ErrorDeallocate(info,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
      deallocate(mat%q, stat=info)
#ifdef Memory
      call ErrorDeallocate(info,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
  end subroutine DestroyMatrixCRS

!> \brief deallocates the space used by a matrix of matrixSparseType CCS flavoue
!> \author Alin M Elena
!> \date 8th of January 2010
!> \param mat matrixSparseType, source matrix
!> \param io ioType, structure that has the I/O info
!> \param substr character, calling routine
  subroutine DestroyMatrixCCS(mat,substr,io)
    character (len=*), parameter :: sMyName = "DestroyMatrixCRS"
    type(matrixSparseType), intent(inout) :: mat
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr

    integer :: info
    integer(kidp2) :: localMem

#ifdef DebugNonStd
    localMem=sizeof(mat%a)
    call MemoryAccountancySub(io,localMem,trim(substr)//"/"//trim(sMyName))
    localMem=sizeof(mat%i)
    call MemoryAccountancySub(io,localMem,trim(substr)//"/"//trim(sMyName))
    localMem=sizeof(mat%p)
    call MemoryAccountancySub(io,localMem,trim(substr)//"/"//trim(sMyName))
#endif
    deallocate(mat%a, stat=info)
#ifdef Memory
     call ErrorDeallocate(info,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
     deallocate(mat%i, stat=info)
#ifdef Memory
      call ErrorDeallocate(info,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
      deallocate(mat%p, stat=info)
#ifdef Memory
      call ErrorDeallocate(info,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
  end subroutine DestroyMatrixCCS

  subroutine DestroyMatrixSparseP(mat,substr,io)
    character (len=*), parameter :: sMyName = "DestroyMatrixSparseP"
    type(matrixSparsePType), intent(inout):: mat
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr

    mat%MPICommunicator = -1
    mat%MPIGroup = -1
    mat%nProcs = -1
    mat%iProc = -1
    mat%iRowsGlobal = -1
    mat%iColsGlobal = -1

    call DestroyMatrixSparse(mat%matSparse,substr,io)

  end subroutine DestroyMatrixSparseP


  subroutine DestroyMatrixSparse(mat,substr,io)
    character (len=*), parameter :: sMyName = "DestroyMatrixSparse"
    type(matrixSparseType), intent(inout):: mat
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr

    if(allocated(mat%a)) then
       call DestroyMatrixCCS(mat,substr,io)
    endif
    if(allocated(mat%b)) then
       call DestroyMatrixCRS(mat,substr,io)
    endif

  end subroutine DestroyMatrixSparse

!> \brief deallocates an array of integers
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 7th of April, 2009
!> \param array integers arrays of integers to be deallocated
!> \param io ioType, structure that has the I/O info
!> \param substr character, calling routine
  subroutine DestroyArrayv1(array,substr,io)
    character (len=*), parameter :: sMyName = "DestroyMatrixArrayv1"
    integer(kidp),intent(inout), allocatable :: array(:)
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr
    integer :: ierror
    integer(kidp2) :: localMem

#ifdef DebugNonStd
    localMem=sizeof(array)
    call MemoryAccountancySub(io,localMem,trim(substr)//"/"//trim(sMyName))
#endif
    deallocate(array,stat=ierror)
#ifdef Memory
    call ErrorDeallocate(ierror,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
  end subroutine DestroyArrayv1

!> \brief deallocates an array of complex
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 7th of April, 2009
!> \param array integers arrays of integers to be deallocated
!> \param io ioType, structure that has the I/O info
!> \param substr character, calling routine
  subroutine DestroyArrayv3(array,substr,io)
    character (len=*), parameter :: sMyName = "DestroyArrayv3"
    complex(kdp),intent(inout), allocatable :: array(:)
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr
    integer :: ierror
    integer(kidp2) :: localMem

#ifdef DebugNonStd
    localMem=sizeof(array)
    call MemoryAccountancySub(io,localMem,trim(substr)//"/"//trim(sMyName))
#endif
    deallocate(array,stat=ierror)
#ifdef Memory
    call ErrorDeallocate(ierror,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
  end subroutine DestroyArrayv3

!> \brief deallocates a 2D array of integers
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 7th of April, 2009
!> \param array integers arrays of integers to be deallocated
!> \param io ioType, structure that has the I/O info
!> \param substr character, calling routine
  subroutine DestroyArrayv2(array,substr,io)
    character (len=*), parameter :: sMyName = "DestroyArrayv2"
    integer(kidp),intent(inout),allocatable :: array(:,:)
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr
    integer :: ierror
    integer(kidp2) :: localMem

#ifdef DebugNonStd
    localMem=sizeof(array)
    call MemoryAccountancySub(io,localMem,trim(substr)//"/"//trim(sMyName))
#endif
    deallocate(array,stat=ierror)
#ifdef Memory
    call ErrorDeallocate(ierror,trim(substr)//"/"//trim(sMyName),io%iout)
#endif
  end subroutine DestroyArrayv2

  real(kdp) function NormL1(a,io)
    character(len=*), parameter :: sMyName="NormL1"
    type(matrixType), intent(in) :: a
    type(ioType),intent(in)  ::  io

    integer :: i,j
    real(kdp) :: aux,anorm

    anorm=0.0_kdp
    do i=1,a%iRows
      aux=0.0_kdp
      do j=1,a%iCols
        aux=aux+abs(a%a(i,j))
      enddo
      if (aux > anorm) then
        anorm=aux
      endif
    enddo
    NormL1=anorm
  end function NormL1

!> \brief increases the space for nonzero entries in a sparse matrix keeping the old ones
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 12th of january, 2010
!> \param A matrixSparseType the matrix to be reshaped
!> \param io ioType, structure that has the I/O info
!> \param substr character, calling routine
!> \param n integer the new number of nonzero elements
  subroutine IncreaseSparse(A,n,substr,io)
    character(len=*), parameter :: sMyName="IncreaseSparse"
    type(matrixSparseType), intent(inout) :: A
    integer, intent(in) :: n
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr
    complex(kdp),allocatable :: w(:)
    integer,allocatable :: i(:)
    integer :: m
    m=A%nnz
    call AllocateArray(m,w,trim(substr)//"/"//trim(sMyName),io)
    call AllocateArray(m,i,trim(substr)//"/"//trim(sMyName),io)
    w=A%a
    i=A%i
    call DestroyArray(A%a,trim(substr)//"/"//trim(sMyName),io)
    call DestroyArray(A%i,trim(substr)//"/"//trim(sMyName),io)
    call AllocateArray(n,A%a,trim(substr)//"/"//trim(sMyName),io)
    call AllocateArray(n,A%i,trim(substr)//"/"//trim(sMyName),io)
    A%a=w
    A%i=i
    A%nnz=n
    call DestroyArray(w,trim(substr)//"/"//trim(sMyName),io)
    call DestroyArray(i,trim(substr)//"/"//trim(sMyName),io)
  end subroutine IncreaseSparse

!> \brief sets the space for a sparse matrix to fit only the nonzero elements
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 12th of january, 2010
!> \param A matrixSparseType the matrix to be reshaped
!> \param io ioType, structure that has the I/O info
!> \param substr character, calling routine
  subroutine DecreaseSparse(A,substr,io)
    character(len=*), parameter :: sMyName="DecreaseSparse"
    type(matrixSparseType), intent(inout) :: A
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr

    integer :: n
    complex(kdp),allocatable :: w(:)
    integer,allocatable :: i(:)
    integer(kidp2) :: localMem

    n=A%p(A%iCols+1)
    call AllocateArray(n,w,trim(substr)//"/"//trim(sMyName),io)
    call AllocateArray(n,i,trim(substr)//"/"//trim(sMyName),io)
    w=A%a(1:n)
    i=A%i(1:n)
    call DestroyArray(A%a,trim(substr)//"/"//trim(sMyName),io)
    call DestroyArray(A%i,trim(substr)//"/"//trim(sMyName),io)

    call AllocateArray(n,A%a,trim(substr)//"/"//trim(sMyName),io)
    call AllocateArray(n,A%i,trim(substr)//"/"//trim(sMyName),io)
    A%a=w
    A%i=i
    A%nnz=n
    call DestroyArray(w,trim(substr)//"/"//trim(sMyName),io)
    call DestroyArray(i,trim(substr)//"/"//trim(sMyName),io)

  end subroutine DecreaseSparse

!> brief computes \f$ x = x + beta * A(:,j) \f$, where x is a dense vctor and A(:,j) is sparse
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 12th of January, 2010
!> \param A matrixSparseType the matrix which column we multiply
!> \param C matrixSparseType re result
!> \param j integer the column that we multiply
!> \param beta complex the element that multiplies column entries
!> \param x complex dense vector keeps the multiplication results
!> \param w integer vector that keeps row entries
!> \param mark integer the marker that indicates if an entry was multiplied or appears for the first time
!> \param nz current number of nonzero entries
!> \remarks result is always of type matrixSparseType
  subroutine Scatter(A,j,beta,w,x,mark,C,nz)
    character(len=*), parameter :: sMyName="Scatter"
    type(matrixSparseType), intent(in) :: A
    type(matrixSparseType), intent(inout) :: C
    complex(kdp), intent(in) :: beta
    integer, intent(inout) :: w(:),nz
    complex(kdp), intent(inout) :: x(:)
    integer, intent(in) :: mark,j

    integer :: i,p
    integer,allocatable :: myCi(:,:)
    integer :: numThreads=1, myNz
    integer :: myThreadId

! !$omp parallel do default(shared) private(i,p) schedule(dynamic)
!     do p=A%p(j),A%p(j+1)-1
!       i=A%i(p)
!       if (w(i)<mark) then
!         w(i)=mark
! !$omp critical
!         C%i(nz)=i
!         nz=nz+1
! !$omp end critical        
!         x(i)=beta*A%a(p)
!       else
!         x(i)=x(i)+beta*A%a(p)
!       endif
!     enddo
! !$omp end parallel do 

!$omp parallel   
!$ numThreads=omp_get_num_threads( )
!$omp end parallel

allocate(myCi(size(C%i),numThreads))
!$omp parallel default(shared) private(i,p,myNz,myThreadId)
   myThreadId=0
!$ myThreadId=omp_get_thread_num()
   myNz = 0   
!$omp do schedule(dynamic)
   do p=A%p(j),A%p(j+1)-1
     i=A%i(p)
     if (w(i)<mark) then
       w(i)=mark
       myNz=myNz+1
       myCi(myNz,myThreadId+1)=i       
       x(i)=beta*A%a(p)
     else
       x(i)=x(i)+beta*A%a(p)
     endif
   enddo
!$omp end do
!$omp critical
   do i=0,myNz-1
     C%i(nz+i)=myCi(i+1,myThreadId+1)
   end do
   nz=nz+myNz
!$omp end critical
!$omp end parallel 
  deallocate(myCi)
  end subroutine Scatter

!> \brief Computes \f$ C=\alpha A + \beta B \f$, where all are sprarse
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 12th of january, 2010
!> \param A,B matrixSparseType the matrices to be added
!> \param C matrixSparseType the matrix holding the answer
!> \param alpha, beta complex(kdp) the coefficients of the addition
!> \param io ioType, structure that has the I/O info
  subroutine MatrixAddSparsev1(C,alpha,A,beta,B,io)
    character(len=*), parameter :: sMyName="MatrixAddSparsev1"
    type(matrixSparseType), intent(inout) :: C
    type(matrixSparseType), intent(in) :: A,B
    complex(kdp), intent(in) :: alpha, beta
    type(ioType), intent(inout) :: io

    integer :: p,j,nz
    integer, allocatable :: w(:)
    complex(kdp), allocatable :: x(:)

    nz=1
    call AllocateArray(A%iRows,x,trim(sMyName),io)
    call AllocateArray(A%iRows,w,trim(sMyName),io)
    w=0
    if (.not.allocated(C%a)) then
!       call AllocateMatrix(A%iRows,B%iCols,C,A%p(A%iCols+1)+B%p(B%iCols+1),sMyName,io)
       call AllocateMatrixCCS(A%iRows,B%iCols,C,A%p(A%iCols+1)+B%p(B%iCols+1),sMyName,io)
    endif
    C%a=0.0_kdp
    do j=1,A%iCols
      C%p(j)=nz
      call Scatter(A,j,alpha,w,x,j+1,C,nz)
      call Scatter(B,j,beta,w,x,j+1,C,nz)
      do p=C%p(j),nz-1
        C%a(p)=x(C%i(p))
      enddo
    enddo
    C%p(B%iCols+1)=nz
    call DecreaseSparse(C,sMyName,io)
    call DestroyArray(w,trim(sMyName),io)
    call DestroyArray(x,trim(sMyName),io)
  end subroutine MatrixAddSparsev1

!> \brief Computes \f$ C=\alpha A + \beta B \f$, where C,A are dense and B sparse
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 27th of january, 2010
!> \param A matrixType the dense matrix to be added
!> \param B matrixSparseType sparse matrix to be added
!> \param C matrixType the dense matrix holding the answer
!> \param alpha, beta complex(kdp) the coefficients of the addition
!> \param io ioType, structure that has the I/O info
  subroutine MatrixAddSparsev2(C,alpha,A,beta,B,io)
    character(len=*), parameter :: sMyName="MatrixAddSparsev2"
    type(matrixType), intent(inout) :: C
    type(matrixSparseType), intent(in) ::B
    type(matrixType), intent(inout) ::A
    complex(kdp), intent(in) :: alpha, beta
    type(ioType), intent(inout) :: io

    integer :: i,j

    if (.not.allocated(C%a)) then
      call AllocateMatrix(A%iRows,A%iCols,C,sMyName,io)
    endif
!$omp parallel do default(shared) private(i,j) schedule(dynamic)
    do j=1,C%iCols
      do i=1,C%iRows
        C%a(i,j)=alpha*A%a(i,j)
      enddo
      do i=B%p(j),B%p(j+1)-1
        C%a(B%i(i),j)=C%a(B%i(i),j)+beta*B%a(i)
      enddo
    enddo
!$omp end parallel do
  end subroutine MatrixAddSparsev2

  !> \brief Computes \f$ A=A + \beta B \f$, where A is dense and B sparse
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 27th of january, 2010
!> \param A matrixType the dense matrix to be added
!> \param B matrixSparseType the sparse matrix to be added
!> \param beta complex(kdp) the coefficients of the addition
!> \param io ioType, structure that has the I/O info
  subroutine MatrixAddSparsev3(A,beta,B,io)
    character(len=*), parameter :: sMyName="MatrixAddSparsev3"
    type(matrixType), intent(inout) :: A
    type(matrixSparseType), intent(in) ::B
    complex(kdp), intent(in) ::  beta
    type(ioType), intent(inout) :: io

    integer :: i,j
!$omp parallel do default(shared) private(i,j) schedule(dynamic)
    do j=1,A%iCols
      do i=B%p(j),B%p(j+1)-1
        A%a(B%i(i),j)=A%a(B%i(i),j)+beta*B%a(i)
      enddo
    enddo
!$omp end parallel do
  end subroutine MatrixAddSparsev3


!> \brief Computes \f$ D=\alpha A + \beta B + \gamma C \f$, where D,C,B are dense and A sparse
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 27th of january, 2010
!> \param B, C matrixType the dense matrix to be added
!> \param A matrixSparseType sparse matrix to be added
!> \param D matrixType the dense matrix holding the answer
!> \param alpha, beta, gamma complex(kdp) the coefficients of the addition
!> \param io ioType, structure that has the I/O info
  subroutine MatrixAddSparsev4(D,alpha,A,beta,B,gamma,C,io)
    character(len=*), parameter :: sMyName="MatrixAddSparsev2"
    type(matrixType), intent(inout) :: D
    type(matrixSparseType), intent(in) ::A
    type(matrixType), intent(in) ::B,C
    complex(kdp), intent(in) :: alpha, beta,gamma
    type(ioType), intent(inout) :: io

    integer :: i,j

    if (.not.allocated(D%a)) then
      call AllocateMatrix(A%iRows,A%iCols,D,sMyName,io)
    endif
!$omp parallel do default(shared) private(i,j) schedule(dynamic)
    do j=1,C%iCols
      do i=1,C%iRows
          D%a(i,j)=beta*B%a(i,j)+gamma*C%a(i,j)
      enddo
      do i=A%p(j),A%p(j+1)-1
        D%a(A%i(i),j)=D%a(A%i(i),j)+alpha*A%a(i)
      enddo
    enddo
!$omp end parallel do
  end subroutine MatrixAddSparsev4

!> \brief computes the product of two sparse matrices \f$ C=\alpha AB \f$
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 12th of January, 2010
!> \param C matrixSparseType result matrix
!> \param A,B matrixSparseType matrices to be multiplied
!> \param alpha complex(kdp) scale factor for multiplication
!> \param io ioType, structure that has the I/O info
!> \remarks result is always of type matrixSparseType
  subroutine ProductCeAxBv2(C,alpha,A,B,io)
    character(len=*), parameter :: sMyName="ProductCeAxBv2"
    type(matrixSparseType),intent(inout) :: C
    type(matrixSparseType),intent(in) :: A,B
    complex(kdp), intent(in) :: alpha
    type(ioType), intent(inout) :: io
    integer :: ierr

    integer :: p,j,nz
    integer, allocatable :: w(:)
    complex(kdp), allocatable :: x(:)

    nz=1
    call AllocateArray(A%iRows,x,trim(sMyName),io)
    call AllocateArray(A%iRows,w,trim(sMyName),io)
    w=0
    if (.not.allocated(C%a)) then
      call AllocateMatrixCCS(A%iRows,B%iCols,C,A%p(A%iCols+1)+B%p(B%iCols+1),sMyName,io)
    endif
    do j=1,B%iCols
    call IncreaseSparse(C,2*C%nnz+A%iRows,sMyName,io)
      C%p(j)=nz
      do p=B%p(j),B%p(j+1)-1
        call Scatter(A,B%i(p),B%a(p),w,x,j+1,C,nz)
      enddo
      do p=C%p(j),nz-1
        C%a(p)=alpha*x(C%i(p))
      enddo
    enddo
    C%p(B%iCols+1)=nz
    call DecreaseSparse(C,sMyName,io)
    call DestroyArray(w,trim(sMyName),io)
    call DestroyArray(x,trim(sMyName),io)
  end subroutine ProductCeAxBv2



!> \brief computes the product \f$ C=\alpha AB \f$, where A, C are dense, B is sparse,
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 27th of January, 2010
!> \param C matrixSparseType result matrix
!> \param A,B matrixSparseType matrices to be multiplied
!> \param alpha complex scaling factor
!> \param io ioType, structure that has the I/O info
  subroutine ProductCeAxBv3(C,alpha,A,B,io)
    character(len=*), parameter :: sMyName="ProductCeAxBv3"
    type(matrixType),intent(inout) :: C
    type(matrixSparseType),intent(in) :: B
    type(matrixType),intent(inout) :: A
    complex(kdp), intent(in) :: alpha
    type(ioType), intent(inout) :: io

    integer :: i,j,k
    if (.not.allocated(C%a)) then
      call AllocateMatrix(A%iRows,A%iCols,C,sMyName,io)
    endif
!$omp parallel do default(shared) private(j,i,k) schedule(dynamic)
    do j=1,C%iCols
       do i=1,C%iRows
          C%a(i,j)=kczero
          do k=B%p(j),B%p(j+1)-1
             C%a(i,j)=C%a(i,j)+alpha*A%a(i,B%i(k))*B%a(k)
          enddo
       enddo
    enddo
 !$omp end parallel do
  end subroutine ProductCeAxBv3



!> \brief computes the product \f$ C=\alpha AB \f$, where B, C are dense, A is sparse,
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 27th of January, 2010
!> \param C matrixSparseType result matrix
!> \param A,B matrixSparseType matrices to be multiplied
!> \param alpha complex scaling factor
!> \param io ioType, structure that has the I/O info
  subroutine ProductCeAxBv4(C,alpha,A,B,io)
    character(len=*), parameter :: sMyName="ProductCeAxBv4"
    type(matrixType),intent(inout) :: C
    type(matrixSparseType),intent(inout) :: A
    type(matrixType),intent(in) :: B
    complex(kdp), intent(in) :: alpha
    type(ioType), intent(inout) :: io

    integer :: i,j,k
    if (.not.allocated(C%a)) then
      call AllocateMatrix(A%iRows,B%iCols,C,sMyName,io)
    endif

    call Col2RowS(A,io)
!$omp parallel do default(shared) private(j,i,k) schedule(dynamic)
    do j=1,C%iCols
       do i=1,C%iRows
          C%a(i,j)=kczero
          do k=A%q(i),A%q(i+1)-1
             C%a(i,j)=C%a(i,j)+alpha*A%b(k)*B%a(A%j(k),j)
          enddo
       enddo
    enddo
!$omp end parallel do
  end subroutine ProductCeAxBv4
 
 
!> \brief converts a sparse matrix A from CRS to CCS scheme
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 27th of January, 2010
!> \param A matrixSparseType matrix that should converted
!> \param io ioType, structure that has the I/O info
!> \remarks the routine practically just fills the space already allocated.
  subroutine Row2ColS(A,io)
    character(len=*), parameter :: sMyName="Row2ColS"
    type(ioType), intent(inout) ::  io
    type(matrixSparseType), intent(inout) :: A
    integer :: i,j,k
    integer, allocatable :: w(:)
    call AllocateArray(A%iCols,w,trim(sMyName),io)
    w=0
    if (.not.allocated(A%p)) then
      call AllocateMatrixCCS(A%iRows,A%icols,A,A%nnz,sMyName,io)
    endif

    do i=1,A%nnz
      w(A%j(i))=w(A%j(i))+1
    enddo
    a%p(1)=1
    do i=2,A%iCols+1
      a%p(i)=a%p(i-1)+w(i-1)
      w(i-1)=a%p(i-1)
    enddo
!!!$omp parallel do default(shared) private(i,j) schedule(dynamic)
    do j=1,A%iRows
      do i=A%q(j),A%q(j+1)-1
!!!$omp critical
           a%i(w(a%j(i)))=j
           a%a(w(a%j(i)))=a%b(i)
           w(a%j(i))=w(a%j(i))+1
!!!$omp end critical        
      enddo
    enddo
!!!$omp end parallel do
!at this point w has the number of nonzero elements in each row
    call DestroyArray(w,trim(sMyName),io)
  end subroutine Row2ColS



!> \brief converts a sparse matrix A from CCS to CRS scheme
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 27th of January, 2010
!> \param A matrixSparseType matrix that should converted
!> \param io ioType, structure that has the I/O info
!> \remarks the routine practically just fills the space already allocated.
  subroutine Col2RowS(A,io)
    character(len=*), parameter :: sMyName="Col2RowS"
    type(ioType), intent(inout) ::  io
    type(matrixSparseType), intent(inout) :: A
    integer :: i,j,k,ind
    integer, allocatable :: w(:)

    if (.not.allocated(A%q)) then
      call AllocateMatrixCRS(A%iRows,A%iCols,A,A%nnz,sMyName,io)
    endif

    call AllocateArray(A%iRows,w,trim(sMyName),io)
    w=0
    do ind=1,a%nnz
      w(A%i(ind))=w(A%i(ind))+1
    enddo

    a%q(1)=1
    do i=2,A%iRows+1
      a%q(i)=a%q(i-1)+w(i-1)
      w(i-1)=a%q(i-1)
    enddo

!!!$omp parallel do default(shared) private(i,j,ind) schedule(dynamic)
    do j=1,A%iCols
      do ind=A%p(j),A%p(j+1)-1
!!!$omp critical
           i=a%i(ind)
           a%j(w(i))=j
           a%b(w(i))=a%a(ind)
           w(i)=w(i)+1
!!!$omp end critical        
      enddo
    enddo
!!!$omp end parallel do
    call DestroyArray(w,trim(sMyName),io)
  end subroutine Col2RowS

!> \brief converts a sparse matrix A from CCS to dense format
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 4th of March, 2010
!> \param B matrixSparseType matrix that should converted
!> \param A matrixType dense matrix that would contain B
!> \param alpha complex scale factor
  subroutine CopySparse2Dense(A,alpha,B)
    character(len=*), parameter :: sMyName="CopySparse2Dense"
    type(matrixSparseType), intent(in) :: B
    type(matrixType), intent(inout) :: A
    complex(kdp), intent(in) :: alpha

    integer :: i,j
!$omp parallel do default(shared) private(i,j) schedule(dynamic)
    do j=1,B%iCols
      A%a(:,j)=kczero
      do i=B%p(j),B%p(j+1)-1
       A%a(B%i(i),j)=alpha*B%a(i)
      enddo
    enddo
!$omp end parallel do
  end subroutine CopySparse2Dense

!> \brief converts a sparse matrix A from CRS to dense format
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 4th of March, 2010
!> \param B matrixSparseType matrix that should converted
!> \param A matrixType dense matrix that would contain B
!> \param alpha complex scale factor
  subroutine CopySparseCRS2Dense(A,alpha,B)
    character(len=*), parameter :: sMyName="CopySparseCRS2Dense"
    type(matrixSparseType), intent(in) :: B
    type(matrixType), intent(inout) :: A
    complex(kdp), intent(in) :: alpha

    integer :: i,j,ind

!$omp parallel do default(shared) private(i,ind) schedule(dynamic)
    do i=1,B%iRows
      A%a(i,:)=kczero
      do ind=B%q(i),B%q(i+1)-1
       A%a(i,B%j(ind))=alpha*B%b(ind)
      enddo
    enddo
!$omp end parallel do
  end subroutine CopySparseCRS2Dense



!> \brief Computes \f$ D=\alpha A + \beta B + \gamma C \f$, where D,C,B,A are dense
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 10th of March, 2010
!> \param B, C, A matrixType the dense matrix to be added
!> \param D matrixType the dense matrix holding the answer
!> \param alpha, beta, gamma complex(kdp) the coefficients of the addition
!> \param io ioType, structure that has the I/O info
!> \remark the routine just adds does not do any check on dimensions or if the matrix space was already allocated
  subroutine MatrixAddv2(D,alpha,A,beta,B,gamma,C,io)
    character(len=*), parameter :: sMyName="MatrixAddv2"
    type(matrixType), intent(inout) :: D    
    type(matrixType), intent(in) ::B,C,A
    complex(kdp), intent(in) :: alpha, beta,gamma
    type(ioType), intent(inout) :: io

    integer :: i,j

!$omp parallel do default(shared) private(j) schedule(dynamic)
    do j=1,C%iCols      
          D%a(:,j)=alpha*A%a(:,j)+beta*B%a(:,j)+gamma*C%a(:,j)      
    enddo
!$omp end parallel do
  end subroutine MatrixAddv2

!> \brief Computes \f$ C=\alpha A + \beta B \f$, where C,B,A are dense
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 10th of March, 2010
!> \param B, A matrixType the dense matrix to be added
!> \param C matrixType the dense matrix holding the answer
!> \param alpha, beta complex(kdp) the scale  coefficients of the addition
!> \param io ioType, structure that has the I/O info
!> \remark the routine just adds does not do any check on dimensions or if the matrix space was already allocated
  subroutine MatrixAddv1(C,alpha,A,beta,B,io)
    character(len=*), parameter :: sMyName="MatrixAddv1"
    type(matrixType), intent(inout) :: C
    type(matrixType), intent(in) ::B,A
    complex(kdp), intent(in) :: alpha, beta
    type(ioType), intent(inout) :: io

    integer :: i,j

!$omp parallel do default(shared) private(j) schedule(dynamic)
    do j=1,C%iCols      
          C%a(:,j)=alpha*A%a(:,j)+beta*B%a(:,j)
    enddo
!$omp end parallel do
  end subroutine MatrixAddv1

!> \brief Computes \f$ B=\alpha A \f$, where B,A are dense
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 10th of March, 2010
!> \param  A matrixType the dense matrix to be scaled ot just copied
!> \param B matrixType the dense matrix holding the answer
!> \param alpha complex(kdp) the scale  coefficients of the addition
!> \remark the routine just adds does not do any check on dimensions or if the matrix space was already allocated
  subroutine CopyMatrix(B,alpha,A)
    character(len=*), parameter :: sMyName="CopyMatrix"
    type(matrixType), intent(inout) :: B,A    
    complex(kdp), intent(in) :: alpha
    
    integer :: j

!$omp parallel do default(shared) private(j) schedule(dynamic)
    do j=1,B%iCols      
          B%a(:,j)=alpha*A%a(:,j)
    enddo
!$omp end parallel do
  end subroutine CopyMatrix

!> \brief converts a dense matrix A from to sparse CRS scheme
!> \author Ivan Rungger, runggeri@tcd.ie
!> \date 18th of August, 2010
!> \param ADense matrixType matrix that should converted
!> \param ASparse matrixSparseType matrix 
!> \param io ioType, structure that has the I/O info
  subroutine MatDenseToMatSparse(ADense,ASparse,io)
    character(len=*), parameter :: sMyName="MatDenseToMatSparse"
    type(ioType), intent(inout) ::  io
    type(matrixType), intent(in) :: ADense
    type(matrixSparseType), intent(inout) :: ASparse
    integer :: i,j,k,ind,nnz

    ASparse%iRows=ADense%iRows
    ASparse%iCols=ADense%iCols
    
    ASparse%q(1)=1
    ind=0
    do i=1,ADense%iCols
      do j=1,ADense%iRows
        if(ADense%a(i,j).ne.0.0D0)then
          ind=ind+1
          ASparse%b(ind)=ADense%a(i,j)
          ASparse%j(ind)=j
        endif
      enddo
      ASparse%q(i+1)=ind+1
    enddo

    ASparse%nnz=ind


  end subroutine MatDenseToMatSparse
 
  subroutine WriteMatrixSparse(mat,label)

  implicit none
    character (len=*), intent(in) :: label
    type(matrixTypeGeneral), intent(in) :: mat

    integer i,ind

    if (outinfo) then
    do i=1,mat%iRows
      do ind=mat%matSparse%q(i),mat%matSparse%q(i+1)-1
         write(12347,*)label,i,mat%matSparse%j(ind),mat%matSparse%b(ind)
      enddo
    enddo
endif

  end subroutine WriteMatrixSparse

  subroutine DuplicateMatrix(A,B,setvalues,io)
    character(len=*), parameter :: sMyName="DuplicateMatrix"
    type(matrixSparseType), intent(in) :: A
    type(matrixSparseType), intent(out) :: B
    type(ioType), intent(inout) :: io
    logical, intent(in) :: setvalues

    logical matcrs,matccs

    matcrs=.false.
    matccs=.false.
    if (allocated(A%a)) matccs=.true.
    if (allocated(A%b)) matcrs=.true.

    if (matccs) then
       if (.not.allocated(B%a)) call AllocateMatrixCCS(A%iRows,A%iCols,B,A%nnz,sMyName,io)
       b%ihorz=a%ihorz
       b%ivert=a%ivert
       b%i=a%i
       b%p=a%p
       if(setvalues) b%a=a%a
    endif

    if (matcrs)then
       if (.not.allocated(B%b)) call AllocateMatrixCRS(A%iRows,A%iCols,B,A%nnz,sMyName,io)
       b%ihorz=a%ihorz
       b%ivert=a%ivert
       b%j=a%j
       b%q=a%q
       if(setvalues)b%b=a%b
    endif

  end subroutine DuplicateMatrix


  subroutine PrintMatrixCCS(mat,substr,io)
    type(matrixSparseType), intent(in):: mat
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: substr

    integer ind,i,mpierror

    do i=1,mat%iCols

      do ind=mat%p(i),mat%p(i+1)-1
        write(12346,*)substr,i,mat%i(ind),DREAL(mat%a(ind)),DIMAG(mat%a(ind))
      enddo

    enddo

  end subroutine PrintMatrixCCS

  subroutine BlockMatToFullMat(iBlocks,As,Bs,Cs,gmat,io)

  integer, intent(in):: iBlocks
  type(matrixSparseType),intent(inout) :: As(iBlocks),Bs(iBlocks-1),Cs(iBlocks-1)
  type(matrixSparseType),intent(out) :: gmat
  type(ioType),intent(inout) :: io

  integer i,j,j2,nnz,irows,icols,ihorz,ivert,ind

  nnz=0
  irows=0
  icols=0
!$omp parallel do default(shared) reduction(+:nnz,irows,icols)  private(i) schedule(dynamic)
  do i=1,iBlocks-1
    nnz=nnz+As(i)%nnz+Bs(i)%nnz+Cs(i)%nnz
    irows=irows+As(i)%irows
    icols=icols+As(i)%icols
  enddo
!$omp end parallel do
  nnz=nnz+As(iBlocks)%nnz
  irows=irows+As(iBlocks)%irows
  icols=icols+As(iBlocks)%icols
  write(12346,*)"nnz_total,irows,icols=",nnz,irows,icols

  call AllocateMatrixCRS(irows,icols,gmat,nnz,"BlockMatToFullMat",io)

  do i=1,iBlocks-1
    call Col2RowS(As(i),io)
    call Col2RowS(Bs(i),io)
    call Col2RowS(Cs(i),io)
  enddo
  call Col2RowS(As(iBlocks),io)

  ivert=1
  ihorz=1
  gmat%q(1)=1
  i=1
  ind=1
  do j=2,As(i)%irows+1
    ivert=ivert+1
    gmat%q(ivert)=gmat%q(ivert-1)+As(i)%q(j)-As(i)%q(j-1)+Bs(i)%q(j)-Bs(i)%q(j-1)
    do j2=As(i)%q(j-1),As(i)%q(j)-1
      gmat%b(ind)=As(i)%b(j2)
      gmat%j(ind)=As(i)%j(j2)+ihorz-1
      ind=ind+1
    enddo
    do j2=Bs(i)%q(j-1),Bs(i)%q(j)-1
      gmat%b(ind)=Bs(i)%b(j2)
      gmat%j(ind)=Bs(i)%j(j2)+ihorz-1+As(i)%iRows
      ind=ind+1
    enddo
  enddo

  do i=2,iBlocks-1
    ihorz=ihorz+As(i-1)%iRows
    do j=2,As(i)%irows+1
      ivert=ivert+1
      gmat%q(ivert)=gmat%q(ivert-1)+Cs(i-1)%q(j)-Cs(i-1)%q(j-1)+As(i)%q(j)-As(i)%q(j-1)+Bs(i)%q(j)-Bs(i)%q(j-1)
      do j2=Cs(i-1)%q(j-1),Cs(i-1)%q(j)-1
        gmat%b(ind)=Cs(i-1)%b(j2)
        gmat%j(ind)=Cs(i-1)%j(j2)+ihorz-1-As(i-1)%irows
        ind=ind+1
      enddo
      do j2=As(i)%q(j-1),As(i)%q(j)-1
        gmat%b(ind)=As(i)%b(j2)
        gmat%j(ind)=As(i)%j(j2)+ihorz-1
        ind=ind+1
      enddo
      do j2=Bs(i)%q(j-1),Bs(i)%q(j)-1
        gmat%b(ind)=Bs(i)%b(j2)
        gmat%j(ind)=Bs(i)%j(j2)+ihorz-1+As(i)%iRows
        ind=ind+1
      enddo
    enddo
  enddo

  i=iBlocks
  ihorz=ihorz+As(i-1)%iRows
  do j=2,As(i)%irows+1
    ivert=ivert+1
    gmat%q(ivert)=gmat%q(ivert-1)+Cs(i-1)%q(j)-Cs(i-1)%q(j-1)+As(i)%q(j)-As(i)%q(j-1)
    do j2=Cs(i-1)%q(j-1),Cs(i-1)%q(j)-1
      gmat%b(ind)=Cs(i-1)%b(j2)
      gmat%j(ind)=Cs(i-1)%j(j2)+ihorz-1-As(i-1)%irows
      ind=ind+1
    enddo
    do j2=As(i)%q(j-1),As(i)%q(j)-1
      gmat%b(ind)=As(i)%b(j2)
      gmat%j(ind)=As(i)%j(j2)+ihorz-1
      ind=ind+1
    enddo
  enddo

!  write(12346,*)"ivert,ihorz,irows,icols,ind=",ivert-1,ihorz-1+As(i)%irows,irows,icols,ind


  end subroutine BlockMatToFullMat


  subroutine RandomSparseMatrix(nnz,n,m,mat,DiagonalElementsFilled)
    character(len=*), parameter :: sMyName="RandomSparseMatrix"
    integer, intent(in) :: nnz,n,m
    type(matrixSparseType), intent(inout) :: mat
    logical, intent(in) :: DiagonalElementsFilled

    integer :: i,j,k,ind
    real(kdp) :: kk,cm(2)
    integer ncol(m),w(n),indlistw(n)

    if(DiagonalElementsFilled)then
      if(m>nnz)then
        write(*,*)"not enough non-zero elements for a matrix with all diagonal elements non-zero."
        stop
      endif
      if(m .ne. n)then
        write(*,*)"n needs to be equal to m for a matrix with all diagonal elements non-zero."
        stop
      endif

      ncol=1
      k=m
    else
      ncol=0
      k=0
    endif

    do while (k<nnz)
      call random_number(kk)
      i=int(kk*m)+1

      if(ncol(i)<n)then
        ncol(i)=ncol(i)+1
        k=k+1
      endif

    enddo

!    write(*,*)"nnz",nnz,n,m
!    mat%p(1)=1
!    do i=1,m
!      write(*,*)"ncol_i",i,ncol(i)
!      mat%p(i+1)
!    enddo
    mat%p(1)=1
    do i=1,m
      mat%p(i+1)=mat%p(i)+ncol(i)
    enddo


    w=0
    do j=1,m


      if(DiagonalElementsFilled)then
        mat%i(mat%p(j))=j
        w(j)=1
        k=1
        indlistw(k)=j
      else
        k=0
      endif

      do while (k<ncol(j))
        call random_number(kk)
        i=int(kk*n)+1
  
        if(w(i)==0)then

          mat%i(mat%p(j)+k)=i

          w(i)=1
          k=k+1
          indlistw(k)=i
        endif
      enddo
      do i=1,k
        w(indlistw(i))=0
      enddo

    enddo

    do k=1,nnz
      call random_number(cm)
      cm=cm-0.5_kdp
      mat%a(k)=cmplx(cm(1),cm(2),kdp)
    enddo

!!    write(*,*)"matnnz=",mat%p(m+1)-1
!!    do j=1,m
!!      write(*,*)"matnnzcol=",j,mat%p(j+1)-mat%p(j)
!!      do k=mat%p(j),mat%p(j+1)-1
!!        write(*,*)"matij=",j,mat%i(k),mat%a(k)
!!      enddo
!!    enddo

  end subroutine RandomSparseMatrix



end module mMatrixUtil
