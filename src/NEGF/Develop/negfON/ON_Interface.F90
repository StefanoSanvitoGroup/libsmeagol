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
!                   INVERTONGENERAL2,
!                   INVERTONGENERAL,
!                   INVERTSPARSEON,
!                   INVERTSPARSEONV3,
!                   INVERTSPARSEONV2,
!                   INVERTDENSEON  
! AND
! THE MODULE
!                   MONINTERFACE  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
module mONInterface
 use mConstants
 use mTypes
 use mMatrixUtil
 use mPartition
 use mGutenberg
 use mInverse

 implicit none
 private

  public :: InvertONGeneral
  public :: InvertONGeneral2

contains

  subroutine InvertONGeneral2(N1,gfmat,gfserial,nl,nr,gfout,opindex,solver)


  use negfmod, only: outinfo
    use mMPI_NEGF

    character(len=*), parameter :: sMyName="InvertONGeneral"
    type(matrixTypeGeneral), intent(inout) :: gfmat
    type(matrixTypeGeneral), intent(inout) :: gfout
    integer, intent(in) :: N1,nl,nr,opindex,solver
    type(matrixTypeGeneral), intent(inout) :: gfserial

    integer opindexInternal
    type(ioType) :: io

    io%isDebug=.false.

    opindexInternal=opindex
    if(gfmat%mattype == 0)then

      if(opindex==4)opindexInternal=2
      if(opindex==5)opindexInternal=3

      call InvertDenseON(n1,gfmat%matdense,nl,nr,opindexInternal)
      if(opindexInternal==2.or.opindexInternal==3)then
        gfout%matdense%a(:,1:nl)=gfmat%matdense%a(:,1:nl)
        gfout%matdense%a(:,nl+1:nl+nr)=gfmat%matdense%a(:,n1-nr+1:n1)
      endif
    elseif(gfmat%mattype == 2)then
    if (outinfo) write(12347,*)"serial inversion"
      if(solver == 1)then
        call InvertSparseONv3(N1,gfmat,nl,nr,gfout,opindexInternal)
      else

        if(opindex==4)opindexInternal=2
        if(opindex==5)opindexInternal=3
        call InvertSparseON(N1,gfmat,nl,nr,gfout,opindexInternal)
      endif
    elseif(gfmat%mattype == 3)then
    if (outinfo) write(12347,*)"parallel inversion"
      call CollectMatrixGeneral(gfmat,gfserial,"gfconvert",io)

      if(mynode_inverse==0)then
        if(solver == 1)then
          call InvertSparseONv3(N1,gfserial,nl,nr,gfout,opindexInternal)
        else

          if(opindex==4)opindexInternal=2
          if(opindex==5)opindexInternal=3
          call InvertSparseON(N1,gfserial,nl,nr,gfout,opindexInternal)
        endif
!        call DestroyMatrixGeneral(gfserial,"gfconvert",io)
      endif
    endif
  end subroutine InvertONGeneral2


  subroutine InvertONGeneral(N1,gfmat,nl,nr,gfout,opindex,solver)
    character(len=*), parameter :: sMyName="InvertONGeneral"
    type(matrixTypeGeneral), intent(inout) :: gfmat
    type(matrixTypeGeneral), intent(inout) :: gfout
    integer, intent(in) :: N1,nl,nr,opindex,solver

    integer opindexInternal

    opindexInternal=opindex
    if(gfmat%mattype == 0)then

      if(opindex==4)opindexInternal=2
      if(opindex==5)opindexInternal=3

! this will be called by main.F90
      call InvertDenseON(n1,gfmat%matdense,nl,nr,opindexInternal)
      if(opindexInternal==2.or.opindexInternal==3)then
        gfout%matdense%a(:,1:nl)=gfmat%matdense%a(:,1:nl)
        gfout%matdense%a(:,nl+1:nl+nr)=gfmat%matdense%a(:,n1-nr+1:n1)
      endif
    elseif(gfmat%mattype == 2)then
      if(solver == 1)then
! this will be called by main_sparse.F90
        call InvertSparseONv3(N1,gfmat,nl,nr,gfout,opindexInternal)
      else

        if(opindex==4)opindexInternal=2
        if(opindex==5)opindexInternal=3
        call InvertSparseON(N1,gfmat,nl,nr,gfout,opindexInternal)
      endif
    endif
  end subroutine InvertONGeneral

  subroutine InvertSparseON(N1,gfsparse,nl,nr,gfout,opindex)
    type(matrixTypeGeneral) :: gfsparse
    type(matrixTypeGeneral) :: gfout
    integer, intent(in) :: N1,nl,nr,opindex
    character(len=*), parameter :: sMyName="keldyshimag"
    type(ioType) :: io
    type(matrixType), allocatable :: h0(:), h1(:), hm1(:)
    integer :: iBlocks,i1,i
    type(matrixType), allocatable :: sigmaL(:),sigmaR(:),g0(:),gi1(:),gin(:),ml(:),mm(:),gm1(:),g1(:)

    io%isDebug=.false.

    call PartitionMatrix(h0,h1,hm1,iBlocks,gfsparse,nl,nr,io)
    call FillBlocksFromMatrixsparse(h0,h1,hm1,iBlocks,gfsparse)

    call AllocateArray(iBlocks,sigmaL,sMyName,io)
    call AllocateArray(iBlocks,sigmaR,sMyName,io)
    call AllocateArray(iBlocks,g0,sMyName,io)
    do i=1,iBlocks
      call AllocateMatrix(h0(i)%iRows,h0(i)%iCols,h0(i)%iHorz,h0(i)%iVert,sigmaL(i),sMyName,io)
      sigmaL(i)%a = kczero
      call AllocateMatrix(h0(i)%iRows,h0(i)%iCols,h0(i)%iHorz,h0(i)%iVert,sigmaR(i),sMyName,io)
      sigmaR(i)%a = kczero
      call AllocateMatrix(h0(i)%iRows,h0(i)%iCols,h0(i)%iHorz,h0(i)%iVert,g0(i),sMyName,io)
      g0(i)%a = kczero
    enddo

    call InverseDiagonalBlocks(h0,h1,hm1,g0,sigmaL,sigmaR,iBlocks,io)

    do i=1,iBlocks
      call DestroyMatrix(sigmaR(i),sMyName,io)
    enddo
    call DestroyArray(sigmaR,sMyName,io)


    if(opindex==1.or.opindex==3)then

      call AllocateArray(iBlocks,g1,sMyName,io)
      call AllocateArray(iBlocks,gm1,sMyName,io)

      do i=1,iBlocks-1
        call AllocateMatrix(h1(i)%iRows,h1(i)%iCols,h1(i)%iHorz,h1(i)%iVert,g1(i),sMyName,io)
        call AllocateMatrix(hm1(i)%iRows,hm1(i)%iCols,hm1(i)%iHorz,hm1(i)%iVert,gm1(i),sMyName,io)
      enddo

      call InverseOffDiagonalBlocks(h0,h1,hm1,g0,g1,gm1,sigmaL,iBlocks,io)

      call CopySparseBlocks(gfsparse,g0,iBlocks)
      call CopySparseBlocks(gfsparse,g1,iBlocks)
      call CopySparseBlocks(gfsparse,gm1,iBlocks)

      do i=1,iBlocks-1
        call DestroyMatrix(g1(i),sMyName,io)
        call DestroyMatrix(gm1(i),sMyName,io)
      enddo
      call DestroyArray(g1,sMyName,io)
      call DestroyArray(gm1,sMyName,io)

    endif

    do i=2,iBlocks-1
      call DestroyMatrix(g0(i),sMyName,io)
    enddo

    do i=1,iBlocks
      call DestroyMatrix(sigmaL(i),sMyName,io)
    enddo
    call DestroyArray(sigmaL,sMyName,io)

    if(opindex==2.or.opindex==3)then


      gfout%matdense%a=0.0_kdp
      call AllocateArray(iBlocks,gin,sMyName,io)
      call AllocateArray(iBlocks-1,mm,sMyName,io)

      do i=1,iBlocks
          call AllocateMatrix(h0(i)%iRows,h0(iBlocks)%iCols,h0(iBlocks)%iHorz,h0(i)%iVert,gin(i),sMyName,io)
          gin(i)%a = kczero
          if (i<=iBlocks-1) then
            call AllocateMatrix(h0(i)%iRows,h0(i+1)%iCols,h0(i)%iHorz,h0(1)%iVert,mm(i),sMyName,io)
            mm(i)%a = kczero
          endif
      enddo

      call InverseLastColumnBlocks(h0,h1,hm1,gin,g0(iBlocks),mm,iBlocks,io)
      call CopyDenseBlocksShift(gfout%matdense%a,n1,nl+nr,gin,iBlocks,0,n1-nr-nl)

      do i=1,iBlocks-1
        call DestroyMatrix(gin(i),sMyName,io)
        call DestroyMatrix(mm(i),sMyName,io)
      enddo
      call DestroyMatrix(gin(iBlocks),sMyName,io)
      call DestroyArray(mm,sMyName,io)
      call DestroyArray(gin,sMyName,io)



      call AllocateArray(iBlocks,gi1,sMyName,io)
      call AllocateArray(iBlocks-1,ml,sMyName,io)
      do i=1,iBlocks
          call AllocateMatrix(h0(i)%iRows,h0(1)%iCols,h0(1)%iHorz,h0(i)%iVert,gi1(i),sMyName,io)
          gi1(i)%a = kczero
          if (i>=2) then
            call AllocateMatrix(h0(i)%iRows,h0(i-1)%iCols,h0(i)%iHorz,h0(1)%iVert,ml(i-1),sMyName,io)
            ml(i-1)%a = kczero
          endif
      enddo

      call InverseFirstColumnBlocks(h0,h1,hm1,gi1,g0(1),ml,iBlocks,io)
      call CopyDenseBlocksShift(gfout%matdense%a,n1,nl+nr,gi1,iBlocks,0,0)

      do i=1,iBlocks-1
        call DestroyMatrix(gi1(i),sMyName,io)
        call DestroyMatrix(ml(i),sMyName,io)
      enddo
      call DestroyMatrix(gi1(iBlocks),sMyName,io)
      call DestroyArray(gi1,sMyName,io)
      call DestroyArray(ml,sMyName,io)

    endif



    call DestroyMatrix(g0(1),sMyName,io)
    call DestroyMatrix(g0(iBlocks),sMyName,io)


    do i=1,iBlocks
      call DestroyMatrix(h0(i),sMyName,io)
    enddo
    do i=1,iBlocks-1
      call DestroyMatrix(h1(i),sMyName,io)
      call DestroyMatrix(hm1(i),sMyName,io)
    enddo
    call DestroyArray(g0,sMyName,io)
    call DestroyArray(h0,sMyName,io)
    call DestroyArray(h1,sMyName,io)
    call DestroyArray(hm1,sMyName,io)

  end subroutine InvertSparseON

subroutine InvertSparseONv3(N1,gfsparse,nl,nr,gfout,opindex)
    character(len=*), parameter :: sMyName="InvertSparseONv3"
    type(matrixTypeGeneral) :: gfsparse
    type(matrixTypeGeneral) :: gfout
    integer, intent(in) :: N1,nl,nr,opindex
    type(ioType) :: io
    type(matrixSparseType), allocatable :: h0(:), h1(:), hm1(:)
    integer :: iBlocks,i1,i
    type(matrixType), allocatable :: sigmaL(:),gi1(:),gin(:),mm(:)
    type(matrixType), allocatable ::  g0(:) 
    type(matrixType) :: g1n

    io%isDebug=.false.

    call Row2ColS(gfsparse%MatSparse,io)
    call PartitionMatrix(h0,h1,hm1,iBlocks,gfsparse,nl,nr,io)
    call FillBlocksFromMatrixSparse2Sparse(h0,h1,hm1,iBlocks,gfsparse,io)
!
    call AllocateArray(iBlocks,sigmaL,sMyName,io)
    call AllocateArray(iBlocks,g0,sMyName,io)
    do i=1,iBlocks
      call AllocateMatrix(h0(i)%iRows,h0(i)%iCols,h0(i)%iHorz,h0(i)%iVert,sigmaL(i),sMyName,io)
      sigmaL(i)%a = kczero
    enddo
!
!
    call InverseDiagonalOffdiagonalBlocksSparse(h0,h1,hm1,g0,sigmaL,iBlocks,gfsparse,opindex,io)
!  Gilles: here we can call return
! ! this does compute the diagonal and offdiagonal blocks in one go.
!
    do i=1,iBlocks
      call DestroyMatrix(sigmaL(i),sMyName,io)
    enddo
    call DestroyArray(sigmaL,sMyName,io)

    if(opindex==2.or.opindex==3)then


      gfout%matdense%a=0.0_kdp
      call AllocateArray(iBlocks,gin,sMyName,io)
      call AllocateArray(iBlocks-1,mm,sMyName,io)

      do i=1,iBlocks
          call AllocateMatrix(h0(i)%iRows,h0(iBlocks)%iCols,h0(iBlocks)%iHorz,h0(i)%iVert,gin(i),sMyName,io)
          gin(i)%a = kczero
          if (i<=iBlocks-1) then
            call AllocateMatrix(h0(i)%iRows,h0(i+1)%iCols,h0(i)%iHorz,h0(1)%iVert,mm(i),sMyName,io)
            mm(i)%a = kczero
          endif
      enddo
!

      call InverseLastColumnBlocksSparse(h0,h1,hm1,gin,g0(iBlocks),mm,iBlocks,io)

      do i=1,iBlocks-1
        call DestroyMatrix(mm(i),sMyName,io)
      enddo
      call DestroyArray(mm,sMyName,io)

      call CopyDenseBlocksShift(gfout%matdense%a,n1,nl+nr,gin,iBlocks,0,n1-nr-nl)
!
      do i=1,iBlocks-1
        call DestroyMatrix(gin(i),sMyName,io)
      enddo
      call DestroyMatrix(gin(iBlocks),sMyName,io)
      call DestroyArray(gin,sMyName,io)



      call AllocateArray(iBlocks,gi1,sMyName,io)
      call AllocateArray(iBlocks-1,mm,sMyName,io)
      do i=1,iBlocks
          call AllocateMatrix(h0(i)%iRows,h0(1)%iCols,h0(1)%iHorz,h0(i)%iVert,gi1(i),sMyName,io)
          gi1(i)%a = kczero
          if (i>=2) then
            call AllocateMatrix(h0(i)%iRows,h0(i-1)%iCols,h0(i)%iHorz,h0(1)%iVert,mm(i-1),sMyName,io)
            mm(i-1)%a = kczero
          endif
      enddo

     call InverseFirstColumnBlocksSparse(h0,h1,hm1,gi1,g0(1),mm,iBlocks,io)
     call CopyDenseBlocksShift(gfout%matdense%a,n1,nl+nr,gi1,iBlocks,0,0)
!
      do i=1,iBlocks-1
        call DestroyMatrix(gi1(i),sMyName,io)
        call DestroyMatrix(mm(i),sMyName,io)
      enddo
      call DestroyMatrix(gi1(iBlocks),sMyName,io)
      call DestroyArray(gi1,sMyName,io)
      call DestroyArray(mm,sMyName,io)

    endif

    call DestroyMatrix(g0(iBlocks),sMyName,io)

    if(opindex==4.or.opindex==5)then

      gfout%matdense%a = kczero

      call AllocateMatrix(h0(iBlocks)%iRows,h0(1)%iCols,h0(1)%iHorz,h0(iBlocks)%iVert,g1n,sMyName,io)
      g1n%a = kczero

      call AllocateArray(iBlocks-1,mm,sMyName,io)
      do i=2,iBlocks
        call AllocateMatrix(h0(i)%iRows,h0(i-1)%iCols,h0(i)%iHorz,h0(1)%iVert,mm(i-1),sMyName,io)
        mm(i-1)%a = kczero
      enddo

      call Inverse1NBlocksSparse(h0,h1,hm1,g1n,g0(1),mm,iBlocks,io)
      gfout%matdense%a=g1n%a

      call DestroyMatrix(g1n,sMyName,io)
      do i=1,iBlocks-1
        call DestroyMatrix(mm(i),sMyName,io)
      enddo
      call DestroyArray(mm,sMyName,io)


    endif

    call DestroyMatrix(g0(1),sMyName,io)


    do i=1,iBlocks
      call DestroyMatrixCCS(h0(i),sMyName,io)
    enddo
    do i=1,iBlocks-1
      call DestroyMatrixCCS(h1(i),sMyName,io)
      call DestroyMatrixCCS(hm1(i),sMyName,io)
    enddo
    call DestroyArray(g0,sMyName,io)
    call DestroyArray(h0,sMyName,io)
    call DestroyArray(h1,sMyName,io)
    call DestroyArray(hm1,sMyName,io)

  end subroutine InvertSparseONv3



subroutine InvertSparseONv2(N1,gfsparse,nl,nr,gfout,opindex)
    character(len=*), parameter :: sMyName="InvertSparseONv2"
    type(matrixTypeGeneral) :: gfsparse
    type(matrixTypeGeneral) :: gfout
    integer, intent(in) :: N1,nl,nr,opindex
    type(ioType) :: io
    type(matrixSparseType), allocatable :: h0(:), h1(:), hm1(:)
    integer :: iBlocks,i1,i
    type(matrixType), allocatable :: sigmaL(:),gi1(:),gin(:),mm(:)
    type(matrixType), allocatable ::  g0(:),gm1(:),g1(:) ! these should become sparse in stage 2
    type(matrixType), allocatable :: sigmaR(:),ml(:) ! these two should be phased out

    io%isDebug=.false.

    call Row2ColS(gfsparse%MatSparse,io)
    call PartitionMatrix(h0,h1,hm1,iBlocks,gfsparse,nl,nr,io)
    call FillBlocksFromMatrixSparse2Sparse(h0,h1,hm1,iBlocks,gfsparse,io)
!
    call AllocateArray(iBlocks,sigmaL,sMyName,io)
    call AllocateArray(iBlocks,sigmaR,sMyName,io)
    call AllocateArray(iBlocks,g0,sMyName,io)
    do i=1,iBlocks
      call AllocateMatrix(h0(i)%iRows,h0(i)%iCols,h0(i)%iHorz,h0(i)%iVert,sigmaL(i),sMyName,io)
      sigmaL(i)%a = kczero
      call AllocateMatrix(h0(i)%iRows,h0(i)%iCols,h0(i)%iHorz,h0(i)%iVert,sigmaR(i),sMyName,io)
      sigmaR(i)%a = kczero
      call AllocateMatrix(h0(i)%iRows,h0(i)%iCols,h0(i)%iHorz,h0(i)%iVert,g0(i),sMyName,io)
      g0(i)%a = kczero
    enddo
!
    call InverseDiagonalBlocksSparse(h0,h1,hm1,g0,sigmaL,sigmaR,iBlocks,io)
! ! this should compute the diagonal and offdiagonal blocks in one go.
!
    do i=1,iBlocks
      call DestroyMatrix(sigmaR(i),sMyName,io)
    enddo
    call DestroyArray(sigmaR,sMyName,io)
!
!
    if(opindex==1.or.opindex==3)then

      call AllocateArray(iBlocks-1,g1,sMyName,io)
      call AllocateArray(iBlocks-1,gm1,sMyName,io)

      do i=1,iBlocks-1
        call AllocateMatrix(h1(i)%iRows,h1(i)%iCols,h1(i)%iHorz,h1(i)%iVert,g1(i),sMyName,io)
        call AllocateMatrix(hm1(i)%iRows,hm1(i)%iCols,hm1(i)%iHorz,hm1(i)%iVert,gm1(i),sMyName,io)
      enddo
!
       call InverseOffDiagonalBlocksSparse(h0,h1,hm1,g0,g1,gm1,sigmaL,iBlocks,io)
!
      call CopySparseBlocks(gfsparse,g0,iBlocks)
      call CopySparseBlocks(gfsparse,g1,iBlocks-1)
      call CopySparseBlocks(gfsparse,gm1,iBlocks-1)
!
      do i=1,iBlocks-1
        call DestroyMatrix(g1(i),sMyName,io)
        call DestroyMatrix(gm1(i),sMyName,io)
      enddo
      call DestroyArray(g1,sMyName,io)
      call DestroyArray(gm1,sMyName,io)

    endif
!
    do i=2,iBlocks-1
      call DestroyMatrix(g0(i),sMyName,io)
    enddo

    do i=1,iBlocks
      call DestroyMatrix(sigmaL(i),sMyName,io)
    enddo
    call DestroyArray(sigmaL,sMyName,io)

    if(opindex==2.or.opindex==3)then


      gfout%matdense%a=0.0_kdp
      call AllocateArray(iBlocks,gin,sMyName,io)
      call AllocateArray(iBlocks-1,mm,sMyName,io)

      do i=1,iBlocks
          call AllocateMatrix(h0(i)%iRows,h0(iBlocks)%iCols,h0(iBlocks)%iHorz,h0(i)%iVert,gin(i),sMyName,io)
          gin(i)%a = kczero
          if (i<=iBlocks-1) then
            call AllocateMatrix(h0(i)%iRows,h0(i+1)%iCols,h0(i)%iHorz,h0(1)%iVert,mm(i),sMyName,io)
            mm(i)%a = kczero
          endif
      enddo

      call InverseLastColumnBlocksSparse(h0,h1,hm1,gin,g0(iBlocks),mm,iBlocks,io)
      call CopyDenseBlocksShift(gfout%matdense%a,n1,nl+nr,gin,iBlocks,0,n1-nr-nl)
!
      do i=1,iBlocks-1
        call DestroyMatrix(gin(i),sMyName,io)
        call DestroyMatrix(mm(i),sMyName,io)
      enddo
      call DestroyMatrix(gin(iBlocks),sMyName,io)
      call DestroyArray(mm,sMyName,io)
      call DestroyArray(gin,sMyName,io)



      call AllocateArray(iBlocks,gi1,sMyName,io)
      call AllocateArray(iBlocks-1,ml,sMyName,io)
      do i=1,iBlocks
          call AllocateMatrix(h0(i)%iRows,h0(1)%iCols,h0(1)%iHorz,h0(i)%iVert,gi1(i),sMyName,io)
          gi1(i)%a = kczero
          if (i>=2) then
            call AllocateMatrix(h0(i)%iRows,h0(i-1)%iCols,h0(i)%iHorz,h0(1)%iVert,ml(i-1),sMyName,io)
            ml(i-1)%a = kczero
          endif
      enddo

     call InverseFirstColumnBlocksSparse(h0,h1,hm1,gi1,g0(1),ml,iBlocks,io)
     call CopyDenseBlocksShift(gfout%matdense%a,n1,nl+nr,gi1,iBlocks,0,0)
!
      do i=1,iBlocks-1
        call DestroyMatrix(gi1(i),sMyName,io)
        call DestroyMatrix(ml(i),sMyName,io)
      enddo
      call DestroyMatrix(gi1(iBlocks),sMyName,io)
      call DestroyArray(gi1,sMyName,io)
      call DestroyArray(ml,sMyName,io)

    endif

    call DestroyMatrix(g0(1),sMyName,io)
    call DestroyMatrix(g0(iBlocks),sMyName,io)


    do i=1,iBlocks
      call DestroyMatrixCCS(h0(i),sMyName,io)
    enddo
    do i=1,iBlocks-1
      call DestroyMatrixCCS(h1(i),sMyName,io)
      call DestroyMatrixCCS(hm1(i),sMyName,io)
    enddo
    call DestroyArray(g0,sMyName,io)
    call DestroyArray(h0,sMyName,io)
    call DestroyArray(h1,sMyName,io)
    call DestroyArray(hm1,sMyName,io)

  end subroutine InvertSparseONv2



  subroutine InvertDenseON(N1,GF_iter,nl,nr,opindex)

    type(matrixType) :: GF_iter
    integer, intent(in) :: N1,nl,nr,opindex
    character(len=*), parameter :: sMyName="keldyshimag"
    type(ioType) :: io
    type(matrixType), allocatable :: h0(:), h1(:), hm1(:)
    integer :: iBlocks
    type(matrixType), allocatable :: sigmaL(:),sigmaR(:),g0(:),gi1(:),gin(:),ml(:),mm(:),gm1(:),g1(:)

    io%isDebug=.false.

    call PartitionMatrix(h0,h1,hm1,iBlocks,gf_iter%a,gf_iter%iRows,gf_iter%iCols,nl,nr,io)
    call FillBlocksFromMatrix(h0,h1,hm1,iBlocks,gf_iter%a)

!notes:
!sparse matrices: h0, h1, hm1,    , sigmaL, sigmaR
!we do not need to allocate g0, g1, gm1: we simply overwrite the original input matrix as we go along
!we do not need to allocate gfout before entering, it can be allocated inside the routine at the end, we can allocate g1i and gin at the end when we need the, note also that we do not need it for the equilibrium case
    call AllocateSpace(iBlocks,sigmaL,sigmaR,g0,g1,gm1,h0,h1,hm1,7,io,gi1=gi1,gin=gin,ml=ml,mm=mm)

    call InverseDiagonalBlocks(h0,h1,hm1,g0,sigmaL,sigmaR,iBlocks,io)

    if(opindex==1.or.opindex==3)then
      call InverseOffDiagonalBlocks(h0,h1,hm1,g0,g1,gm1,sigmaL,iBlocks,io)
      call CopyDenseBlocks(gf_iter%a,n1,g0,iBlocks)
      call CopyDenseBlocks(gf_iter%a,n1,g1,iBlocks-1)
      call CopyDenseBlocks(gf_iter%a,n1,gm1,iBlocks-1)
    endif

!here we can deallocate sigmal, sigmar, and then allocate mm and ml
!out of equilibrium we do not even need g0, and for sure not g1 and gm1

! for the transmission calculation we need only a small portion of G

    if(opindex==2.or.opindex==3)then
      call InverseLastColumnBlocks(h0,h1,hm1,gin,g0(iBlocks),mm,iBlocks,io)
      call InverseFirstColumnBlocks(h0,h1,hm1,gi1,g0(1),ml,iBlocks,io)
      call CopyDenseBlocks(gf_iter%a,n1,gi1,iBlocks)
      call CopyDenseBlocks(gf_iter%a,n1,gin,iBlocks)
    endif

    call DeallocateSpace(iBlocks,sigmaL,sigmaR,g0,g1,gm1,h0,h1,hm1,7,io,gi1=gi1,gin=gin,ml=ml,mm=mm)

  end subroutine InvertDenseON

end module mONInterface
