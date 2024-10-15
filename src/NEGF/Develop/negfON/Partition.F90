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
!                   GETHORZSS2,
!                   GETVERTSS2,
!                   GETHORZSS,
!                   GETVERTSS,
!                   GETHORZSSVERTSSSPARSEROWSTORED,
!                   COUNTBLOCKS,
!                   PARTITION,
!                   PARTITIONSPARSE,
!                   PARTITIONANDREADSPARSEROWSTOREDASCII,
!                   GETHORZSSVERTSSSPARSEROWSTOREDBINARY,
!                   PARTITIONANDREADSPARSEROWSTOREDBINARY,
!                   PARTITIONGENERALMATRIX,
!                   PARTITIONMATRIXSPARSEBLOCKS,
!                   PARTITIONDENSEMATRIX,
!                   FILLBLOCKSFROMMATRIX,
!                   FILLBLOCKSFROMMATRIXSPARSE,
!                   FILLBLOCKSFROMMATRIXSPARSE2SPARSE  
! AND
! THE MODULES
!                   MPARTITION  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
 !> \brief partitions a sparse matrix in tridiagonal form
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 1st of April 2009, alin.elena@ichec.ie
!> \remarks
!> \todo
!> \bug
module mPartition
  use mConstants
  use mTypes
  use mMatrixUtil
  use mGutenberg
  use mReadData
  use mUseful
  use negfmod, only: outinfo
  implicit none
  private
!
  public :: PartitionMatrix
  public :: FillBlocksFromMatrix
  public :: FillBlocksFromMatrixSparse
  public :: FillBlocksFromMatrixSparse2Sparse

  interface PartitionMatrix
    module procedure PartitionAndReadSparseRowStoredASCII, PartitionAndReadSparseRowStoredBinary, &
            PartitionDenseMatrix, PartitionGeneralMatrix, PartitionMatrixSparseBlocks
  end interface

!
 contains


!> \brief  This subroutine gets the index of the first and last non zero elements
!> in A along each horizontal row
!> \author Michael Browne and Ivan Rungger
!> \date 2009
!> \param mat matrixType, the matrx that contains the data
!> \param horzss integer Nx2 matrix where the first value is the start and second is the stop index for each row.
!> \param io ioType, structure that has the I/O info
!> \remarks in practice this information will be available from the structure
!>  describing sparse arrays but this function allows one to work with dense
!>  arrays
  subroutine GetHorzss2(mat,horzss,io)
    character (len=*), parameter :: sMyName = "GetHorzss2"
    type(matrixTypeGeneral), intent(in)  :: mat
    integer, intent(inout) :: horzss(:,:)
    type(ioType), intent(in) :: io
    integer:: iRows
    integer :: startidx, stopidx, i, j, line

    startidx = 1
    stopidx = 2

    iRows=mat%iRows

    horzss = 0


    do line = 1, iRows
      horzss(line,startidx) = mat%matSparse%j(mat%matSparse%q(line))
      if(line < iRows)then
        horzss(line,stopidx) = mat%matSparse%j(mat%matSparse%q(line+1)-1)
      else
        horzss(line,stopidx) = mat%matSparse%j(mat%matSparse%nnz)
      endif
    enddo

  end subroutine GetHorzss2
!

!> \brief  This subroutine gets the index of the first and last non zero elements
!> in A along each horizontal row
!> \author Michael Browne and Ivan Rungger
!> \date 2009
!> \param mat matrixType, the matrx that contains the data
!> \param vertss integer Nx2 matrix where the first value is the start and second is the stop index for each row.
!> \param io ioType, structure that has the I/O info
!> \remarks in practice this information will be available from the structure
!>  describing sparse arrays but this function allows one to work with dense
!>  arrays
  subroutine GetVertss2(mat,vertss,io)
    character (len=*), parameter :: sMyName = "GetVertss2"
    type(matrixTypeGeneral), intent(in)  :: mat
    integer, intent(inout) :: vertss(:,:)
    type(ioType), intent(in) :: io
    integer:: iCols,iRows
    integer :: startidx, stopidx, i, j, line,jj,ind

    startidx = 1
    stopidx = 2

    iRows=mat%iRows
    iCols=mat%iCols



    vertss(:,startidx)=iRows
    vertss(:,stopidx)=1

    do line = 1, iRows
      do ind = mat%matSparse%q(line),mat%matSparse%q(line+1)-1
        if(line < vertss(mat%matSparse%j(ind),startidx))then
          vertss(mat%matSparse%j(ind),startidx)=line
        endif
        if(line > vertss(mat%matSparse%j(ind),stopidx))then
          vertss(mat%matSparse%j(ind),stopidx)=line
        endif
      enddo
    enddo

  end subroutine GetVertss2
!


!> \brief  This subroutine gets the index of the first and last non zero elements
!> in A along each horizontal row
!> \author Michael Browne/Alin M Elena/Ivan Rungger
!> \date 2009
!> \param mat matrixType, the matrx that contains the data
!> \param horzss integer Nx2 matrix where the first value is the start and second is the stop index for each row.
!> \param io ioType, structure that has the I/O info
!> \param iRows, iCols integer dimensions of the matrix
!> \remarks in practice this information will be available from the structure
!>  describing sparse arrays but this function allows one to work with dense
!>  arrays
  subroutine GetHorzss(mat,iRows,iCols,horzss,io)
    character (len=*), parameter :: sMyName = "GetHorzss"
    complex(kdp), intent(in) :: mat(:,:)
    integer, intent(in) :: iRows,iCols
    integer, intent(inout) :: horzss(:,:)
    type(ioType), intent(in) :: io
    integer :: startidx, stopidx, i, j, line

    startidx = 1
    stopidx = 2

    horzss = 0
    if (io%isDebug) then
      write(io%iout,*) 'GetHorzss A:'
    endif
    do line = 1, iRows
      j = 1
      do
        if ((abs(mat(line,j)) > (0.0_kdp + keps)) .or. (j>=iCols)) exit
        j = j + 1
      end do
      horzss(line,startidx) = j
      j=iRows
      do
        if ((abs(mat(line,j)) > (0.0_kdp + keps)) .or. (j<=1)) exit
        j = j - 1
      end do
      horzss(line,stopidx) = j
      if (io%isDebug) then
        write(io%iout,*) 'horz: start:', horzss(line,startidx) ,'stop', horzss(line,stopidx)
      endif
    end do
  end subroutine GetHorzss
!
!> \brief This subroutine gets the index of the first and last non zero elements
!> in A along each vertical column
!> \author Michael Browne
!> \date 2009
!> \param mat matrixType, the matrx that contains the data
!> \param vertss integer Nx2 matrix where the first value is the start and second is the stop index for each row.
!> \param iRows, iCols integer dimension of the matrix
!> \param io ioType, structure that has the I/O info
!>\remarks in practice this information will be available from the structure
!> describing sparse arrays but this function allows one to work with dense
!> arrays
  subroutine GetVertss(mat,iRows,iCols,vertss,io)
    character (len=*), parameter :: sMyName = "GetVertss"
    complex(kdp), intent(in) :: mat(:,:)
    integer, intent(in) :: iRows,iCols
    integer, intent(inout) ::vertss(:,:)
    type(ioType), intent(in) :: io
    integer :: startidx, stopidx, j, col

    startidx = 1
    stopidx = 2
    vertss = 0
    if (io%isDebug) then
      write(io%iout,*) 'GetVertss A:'
    endif
    do col = 1, iCols
      j = 1
      do
        if ((abs(mat(j,col)) > 0.0_kdp + keps) .or. (j>=iRows)) exit
        j = j + 1
      end do
      vertss(col,startidx) = j
      j=iRows
      do
        if ((abs(mat(j,col)) > 0.0_kdp + keps) .or. (j<=1)) exit
        j = j - 1
      end do
      vertss(col,stopidx) = j
      if (io%isDebug) then
        write(io%iout,*) 'vert: start:',vertss(col,startidx) ,'stop',vertss(col,stopidx)
      endif
    end do
  end subroutine GetVertss
!

  subroutine GetHorzssVertssSparseRowStored(sFile,horzss,vertss,iRows,iCols,io)
    character (len=*), parameter :: sMyName = "GetHorzssVertssSparseRowStored"
    character(len=klw), intent(in) :: sFile
    type(ioType), intent(inout) :: io
    integer, intent(inout) :: iRows,iCols
    integer, intent(inout),allocatable :: horzss(:,:),vertss(:,:)
    integer :: ierror
    character(len=1) :: dummy1, dummy3
    character(len=klw) :: dummy2
    integer :: iLine(0:2,2)
    integer :: i,j,k

    iLine=0
    open(unit=8, file=trim(sFile),status='OLD',action='read',iostat=ierror)
#ifdef Reading
    call ErrorRead(ierror,sMyName,sfile,io%iout)
#endif
    read(8,*,iostat=ierror)dummy1,dummy2,dummy3,iRows,iCols
#ifdef Reading
    call ErrorRead(ierror,sMyName,sfile,io%iout)
#endif
    if (io%isDebug) then
      write(io%iout,'(a,a,x,a,x,a,x,a,x,i0,x,i0)')"reading: ", trim(sFile), trim(dummy1), trim(dummy2), &
      trim(dummy3), iRows,iCols
    endif
    call AllocateArray(iCols,horzss,sMyName,io)
    call AllocateArray(iRows,vertss,sMyName,io)

    vertss(:,1)=iCols+1
    vertss(:,2)=-1
    ! this line should read the nonzero elements
    read(8,*)
    ierror = 0
    k=0
    read(8,*,iostat=ierror)iLine(k,1),iLine(k,2)
#ifdef Reading
    call ErrorRead(ierror,sMyName,sfile,io%iout)
#endif
    vertss(iLine(k,2),1)=iLine(k,1)
    vertss(iLine(k,2),2)=iLine(k,1)
    k=2
    do while (ierror == 0)
      do
        read(8,*,iostat=ierror)iLine(k,1),iLine(k,2)
        if (ierror==0) then
          if (iLine(k,1)<vertss(iLine(k,2),1)) then
            vertss(iLine(k,2),1)=iLine(k,1)
          endif
          if (iLine(k,1)>vertss(iLine(k,2),2)) then
            vertss(iLine(k,2),2)=iLine(k,1)
          endif
        endif
        if ((ierror==0) .and.(iLine(k,1)==iLine(0,1))) then
          iLine(k-1,1)=iLine(k,1)
          iLine(k-1,2)=iLine(k,2)
        endif
        if ((ierror==0) .and.(iLine(k,1)/=iLine(0,1))) then
          horzss(iLine(0,1),1)=iLine(0,2)
          horzss(iLine(0,1),2)=iLine(k-1,2)
          iLine(0,1)=iLine(k,1)
          iLine(0,2)=iLine(k,2)
          exit
         elseif (ierror /=0 ) then
          horzss(iLine(0,1),1)=iLine(0,2)
          horzss(iLine(0,1),2)=iLine(k-1,2)
          exit
        endif
      enddo
    enddo
    if (.false.) then
      do k=1,iRows
        write(*,'(a,i0,x,i0,x,i0)')"  line: ",k,horzss(k,1),horzss(k,2)
        write(*,'(a,i0,x,i0,x,i0)')"column: ",k,vertss(k,1),vertss(k,2)
      enddo
    endif
    close(8)
  end subroutine GetHorzssVertssSparseRowStored

  subroutine CountBlocks(vertss,horzss,iBlocks,iRows,iCols,nl,nr,io)
    character (len=*), parameter :: sMyName = "CountBlocks"
    integer, intent(inout) :: vertss(:,:),horzss(:,:)
    integer, intent(inout) :: iBlocks
    integer, intent(in) :: iRows,iCols,nl,nr
    type(ioType), intent(inout) :: io

    integer :: itmpH1Horz, itmpH1Vert, itmpHm1Horz, itmpHm1Vert, itmpHorz, itmpVert,itmpH1HorzEnd,itmpHm1VertEnd
    integer :: h0ncols, h0nrows, hm1ncols,hm1nrows,h1ncols,h1nrows,h0size,h0nrowslast
    integer :: startidx,stopidx
    integer :: j,i

!!!!!!!!!!!!!!!!!!!!!!!
!!! compute the number of blocks
! initialise the origins and block counter to 1
    startidx = 1
    stopidx = 2
    iBlocks = 1

    ! this block finds the size of h0(1) which is the starting point
    ! for the rest of the process
    if (io%isDebug) then
      write(io%iout,'(a)') 'Finding h0(1) '
    endif
    ! set h0(1) origin, this is known
    itmpHorz = 1
    itmpVert = 1
    itmpH1Vert=1
    ! find point where horz start of data is no longer in the
    ! origin column i.e an indent has occurred the result h0nrows is
    ! the number of rows in the h0(1) block
    if(nl.eq.0)then
      j = itmpVert
      do
        if (horzss(j,startidx) > itmpHorz) exit
        j = j + 1
        if (j>iRows) exit
      end do
      h0nrows = j -  itmpVert

      ! find point where vert start of data is no longer in the
      ! origin row i.e an indent has occurred the result h0ncols is
      ! the number of columns in the h0(1) block
      j = itmpHorz
      do
        if (vertss(j,startidx) > itmpVert) exit
        j = j + 1
        if (j>iCols) exit
      end do
      h0ncols = j - itmphorz

    else
      h0nrows = nl
      h0ncols = nl
    endif


    ! now work throught the rest of the blocks in a loop stopping when the origin
    ! of the next block goes outside the bounds of the source array
    do
!       if ((itmpHorz >= iRows) .or. ( itmpVert >= iCols)) exit
!        write(*,*) 'Starting pass', iBlocks,h0ncols,h0nrows

        ! size used to represent a square with sides equal to the lognest side
      h0size = max(h0nrows,h0ncols)

        ! temporary print some debug info
        ! Should this be a failure case??? - yes
      if (h0nrows /= h0ncols) then
         write(io%iout,'(a,i0,1x,i0)') 'Error: h0nrows /= h0ncols', h0nrows,h0ncols
      end if
      if (.false.) then
         write (*,'(a,i0)') 'h0 ivert:',itmpVert
         write (*,'(a,i0)') 'h0 ihorz:',itmpHorz
         write (*,'(a,i0)') 'h0nrows:',h0nrows
         write (*,'(a,i0)') 'h0ncols:',h0ncols
         write (*,'(a,i0)') 'h0size:',h0size
      endif

      ! set known h1 hm1 origins based on h0 origin
      itmpH1Vert = itmpVert
      itmpHm1Horz = itmpHorz

        ! set known h1 hm1 origins based on h0 origin and sizes
      h1nrows = h0size
        ! do we need a loop here to determine which is the exterme case?
      itmpH1Horz = itmpHorz + h0size
        !write (*,*) 'prep:',horzss(h0(iBlocks)%verts + h0size-1, stopidx),h1(iBlocks)%ihorz,'+1'

      itmpH1HorzEnd=0
      do i=itmpVert,itmpVert + h0size-1
        if(itmpH1HorzEnd.lt.horzss(i, stopidx))then
          itmpH1HorzEnd=horzss(i, stopidx)
        endif
      enddo
      h1ncols = itmpH1HorzEnd - itmpH1Horz+1
      if(h1ncols==0.and.itmpH1Horz<=iCols)then
        h1ncols=1
      endif

!      write(*,*)"ih1ncols=",h1ncols,iBlocks,itmpH1Horz,h0size,horzss(itmpVert + h0size-1, stopidx)  ! = horzss(itmpVert + h0size-1, stopidx) - itmpH1Horz+1

      hm1ncols = h0size
      itmpHm1Vert = itmpVert + h0size

      itmpHm1VertEnd=0
      do i=itmpHorz,itmpHorz + h0size-1
        if(itmpHm1VertEnd.lt.vertss(i, stopidx))then
          itmpHm1VertEnd=vertss(i, stopidx)
        endif
      enddo
      hm1nrows = itmpHm1VertEnd - itmpHm1Vert+1
      if(hm1nrows==0.and.itmpHm1Vert<=iRows)then
        hm1nrows=1
      endif
!      hm1nrows =vertss(itmpHorz + h0size-1, stopidx) - itmpHm1Vert + 1



!      write(*,*)"ihm1nrows=",hm1nrows,iBlocks
      h1nrows = max(h1nrows, hm1ncols)
      h1ncols = max(h1ncols, hm1nrows)
      hm1nrows = max(hm1nrows, h1ncols)
      hm1ncols = max(hm1ncols, h1nrows)
        ! debug info
      if (.false.) then
        write (*,'(a,i0)') 'h1ivert:',itmpH1Vert
        write (*,'(a,i0)') 'h1 ihorz:',itmpH1Horz
        write (*,'(a,i0)') 'h1nrows:',h1nrows
        write (*,'(a,i0)') 'h1ncols:',h1ncols
      endif

        ! debug info
      if (.false.) then
        write (*,'(a,i0)') 'hm1ivert:',itmpHm1Vert
        write (*,'(a,i0)') 'hm1 ihorz:',itmpHm1Horz
        write (*,'(a,i0)') 'hm1nrows:',hm1nrows
        write (*,'(a,i0)') 'hm1ncols:',hm1ncols
      endif
!      write(*,*)"blocksize=",iBlocks,itmpHorz,itmpVert,h0nrows,h0ncols
      h0nrowslast=h0nrows

      iBlocks = iBlocks + 1
      if ((h1ncols == 0) .or. (hm1nrows==0)) then
!          write (io%iout,'(a,i0)') "++This matrix is not tridiagonal, please check it!!!"
!         write(*,*)"exitingloop",h1ncols,hm1nrows
         exit
      else
!         write(*,*)"notexitingloop",h1ncols,hm1nrows
      endif
        ! increment the block counter

        ! set the origin of the next h0 block based on that and the size of the last
      itmpHorz = itmpHorz + h0ncols
      itmpVert = itmpVert + h0nrows
        ! set the known sizes of the next h0 block
      h0ncols = h1ncols
      h0nrows = hm1nrows
        ! debug info
      if (.false.) then
        write(*,'(a,i0,a,i0)') 'next origin: V: ',  itmpVert, 'H: ', itmpHorz
      endif
    end do

    iBlocks = iBlocks - 1
!    write(*,'(a,i0)') "basic Blocks counted! found ",iBlocks
    if(nr > h0nrowslast)then
      iBlocks=iBlocks-1
    endif
!    write(*,'(a,i0)') "Blocks counted! found ",iBlocks
  end subroutine CountBlocks

!> \brief partition the dense matrix A
!> \author Michael Browne/Ain M Elena/Ivan Rungger
!> \date 2009
!> \param h0 matrixSparseType array containing the diagonal blocks
!> \param h1 matrixSparseType array containing the blocks above diagonal
!> \param hm1 matrixSparseType array containing the blocks below diagonal
!> \param horzss, vertss integer, artays containing information about indeces of the blocks
!> \param iRows, iCols integer dimension of the  matrix
!> \param nl,nr integer dimensions of the first and last block
!> \param iBlocks integer no of diagonal blocks
!> \param io ioType, structure that has the I/O info
!> \remarks
  subroutine Partition(h0,h1,hm1,horzss,vertss,iRows,iCols,nl,nr,iBlocks,io)
  implicit none
    character(len=*),parameter :: sMyName="Partition"
    type(matrixType), intent(inout) :: h0(:), h1(:), hm1(:)
    type(ioType),intent(inout) :: io
    integer, intent(in) :: iRows,iCols,nl,nr,iBlocks
    integer, intent(inout) :: horzss(:,:),vertss(:,:)
    integer :: iB
    integer :: i,j,N
    integer :: vertorig, horzorig, h0nrows, h0ncols, h1nrows, h1ncols, hm1nrows, hm1ncols
    integer :: startidx, stopidx, h0size, ierror
    integer :: itmpHorz, itmpVert,itmpH1HorzEnd,itmpHm1VertEnd
    character(len=klw) :: str

    startidx = 1
    stopidx = 2
    ! initialise the origins and block counter to 1
    iB = 1

    ! this block finds the size of h0(1) which is the starting point
    ! for the rest of the process
    if (io%isDebug) then
      write(io%iout,'(a,i0)') 'Finding h0(1)', iB
    endif
    ! set h0(1) origin, this is known
    h0(iB)%iHorz = 1
    h0(iB)%iVert = 1



    if(nl == 0)then
      ! find point where horz start of data is no longer in the
      ! origin column i.e an indent has occurred the result h0nrows is
      ! the number of rows in the h0(1) block
      j = h0(iB)%iVert
      do
        if (horzss(j,startidx) > h0(iB)%iHorz) exit
        j = j + 1
      end do
      h0nrows = j -  h0(iB)%iVert

      ! find point where vert start of data is no longer in the
      ! origin row i.e an indent has occurred the result h0ncols is
      ! the number of columns in the h0(1) block
      j = h0(iB)%iHorz
      do
        if (vertss(j,startidx) > h0(iB)%iVert) exit
        j = j + 1
      end do
      h0ncols = j - h0(iB)%iHorz
    else
      h0nrows = nl
      h0ncols = nl
    endif
    itmpHorz = h0(iB)%iHorz
    itmpVert = h0(iB)%iVert
!    write(*,*)"h0nrows",h0ncols,h0nrows
    ! now work throught the rest of the blocks in a loop stopping when the origin
    ! of the next block goes outside the bounds of the source array
    do
!       if ((itmpHorz >= iRows) .or. (itmpVert >= iCols)) exit
      if (outinfo) write(12346,*)"ibp0=",ib,h0ncols,h0nrows,iBlocks
      h0(iB)%iHorz = itmpHorz
      h0(iB)%iVert = itmpVert
               ! debug info

!        write(*,'(a,i0,a,i0)') 'next origin: V: ',  h0(iB)%iVert, 'H: ', h0(iB)%iHorz

!        write(*,'(a,i0)') 'Starting pass', iB

        ! size used to represent a square with sides equal to the lognest side
      h0size = max(h0nrows,h0ncols)

        ! temporary print some debug info
        ! Should this be a failure case??? - yes
      if (h0nrows /= h0ncols) then
!         write(io%iout,'(a,i0,1x,i0)') 'Error: h0nrows /= h0ncols', h0nrows,h0ncols
      end if
      if (io%isDebug) then
         write (io%iout,'(a,i0)') 'h0 ivert:',h0(iB)%iVert
         write (io%iout,'(a,i0)') 'h0 ihorz:',h0(iB)%iHorz
         write (io%iout,'(a,i0)') 'h0nrows:',h0nrows
         write (io%iout,'(a,i0)') 'h0ncols:',h0ncols
         write (io%iout,'(a,i0)') 'h0size:',h0size
      endif
        ! allocate h0() based on the rows and cols already determined
      call AllocateMatrix(h0nrows,h0ncols,h0(iB)%iHorz,h0(iB)%iVert,h0(iB),sMyName,io)
      h0(iB)%a = cmplx(0.0_kdp,0.0_kdp,kdp)

      ! set known h1 hm1 origins based on h0 origin
      h1(iB)%iVert = h0(iB)%iVert
      hm1(iB)%iHorz = h0(iB)%iHorz

        ! set known h1 hm1 origins based on h0 origin and sizes
      h1nrows = h0size
        ! do we need a loop here to determine which is the exterme case?
      h1(iB)%iHorz = h0(iB)%iHorz + h0size
        !write (*,*) 'prep:',horzss(h0(iB)%verts + h0size-1, stopidx),h1(iB)%ihorz,'+1'

      itmpH1HorzEnd=0
      do i=h0(iB)%iVert,h0(iB)%iVert + h0size-1
        if(itmpH1HorzEnd.lt.horzss(i, stopidx))then
          itmpH1HorzEnd=horzss(i, stopidx)
        endif
      enddo
      h1ncols = itmpH1HorzEnd - h1(iB)%iHorz +1


      hm1ncols = h0size
      hm1(iB)%iVert = h0(iB)%iVert + h0size

      itmpHm1VertEnd=0
      do i=h0(iB)%iHorz,h0(iB)%iHorz + h0size-1
        if(itmpHm1VertEnd<vertss(i, stopidx))then
          itmpHm1VertEnd=vertss(i, stopidx)
        endif
      enddo
      hm1nrows = itmpHm1VertEnd - hm1(iB)%iVert +1

      h1nrows = max(h1nrows, hm1ncols)
      h1ncols = max(h1ncols, hm1nrows)
      hm1nrows = max(hm1nrows, h1ncols)
      hm1ncols = max(hm1ncols, h1nrows)

      if(iB==iBlocks-2)then
        h1ncols = iCols - h1(iB)%iHorz - nr + 1
        hm1nrows = h1ncols
!        write(*,*)"h1ncolsshift=",h1ncols,iCols , h1(iB)%iHorz , nr
      elseif(iB==iBlocks-1)then
        h1ncols = nr
        hm1nrows = h1ncols
!        write(*,*)"h1ncolsshift2=",h1ncols,iCols , h1(iB)%iHorz , nr
      endif
!      write(*,*)"h1ncols=",h1ncols

      if ((h1ncols == 0) .or. (hm1nrows==0)) then
         exit
      endif
        ! debug info
      if (io%isDebug) then
        write (io%iout,'(a,i0)') 'h1ivert:',h1(iB)%iVert
        write (io%iout,'(a,i0)') 'h1 ihorz:',h1(iB)%iHorz
        write (io%iout,'(a,i0)') 'h1nrows:',h1nrows
        write (io%iout,'(a,i0)') 'h1ncols:',h1ncols
      endif
        ! allocate space for the h1() block
      call AllocateMatrix(h1nrows, h1ncols,h1(iB)%iHorz,h1(iB)%iVert,h1(iB),sMyName,io)
      h1(iB)%a = cmplx(0.0_kdp,0.0_kdp,kdp)

        ! debug info
      if (io%isDebug) then
        write (io%iout,'(a,i0)') 'hm1ivert:',hm1(iB)%iVert
        write (io%iout,'(a,i0)') 'hm1 ihorz:',hm1(iB)%iHorz
        write (io%iout,'(a,i0)') 'hm1nrows:',hm1nrows
        write (io%iout,'(a,i0)') 'hm1ncols:',hm1ncols
      endif
        ! allocate space for the hm1() block
      call AllocateMatrix(hm1nrows,hm1ncols,hm1(iB)%iHorz,hm1(iB)%iVert,hm1(iB),sMyName,io)
      hm1(iB)%a = cmplx(0.0_kdp,0.0_kdp,kdp)

        ! increment the block counter
      iB = iB + 1
        ! set the origin of the next h0 block based on that and the size of the last
     itmpHorz = h0(iB-1)%iHorz + h0ncols
     itmpVert = h0(iB-1)%iVert + h0nrows

        ! set the known sizes of the next h0 block
      h0ncols = h1ncols
      h0nrows = hm1nrows
    end do
  end subroutine Partition

!> \brief partition the dense matrix A
!> \author Michael Browne/Alin M Elena/Ivan Rungger
!> \date 2009, 11/02/2010
!> \param h0 matrixSparseType array containing the diagonal blocks
!> \param h1 matrixSparseType array containing the blocks above diagonal
!> \param hm1 matrixSparseType array containing the blocks below diagonal
!> \param horzss, vertss integer, artays containing information about indeces of the blocks
!> \param iRows, iCols integer dimension of the  matrix
!> \param nl,nr integer dimensions of the first and last block
!> \param iBlocks integer no of diagonal blocks
!> \param io ioType, structure that has the I/O info
!> \remarks
  subroutine PartitionSparse(h0,h1,hm1,horzss,vertss,iRows,iCols,nl,nr,iBlocks,io)
    character(len=*),parameter :: sMyName="Partition"
    type(matrixSparseType), intent(inout) :: h0(:), h1(:), hm1(:)
    type(ioType),intent(inout) :: io
    integer, intent(in) :: iRows,iCols,nl,nr,iBlocks
    integer, intent(inout) :: horzss(:,:),vertss(:,:)
    integer :: iB
    integer :: i,j,N
    integer :: vertorig, horzorig, h0nrows, h0ncols, h1nrows, h1ncols, hm1nrows, hm1ncols
    integer :: startidx, stopidx, h0size, ierror
    integer :: itmpHorz, itmpVert,itmpH1HorzEnd,itmpHm1VertEnd
    character(len=klw) :: str

    startidx = 1
    stopidx = 2
    ! initialise the origins and block counter to 1
    iB = 1

    ! this block finds the size of h0(1) which is the starting point
    ! for the rest of the process
    if (io%isDebug) then
      write(io%iout,'(a,i0)') 'Finding h0(1)', iB
    endif
    ! set h0(1) origin, this is known
    h0(iB)%iHorz = 1
    h0(iB)%iVert = 1



    if(nl==0)then
      ! find point where horz start of data is no longer in the
      ! origin column i.e an indent has occurred the result h0nrows is
      ! the number of rows in the h0(1) block
      j = h0(iB)%iVert
      do
        if (horzss(j,startidx) > h0(iB)%iHorz) exit
        j = j + 1
      end do
      h0nrows = j -  h0(iB)%iVert

      ! find point where vert start of data is no longer in the
      ! origin row i.e an indent has occurred the result h0ncols is
      ! the number of columns in the h0(1) block
      j = h0(iB)%iHorz
      do
        if (vertss(j,startidx) > h0(iB)%iVert) exit
        j = j + 1
      end do
      h0ncols = j - h0(iB)%iHorz
    else
      h0nrows = nl
      h0ncols = nl
    endif
    itmpHorz = h0(iB)%iHorz
    itmpVert = h0(iB)%iVert
!    write(*,*)"h0nrows",h0ncols,h0nrows
    ! now work throught the rest of the blocks in a loop stopping when the origin
    ! of the next block goes outside the bounds of the source array
    do
!       if ((itmpHorz >= iRows) .or. (itmpVert >= iCols)) exit
      if (outinfo) write(12346,*)"ibp=",ib,h0ncols,h0nrows,iBlocks
      h0(iB)%iHorz = itmpHorz
      h0(iB)%iVert = itmpVert
               ! debug info

!        write(*,'(a,i0,a,i0)') 'next origin: V: ',  h0(iB)%iVert, 'H: ', h0(iB)%iHorz

!        write(*,'(a,i0)') 'Starting pass', iB

        ! size used to represent a square with sides equal to the lognest side
      h0size = max(h0nrows,h0ncols)

        ! temporary print some debug info
        ! Should this be a failure case??? - yes
      if (h0nrows /= h0ncols) then
!         write(io%iout,'(a,i0,1x,i0)') 'Error: h0nrows /= h0ncols', h0nrows,h0ncols
      end if
      if (io%isDebug) then
         write (io%iout,'(a,i0)') 'h0 ivert:',h0(iB)%iVert
         write (io%iout,'(a,i0)') 'h0 ihorz:',h0(iB)%iHorz
         write (io%iout,'(a,i0)') 'h0nrows:',h0nrows
         write (io%iout,'(a,i0)') 'h0ncols:',h0ncols
         write (io%iout,'(a,i0)') 'h0size:',h0size
      endif
        ! allocate h0() based on the rows and cols already determined
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Now only space for non zero elements should be allocated
!       call AllocateMatrix(h0nrows,h0ncols,h0(iB)%iHorz,h0(iB)%iVert,h0(iB),sMyName,io)
!       h0(iB)%a = cmplx(0.0_kdp,0.0_kdp,kdp)
      ! for the moment I set up only the dimension of the matrix, the allocation should take place when we know the
! number of non-zero entries
        h0(iB)%iRows=h0nrows;h0(iB)%iCols=h0ncols

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! set known h1 hm1 origins based on h0 origin
      h1(iB)%iVert = h0(iB)%iVert
      hm1(iB)%iHorz = h0(iB)%iHorz

        ! set known h1 hm1 origins based on h0 origin and sizes
      h1nrows = h0size
        ! do we need a loop here to determine which is the exterme case?
      h1(iB)%iHorz = h0(iB)%iHorz + h0size
!      write (12347,*) 'prep:',h0(iB)%iHorz, h1(iB)%ihorz,h0size

      itmpH1HorzEnd=0
      do i=h0(iB)%iVert,h0(iB)%iVert + h0size-1
        if(itmpH1HorzEnd<horzss(i, stopidx))then
          itmpH1HorzEnd=horzss(i, stopidx)
        endif
      enddo
      h1ncols = itmpH1HorzEnd - h1(iB)%iHorz +1

      if(h1ncols==0.and.h1(iB)%iHorz<=iCols)then
        h1ncols=1
      endif


      hm1ncols = h0size
      hm1(iB)%iVert = h0(iB)%iVert + h0size

      itmpHm1VertEnd=0
      do i=h0(iB)%iHorz,h0(iB)%iHorz + h0size-1
        if(itmpHm1VertEnd<vertss(i, stopidx))then
          itmpHm1VertEnd=vertss(i, stopidx)
        endif
      enddo
      hm1nrows = itmpHm1VertEnd - hm1(iB)%iVert +1

      if(hm1nrows==0.and.hm1(iB)%iVert<=iRows)then
        hm1nrows=1
      endif

      h1nrows = max(h1nrows, hm1ncols)
      h1ncols = max(h1ncols, hm1nrows)
      hm1nrows = max(hm1nrows, h1ncols)
      hm1ncols = max(hm1ncols, h1nrows)

      if(iB==iBlocks-2)then
        h1ncols = iCols - h1(iB)%iHorz - nr + 1
        hm1nrows = h1ncols
!        write(*,*)"h1ncolsshift=",h1ncols,iCols , h1(iB)%iHorz , nr
      elseif(iB==iBlocks-1)then
        h1ncols = nr
        hm1nrows = h1ncols
!        write(*,*)"h1ncolsshift2=",h1ncols,iCols , h1(iB)%iHorz , nr
      endif
!      write(*,*)"h1ncols=",h1ncols

      if ((h1ncols == 0) .or. (hm1nrows==0)) then
         exit
      endif
        ! debug info
      if (io%isDebug) then
        write (io%iout,'(a,i0)') 'h1ivert:',h1(iB)%iVert
        write (io%iout,'(a,i0)') 'h1 ihorz:',h1(iB)%iHorz
        write (io%iout,'(a,i0)') 'h1nrows:',h1nrows
        write (io%iout,'(a,i0)') 'h1ncols:',h1ncols
      endif
        ! allocate space for the h1() block
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Now only space for non zero elements should be allocated
!       call AllocateMatrix(h1nrows, h1ncols,h1(iB)%iHorz,h1(iB)%iVert,h1(iB),sMyName,io)
!       h1(iB)%a = cmplx(0.0_kdp,0.0_kdp,kdp)
      ! for the moment I set up only the dimension of the matrix, the allocation should take place when we know the
! number of non-zero entries
        h1(iB)%iRows=h1nrows; h1(iB)%iCols=h1ncols

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! debug info
      if (io%isDebug) then
        write (io%iout,'(a,i0)') 'hm1ivert:',hm1(iB)%iVert
        write (io%iout,'(a,i0)') 'hm1 ihorz:',hm1(iB)%iHorz
        write (io%iout,'(a,i0)') 'hm1nrows:',hm1nrows
        write (io%iout,'(a,i0)') 'hm1ncols:',hm1ncols
      endif
        ! allocate space for the hm1() block
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Now only space for non zero elements should be allocated
!       call AllocateMatrix(hm1nrows,hm1ncols,hm1(iB)%iHorz,hm1(iB)%iVert,hm1(iB),sMyName,io)
!       hm1(iB)%a = cmplx(0.0_kdp,0.0_kdp,kdp)
! for the moment I set up only the dimension of the matrix, the allocation should take place when we know the
! number of non-zero entries
        hm1(iB)%iRows=hm1nrows;hm1(iB)%iCols=hm1ncols
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! increment the block counter
      iB = iB + 1
        ! set the origin of the next h0 block based on that and the size of the last
     itmpHorz = h0(iB-1)%iHorz + h0ncols
     itmpVert = h0(iB-1)%iVert + h0nrows

        ! set the known sizes of the next h0 block
      h0ncols = h1ncols
      h0nrows = hm1nrows
    end do
  end subroutine PartitionSparse

  subroutine PartitionAndReadSparseRowStoredASCII(sFile,h0,h1,hm1,iBlocks,iRows,iCols,nl,nr,io)
    character (len=*), parameter :: sMyName = "PartitionAndReadSparseRowStoredASCII"
    character(len=klw), intent(inout) :: sfile
    type(ioType), intent(inout) :: io
    integer, intent(inout) :: iBlocks,iRows,iCols
    integer, intent(in) :: nl,nr
    type(matrixType), intent(inout),allocatable :: h0(:),h1(:),hm1(:)
    integer, allocatable ::  horzss(:,:), vertss(:,:)
    integer :: ierror

    type(matrixType) :: tmpA
    integer :: i

!    write(*,*)"ascii",irows
    call GetHorzssVertssSparseRowStored(sFile,horzss,vertss,iRows,iCols,io)
!    write(*,*)"ascii2",irows,icols
    call CountBlocks(vertss,horzss,iBlocks,iRows,iCols,nl,nr,io)
!    write(*,*)"ascii3",iRows,iBlocks,iCols

    call AllocateArray(iBlocks,h0,sMyName,io)
! apparently we allocate more memory but is not the case. check the partition algorithm
    call AllocateArray(iBlocks,h1,sMyName,io)
    call AllocateArray(iBlocks,hm1,sMyName,io)
    call Partition(h0,h1,hm1,horzss,vertss,iRows,iCols,nl,nr,iBlocks,io)
!    write(*,*)"ascii4"

    call ReadSparseRowStored(sFile,tmpA,io)
!    write(*,*)"ascii5"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,iBlocks
      call FillBlock(tmpA, h0(i))
    end do
    do i=1,iBlocks-1
      call FillBlock(tmpA, hm1(i))
      call FillBlock(tmpA, h1(i))
    end do

    call DestroyMatrix(tmpA,sMyName,io)
    call DestroyArray(horzss,sMyName,io)
    call DestroyArray(vertss,sMyName,io)
!    write(*,*)"ascii6"

  end subroutine PartitionAndReadSparseRowStoredASCII

  subroutine GetHorzssVertssSparseRowStoredBinary(sFile,horzss,vertss,iRows,iCols,magic,endian,io)
    character (len=*), parameter :: sMyName = "GetHorzssVertssSparseRowStoredBinary"
    character(len=klw), intent(in) :: sFile
    type(ioType), intent(inout) :: io
    integer, intent(inout) :: iRows,iCols
    integer, intent(inout),allocatable :: horzss(:,:),vertss(:,:)
    integer, intent(in) :: magic
    character(len=*), intent(in) :: endian

    integer :: ierror,iNonZero
    character(len=1) :: dummy1, dummy3
    character(len=klw) :: dummy2
    integer :: iLine(0:2,2)
    integer :: i,j,k,imagic, maxNzRows
    integer, allocatable:: nzRows(:), tmpRow(:)

  open(unit=13, file=trim(sfile),  status='old', access='stream',  form='unformatted',convert=trim(endian),action='read',iostat=ierror)
#ifdef Reading
    call ErrorRead(ierror,sMyName,sfile,io%iout)
#endif
    read(13,iostat=ierror) imagic
#ifdef Reading
    call ErrorRead(ierror,sMyName,sfile,io%iout)
#endif

    if (io%isDebug) then
      write(io%iout,'(a,i0)')trim(sMyName)//"cookie: ", imagic
    endif
    if (magic /= imagic) then
      write(io%iout,'(a,i0,a,i0)')"Wrong cookie read: ", imagic, "expected ", magic
      stop
    endif
    read(13,iostat=ierror)iRows,iCols,iNonZero
#ifdef Reading
    call ErrorRead(ierror,sMyName,sfile,io%iout)
#endif
    if (io%isDebug) then
      write(io%iout,'(a,i0,x,i0,x,i0)')trim(sMyName)//"rows, cols, nonzero: ", iRows, iCols, iNonZero
    endif
    call AllocateArray(iCols,horzss,sMyName,io)
    call AllocateArray(iRows,vertss,sMyName,io)

    vertss(:,1)=iCols+1
    vertss(:,2)=-1
    call AllocateArray(iRows,nzrows,sMyName,io)
    read(13,iostat=ierror) nzRows
#ifdef Reading
    call ErrorRead(ierror,sMyName,sfile,io%iout)
#endif
    if (io%isDebug) then
      write(io%iout,'(a)', advance="no")"non zero: "
      do i=1,iRows
        write(io%iout,'(i0,x)', advance="no")nzRows(i)
      enddo
      write(io%iout,*)
    endif
    maxNzRows=maxval(nzRows)
    call AllocateArray(maxNzRows,tmpRow,sMyName,io)
    do i=1,iRows
      read(13,iostat=ierror)tmpRow(1:nzRows(i))
#ifdef Reading
      call ErrorRead(ierror,sMyName,sfile,io%iout)
#endif

     if (io%isDebug) then
      write(io%iout,'(a,i0,a)', advance="no")"row ",i,": "
      do k=1,nzRows(i)
        write(io%iout,'(i0,x)', advance="no")tmpRow(k)
      enddo
      write(io%iout,*)
    endif
     tmpRow=tmpRow+1 ! in file they are written in c style
      horzss(i,1)=tmpRow(1)
      horzss(i,2)=tmpRow(nzRows(i))
      do j=1,nzRows(i)
        if(vertss(tmpRow(j),1) > i) then
          vertss(tmpRow(j),1) = i
        endif
        if(vertss(tmpRow(j),2) < i) then
          vertss(tmpRow(j),2) = i
        endif
      enddo
    enddo
    if (io%isDebug) then
      do k=1,iRows
        write(io%iout,'(a,i0,x,i0,x,i0)')"  line: ",k,horzss(k,1),horzss(k,2)
        write(io%iout,'(a,i0,x,i0,x,i0)')"column: ",k,vertss(k,1),vertss(k,2)
      enddo
    endif
    call DestroyArray(tmpRow,sMyName,io)
    call DestroyArray(nzRows,sMyName,io)
    close(13)
  end subroutine GetHorzssVertssSparseRowStoredBinary


  subroutine PartitionAndReadSparseRowStoredBinary(sFile,h0,h1,hm1,iBlocks,iRows,iCols,nl,nr,magic,endian,io)
    character (len=*), parameter :: sMyName = "PartitionAndReadSparseRowStoredBinary"
    character(len=klw), intent(inout) :: sfile
    type(ioType), intent(inout) :: io
    integer, intent(inout) :: iBlocks,iRows,iCols
    integer, intent(in) :: nl,nr
    type(matrixType), intent(inout),allocatable :: h0(:),h1(:),hm1(:)
    character(len=*), intent(in) :: endian
    integer :: magic
    integer, allocatable ::  horzss(:,:), vertss(:,:)
    integer :: ierror

    type(matrixType) :: tmpA
    integer :: i

!    write(*,*)"binary"
    stop
    call GetHorzssVertssSparseRowStoredBinary(sFile,horzss,vertss,iRows,iCols,magic,endian,io)
    call CountBlocks(vertss,horzss,iBlocks,iRows,iCols,nl,nr,io)

    call AllocateArray(iBlocks,h0,sMyName,io)
! apparently we allocate more memory but is not the case. check the partition algorithm
    call AllocateArray(iBlocks,h1,sMyName,io)
    call AllocateArray(iBlocks,hm1,sMyName,io)
    call Partition(h0,h1,hm1,horzss,vertss,iRows,iCols,nl,nr,iBlocks,io)

    call ReadBinaryPestsc(sfile,tmpA,magic,endian,io)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1,iBlocks
      call FillBlock(tmpA, h0(i))
    end do
    do i=1,iBlocks-1
      call FillBlock(tmpA, hm1(i))
      call FillBlock(tmpA, h1(i))
    end do

    call DestroyMatrix(tmpA,sMyName,io)
    call DestroyArray(horzss,sMyName,io)
    call DestroyArray(vertss,sMyName,io)
  end subroutine PartitionAndReadSparseRowStoredBinary


  subroutine PartitionGeneralMatrix(h0,h1,hm1,iBlocks,matA,nl,nr,io)
    character (len=*), parameter :: sMyName = "PartitionGeneralMatrix"
    type(ioType), intent(inout) :: io
    type(matrixTypeGeneral) :: matA
    integer, intent(inout) :: iBlocks
    integer, intent(in) :: nl,nr
    type(matrixType), intent(inout),allocatable :: h0(:),h1(:),hm1(:)

    integer, allocatable ::  horzss(:,:), vertss(:,:)
    integer :: ierror,iCols,iRows

    integer :: i

    iCols=matA%iCols
    iRows=matA%iRows
    call AllocateArray(iCols,horzss,sMyName,io)
    call AllocateArray(iRows,vertss,sMyName,io)
    call GetHorzss2(matA,horzss,io)
!    do i=1,irows
!      write(12347,*)"deltah=",i,horzss(i,1),horzss(i,2)
!    enddo
    call GetVertss2(matA,vertss,io)
    call CountBlocks(vertss,horzss,iBlocks,iRows,iCols,nl,nr,io)

    call AllocateArray(iBlocks,h0,sMyName,io)
    call AllocateArray(iBlocks,h1,sMyName,io)
    call AllocateArray(iBlocks,hm1,sMyName,io)
    call Partition(h0,h1,hm1,horzss,vertss,iRows,iCols,nl,nr,iBlocks,io)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call DestroyArray(horzss,sMyName,io)
    call DestroyArray(vertss,sMyName,io)
  end subroutine PartitionGeneralMatrix

  subroutine PartitionMatrixSparseBlocks(h0,h1,hm1,iBlocks,matA,nl,nr,io)
    character (len=*), parameter :: sMyName = "PartitionGeneralMatrix"
    type(ioType), intent(inout) :: io
    type(matrixTypeGeneral) :: matA
    integer, intent(inout) :: iBlocks
    integer, intent(in) :: nl,nr
    type(matrixSparseType), intent(inout),allocatable :: h0(:),h1(:),hm1(:)

    integer, allocatable ::  horzss(:,:), vertss(:,:)
    integer :: ierror,iCols,iRows

    integer :: i

    iCols=matA%iCols
    iRows=matA%iRows
    call AllocateArray(iCols,horzss,sMyName,io)
    call AllocateArray(iRows,vertss,sMyName,io)
    call GetHorzss2(matA,horzss,io)
!    do i=1,irows
!      write(12347,*)"deltah=",i,horzss(i,1),horzss(i,2)
!    enddo
    call GetVertss2(matA,vertss,io)
    call CountBlocks(vertss,horzss,iBlocks,iRows,iCols,nl,nr,io)

    call AllocateArray(iBlocks,h0,sMyName,io)
    call AllocateArray(iBlocks,h1,sMyName,io)
    call AllocateArray(iBlocks,hm1,sMyName,io)
!!!!!!!!!! no space is allocated by the PartitionSparse routine...
    call PartitionSparse(h0,h1,hm1,horzss,vertss,iRows,iCols,nl,nr,iBlocks,io)
    call DestroyArray(horzss,sMyName,io)
    call DestroyArray(vertss,sMyName,io)
  end subroutine PartitionMatrixSparseBlocks


  subroutine PartitionDenseMatrix(h0,h1,hm1,iBlocks,matA,iRows,iCols,nl,nr,io)
    character (len=*), parameter :: sMyName = "PartitionDenseMatrix"
    type(ioType), intent(inout) :: io
    integer, intent(inout) :: iBlocks
    complex(kdp), intent(inout) :: matA(:,:)
    integer, intent(in) :: nl,nr
    integer, intent(inout) :: iRows,iCols
    type(matrixType), intent(inout),allocatable :: h0(:),h1(:),hm1(:)

    integer, allocatable ::  horzss(:,:), vertss(:,:)
    integer :: ierror

    integer :: i
!    write(*,*)"dense"

    call AllocateArray(iCols,horzss,sMyName,io)
    call AllocateArray(iRows,vertss,sMyName,io)
    call GetHorzss(matA,iRows,iCols,horzss,io)
    call GetVertss(matA,iRows,iCols,vertss,io)
    call CountBlocks(vertss,horzss,iBlocks,iRows,iCols,nl,nr,io)

    call AllocateArray(iBlocks,h0,sMyName,io)
    call AllocateArray(iBlocks,h1,sMyName,io)
    call AllocateArray(iBlocks,hm1,sMyName,io)
    call Partition(h0,h1,hm1,horzss,vertss,iRows,iCols,nl,nr,iBlocks,io)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call DestroyArray(horzss,sMyName,io)
    call DestroyArray(vertss,sMyName,io)
  end subroutine PartitionDenseMatrix

   subroutine FillBlocksFromMatrix(h0,h1,hm1,iBlocks,matA)
    character(len=*), parameter :: sMyName="FillBlocksFromMatrix"
    integer, intent(inout) :: iBlocks
    complex(kdp), intent(inout) :: matA(:,:)
    type(matrixType), intent(inout) :: h0(:),h1(:),hm1(:)

    integer :: i

!    write(*,*)"blocks=",iblocks
    do i=1,iBlocks
!      write(*,*)"filling blockb=",i
      call FillBlock(matA, h0(i))
!      write(*,*)"filling blocka=",i
    end do
    do i=1,iBlocks-1
      call FillBlock(matA, hm1(i))
      call FillBlock(matA, h1(i))
    end do

   end subroutine FillBlocksFromMatrix


   subroutine FillBlocksFromMatrixSparse(h0,h1,hm1,iBlocks,matB)
    character(len=*), parameter :: sMyName="FillBlocksFromMatrixSparse"
    integer, intent(inout) :: iBlocks
    type(matrixTypeGeneral), intent(in)  :: matB
    type(matrixType), intent(inout) :: h0(:),h1(:),hm1(:)

    integer :: i

!    write(*,*)"blocks=",iblocks
    do i=1,iBlocks
!      write(*,*)"filling blockb=",i
!      call FillBlock(matA, h0(i))
      call FillBlock(matB, h0(i))
!      write(*,*)"filling blocka=",i
    end do
    do i=1,iBlocks-1
      call FillBlock(matB, hm1(i))
      call FillBlock(matB, h1(i))
!      call FillBlock(matA, hm1(i))
!      call FillBlock(matA, h1(i))
    end do

   end subroutine FillBlocksFromMatrixSparse

  subroutine FillBlocksFromMatrixSparse2Sparse(h0,h1,hm1,iBlocks,matB,io)
    character(len=*), parameter :: sMyName="FillBlocksFromMatrixSparse2Sparse"
    integer, intent(inout) :: iBlocks
    type(matrixTypeGeneral), intent(in)  :: matB
    type(matrixSparseType), intent(inout) :: h0(:),h1(:),hm1(:)
    type(ioType), intent(inout) :: io

    integer :: i

!    write(*,*)"blocks=",iblocks
    do i=1,iBlocks
!      write(*,*)"filling blockb=",i
!      call FillBlock(matA, h0(i))
      call FillBlock(matB%matSparse, h0(i),io)
!      write(*,*)"filling blocka=",i
    end do
    do i=1,iBlocks-1
      call FillBlock(matB%matSparse, hm1(i),io)
      call FillBlock(matB%matSparse, h1(i),io)
    end do

   end subroutine FillBlocksFromMatrixSparse2Sparse


end module mPartition
