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
!                   READDENSEMATRIX,
!                   READSPARSEROWSTORED,
!                   READBINARYPESTSC  
! AND
! THE MODULE
!                   MREADDATA  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
!> \brief deals with the input of data
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 1st of April 2009
!> \remarks
!>
!> \todo
!> \bug
module mReadData
  use mConstants
  use mTypes
  use mMatrixUtil
  use mUseful
  implicit none
  private

  public :: ReadDenseMatrix
  public :: ReadSparseRowStored
  public :: ReadBinaryPestsc
 contains

!> \brief reads a dense matrix from a file
!> \details this subroutine reads in a dense matrix in ascii format
!> into a matrix A which is also dense. The matrix is N by N
!> it also writes it our with implied loops
!> to allow for fortran style storage it transposes the matrix
!> after reading it
!> \author Michael Browne
!> \date 2009
!> \param mat matrixType, the matrix in which we read the data
!> \param sFile charcter the file from which we read the data
!> \param io ioType, structure that has the I/O info
!> \remarks Removed the write statement
  subroutine ReadDenseMatrix(sFile,mat,io)
    character (len=*), parameter :: sMyName = "ReadDenseMatrix"
    type(matrixType), intent(inout) :: mat
    character(len=klw), intent(in) :: sFile
    type(ioType), intent(in) :: io
    integer :: ierror,i,j
    real(kdp) :: aux


    mat%a = cmplx(0.0_kdp,0.0_kdp,kdp)
    open(unit=8, file=trim(sFile),status='OLD',action='read',iostat=ierror)
#ifdef Reading
    call ErrorRead(ierror,sMyName,sfile,io%iout)
#endif
    read(8,*,iostat=ierror)mat%a
#ifdef Reading
    call ErrorRead(ierror,sMyName,sfile,io%iout)
#endif
    close(unit=8)

!     mat%a = transpose(mat%a)

  end subroutine ReadDenseMatrix


!> \brief reads a sparse matrix
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 29th of April, 2009
!> \param mat matrixType, the matrix in which we read the data
!> \param sFile charcter the file from which we read the data
!> \param io ioType, structure that has the I/O info

  subroutine ReadSparseRowStored(sFile,mat,io)
    character(len=*), parameter :: sMyName="ReadSparseRowStored"
    character(len=klw), intent(in) :: sFile
    type(ioType), intent(inout) :: io
    type(matrixType), intent(inout) :: mat

    integer :: ierror
    integer :: iNonZero
    character(len=1) :: dummy1, dummy3
    character(len=klw) :: dummy2
    integer :: k,l,i,iRows,iCols
    real(kdp) :: rex,imx

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
    read(8,*,iostat=ierror)dummy1,dummy2,dummy3,iNonZero
#ifdef Reading
    call ErrorRead(ierror,sMyName,sfile,io%iout)
#endif
    if (io%isDebug) then
      write(io%iout,'(a,a,x,a,x,a,x,a,x,i0)')"reading: ", trim(sFile), trim(dummy1), trim(dummy2), &
      trim(dummy3), iNonZero
    endif
    call AllocateMatrix(iRows,iCols,mat,sMyName,io)
    do i=1,iNonZero
      read(8,*,iostat=ierror)k,l,rex,imx
      if (ierror /= 0) then
        write(io%iout,'(a,i0,a,i0)') 'Error: read() ierror = ', ierror, "at line ", i+2
        stop
      end if
      mat%a(k,l)=cmplx(rex,imx,kdp)
    enddo
    close(8)
  end subroutine ReadSparseRowStored

  subroutine ReadBinaryPestsc(sfile,mat,magic,endian,io)
    character(len=*), parameter :: sMyName="ReadBinaryPestsc"
    type(ioType), intent(inout) :: io
    character(len=*), intent(in) :: sfile,endian
    type(matrixType), intent(inout) :: mat
    integer, intent(in) :: magic

    integer :: imagic, ierror,i,j,k
    integer :: iRows, iCols, iNonzero,maxNzRow
    integer, allocatable :: nzCols(:),nzRows(:)
    complex(kdp), allocatable :: tmpRow(:)

    open(unit=10, file=trim(sfile),  status='old', access='stream',  form='unformatted',convert=trim(endian),action='read')

    read(10,iostat=ierror) imagic
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
    read(10,iostat=ierror)iRows,iCols,iNonZero
#ifdef Reading
    call ErrorRead(ierror,sMyName,sfile,io%iout)
#endif
    if (io%isDebug) then
      write(io%iout,'(a,i0,x,i0,x,i0)')trim(sMyName)//"rows, cols, nonzero: ", iRows, iCols, iNonZero
    endif

    call AllocateMatrix(iRows,iCols,mat,sMyName,io)
    call AllocateArray(iRows,nzRows,sMyName,io)
    call AllocateArray(iNonZero,nzCols,sMyName,io)
    read(10,iostat=ierror) nzRows
#ifdef Reading
    call ErrorRead(ierror,sMyName,sfile,io%iout)
#endif
    if (io%isDebug) then
      do i=1,iRows
        write(io%iout,'(i0,x)', advance="no")nzRows(i)
      enddo
      write(io%iout,*)
    endif
    read(10,iostat=ierror) nzCols
#ifdef Reading
    call ErrorRead(ierror,sMyName,sfile,io%iout)
#endif
    nzCols = nzCols+1 ! are they sure C style???? then I should take care of indeces to avoid silly increments
    if (io%isDebug) then
      j=0
      do k=1,iRows
        write(io%iout,'(a,i0,a)', advance="no")"Col: ", k," : "
        do i=1,nzRows(k)
          write(io%iout,'(i0,x)',advance="no") nzCols(j+i)
        enddo
        j=j+nzRows(k)
        write(io%iout,*)
      enddo
    endif
    maxNzRow = maxval(nzRows)
    call AllocateArray(maxNzRow,tmpRow,sMyName,io)
    j=0
    do k=1,iRows
      read(10,iostat=ierror)tmpRow(1:nzRows(k))
#ifdef Reading
      call ErrorRead(ierror,sMyName,sfile,io%iout)
#endif
      do i=1,nzRows(k)
        mat%a(nzCols(j+i),k)=tmpRow(i)
      enddo
      j=j+nzRows(k)
    enddo
    call DestroyArray(tmpRow,sMyName,io)
    call DestroyArray(nzCols,sMyName,io)
    call DestroyArray(nzRows,sMyName,io)
    close(10)
  end subroutine ReadBinaryPestsc

end module mReadData
