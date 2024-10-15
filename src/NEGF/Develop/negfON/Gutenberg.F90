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
!                   WRITEMATRIX,
!                   WRITEBLOCK,
!                   WRITEBLOCKS  
! AND
! THE MODULE
!                   MGUTENBERG  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
!> \brief controls all the output
!> \author Alin M Elena
!> \date 1st of April 2009
!> \remarks
!> \todo
!> \bug
module mGutenberg
  use mConstants
  use mTypes
  implicit none
  private
!
  public :: WriteMatrix
  public :: WriteBlock
  public :: WriteBlocks
!
 contains

!> \brief write out a matrix with a label
!> \author Michael Browne
!> \date 2009
!> \param mat matrixType, the matrix to be printed
!> \param str character, the label to be printed
!> \param io ioType, structure that has the I/O info
!> \remarks
  subroutine WriteMatrix(mat,str,io)
    character (len=*), parameter :: sMyName = "WriteMatrix"
    type(matrixType), intent(in) :: mat
    character(len=*), intent(in) :: str
    type(ioType), intent(in) :: io
    integer :: i,j

    write(io%iout,'(a,a,i0,a,i0)') 'writing real part of matrix: ', trim(str), mat%iRows,"x",mat%iCols
    do i=1,mat%iRows
      do j=1,mat%iCols
        write(io%iout,'(G15.4)',advance="no") real(mat%a(i,j),kdp)
      enddo
      write(io%iout,*)
    end do
    write(io%iout,'(a,a,i0,a,i0)') 'writing complex part of matrix: ', trim(str), mat%iRows,"x",mat%iCols
    do i=1,mat%iRows
      do j=1,mat%iCols
        write(io%iout,'(G15.4)',advance="no") aimag(mat%a(i,j))
      enddo
      write(io%iout,*)
    end do
  end subroutine WriteMatrix

!> \brief write out a subblock of a matrix
!> \author Michael Browne
!> \date 2009
!> \param blk blockType, the matrix to be printed
!> \param str character, the label to be printed
!> \param io ioType, structure that has the I/O info
  subroutine WriteBlock(blk, str,io)
    character (len=*), parameter :: sMyName = "WriteBlock"
    type(matrixType), intent(in) :: blk
    character(len = *), intent(in) :: str
    type(ioType), intent(in) :: io
    integer :: i,j

    write(io%iout,'(a,a,i0,a,i0,a,i0,a,i0,a)') 'writing real part of subblock: ', &
      trim(str),blk%iRows,"x",blk%iCols,"(",blk%iVert,",",blk%iHorz,")"
    do i=1,blk%iRows
      do j=1,blk%iCols
        write(io%iout,'(G10.4)',advance="no") real(blk%a(i,j),kdp)
      enddo
      write(io%iout,*)
    end do
    write(io%iout,'(a,a,i0,a,i0,a,i0,a,i0,a)') 'writing complex part of subblock: ', &
      trim(str),blk%iRows,"x",blk%iCols,"(",blk%iVert,",",blk%iHorz,")"
    do i=1,blk%iRows
      do j=1,blk%iCols
        write(io%iout,'(G10.4)',advance="no") aimag(blk%a(i,j))
      enddo
      write(io%iout,*)
    end do
  end subroutine WriteBlock

!> \brief write out a blocks of the matrix
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 1st of April, 2009
!> \param h0 blockType array containing the blocks
!> \param iBlocks integer no of diagonal blocks
!> \param str character, the label to be printed
!> \param io ioType, structure that has the I/O info
  subroutine WriteBlocks(h0,iBlocks,str,io)
    character (len=*), parameter :: sMyName = "WriteBlocks"
    type(matrixType), intent(inout) :: h0(:)
    type(ioType),intent(inout) :: io
    integer, intent(in) :: iBlocks
    character(len=*), intent(in) :: str
    integer :: i
    character(len=klw) :: str1

    do i=1,iBlocks
      write(str1,'(a,a,i0,a)') trim(str),'(',i,')'
      call WriteBlock(h0(i), str1,io)
    end do
  end subroutine WriteBlocks
end module mGutenberg
