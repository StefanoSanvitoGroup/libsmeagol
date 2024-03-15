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
!                   ERRORREAD,
!                   ERRORALLOCATE,
!                   ERRORDEALLOCATE  
! AND
! THE MODULE
!                   MUSEFUL  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
!> \brief general useful routines
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 1st of April 2009, alin.elena@ichec.ie
!> \remarks
!>
!> \todo
!> \bug
module mUseful
  use mConstants
  implicit none
  private
!
  public :: ErrorRead
  public :: ErrorAllocate
  public :: ErrorDeallocate

!
 contains

  subroutine ErrorRead(ierror,caller,sfile,iout)
   integer, intent(in) :: ierror,iout
   character(len=*), intent(in) :: caller, sfile

   if (ierror/=0) then
      write(iout,'(a,a,i0)')trim(caller),": Error reading file "//trim(sfile)//" code: ",ierror
      stop
    endif

  end subroutine ErrorRead

  subroutine ErrorAllocate(ierror,caller,iout)
    integer, intent(in) :: ierror,iout
    character(len=*), intent(in) :: caller
    if (ierror /= 0) then
      write(iout,'(a,i0)')trim(caller)//": Error allocating memory, probably you are out of it, check the error code: ", ierror
      stop
    endif
  end subroutine ErrorAllocate

  subroutine ErrorDeallocate(ierror,caller,iout)
    integer, intent(in) :: ierror,iout
    character(len=*), intent(in) :: caller
    if (ierror /= 0) then
      write(iout,'(a,i0)')trim(caller)//": Error deallocating memory, check the error code: ", ierror
      stop
    endif
  end subroutine ErrorDeallocate
end module mUseful
