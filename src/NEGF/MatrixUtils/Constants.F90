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
! THE MODULE
!                   MCONSTANTS  
! IN THIS FILE IS LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
!> \brief contains all the constants used through the programme
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 1st of April 2009
!> \remarks
!> \todo
!> \bug
module mConstants
  implicit none
  private

  integer, parameter, public :: kdp = kind (1.0d0)!< defined precision for reals
  integer, parameter, public :: kidp = kind (1)!< defined precision for integers
  integer, parameter, public :: kidp2 = 2*kidp!< defined precision for integers
  integer, parameter, public :: klw = 25 !< the length of a string (word)
  integer, parameter, public :: kll = 255 !< the length of a line
  real(kdp), parameter, public :: keps=10*tiny(1.0_kdp) !< how small is zero
  real(kdp), parameter, public :: kzero=0.0_kdp !< zero
  real(kdp), parameter, public :: kone=1.0_kdp !< one
  character(len=1), parameter, public :: kn="N" !< letter N
  complex(kdp), parameter, public :: kcone=cmplx(kone,kzero,kdp) !< 1+0i
  complex(kdp), parameter, public :: kczero=cmplx(kzero,kzero,kdp) !< 0+0i

end module mConstants
