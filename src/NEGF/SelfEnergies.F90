module mSelfenergies

  use mConstants
  use mTypes

  implicit none
  private
  
  public :: SelfEnergyGeneral
  public :: SetOptionsSelfEnergies
  
  integer, allocatable, save :: ndivxy(:,:)

  contains

  subroutine SetOptionsSelfEnergies(ndivxy_in,nleads)

  integer, intent(in) :: nleads
  integer, intent(in) :: ndivxy_in(nleads,2)

  allocate(ndivxy(nleads,2))
  ndivxy=ndivxy_in

  end subroutine SetOptionsSelfEnergies
      
  subroutine SelfEnergyGeneral(side,n,e,h0,h1,s0,s1,sigma,nchan,delta,DoFourier)

! **********************************************************************
! Calculates the self-energies, based on the singularity-free scheme,
! Ref: I. Rungger and S. Sanvito, Phys Rev B 78, 035407 (2008)
!
! Written by Ivan Rungger, October 2013
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

  use negfmod,       only : TransmissionMatrix,em_Last_SCF_Step
  use mSigmaMethod1, only : GetSelfEnergy
  use mSigmaFourier, only : selfenergy_k
  use mComputeULR,   only : PhiS,InitPhiS

  character(len=1), intent(in) :: side
  integer, intent(in)          :: n
  complex(kdp), intent(in)     :: e
  real(kdp), intent(in)        :: delta
  complex(kdp), intent(inout)  :: h0(n,n),h1(n,n),s0(n,n),s1(n,n) ! within selfenergy_k these matrices can potentially be changed, therefore intent(inout)
  integer, intent(out)         :: nchan
  complex(kdp), intent(out)    :: sigma(n,n)
  logical, intent(in)          :: DoFourier 
   
  integer il
 
  if(side.eq.'L')then
    il=1
  else
    il=2
  endif

  if(TransmissionMatrix.and.em_Last_SCF_Step) call InitPhiS(PhiS(il),side,ndivxy(il,:),n,e)

  if(maxval(ndivxy(il,:))==1)then
    call GetSelfEnergy(side,n,e,h0,h1,s0,s1,sigma,nchan,delta)
  else
    call selfenergy_k(side,ndivxy(il,:),n ,e,h0,h1,s0,s1,sigma ,nchan,delta,DoFourier)
  endif

  end subroutine SelfEnergyGeneral


end module mSelfenergies
