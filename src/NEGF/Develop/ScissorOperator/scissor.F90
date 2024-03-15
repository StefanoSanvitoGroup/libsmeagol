module ScissorOperator
use mTypes, only: matrixTypeGeneral
implicit none
private

  ! A dummy ScissorOperator module. The original module
  ! relies upon the SIESTA-specific fdf library to read
  ! input files.

  ! global-scope parameters that enable ScissorOperator code in SMEAGOL.
  ! Disable it as there is no actual ScissorOperator code here
  logical, public :: sco=.false.
  logical, public :: SCO_leads=.false.
  logical, public :: SCOAndLast=.false.
  logical, public :: scoSCF=.false.

  integer, public :: SCO_nob
  integer, public :: SCO_istart
  logical, public :: SCOSetHamiltonianBlock=.false.
  double complex, allocatable, public :: SCO_Hblock(:,:,:) ! Contains the subspace

  public :: SCOApplyK, SCOApplyK_nc, SCOLeads

  contains

  subroutine SCOApplyK(nspin, maxnh, H, S, xij, kpoint, node, nodes, nuo, no, numh, listhptr, listh, indxuo,&
                       hgeneralp, sgeneralp, NspinComplexMatrix)
  integer, intent(in) :: nspin, maxnh, node, nodes, nuo, no, numh(nuo), listhptr(nuo), listh(maxnh), indxuo(no)
  double precision, intent(in) :: H(maxnh,nspin), S(maxnh)
  double precision, intent(in) :: xij(3,maxnh)
  double precision, intent(in):: kpoint(3)
  integer, intent(in) :: NspinComplexMatrix
  type(matrixTypeGeneral), intent(in) :: hgeneralp(NspinComplexMatrix), sgeneralp

  CALL stop_not_implemented("[STOP] SCOApplyK is not implemented")
  end subroutine SCOApplyK

  subroutine SCOApplyK_nc(hgeneralp, sgeneralp, nuo, NspinComplexMatrix, nl, nr)
  integer, intent(in) :: nuo, NspinComplexMatrix, nl, nr
  type(matrixTypeGeneral), intent(in) :: hgeneralp(NspinComplexMatrix), sgeneralp

  CALL stop_not_implemented("[STOP] SCOApplyK_nc is not implemented")
  end subroutine SCOApplyK_nc


  subroutine SCOLeads(s0, h0, n, nspin, nl, nr)
  integer, intent(in) :: n, nspin, nl, nr
  double complex, intent(in) :: s0(n,n)
  double complex, intent(inout) :: h0(n,n,nspin)

  CALL stop_not_implemented("[STOP] SCOLeads is not implemented")
  end subroutine SCOLeads

end module ScissorOperator

