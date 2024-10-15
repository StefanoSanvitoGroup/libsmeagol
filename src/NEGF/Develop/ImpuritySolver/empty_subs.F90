!modules
  module mImpuritySolver

  use mConstants
  use mTypes

  public PrintGFHSMatsubaraGeneral
  public AddRhoTilde
  public CalculateSinv
  public RegularizeGFPolesGeneral
  public CTQMCHyb_ImpuritySolverInterface


  private 

  contains 

  subroutine PrintGFHSMatsubaraGeneral(temp,e,ef,gf,hgeneralp,sgeneralp,sigmal,nl,sigmar,nr,writeHeader,node,PrintImpurityGfMatsubara,CallImpuritySolver,ie,ne,ispin,nspin)

  integer, intent(in)                 :: node,nl,nr,ie,ne,ispin,nspin
  complex(kdp), intent(in)            :: e
  real(kdp), intent(in)               :: ef,temp
  type(matrixTypeGeneral), intent(in) :: gf
  type(matrixTypeGeneral), intent(in) :: hgeneralp,sgeneralp
  complex(kdp), intent(in)            :: sigmal(nl,nl),sigmar(nr,nr)
  logical, intent(in)                 :: writeHeader,PrintImpurityGfMatsubara,CallImpuritySolver

  end subroutine PrintGFHSMatsubaraGeneral

  subroutine AddRhoTilde(rhogeneral,gfmattype,nl,nr,set_rho_boundary)

  integer, intent(in)                    :: gfmattype
  type(matrixTypeGeneral), intent(inout) :: rhogeneral
  integer, intent(in) :: nl,nr
  logical, intent(in):: set_rho_boundary
  call stop_not_implemented("Input options related to Boundstates is not implemented yet.")

  end subroutine AddRhoTilde

  subroutine CalculateSinv(sgeneralp,gfmattype)

  integer, intent(in)                    :: gfmattype
  type(matrixTypeGeneral), intent(inout) :: sgeneralp
  call stop_not_implemented("Input options related to Boundstates is not implemented yet.")

  end subroutine CalculateSinv

  subroutine RegularizeGFPolesGeneral(e,ef,gf)

  complex(kdp), intent(in)            :: e
  real(kdp), intent(in)               :: ef
  type(matrixTypeGeneral), intent(inout) :: gf
  call stop_not_implemented("Input options related to Boundstates is not implemented yet.")

  end subroutine RegularizeGFPolesGeneral

  subroutine CTQMCHyb_ImpuritySolverInterface(temp,node,nnodes,comm)

  integer, intent(in)     :: node,nnodes,comm
  real(kdp), intent(in)   :: temp

  call stop_not_implemented("Input options related to Boundstates is not implemented yet.")

  end subroutine CTQMCHyb_ImpuritySolverInterface



  end module mImpuritySolver

  module mImpuritySolverParameters

  private

  public :: read_options_ImpuritySolver

  contains

  subroutine read_options_ImpuritySolver(CallImpuritySolver,norb,mynode)

  logical, intent(in) :: CallImpuritySolver
  integer, intent(in) :: norb
  integer, intent(in) :: mynode

  end subroutine read_options_ImpuritySolver



  end module mImpuritySolverParameters

