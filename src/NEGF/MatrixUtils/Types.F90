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
!                   MTYPES  
! IN THIS FILE IS LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
!> \brief definitions of the types used
!> \author Alin M Elena, alin.elena@ichec.ie
!> \date 1st of April 2009, alin.elena@ichec.ie
!> \remarks
!>
!> \todo
!> \bug
module mTypes
  use mConstants
  implicit none
  private
!> \brief This structure holds: a dense matrix or a block matrix if iHorz, iCols are
  type, public :: matrixType
    integer :: iRows, iCols !< no of rows and columns of the matrix
    integer :: iHorz=0, iVert=0 !< are the indices of the the top left corner of the data in the case we want to store  a block
    complex(kdp), allocatable:: a(:,:) !< block in the original array
  end type matrixType


!> \brief This structure holds: a matrix/block matrix CCS or CRS sparse format
!> \remark we use here the CCS matrix format. Elements are stored compressed by column, row by row in array a
!> the rows for column j are stored in array i from i(p(j)) to i(p(j+1)-1)
!> CRS Elements are stored compressed by row, column by  column in array b
!> the columns for row i are stored in array j from j(q(i)) to j(q(i+1)-1)
  type, public :: matrixSparseType
    integer :: iRows, iCols !< no of rows and columns of the matrix
    integer :: iHorz=0, iVert=0 !< are the indices of the the top left corner of the data
    integer :: nnz !< no of non-zero elements
!for crs we use i, p, a 
    integer, allocatable:: i(:) !< array containing the row index of the elements, length=nnz
    integer, allocatable:: p(:) !< array of dimension iCols+1 containing the start index of each col, last element is nnz+1, first element is 1
    complex(kdp), allocatable:: a(:) !< array containing the matrix elements, length=nnz, in CCS
!for ccs we use j, q, b
    integer, allocatable:: j(:) !< array containing the column index of the elements, length=nnz
    integer, allocatable:: q(:) !< array of dimension iRows+1 containing the start index of each row, last element is nnz+1, first element is 1
    complex(kdp), allocatable:: b(:) !< array containing the matrix elements, length=nnz in CRS
  end type matrixSparseType

!> \brief This structure holds: a matrix/block matrix CCS or CRS sparse parallel format
!> \remark we use here the CCS matrix format. Elements are stored compressed by column, row by row in array a
!> the rows for column j are stored in array i from i(p(j)) to i(p(j+1)-1)
!> CRS Elements are stored compressed by row, column by  column in array b
!> the columns for row i are stored in array j from j(q(i)) to j(q(i+1)-1)
  type, public :: matrixSparsePType
    integer :: MPICommunicator !< MPI communicator over which the matrix is distributed
    integer :: MPIGroup !< MPI group over which the matrix is distributed
    integer :: nProcs !< no of processors on which the matrix is stored
    integer :: iProc !< processors rank within the MPI communicator or group
    integer :: iRowsGlobal, iColsGlobal !< no of rows and columns of the full matrix
    integer, allocatable :: iHorzGlobal(:), iVertGlobal(:) !< indices of top left corners of the full matrix for all nProcs processors
    type(matrixSparseType)  :: matSparse !< locally stored matrix
  end type matrixSparsePType


!> \brief This structure holds: a general matrix, can be spares or dense
  type, public :: matrixTypeGeneral
    integer :: mattype !< the flag mattype sets the type of the matrix (for mattype==0 the matrix is a dense matrix, default; for mattype==2 the matrix is a sparse matrix; 3 for sparse parallel);
    integer :: iRows, iCols !< no of rows and columns of the matrix
    integer :: iHorz=0, iVert=0 !< are the indices of the the top left corner of the data in the case we want to store  a block
    type(matrixType) :: matDense !< allocated if the matrix is stored in dense format (mattype==0)
    type(matrixSparseType)  :: matSparse !< allocated if the matrix is stored in sparse format (mattype==1)
    type(matrixSparsePType)  :: matSparseP !< allocated if the matrix is stored in sparse parallel format (mattype==2)
!    type(matrixTypeGeneral),  pointer  :: matGeneral !< a matrix of type general, allocated if needed, e.g. for matrix conversion
  end type matrixTypeGeneral

!> \brief the type that holds the info about I/O units and files
  type, public :: ioType
    integer :: iout !< the unit no for the standard output
    character(len=klw) :: sFileOut !< filename to which the standard output is writen
    logical :: isDebug !< toggles debug info printing
    real(kdp) :: memCount !< amount of memory dynamically allocated in kiB
  end type ioType

!> \brief This structure holds the informations related to a single self-energy
  type, public :: SelfEnergyType
    complex(kdp), allocatable:: sigma(:,:) !< the self-energy matrix
    integer :: n !< the dimension of the self-energy
    CHARACTER(LEN=1) :: Side !< the side of the self-energy (either 'L' or 'R')
    complex(kdp) :: e !< energy at which the self-energy is calculated
    integer :: nchannels !< the number of open channels
    integer :: InfoSigma !< informations on how the self-energy is stored
    integer :: node !< the processor on which the self-energy is stored

!    integer :: iSide !< the side of the self-energy (-1 for left side, 1 for right side)
!    integer, allocatable:: gpSigma(:,:) !< group of processors where the selfenergies are stored (dimensions nLeads,nSpin)
  end type SelfEnergyType

!> \brief This structure holds the informations related to an energy grid
  type, public :: EnergyGridType
    complex(kdp), allocatable :: e(:) !< energies at which the self-energy is calculated
    complex(kdp), allocatable :: w(:) !< integral weight of the point
    integer,      allocatable :: ig(:) !< global index of the energy point
    complex(kdp), allocatable :: eGlobal(:) !< energies at which the self-energy is calculated, for all the processors
    type(SelfEnergyType), allocatable :: sigma(:,:,:,:) !< self-energies
    integer :: nEnergies !< number of energy points on the node
    integer :: nEnergiesGlobal !< number of global energy points
    integer :: nSpin !< number of spins
    integer :: nLeads !< number of leads
    integer :: nK !< number of k-points

    real(kdp) :: v !< the bias voltage
    real(kdp) :: deltasigma !< the imaginary  part added when calculating the self-energies
    integer, allocatable :: leadsTotalDim(:)
    integer :: GridType !< type of grid, 1 for adaptive, 0 for standard
    integer :: InfoSigma !< informations on how the self-energy is stored
    character(len=klw) :: sLabel !< filename to which the standard output is writen
    character(len=klw) :: SigmaSuffix !< filename to which the standard output is writen

!    integer :: nEnergiesSigma !< number of energy points on the node
!    integer :: nSpinSigma !< number of spins
!    integer :: nLeadsSigma !< number of leads
!    integer :: nKSigma !< number of k-points

  end type EnergyGridType

!!  Not set variables
!!
!!  ERealGrid, EImagGrid:
!!    complex(kdp), allocatable :: eGlobal(:) !< energies at which the self-energy is calculated, for all the processors
!!
!!  ETransmGrid:
!!    complex(kdp), allocatable :: w(:) !< integral weight of the point

!> \brief This structure holds the informations related to a set of scattering states
  type, public :: ScatteringStates
    integer :: nTstates(2) !< 1 is for incoming states into the scattering region, 2 is for outgoing states out of the scattering region
!xxx make ntstates_in and ntstates_out
    integer :: n !< the dimension of the lead
    complex(kdp), allocatable:: phi_in(:,:) !< the incoming scattering vectors
    complex(kdp), allocatable:: phit_in(:,:) !< the duals of phi_in
    complex(kdp), allocatable:: k_in(:) !< the k vector of phi_in
    complex(kdp), allocatable:: z_in(:) !< the e^(i k) vector of phi_in
    complex(kdp), allocatable:: v_in(:) !< the group velocity of phi_in
    complex(kdp), allocatable:: mu_in(:,:) !< the spin vector of phi_in
    complex(kdp), allocatable:: phiAll_in(:,:) !< the incoming scattering vectors
    complex(kdp), allocatable:: phitAll_in(:,:) !< the duals of phi_in

    complex(kdp), allocatable:: phi_out(:,:) !< the outgoing scattering vectors
    complex(kdp), allocatable:: phit_out(:,:) !< the duals of the outgoing scattering vectors
    complex(kdp), allocatable:: k_out(:) !< the k vector of phi_out
    complex(kdp), allocatable:: z_out(:) !< the e^(i k) vector of phi_out
    complex(kdp), allocatable:: v_out(:) !< the group velocity of phi_out
    complex(kdp), allocatable:: mu_out(:,:) !< the spin vector of phi_out
    complex(kdp), allocatable:: phiAll_out(:,:) !< the outgoing scattering vectors
    complex(kdp), allocatable:: phitAll_out(:,:) !< the duals of the outgoing scattering vectors

    complex(kdp) :: e !< energy at which the scattering states are calculated
    CHARACTER(LEN=1) :: Side !< the side of the lead
!!    real(kdp), allocatable :: k(:) !< kpoint in the x-y plane of the Fourier transform; if allocated with size 2 then it implies that this scattering state is a contribution to a full set of Fourier scattering states
!optional : add kx,ky point at which scattering states are calculated

  end type ScatteringStates


!> \brief This structure holds the informations for the Fourier components of all scattering states
  type, public :: FourierScatteringStates
    integer :: n !< the dimension of the lead
    integer :: ns !< the dimension of the Fourier parts of the lead
    integer :: ndivxy(2) !< number of sigma.nx and sigma.ny for this lead; if sigma.nx and ny is 1, then it is not fourier decomposed
    integer :: ik(2) !< ikx,iky of the element of FourierSstates that needs to be edited

    real(kdp), allocatable :: kxy(:,:,:) !< list of kpoints in the x-y plane of the Fourier transform
    type(ScatteringStates), allocatable :: FourierSstates(:,:) !< fourier Sstates (small system)
    type(ScatteringStates), allocatable :: Sstates(:,:) !< Full s states (large system)
!!    type(ScatteringStates), allocatable :: SingleSstates(:) !< Full s states (large system) for only one kx-ky-point

    complex(kdp) :: e !< energy at which the scattering states are calculated
    CHARACTER(LEN=1) :: Side !< the side of the lead
!optional : add kx,ky point at which scattering states are calculated
  end type FourierScatteringStates



end module mTypes
