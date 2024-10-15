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
!                   NEGFMOD  
! IN THIS FILE IS LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
      MODULE negfmod
 
! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************
     
      IMPLICIT NONE
      INTEGER     ikpmod,itermod,ivmod,bs_nmid1,bs_nmid2,nbss,ndivxy(2)
      integer     ndivxyNL(2,2) !first index is the number of leads, second index is 2 for kx and ky
      INTEGER     ikpmodK
      INTEGER, ALLOCATABLE ::    nebss_bs(:,:)
      double precision, ALLOCATABLE ::    deltabss_bs(:)
      DOUBLE PRECISION, ALLOCATABLE :: lkpointsmod(:,:)
      DOUBLE PRECISION :: kpointmod(3),lucmod(3),vmod,dosk,doskvv,wkmod,
     .    deltamin,setemax,setemin,critam,bs_tol,
     .    bs_min,deltaimag,deltaini,deltatode,SigmaWideBand
      double precision :: weightrho

      DOUBLE PRECISION :: vgate,zgatel,zgater,alphag
      double precision tolab
      integer em_nau,em_nuo,em_nas,em_nso,em_iscf,storesigma
      integer bssc,bseskip,bs_nmid,bs_method,
     .    bsrun,bsskip,maxdepth,ndivisions
      integer, pointer ::  em_isa(:),em_iaorb(:),
     .    em_iphorb(:)
      logical  dsigmade,bs_add,setene,cbspart,evec,
     .    bs_writetrc,writetk,negfon,emldos,emldos2,emSTTk,emSTT,
     .    emSTTLin
      logical SkipEquilibriumRho,SkipNonEquilibriumRho
      double precision emSTTkScale
      DOUBLE COMPLEX, ALLOCATABLE, DIMENSION (:,:,:) :: rhobs_delta,
     .    aux_par
      CHARACTER(LEN=15) :: integraltype
      CHARACTER(LEN=12) :: gridmethod
      LOGICAL :: pdosgs,skipright,emtimings,WriteNk, outinfo

      LOGICAL m_usehinv,m_usevinv,m_callsvd,m_complexbands,m_dosleads,
     .    leadsdos,leadspdos,empdos,emdos,empdosk,empdoskso
      DOUBLE PRECISION m_tolki,m_svdtolmax,m_svdtolmin,m_dsigmamax,
     .    m_rnoise,m_skipsvd,m_svdtolzi

      integer nenet,inversion_solver,nprocs_hs
      logical TransmOverk,transmembed,writeevsigma,skiptransm,
     &    deauto,trcef,geterrsigma,getsigma,different_eflr,emforces
      logical WriteIkTrcK
      DOUBLE PRECISION :: tolsigma,deltaimagtrc,trcde
      DOUBLE PRECISION :: TEnergI,TEnergF,maxkappa
      double precision, allocatable :: TrcAverage(:)
      double PRECISION    ef_em, ef_l,ef_r
      logical    gamma_negf,sigmatodisk,sigmafromdisk,overwritehs
      DOUBLE COMPLEX, Allocatable :: matbuf1(:,:,:),matbuf2(:,:,:),
     .    matbuf3(:,:,:)
      integer negf_base_comm
      DOUBLE PRECISION :: ThetaLeadsL,PhiLeadsL,ThetaLeadsR,PhiLeadsR
      logical timereversal,curr_dist,curr_distK, curr_distKEne,
     .    curr_distTRC
      logical TransmissionChannels,TransmissionChannelsWFS
      integer TransmissionChannelsWFSSkipEne
      integer TransmissionChannelsWFSSkipKP
      logical TransmissionMatrix,TransmissionMatrixWFS,
     .    TransmissionMatrixPDOS
      integer TransmissionMatrixPDOSNWrite
      integer TransmissionMatrixiSetPhase
      integer GetRhoSingleLead,MinChannelIndex,MaxChannelIndex
      logical TransmissionRL
      logical GetT0S
      integer nSplit
      integer NParallelK
      double precision, parameter :: RyToeV=13.605698066D0
      double precision, parameter :: BohrToAng=0.529177249D0
      integer onekp
      logical Add_Rho_Gate,debug_gate,
     .    AddVgCompensatingIsolatedLocalCharges
      double precision Q_gate,RGlxmin,RGlxmax,RGlymin,RGlymax,RGlzmin,
     .    RGlzmax
      logical Set_RhoBoundary_Leads,Set_RhoBoundaryOverlap_Leads
      logical Set_HBoundary_Leads,Set_HLR_Zero
      integer n_Replace_L,n_Replace_R
      logical TRCAddVCDFT
      integer ik_select

      integer N_LocalCharges
      double precision, allocatable :: x_LocalCharges(:,:)
      double precision, allocatable :: r_LocalCharges(:)
      double precision, allocatable :: q_LocalCharges(:)

      integer EM_NSPINBLocks
      logical em_Last_SCF_Step,em_endcode1
      double precision nlSetZeroRatio, nrSetZeroRatio
      double precision HSTol
      logical RhoSetZeroIfHZero
      integer NRunSigmaMax

!current distribution

      integer curr_nb1,curr_nb2
      double precision curr_fl_L, curr_fr_L, curr_fl_R, curr_fr_R
      integer em_nflux
      integer, allocatable :: em_nbfluxStart(:)

!impurity solver parameters
      logical PrintImpurityGfMatsubara, ComputeImpurityGfMatsubara, 
     .    CallImpuritySolver

****BEGIN electron-phonon coupling********
      logical calcmpsh
      integer nmode, nfine_Vb, rmoinit, rmolast, rorbinit, rorblast
      integer maxnh_eph, ihead_eph, itail_eph, na_fc_tot
      integer, allocatable :: naorbs(:)
      integer,  dimension(:), allocatable :: mode_index
      double precision, dimension(:), allocatable :: mode_freq
      double precision, dimension(:), allocatable :: mode_damp
      double precision  temp_ph
      character*8 LOE_Damp_Type
      character*25 fildhdq, filiets

      integer, dimension(:), allocatable :: 
     .              listh_eph, listhptr_eph, numh_eph, indxuo_eph
      double precision, dimension(:,:), allocatable :: xij_eph

      double precision, dimension(:,:,:), allocatable :: t_ec_Vb
      double precision, dimension(:,:,:), allocatable :: t_ecL_Vb
      double precision, dimension(:,:,:), allocatable :: t_ecR_Vb
      double precision, dimension(:,:,:), allocatable :: t_AsymL_Vb
      double precision, dimension(:,:,:), allocatable :: t_AsymR_Vb
      double precision, dimension(:,:,:), allocatable :: t_inel_Vb

      double precision, dimension(:,:,:), allocatable :: t_in_vib
      double precision, dimension(:,:,:), allocatable :: t_inL_vib
      double precision, dimension(:,:,:), allocatable :: t_inR_vib
      double precision, dimension(:,:,:), allocatable :: t_eh_vib
**** END  electron-phonon coupling********

      logical em_NonSCF

!      integer em_nau,em_nuo,em_isa(em_nau),em_iaorb(em_nuo),
!    .    em_iphorb(em_nuo)
      END MODULE negfmod

