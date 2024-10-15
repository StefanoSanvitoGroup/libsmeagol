!
! Copyright (c) Smeagol Authors:
! A. R. Rocha, V. Garcia-Suarez, S. Bailey, C. J. Lambert, J. Ferrer and
! S. Sanvito 2003-2005
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
! SMEAGOL IS DISTRIBUTED ONLY THROUGH THE OFICIAL WEBSITE (www.smeagol.tcd.ie)
! UPON COMPLETION OF THE "SMEAGOL ACADEMIC LICENSE" .
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
  SUBROUTINE TRANSM(ik,nk,N1,NL,NR,NSPIN,V,IV, NspinBlocks, &
             NspinComplexMatrix,slabel,kpoint,wk, &
             H0L,H0R,S0L,S0R,S1Li,S1Ri,H1Li,H1Ri,ef,T, &
             hgeneral,sgeneral,rhogeneral,ematgeneral, &
             istep, inicoor, idyn)

! *****************************************************************
! Calculates the transmission coefficients after a self consistent
! DFT calculation of a nanoscopic system coupled to
! charge reservoirs.
!
! Written by Alexandre Reily Rocha, June 2003
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: rochaa@tcd.ie
! ***************************** HISTORY ***********************************
! Original version: June 2003
! Altered by Ivan Rungger, October 2008
!   Added the output of different quantities (e.g. such as the PDOS),and
!   include the use of the new method for the calculation of the
!   self-energies. The format of the output has been changed.
!! -----------------------------------------------------------------
!! Add arguments idyn, inicoor, istep for MD
!! Meilin Bai, Dec 2012
!!
! ***************************** INPUT *****************************
! integer ik                    : k-point index
! integer nk                    : Number of k-points (time rev. symm.)
! integer N1                    : Dimension of the basis orbitals
! integer NL                    : Number of basis orbitals in the left
!                            lead, including spin components
! integer NR                    : Number of basis orbitals in the right
! integer NSPIN                 : Number of spin components
! integer IV                    : Voltage step
! complex*8 H0L(NL,NL,NSPIN)    : Hamiltonian left lead
! complex*8 H1L(NL,NL,NSPIN)    : Coupling Matrix left lead
! complex*8 S0L(NL,NL,NSPIN)    : on-site Overlap Matrix left lead
! complex*8 S1L(NL,NL,NSPIN)    : neighbour Overlap Matrix left lead
! complex*8 H0R(N2,N2,NSPIN)    : Hamiltonian right lead
! complex*8 H1R(N2,N2,NSPIN)    : Coupling M. right lead
! complex*8 S0R(N2,N2,NSPIN)    : on-site Overlap M. right lead
! complex*8 S1R(N2,N2,NSPIN)    : neighbour Overlap M. right lead
! complex*8 QL(NL,NL,NSPIN)     : transformation Matrix for decimation
! complex*8 QR(NR,NR,NSPIN)     : transformation Matrix for decimation
! CHARACTER side_rankL(NSPIN)   : Side to which transformation must be
!                                 applied (left lead)
! CHARACTER side_rankR(NSPIN)   : Side to which transformaiont must be
!                                 applied (right lead)
! integer  numberL(NSPIN)       : Number of decoupled states on the left
!                                 lead
! integer  numberR(NSPIN)       : Number of decoupled states on the right
!                                 lead
! character*20  slabel          : System label
! real*8   wk                   : Weight for k points
! *********************** SAVED VARIABLES ***************************
! integer  NeneT              : Number of energy points for Transmission
! double precision TEnergI        : Initial Energy for Transmission
! double precision TEnergF        : Final Energy for Transmission
! ************************ OUTPUT/SAVED ******************************
! TE(Nenerg,NSPIN)              : Transmission Coefficients
! ********************************************************************

      USE precision
      use sigma, only: sigma_method
      use negfmod
      use mTypes
      use mMatrixUtil
      use mONInterface
      use mSigmaMethod1
      use mMPI_NEGF
      use mEnergyGrid
      use mTransmissionDecomposition, only : Tmatrix,FindnrN,TransmissionChannelsDecomposition,ShiftRightSelfenergy,TransmissionSpinComponents_general
      use negfcoop ! COOP/COHP
      use mImpuritySolver,only: PrintGFHSMatsubaraGeneral,CalculateSinv
      use mCurrentDistribution, only : CurrentDistributionGeneral

      IMPLICIT NONE

#ifdef MPI
      INTEGER, DIMENSION (MPI_STATUS_SIZE) :: istatus
      INTEGER :: MPIerror
#endif

      INTEGER :: Nnodes
      INTEGER, intent(in) :: N1,NL,NR,NSPIN,IV,NspinComplexMatrix,NspinBlocks, idyn, istep, inicoor
      INTEGER, SAVE :: N1Half,NlHalf,NrHalf
      CHARACTER :: slabel*20, paste*35, filetr*35, filetrk*35,chivv*8
      CHARACTER :: fileempdos*100

      character pasbias*25


      DOUBLE PRECISION, ALLOCATABLE, SAVE :: TE(:,:),TEk(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: TERL(:,:),TEkRL(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE :: TELL(:,:),TEL(:,:)
      DOUBLE PRECISION, DIMENSION (2) :: kpoint
      INTEGER, ALLOCATABLE, SAVE :: NCHANk(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE, SAVE ::  dos(:,:),dosvv(:,:),dos2(:,:),NCHANL(:,:),  NCHANR(:,:),emdostot(:,:),empdostot(:,:,:)
      ! COOP/COPH
      double precision, allocatable, save :: coop1(:,:,:), coop2(:,:,:),cohp1(:,:,:), cohp2(:,:,:)
      integer ibond

      DOUBLE PRECISION, ALLOCATABLE ::  empdosbuf(:,:),empdostotkpt(:,:,:)
      double complex, allocatable, save:: Tr4G(:,:)
      double precision, allocatable, save :: tchannels(:,:,:),tchannelsK(:,:,:,:)
      integer, allocatable, save :: nchannelSLK(:,:,:),nchannelsRK(:,:,:)
      double precision, allocatable :: tchannelsEne(:)
      double precision emdostotk(NspinBlocks), empdostotk(em_nuo,NspinBlocks)
      double complex, allocatable :: Tr4GK(:)
      integer    nTr4G
      DOUBLE COMPLEX, PARAMETER :: zi=(0.D0,1.D0)
      DOUBLE COMPLEX, DIMENSION (NL,NL) :: S0L
      DOUBLE COMPLEX, DIMENSION (NL,NL) :: S1Li
      DOUBLE COMPLEX, DIMENSION (NL,NL,NSPIN) :: H0L,H1L,H1Li,QL,S1L
      DOUBLE COMPLEX, DIMENSION (NR,NR) :: S0R
      DOUBLE COMPLEX, DIMENSION (NR,NR) :: S1Ri
      DOUBLE COMPLEX, DIMENSION (NR,NR,NSPIN) :: H0R,H1R,H1Ri,QR,S1R

      DOUBLE COMPLEX, ALLOCATABLE :: GF1_aux(:,:), GF2_aux(:,:)
      DOUBLE COMPLEX, allocatable :: Gamma1_aux(:,:)
      DOUBLE COMPLEX, allocatable :: Gamma2_aux(:,:)
      DOUBLE COMPLEX, allocatable :: GammaGf1_RL(:,:)
      DOUBLE COMPLEX, allocatable :: GammaGf2_RL(:,:)
      DOUBLE COMPLEX, allocatable :: GammaGf1_LL(:,:)
      DOUBLE COMPLEX, allocatable :: GammaGf2_LL(:,:)
      DOUBLE COMPLEX, allocatable :: GammaGf3_LL(:,:)

      DOUBLE COMPLEX, allocatable :: GammaL(:,:)
      DOUBLE COMPLEX, allocatable :: GammaR(:,:)


      type(matrixTypeGeneral) :: gfgeneral,gfout

      type(matrixTypeGeneral),intent(inout) :: ematgeneral(NspinComplexMatrix)
      type(matrixTypeGeneral) :: hgeneral(NspinComplexMatrix),rhogeneral(NspinComplexMatrix),sgeneral
      integer  nnzrow(n1),nnz
      type(ioType) :: io
      integer gfmattype

      INTEGER :: nk,ik,iufiletr,nrchan,iukfiletr

      INTEGER, allocatable :: IPIV(:)
      INTEGER     indt,NM
      DOUBLE COMPLEX, ALLOCATABLE :: work(:),GF_iter_dag(:,:)
      DOUBLE COMPLEX,ALLOCATABLE,SAVE :: fHL(:,:),fHR(:,:),fH0L(:,:), fHLmult(:,:),fHRmult(:,:), GfLR(:,:), GfRL(:,:),GfLL(:,:),GfRR(:,:),fHmult2(:,:),fHmult3(:,:),   fH0R(:,:),fH0(:,:)
      DOUBLE COMPLEX, ALLOCATABLE :: GF_iter1l(:,:),GF_iter1r(:,:)
      DOUBLE COMPLEX, ALLOCATABLE :: GF_iterLR(:,:),GF_iterRL(:,:)
      DOUBLE COMPLEX, ALLOCATABLE :: gamma1op(:,:),gamma2op(:,:)
      DOUBLE COMPLEX, ALLOCATABLE, SAVE :: WORK3(:)
      DOUBLE COMPLEX, ALLOCATABLE, SAVE :: WORK2L(:), WORK2R(:)
      INTEGER,ALLOCATABLE,SAVE ::  IPIV2L(:),IPIV2R(:),IPIV3(:,:)
      DOUBLE COMPLEX,PARAMETER :: alpha=(1.D0,0D0)
      DOUBLE PRECISION :: V,ef,T,dosk2
      double complex ei

      INTEGER, DIMENSION (NSPIN) :: numberL,numberR

      DOUBLE COMPLEX :: ei0,ener_sigma
      double precision dsigma
      double complex weightc,signweight,weightcurrent,flrbuf
      double complex, allocatable :: kmat(:,:),GfLN(:,:,:),SigmaRN(:,:,:),SigmaLN(:,:,:)
      double complex, allocatable :: KC(:,:,:),KmC(:,:,:),KTip(:,:,:)
      integer nrN

      INTEGER :: I,II,JJ,ISPIN,L,J,k,MyNode,ind,opindex,nt,info,ispin2,i1,i2

      DOUBLE PRECISION :: wk,TEinterm, pi
      double PRECISION, save :: deltaene,deltaenes,deltaenebig
      double precision kbtoltrans,fl,fr,efl,efr
      double complex flc, frc, flbuf, frbuf

      external pasbias
      integer*4:: sc_0,sc_1,sc_r,sc_m
      integer*4:: sc_0b,sc_1b,sc_rb,sc_mb
      integer :: nleadslr=2
      integer, allocatable :: leadsdim(:)
      integer ntrcchannels
      logical outputwfs
      double precision, allocatable :: LeadsVoltageShift(:)
      integer nb
      integer ie_global
      logical writeGFHSheader

      logical frstme  !!
      save frstme     !!
      data frstme /.true./ !! Meilin Bai

!      CALL TIMER('TRANSI',1)
!      CALL TIMER('TRANA1',1)
#ifdef MPI
      CALL MPI_COMM_SIZE(negf_comm,Nnodes,MPIerror)
      CALL MPI_COMM_RANK(negf_comm,mynode,MPIerror)
#else
      MyNode=0
      Nnodes =1
#endif

      if(.not.negfon)then
        gfmattype=0 !for dense GF matrix
      else
        gfmattype=2 !for sparse GF matrix
      endif

      if(v.gt.1d-3)then
        leadsdos=.false.
      endif

      pi=3.141592653589793D0

      if(ComputeImpurityGfMatsubara.and.ik.eq.1)then
!if we use k-points this needs to be re-run for each kpoints, or else
!stored in memory for all k-points
        call CalculateSinv(sgeneral,gfmattype)
      endif


      IF (ik.eq.1) THEN

        if(mynode_inverse.eq.0)then
          allocate(leadsdim(nleadslr),LeadsVoltageShift(nleadslr))
          leadsdim(1)=nl
          leadsdim(2)=nr
          LeadsVoltageShift(1)=0.5D0 * V
          LeadsVoltageShift(2)=-0.5D0 * V

          call energygrid_transm(slabel,nspin,V,T,Ef, &
          nleadslr,leadsdim,nk, LeadsVoltageShift, &
          deltaene,deltaenes, deltaenebig)
          deallocate(leadsdim,LeadsVoltageShift)
        endif

        N1Half=N1/2
        NlHalf=NL/2
        NrHalf=NR/2

      ENDIF

      nTr4G=74 ! 4*2*4*2 + 2 + 4*2 for Gamma_i * G_j^dagger * Gamma_k * G_l + Gamma_i * G_i + Gamma_Lup * G_up^dagger * Gamma_Rup * G_down +Gamma_Ldown * G_down^dagger * Gamma_Rdown * G_down
      ntrcchannels=MaxChannelIndex-MinChannelIndex+1
      IF (mynode_inverse.eq.0.and.ik .EQ. 1) Then

        if(MyNode.EQ.0) THEN
!          write(*,*)"deltaimagtrc=",deltaimagtrc
          ALLOCATE(NCHANL(ETransmGrid%nEnergiesGlobal,NSPIN), NCHANR(ETransmGrid%nEnergiesGlobal,NSPIN))
          ALLOCATE(TE(ETransmGrid%nEnergiesGlobal,NSPIN))
          if(TransmissionRL)then
            ALLOCATE(TERL(ETransmGrid%nEnergiesGlobal,NSPIN))
            ALLOCATE(TELL(ETransmGrid%nEnergiesGlobal,NSPIN))
            ALLOCATE(TEL(ETransmGrid%nEnergiesGlobal,NSPIN))
          endif
          if(m_dosleads)then
            ALLOCATE(dos(ETransmGrid%nEnergiesGlobal,NSPIN))
            ALLOCATE(dosvv(ETransmGrid%nEnergiesGlobal,NSPIN))
          endif
          if(leadsdos)then
            ALLOCATE(dos2(ETransmGrid%nEnergiesGlobal,NSPIN))
          endif
          if(emdos)then
            ALLOCATE(emdostot(ETransmGrid%nEnergiesGlobal,NspinBlocks))
          endif
          if(emdos)then
            if(empdos)then
              ALLOCATE(empdostot(em_nuo,ETransmGrid%nEnergiesGlobal,NspinBlocks))
            endif
          endif
          ! COOP
          if (coopinfo%ccoop) then
              allocate(coop1(ETransmGrid%nEnergiesGlobal,nspin,coopinfo%nbond))
              allocate(coop2(ETransmGrid%nEnergiesGlobal,nspin,coopinfo%nbond))
              allocate(cohp1(ETransmGrid%nEnergiesGlobal,nspin,coopinfo%nbond))
              allocate(cohp2(ETransmGrid%nEnergiesGlobal,nspin,coopinfo%nbond))
          endif
          !
          if(TransmissionChannels)then
            ALLOCATE(tchannels(ntrcchannels,ETransmGrid%nEnergiesGlobal,NSPIN))
            if(TransmOverK)then
              ALLOCATE(tchannelsK(ntrcchannels,ETransmGrid%nEnergiesGlobal,NSPIN,nk))
            endif
          endif
          if(GetT0S)then
            allocate(Tr4G(nTr4G,ETransmGrid%nEnergiesGlobal))
          endif
          if (TransmOverK) then
            ALLOCATE(TEk(ETransmGrid%nEnergiesGlobal,NSPIN,nk))
            ALLOCATE(nchannelsLK(ETransmGrid%nEnergiesGlobal,NSPIN,nk))
            ALLOCATE(nchannelsRK(ETransmGrid%nEnergiesGlobal,NSPIN,nk))
            TEk = 0.d0
            nchannelsLK = 0
            nchannelsRK = 0
          endif
        else
          ALLOCATE(TE(ETransmGrid%nEnergies,NSPIN))
          if(TransmissionRL)then
            ALLOCATE(TERL(ETransmGrid%nEnergies,NSPIN))
            ALLOCATE(TELL(ETransmGrid%nEnergies,NSPIN))
            ALLOCATE(TEL(ETransmGrid%nEnergies,NSPIN))
          endif
          ALLOCATE(NCHANL(ETransmGrid%nEnergies,NSPIN),  NCHANR(ETransmGrid%nEnergies,NSPIN))
          if(m_dosleads)then
            ALLOCATE(dos(ETransmGrid%nEnergies,NSPIN))
            ALLOCATE(dosvv(ETransmGrid%nEnergies,NSPIN))
          endif
          if(leadsdos)then
            ALLOCATE(dos2(ETransmGrid%nEnergies,NSPIN))
          endif
          if(emdos)then
            ALLOCATE(emdostot(ETransmGrid%nEnergies,NspinBlocks))
            if(empdos)then
              ALLOCATE(empdostot(em_nuo,ETransmGrid%nEnergies,NspinBlocks))
            endif
          endif

          ! COOP
          if (coopinfo%ccoop) then
              allocate(coop1(ETransmGrid%nEnergies,nspin,coopinfo%nbond))
              allocate(coop2(ETransmGrid%nEnergies,nspin,coopinfo%nbond))
              allocate(cohp1(ETransmGrid%nEnergies,nspin,coopinfo%nbond))
              allocate(cohp2(ETransmGrid%nEnergies,nspin,coopinfo%nbond))
          endif
          !
          if(TransmissionChannels)then
            ALLOCATE(tchannels(ntrcchannels,ETransmGrid%nEnergies,NSPIN))
            if(TransmOverK)then
              ALLOCATE(tchannelsK(ntrcchannels,ETransmGrid%nEnergies,NSPIN,nk))
            endif
          endif
          if(GetT0S)then
            allocate(Tr4G(nTr4G,ETransmGrid%nEnergies))
          endif
          if (TransmOverK) then
            ALLOCATE(TEk(ETransmGrid%nEnergies,NSPIN,nk))
            ALLOCATE(nchannelsLK(ETransmGrid%nEnergies,NSPIN,nk))
            ALLOCATE(nchannelsRK(ETransmGrid%nEnergies,NSPIN,nk))
            TEk = 0.d0
            nchannelsLK = 0
            nchannelsRK = 0
          endif
        endif ! Mynode =0

        TE=0.D0
        if(TransmissionRL) then
          TERL=0.D0
          TELL=0.D0
          TEL=0.D0
        endif
        NCHANL = 0D0
        NCHANR = 0D0

        if(m_dosleads)then
          dos = 0D0
          dosvv = 0D0
        endif
        if(leadsdos)then
          dos2 = 0D0
        endif
        if(emdos)then
          emdostot=0D0
          if(empdos)then
            empdostot=0D0
          endif
        endif
        ! COOP
        if (coopinfo%ccoop) then
            coop1 =0.d0
            coop2 =0.d0
            cohp1 =0.d0
            cohp2 =0.d0
        endif
        !
        if(TransmissionChannels)then
          tchannels=0D0
          if(TransmOverK)then
            tchannelsK=0D0
          endif
        endif
        if(GetT0S)then
          Tr4G=0.0D0
        endif
      ENDIF ! mynode_inverse .eq.0 .and. ik .eq.1

      IF (mynode_inverse.eq.0) Then
        if(MyNode.EQ.0) THEN
          if(empdosk) ALLOCATE(empdostotkpt(em_nuo,ETransmGrid%nEnergiesGlobal,NspinBlocks))
        else
          if(empdosk) ALLOCATE(empdostotkpt(em_nuo,ETransmGrid%nEnergies,NspinBlocks))
        endif
        if(empdosk)empdostotkpt=0.0D0
      endif

!      CALL TIMER('TRANA1',2)
!      CALL TIMER('TRANA2',1)

      if((MyNode.eq.0).and.getsigma.and.(ik .EQ. 1)) &
      call writesigmainfo(slabel,n1,nl,nr, Nnodes,&
      ETransmGrid%nEnergies,NSPIN,nk,ETransmGrid%e(1), &
      tenergi,deltaenebig,deltaenes, v,ef)

      if(mynode_inverse.eq.0.and..not.(negfon))then
        allocate(work(N1**2),GF_iter_dag(N1,N1))
      endif

      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
      endif

      weightc=1.D0/(1.0D0 * ETransmGrid%nEnergiesGlobal)


      eneloop: DO I=1,ETransmGrid%nEnergies

        call fermi_distribution_real(dreal(ETransmGrid%e(I)-ef-v*0.5D0),T,fl)
        call fermi_distribution_real(dreal(ETransmGrid%e(I)-ef+v*0.5D0),T,fr)
!        write(12347,*)"fl,fr=",dreal((ETransmGrid%e(I)-ef)*13.6056981D0),fl,fr
        
        spinloop: DO ISPIN=1,ETransmGrid%nSpin

          if(mynode_inverse.eq.0)then

            if(emtimings)then
              CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')  'tbeforebarrier',(sc_1-sc_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
            endif
#ifdef MPI
            CALL MPI_BARRIER(inverseheads_comm, MPIerror)
#endif
            if(emtimings)then
              CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')  'tafterbarrier',(sc_1-sc_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
            endif

            allocate(LeadsVoltageShift(nleadslr))
            LeadsVoltageShift(1)=0.5D0 * V
            LeadsVoltageShift(2)=-0.5D0 * V

            
            if(m_complexbands)then
              write(12346,*)"start cb output for ispin,ik=",ispin,ikpmod
            endif
            call get_selfenergies(i,ispin,ik,ETransmGrid,v,0,LeadsVoltageShift,sigmatodisk)
            if(m_complexbands)then
              write(12346,*)"end cb output for ispin,ik=",ispin,ikpmod
            endif

            deallocate(LeadsVoltageShift)

            if(emtimings)then
              CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')  'senetrc',(sc_1-sc_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
            endif
!              CALL TIMER('TRANc2',1)


!            CALL TIMER('TRANc1',1)

              NCHANL(I,ISPIN) = NCHANL(I,ISPIN) + wk * ETransmGrid%sigma(1,1,1,1)%nchannels
              if(m_dosleads)then
                dos(I,ISPIN) = dos(I,ISPIN) + dosk * wk 
                dosvv(I,ISPIN) = dosvv(I,ISPIN) + doskvv * wk 
              endif
              NCHANR(I,ISPIN) = NCHANR(I,ISPIN) + wk *  ETransmGrid%sigma(2,1,1,1)%nchannels
              if(TransmOverK)then
                nchannelsLK(i,ispin,ik)=ETransmGrid%sigma(1,1,1,1)%nchannels
                nchannelsRK(i,ispin,ik)=ETransmGrid%sigma(2,1,1,1)%nchannels
              endif

            if (geterrsigma)then
              write(12347,*)"dsigmainfot=",ETransmGrid%sigma(1,1,1,1)%e,ispin,i,ikpmod
              call check_error_sigma(dsigma,.true.,ETransmGrid%sigma(1,1,1,1)%e,'L',nl, &
                h0L(:,:,ispin),H1Li(:,:,ispin),S0L,S1Li,  ETransmGrid%sigma(1,1,1,1)%sigma)
              call check_error_sigma(dsigma,.true.,ETransmGrid%sigma(2,1,1,1)%e,'R',nr, &
               h0r(:,:,ispin),H1ri(:,:,ispin),S0r,S1ri,   ETransmGrid%sigma(2,1,1,1)%sigma)
            endif

            if(emtimings)then
              CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')  'senechecksigma',(sc_1-sc_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
            endif
! begin get leads dos
          if(leadsdos)then

            call dostrc(NL,ETransmGrid%sigma(1,1,1,1)%sigma, ETransmGrid%sigma(2,1,1,1)%sigma,H0L(:,:,ISPIN), &
              H1Li(:,:,ISPIN),S0L,S1Li ,ETransmGrid%sigma(1,1,1,1)%e, ispin,dosk2,ef,leadspdos)

            dos2(i,ispin)=dos2(i,ispin)+dosk2 * wk

          endif
! end get leads dos

! begin get self-energy eigenvalues
          if(writeevsigma)then
            ei0=ETransmGrid%e(I)-v*0.5D0
            call evsigma(nl,ETransmGrid%sigma(1,1,1,1)%sigma,nr, ETransmGrid%sigma(2,1,1,1)%sigma,ispin,ik,ei0)
          endif
! end get self-energy eigenvalues

            if(skiptransm)then
              if(TransmissionMatrix)then

                if(.not.allocated(GfRL))allocate(GfRL(NR,NL))
                if(.not.allocated(GfLL))allocate(GfLL(NL,NL))
                if(.not.allocated(GfLR))allocate(GfLR(NL,NR))
                if(.not.allocated(GfRR))allocate(GfRR(NR,NR))
                GfRL=0.0D0
                GfLL=0.0D0
                GfLR=0.0D0
                GfRR=0.0D0

                call Tmatrix(ETransmGrid%e(i),ef,nl,nr,GfRL,GfLL,GfLR,GfRR,ETransmGrid%sigma(1,1,1,1)%sigma,ETransmGrid%sigma(2,1,1,1)%sigma,ispin,ETransmGrid%nspin)

                deallocate(GfRL)
                deallocate(GfLL)
                deallocate(GfRR)
                deallocate(GfLR)
              endif

              call deallocate_selfenergies(i,ispin,ik,ETransmGrid)

              if(ISPIN.lt.NSPIN)then
                cycle spinloop
              endif
              if(ISPIN.eq.NSPIN)then
                cycle eneloop
              endif
            endif

!          CALL TIMER('TRANc1',2)

            if(emtimings)then
              CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')  'dosltrc',(sc_1-sc_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
            endif
!              CALL TIMER('TRANc2',1)
            if(gfmattype.eq.0)then
              nnz=n1*n1
            elseif(gfmattype.eq.2)then
              call findnnzgf2(nnz,nnzrow,n1,n1,nl,nr, ETransmGrid%sigma(1,1,1,1)%sigma, ETransmGrid%sigma(2,1,1,1)%sigma,hgeneral(ispin))
            endif

            if(emtimings)then
              CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')'findnnzt',(sc_1-sc_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
            endif

            call AllocateMatrixGeneral(n1,n1,nnz,gfmattype,gfgeneral, "keldyshreal", io)

            ei0=ETransmGrid%e(i)+zi*deltaimagtrc
            call setgfelementsgeneral_nc(ei0,NspinComplexMatrix,ispin, &
            gfgeneral,nnz,n1,nl,nr,  ETransmGrid%sigma(1,1,1,1)%sigma, &
             ETransmGrid%sigma(2,1,1,1)%sigma, hgeneral,sgeneral)
!end allocating sparse matrix

            if(emtimings)then
              CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')'setgft',(sc_1-sc_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
            endif



            if(mynode_inverse.eq.0)then
!              CALL TIMER('gfinvt',1)
              If(.not.negfon)then

                if(GetT0S)then

                  if(nsplit.gt.0)then
                    nb=nsplit-1
                    allocate(kmat(n1-nb,nb))
                    kmat=gfgeneral%matdense%a(nsplit:n1,1:nsplit-1)
                    call FindnrN(kmat,n1-nb,nb,nrN)
                    deallocate(kmat)
                  else
                    nb=n1
                    nrN=nr
                  endif
                  if (outinfo) write(12347,*)"nrN,nb=",nrN,nb

                  if(.not.allocated(KC))then
                    allocate(KC(nrN,n1-nb,nspin))
                    allocate(KmC(n1-nb,nrN,nspin))
                    allocate(KTip(n1-nb,n1-nb,nspin))
                  endif
                  KC(:,:,ispin)=gfgeneral%matdense%a(nb-nrN+1:nb,nb+1:n1)
                  KmC(:,:,ispin)=gfgeneral%matdense%a(nb+1:n1,nb-nrN+1:nb)
                  KTip(:,:,ispin)=gfgeneral%matdense%a(nb+1:n1,nb+1:n1)

                endif

                allocate(IPIV(n1))
                CALL ZGETRF(N1,N1,gfgeneral%matdense%a,N1,IPIV,   INFO)
                CALL ZGETRI(N1,gfgeneral%matdense%a,  N1,IPIV,WORK,N1**2,INFO)
                deallocate(ipiv)

                DO II=1,N1
                  DO JJ=1,N1
                    GF_iter_dag(II,JJ)=  DCONJG(gfgeneral%matdense%a(JJ,II))
                  ENDDO
                ENDDO

              else

                if((emdos).and.(GetRhoSingleLead.ne.0.or.TransmissionMatrix.or.TransmissionRL.or.TransmissionChannels))then
                  opindex=3
                elseif(emdos.and.GetRhoSingleLead.eq.0)then
                  opindex=5
                elseif(TransmissionMatrix.or.TransmissionRL.or.TransmissionChannels)then
                  opindex=2
                else
                  opindex=4
                endif

                if((opindex==2).or.(opindex==3))then
                  call AllocateMatrixGeneral(n1,nl+nr,n1*(nl+nr),0,gfout,"transm", io)
                else
                  call AllocateMatrixGeneral(nr,nl,nl*nr,  0,gfout,"transm", io)
                endif

!                if(mynode==0)then
!                  call PrintMatrixCRS(gfgeneral%matSparse,"bf",io)
!                endif
                call InvertONGeneral(N1,gfgeneral,nl,nr,gfout,opindex,  inversion_solver)
!                if(mynode==0)then
!                  call PrintMatrixCRS(gfgeneral%matSparse,"gf",io)
!                endif
!#ifdef MPI
!                CALL MPI_BARRIER(inverseheads_comm, MPIerror)
!#endif
!                call stopnegf


                if(.not.allocated(GfRL))allocate(GfRL(NR,NL))
                if(TransmissionRL.or.TransmissionMatrix)then
                  if(.not.allocated(GfLR))allocate(GfLR(NL,NR))
                  if(.not.allocated(GfLL))allocate(GfLL(NL,NL))
                endif

                if((opindex==2).or.(opindex==3))then
                  if((emldos2.or.emdos).and.(GetRhoSingleLead.ne.0))then
                    call AllocateMatrix(n1,nl+nr,gfgeneral%matdense,"transm",io)
                    gfgeneral%matdense%a=gfout%matdense%a
                  endif
                  if(TransmissionMatrix.or.TransmissionRL)then
                    GfLR=gfout%matdense%a(1:nl,nl+1:nl+nr)
                    GfLL=gfout%matdense%a(1:nl,1:nl)

                    if(TransmissionMatrix)then
                      allocate(gfRR(NR,NR))
                      GfRR=gfout%matdense%a(n1-nr+1:n1,nl+1:nl+nr)
                    endif
  
                  endif
                  GfRL=gfout%matdense%a(n1-nr+1:n1,1:nl)

                  if(TransmissionChannels)then
                    allocate(GF_iter1l(n1,nl))
                    allocate(GF_iter1r(n1,nr))
                    allocate(GF_iterLR(nr,nl))
                    allocate(GF_iterRL(nl,nr))
                    GF_iter1l(:,1:nl)=gfout%matdense%a(:,1:nl)
                    GF_iter1r(:,1:nr)=gfout%matdense%a(:,nl+1:nl+nr)
                    GF_iterLR(1:nr,1:nl)=gfout%matdense%a(n1-nr+1:n1,1:nl)
                    GF_iterRL(1:nl,1:nr)=gfout%matdense%a(1:nl,nl+1:nl+nr)
                  endif

                else
                  GfRL=gfout%matdense%a(1:nr,1:nl)
                endif

                call DestroyMatrixGeneral(gfout,"transm",io)

              EndIf
            else
                if (outinfo) write(12347,*)"not calculating T"
            endif

            if(ComputeImpurityGfMatsubara)then
              if(i==1.and.mynode_inverse==0)then
                writeGFHSheader=.true.
!                write(12347,*)"setting writeheader T"
              else
                writeGFHSheader=.false.
!                write(12347,*)"setting writeheader F"
              endif

              call PrintGFHSMatsubaraGeneral(T,ei0,ef,gfgeneral,hgeneral(ispin),sgeneral, ETransmGrid%sigma(1,1,1,1)%sigma,nl, ETransmGrid%sigma(2,1,1,1)%sigma,nr,writeGFHSheader, myhead,PrintImpurityGfMatsubara,CallImpuritySolver,ETransmGrid%ig(i),ETransmGrid%nEnergies,ispin,ETransmGrid%nSpin)
            endif

!              CALL TIMER('gfinvt',2)

            if(emtimings)then
              CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
              write(12347,'(A,f12.6)') 'invertt',(sc_1-sc_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
            endif

            if(mynode_inverse.eq.0)then
              if(getsigma)then
!                write(*,*)"writesigma, node=",mynode
                call writesigma(slabel,ispin,mynode,ik,i,n1,nl,nr, ETransmGrid%e(i),ef,   ETransmGrid%sigma(1,1,1,1)%sigma,ETransmGrid%sigma(2,1,1,1)%sigma, gfgeneral%matdense%a)
!              else
!                write(*,*)"writesigma, not ouputting, node=",mynode
              endif


!              CALL TIMER('TRANc2',2)
!              CALL TIMER('TRANc3',1)
! begin get em dos
              if(TransmissionMatrix)then
  
                if(.not.negfon)then
                  allocate(GfRL(NR,NL))
                  GfRL=gfgeneral%matdense%a(n1-nr+1:n1,1:nl)
                  allocate(gfLL(NL,NL))
                  gfLL=gfgeneral%matdense%a(1:nl,1:nl)

                  allocate(GfLR(NL,NR))
                  GfLR=gfgeneral%matdense%a(1:nl,n1-nr+1:n1)
                  allocate(gfRR(NR,NR))
                  gfRR=gfgeneral%matdense%a(n1-nr+1:n1,n1-nr+1:n1)

                endif

                call Tmatrix(ETransmGrid%e(i),ef,nl,nr,GfRL,gfLL,gfLR,gfRR,ETransmGrid%sigma(1,1,1,1)%sigma,ETransmGrid%sigma(2,1,1,1)%sigma,ispin,ETransmGrid%nspin)
                if(.not.negfon)then
                  deallocate(GfRL)
                  deallocate(GfLL)
                  deallocate(GfLR)
                endif
                deallocate(GfRR)
              endif
              if(emdos)then

                if(mynode_inverse.eq.0.and..not.emldos2.and.GetRhoSingleLead.ne.0)then
!                  write(12347,*)"setting rho and omega to zero"
                  call SetRhoToZero(NspinComplexMatrix,rhogeneral)
                  call SetRhoToZero(NspinComplexMatrix,ematgeneral)
                endif

                if(mynode_inverse.eq.0.and.(emldos2.or.GetRhoSingleLead.ne.0))then

                  if(emtimings)then
                    CALL SYSTEM_CLOCK(sc_0b,sc_rb,sc_mb)
                  endif

                  if(GetRhoSingleLead==0)then

                    if(abs(curr_fl_R) > 1.000000001D0)then
                      flbuf=fl
                      if(curr_fl_R<0.D0)flbuf=-flbuf
                    else
                      flbuf=curr_fl_R
                    endif
                    if(abs(curr_fr_R) > 1.000000001D0)then
                      frbuf=fr
                      if(curr_fr_R<0.0D0)frbuf=-frbuf
                    else
                      frbuf=curr_fr_R
                    endif

                    frc=flbuf+frbuf
                    flc=0.0D0

!                    write(12347,*)"gr0:fl,fr=",dreal(ETransmGrid%e(I)),flc,frc
                    call updaterho_nc(rhogeneral,ematgeneral,emforces,&
                    ispin,NspinComplexMatrix,gfgeneral,nl,nr,gfmattype,&
                    -ETransmGrid%w(i),flc,frc,1.0D0,ei0,.true.)
                  else
                    if(abs(GetRhoSingleLead)==1.or.abs(GetRhoSingleLead)==3)then

                      if(abs(curr_fl_L) > 1.000000001D0)then
                        flbuf=fl
                        if(curr_fl_L<0.D0)flbuf=-flbuf
                      else
                        flbuf=curr_fl_L
                      endif
                      if(abs(curr_fr_L) > 1.000000001D0)then
                        frbuf=fr
                        if(curr_fr_L<0.0D0)frbuf=-frbuf
                      else
                        frbuf=curr_fr_L
                      endif
                      signweight= ((1.D0,0D0)/(2.D0*PI))* ETransmGrid%w(I) * (flbuf+frbuf)
!xxx: check units of linear response STT, and check why we set signweight to 1 if not emLDOS2
!xxx: check units when doing no emldos2 and no curr_dist
                    
!                      write(12347,*)"signweight_L=",signweight,flbuf,frbuf,flbuf+frbuf
!                      write(12347,*)"signweight_Lb=",dreal(ETransmGrid%e(i))-ef,dreal(flbuf+frbuf)

                      if(NspinComplexMatrix<=2)then
                        call RhoSingleLead(ETransmGrid%e(i),gfgeneral%matdense%a(1:n1,1:nl),n1,nl,nr,nl,ETransmGrid%sigma(1,1,1,1)%sigma,rhogeneral,ematgeneral,NspinComplexMatrix,signweight,ispin)
                      else
                        allocate(GF_iter1l(n1,nl))
                        GF_iter1l(:,1:NlHalf)=gfgeneral%matdense%a(1:n1,1:NlHalf)
                        GF_iter1l(:,NlHalf+1:nl)=gfgeneral%matdense%a(1:n1,N1Half+1:N1Half+NlHalf)
                        call RhoSingleLead(ETransmGrid%e(i),gf_iter1l,n1,nl,nr,nl,ETransmGrid%sigma(1,1,1,1)%sigma,rhogeneral,ematgeneral,NspinComplexMatrix,signweight,ispin)
                        deallocate(gf_iter1l)
                      endif

                    endif
                    if(abs(GetRhoSingleLead)==2.or.abs(GetRhoSingleLead)==3)then

                      if(abs(curr_fl_R) > 1.000000001D0)then
                        flbuf=fl
                        if(curr_fl_R<0.D0)flbuf=-flbuf
                      else
                        flbuf=curr_fl_R
                      endif
                      if(abs(curr_fr_R) > 1.000000001D0)then
                        frbuf=fr
                        if(curr_fr_R<0.0D0)frbuf=-frbuf
                      else
                        frbuf=curr_fr_R * GetRhoSingleLead/abs(GetRhoSingleLead)
                      endif
                      signweight= ((1.D0,0D0)/(2.D0*PI))* ETransmGrid%w(I) * (flbuf+frbuf)

!                      write(12347,*)"signweight_R=",signweight,flbuf,frbuf,flbuf+frbuf
                      if(negfon)then
                        call RhoSingleLead(ETransmGrid%e(i),gfgeneral%matdense%a(1:n1,nl+1:nl+nr),n1,nl,nr,nr,ETransmGrid%sigma(2,1,1,1)%sigma,rhogeneral,ematgeneral,NspinComplexMatrix,signweight,ispin)
                      else
                        allocate(GF_iter1r(n1,nr))
                        if(NspinComplexMatrix<=2)then
                          GF_iter1r(:,1:nr)=gfgeneral%matdense%a(1:n1,n1-nr+1:n1)
                        else
                          GF_iter1r(:,1:NrHalf)=gfgeneral%matdense%a(1:n1,N1Half-NrHalf+1:N1Half)
                          GF_iter1r(:,NrHalf+1:nr)=gfgeneral%matdense%a(1:n1,n1-NrHalf+1:n1)
                        endif
  
                        call RhoSingleLead(ETransmGrid%e(i),gf_iter1r,n1,nl,nr,nr,ETransmGrid%sigma(2,1,1,1)%sigma,rhogeneral,ematgeneral,NspinComplexMatrix,signweight,ispin)
                        deallocate(gf_iter1r)

                      endif
                    endif

                  endif

                  if(emtimings)then
                    CALL SYSTEM_CLOCK(sc_1b,sc_rb,sc_mb)
                    write(12347,'(A,f12.6)') 'RhoSingleLead',(sc_1b-sc_0b)*1.0d0/sc_rb
                  endif


                endif

                if(GetRhoSingleLead==0)then
                  call em_dos_general(n1,em_nuo,NspinBlocks,NspinComplexMatrix,gfgeneral,empdos,emdostotk,empdostotk,sgeneral)
                else
                  if(mynode_inverse.eq.0.and..not.emldos2)then
                    call em_dos_SingleLead_general(n1,em_nuo,NspinBlocks,NspinComplexMatrix,empdos,emdostotk,empdostotk,sgeneral,rhogeneral)
                  endif
                endif

                if(NspinBlocks<=2)then
                  emdostot(i,ispin)=emdostot(i,ispin)+emdostotk(1) * wk
                  if(empdos)then
                    empdostot(:,i,ispin)=empdostot(:,i,ispin)+   empdostotk(:,1) * wk
                  endif
                  if(empdosk) empdostotkpt(:,i,ispin)= empdostotk(:,1)
                else
                  do i1=1,NspinBlocks
                    emdostot(i,i1)=emdostot(i,i1)+emdostotk(i1) * wk
                    if(empdos)then
                      empdostot(:,i,i1)=empdostot(:,i,i1)+   empdostotk(:,i1) * wk
                    endif
                    if(empdosk) empdostotkpt(:,i,i1)= empdostotk(:,i1) 
                  enddo
                endif

              endif
               
              if(emtimings)then
                CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
                write(12347,'(A,f12.6)')  'emdost',(sc_1-sc_0)*1.0d0/sc_r
                CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
              endif


              if(TransmissionChannels)then
!                 here we do Trace[GammaL Gdagger GammaR G]
                allocate(gamma1op(nl,nl))
                allocate(gamma2op(nr,nr))
                if(.not.negfon)then
                  allocate(GF_iter1l(n1,nl))
                  allocate(GF_iter1r(n1,nr))
                  allocate(GF_iterLR(nr,nl))
                  allocate(GF_iterRL(nl,nr))
                  if (NSpinBlocks .le. 2)then
                      GF_iter1l(:,1:nl)=gfgeneral%matdense%a(1:n1,1:nl)
                      GF_iter1r(:,1:nr)=gfgeneral%matdense%a(1:n1,n1-nr+1:n1)
                      GF_iterLR(1:nr,1:nl)=gfgeneral%matdense%a(n1-nr+1:n1,1:nl)
                      GF_iterRL(1:nl,1:nr)=gfgeneral%matdense%a(1:nl,n1-nr+1:n1)
                  else
                      do j=1,2
                          GF_iter1l(:,(j-1)*NLhalf+1:j*NLHalf)=gfgeneral%matdense%a(:,(j-1)*N1Half+1:(j-1)*N1Half+NLHalf)
                          GF_iter1r(:,(j-1)*NRhalf+1:j*NRHalf)=gfgeneral%matdense%a(:,j*N1Half-NRHalf+1:j*N1Half)
                          do k=1,2
                             GF_iterLR((j-1)*NRHalf+1:j*NRHalf,(k-1)*NLhalf+1:k*NLHalf)=gfgeneral%matdense%a(j*N1Half-NRHalf+1:j*N1Half,(k-1)*N1Half+1:(k-1)*N1Half+NLHalf)
                             GF_iterRL((j-1)*NLHalf+1:j*NLHalf,(k-1)*NRhalf+1:k*NRHalf)=gfgeneral%matdense%a((j-1)*N1Half+1:(j-1)*N1Half+NLHalf,k*N1Half-NRHalf+1:k*N1Half)
                          enddo
                      enddo
                  endif
                endif

                allocate(tchannelsEne(ntrcchannels))

                ie_global=myhead+1+(i-1)*nheads
                if(TransmissionChannelsWFS.and.mod(ie_global-1,TransmissionChannelsWFSSkipEne)==0.and.mod(ik-1,TransmissionChannelsWFSSkipKP)==0)then
!                  write(12347,*)"true:",ik,i,ie_global,dreal(ei0),ispin
                  outputwfs=.true.
                else
!                  write(12347,*)"false:",ik,i,ie_global,dreal(ei0),ispin
                  outputwfs=.false.
                endif
                call TransmissionChannelsDecomposition(GF_iter1l,GF_iter1r,GF_iterLR,GF_iterRL,n1,nl,nr,ETransmGrid%sigma(1,1,1,1)%sigma,ETransmGrid%sigma(2,1,1,1)%sigma,nspin,ispin,dreal(ei0),ef,gamma1op,gamma2op,MinChannelIndex,MaxChannelIndex,ntrcchannels,tchannelsEne,ik,ie_global,nk,iv,slabel,kpoint,outputwfs)

                tchannels(:,i,ispin)=tchannels(:,i,ispin)+ tchannelsEne(:) * wk
                if(TransmOverK) tchannelsK(:,i,ispin,ik)=tchannelsK(:,i,ispin,ik)+ tchannelsEne(:) 

                deallocate(tchannelsEne)

                deallocate(gf_iter1l)
                deallocate(gf_iter1r)
                deallocate(GF_iterLR)
                deallocate(GF_iterRL)
                deallocate(gamma1op)
                deallocate(gamma2op)

                if(emtimings)then
                  CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
                  write(12347,'(A,f12.6)')  'tchannels',(sc_1-sc_0)*1.0d0/sc_r
                  CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
                endif
  

              endif


              if(GetT0S)then


                if(emtimings)then
                  CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
                  write(12347,'(A,f12.6)')  't_findnrN',(sc_1-sc_0)*1.0d0/sc_r
                  CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
                endif
  


                if(ispin==1)then
                  allocate(SigmaRN(nrN,nrN,nspin))
                  allocate(SigmaLN(nl,nl,nspin))
                  allocate(GfLN(nb,nb,nspin))
                endif

                if(nsplit.gt.0)then
                  GfLN(:,:,ispin)=gfgeneral%matdense%a(1:nb,1:nb)
                  call ShiftRightSelfenergy(n1,nr,nl,nsplit,KC,KmC,KTip,SigmaRN,nrN,work,ispin)
                else
                   GfLN(:,:,ispin)=gfgeneral%matdense%a
                   SigmaRN(:,:,ispin)=ETransmGrid%sigma(2,1,1,1)%sigma
                endif
                SigmaLN(:,:,ispin)=ETransmGrid%sigma(1,1,1,1)%sigma

                if(emtimings)then
                  CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
                  write(12347,'(A,f12.6)')  't_ShiftSENER',(sc_1-sc_0)*1.0d0/sc_r
                  CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
                endif
  


                if(ispin==nspin) then
                  allocate(Tr4GK(nTr4G))
                  call TransmissionSpinComponents_general(GfLN,nb,nl,nrN,SigmaLN,SigmaRN,nspin,Tr4Gk,nTr4G,dreal(ei0))
                  Tr4G(:,i)=Tr4G(:,i)+ Tr4Gk(:) * wk

                  deallocate(GfLN)
                  deallocate(SigmaRN)
                  deallocate(SigmaLN)
                  deallocate(Tr4Gk)
                  deallocate(KC)
                  deallocate(KmC)
                  deallocate(KTip)
                endif


!                 here we do Trace[GammaR Gdagger GammaL G]
!                 result should be identical to : Trace[GammaR Gdagger GammaL G]

!                call TransmissionSpinComponents(gfgeneral%matdense%a(1:nl,n1-nr+1:n1),nr,nl,ETransmGrid%sigma(2,1,1,1)%sigma,ETransmGrid%sigma(1,1,1,1)%sigma,nspin,ispin,TudTotK,dreal(ei0))

!                 here we do Trace[GammaL Gdagger GammaR G]
!                call TransmissionSpinComponents(gfgeneral%matdense%a(n1-nr+1:n1,1:nl),nl,nr,ETransmGrid%sigma(1,1,1,1)%sigma,ETransmGrid%sigma(2,1,1,1)%sigma,nspin,ispin,TudTotK,dreal(ei0))

!               cheeck: can we really just perform the k-point average like this?

                if(emtimings)then
                  CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
                  write(12347,'(A,f12.6)')  't_T0S',(sc_1-sc_0)*1.0d0/sc_r
                  CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
                endif
 

              endif
! end get em dos
!              CALL TIMER('TRANc3',2)
!              CALL TIMER('TRANc4',1)
      
            if(curr_distKEne.and..not.negfon.and.mynode_inverse.eq.0)then

              call CurrentDistributionGeneral(gfgeneral,n1,nl,nr,ETransmGrid%sigma(1,1,1,1)%sigma,ETransmGrid%sigma(2,1,1,1)%sigma,hgeneral,sgeneral,NspinBlocks,NspinComplexMatrix,weightc,ispin,ei0,ef,v,T,ik,wk)

              if(emtimings)then
                CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
                write(12347,'(A,f12.6)')  't_currdist',(sc_1-sc_0)*1.0d0/sc_r
                CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
              endif

            endif


            allocate(gammal(nl,nl))
            gammal=zi*(ETransmGrid%sigma(1,1,1,1)%sigma-DCONJG(TRANSPOSE(ETransmGrid%sigma(1,1,1,1)%sigma)))
            allocate(gammar(nr,nr))
            gammar=zi*(ETransmGrid%sigma(2,1,1,1)%sigma-DCONJG(TRANSPOSE(ETransmGrid%sigma(2,1,1,1)%sigma)))

            allocate(Gamma1_aux(nl,nr))
            allocate(Gamma2_aux(nr,nl))
!coop
            if (coopinfo%ccoop) then
               if (negfon .and. myhead .eq. 0) then
                  write(6,'(a)') 'Warning: EM.OrderN for COOP/COHP not implemented yet!'
               else
                  call calc_coop(N1, NL, NR, N1Half, NLHalf, NRHalf, NspinBlocks, NspinComplexMatrix, wk, GammaL, GammaR, gfgeneral, hgeneral,sgeneral, coopinfo%nbond, &
                                 coop1(i,ispin,:),coop2(i,ispin,:),cohp1(i,ispin,:),cohp2(i,ispin,:))
                endif
            endif
!end coop

              IF (NspinBlocks.le.2) THEN

                if(negfon)then
                  CALL ZHEMM('L','U',NL,NR,(1.D0,0.D0),gammal,NL,  transpose(DCONJG(GfRL)),NL,(0.D0,0.D0),  Gamma1_aux,NL)
                  CALL ZHEMM('L','U',NR,NL,(1.D0,0.D0),gammar,NR,  GfRL,NR,(0.D0,0.D0),Gamma2_aux,NR)
                  deallocate(GfRL)
                  if(TransmissionRL)deallocate(GfLR,GfLL)
                else
                  CALL ZHEMM('L','U',NL,NR,(1.D0,0.D0),gammal,NL,  GF_iter_dag(1:NL,N1-NR+1:N1),NL,(0.D0,0.D0),   Gamma1_aux,  NL)
                  CALL ZHEMM('L','U',NR,NL,(1.D0,0.D0),gammar,NR,  gfgeneral%matdense%a(N1-NR+1:N1,1:NL),NR,(0.D0,0.D0), Gamma2_aux, NR)
                endif
              ELSE
                if(negfon)then
                  CALL ZHEMM('L','U',NL,NR,(1.D0,0.D0),gammal,NL,  transpose(DCONJG(GfRL)),NL,(0.D0,0.D0),  Gamma1_aux,NL)
                  CALL ZHEMM('L','U',NR,NL,(1.D0,0.D0),gammar,NR,  GfRL,NR,(0.D0,0.D0),Gamma2_aux,NR)
                  deallocate(GfRL)
                  if(TransmissionRL)then

                    allocate(GammaGf1_RL(nl,nr),GammaGf2_RL(nr,nl))
                    CALL ZHEMM('L','U',NL,NR,(1.D0,0.D0),gammal,NL,  GfLR,NL,(0.D0,0.D0),  GammaGf1_RL,NL)
                    CALL ZHEMM('L','U',NR,NL,(1.D0,0.D0),gammar,NR,  transpose(DCONJG(GfLR)),NR,(0.D0,0.D0),GammaGf2_RL,NR)
                    TEinterm = 0.d0
                    DO J=1,NL
                      DO L=1,NR
                        TEinterm = TEinterm +DREAL(GammaGf1_RL(J,L)*  GammaGf2_RL(L,J))
                      ENDDO
                    ENDDO
                    TERL(I,ISPIN) = TERL(I,ISPIN) + wk*TEinterm
                    deallocate(GammaGf1_RL,GammaGf2_RL)


                    allocate(GammaGf1_LL(nl,nl),GammaGf2_LL(nl,nl))
                    CALL ZHEMM('L','U',NL,NL,(1.D0,0.D0),gammal,NL,  GfLL,NL,(0.D0,0.D0),  GammaGf1_LL,NL)
                    CALL ZHEMM('L','U',NL,NL,(1.D0,0.D0),gammal,NL,  transpose(DCONJG(GfLL)),NL,(0.D0,0.D0),GammaGf2_LL,NL)
                    TEinterm = 0.d0
                    DO J=1,NL
                      DO L=1,NL
                        TEinterm = TEinterm +DREAL(GammaGf1_LL(J,L)*  GammaGf2_LL(L,J))
                      ENDDO
                    ENDDO
                    TELL(I,ISPIN) = TELL(I,ISPIN) + wk*TEinterm
                    deallocate(GammaGf1_LL,GammaGf2_LL)

                    allocate(GammaGf3_LL(nl,nl))
                    CALL ZHEMM('L','U',NL,NL,(1.D0,0.D0),gammal,NL, (0.0D0,1.0D0) * (GfLL-transpose(DCONJG(GfLL))),NL,(0.D0,0.D0),  GammaGf3_LL,NL)
                    TEinterm = 0.d0
                    DO J=1,NL
                      TEinterm = TEinterm +DREAL(GammaGf3_LL(J,J))
                    ENDDO
                    TEL(I,ISPIN) = TEL(I,ISPIN) + wk*TEinterm
                    deallocate(GammaGf3_LL)

                  endif

                  if(TransmissionRL.or.TransmissionMatrix)deallocate(gfLL,GfLR)
                else !negfon
                  ALLOCATE(GF1_aux(NL,NR),GF2_aux(NR,NL))
                  DO j=1,2
                    DO k=1,2
                      GF1_aux((j-1)*NlHalf+1:j*NlHalf,(k-1)*NrHalf+1:k*NrHalf) =  GF_iter_dag((j-1)*N1Half+1:(j-1)*N1Half+NlHalf, k*N1Half-NrHalf+1:k*N1Half)
                      GF2_aux((j-1)*NrHalf+1:j*NrHalf,(k-1)*NlHalf+1:k*NlHalf) =  gfgeneral%matdense%a(j*N1Half-NrHalf+1:j*N1Half,(k-1)*N1Half+1:(k-1)*N1Half+NlHalf)
                    ENDDO
                  ENDDO
                  CALL ZHEMM('L','U',NL,NR,(1.D0,0.D0),gammal,NL,GF1_aux,NL,(0.D0,0.D0),Gamma1_aux,NL)
                  CALL ZHEMM('L','U',NR,NL,(1.D0,0.D0),gammar,NR,GF2_aux,NR,(0.D0,0.D0),Gamma2_aux,NR)
                  DEALLOCATE(GF1_aux,GF2_aux)
                endif !negfon
              ENDIF !nspinblock

              deallocate(gammal,gammar)

              TEinterm = 0.d0
              DO J=1,NL
                DO L=1,NR
                  TEinterm = TEinterm +DREAL(Gamma1_aux(J,L)*  Gamma2_aux(L,J))
                ENDDO
              ENDDO

              deallocate(Gamma1_aux,Gamma2_aux)
              if (TEinterm .GE. 0.0d0) then
                TE(I,ISPIN) = TE(I,ISPIN) + wk*TEinterm
                if (TransmOverk) TEk(I,ISPIN,ik) = TEk(I,ISPIN,ik)+ TEinterm
              endif


!              CALL TIMER('TRANc4',2)
              call DestroyMatrixGeneral(gfgeneral,"transm",io)
              if(negfon.and.(opindex==2.or.opindex==3))call DestroyMatrix(gfgeneral%matdense,"transm",io)
            ENDIF

            if(emtimings)then
              CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
              write(12347,'(A,f12.6)') 'multt',(sc_1-sc_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
            endif
            endif


            call deallocate_selfenergies(i,ispin,ik,ETransmGrid)


        ENDDO spinloop

      ENDDO eneloop


!      CALL TIMER('TRANA2',2)
!      CALL TIMER('TRANA3',1)


      if(mynode_inverse.ne.0)return

      nt=ETransmGrid%nEnergies
      if(emtimings)CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
#ifdef MPI
      if (myhead.NE.0) THEN
        do i=1,nspin
          CALL MPI_SEND(TE(1,i),nt,DAT_double,0,6, inverseheads_comm,MPIerror)
          ! COOP
          if (coopinfo%ccoop) then
            do ibond=1,coopinfo%nbond
              CALL MPI_SEND(coop1(1,i,ibond),nt,DAT_double,0,6, inverseheads_comm,MPIerror)
              CALL MPI_SEND(coop2(1,i,ibond),nt,DAT_double,0,6, inverseheads_comm,MPIerror)
              CALL MPI_SEND(cohp1(1,i,ibond),nt,DAT_double,0,6, inverseheads_comm,MPIerror)
              CALL MPI_SEND(cohp2(1,i,ibond),nt,DAT_double,0,6, inverseheads_comm,MPIerror)
            enddo
          endif
          !
          if(TransmissionRL) then
            CALL MPI_SEND(TERL(1,i),nt,DAT_double,0,6, inverseheads_comm,MPIerror)
            CALL MPI_SEND(TELL(1,i),nt,DAT_double,0,6, inverseheads_comm,MPIerror)
            CALL MPI_SEND(TEL(1,i),nt,DAT_double,0,6, inverseheads_comm,MPIerror)
          endif
          CALL MPI_SEND(NCHANL(1,i),nt,DAT_double,0,6, inverseheads_comm,MPIerror)
          CALL MPI_SEND(NCHANR(1,i),nt,DAT_double,0,6, inverseheads_comm,MPIerror)
          if(TransmOverk)then
            CALL MPI_SEND(nchannelsLK(1,i,ik),nt, mpi_integer,0,6, inverseheads_comm,MPIerror)
            CALL MPI_SEND(nchannelsRK(1,i,ik),nt, mpi_integer,0,6, inverseheads_comm,MPIerror)
          endif
          if(m_dosleads)then
            CALL MPI_SEND(dos(1,i),nt,DAT_double,0,6, inverseheads_comm,MPIerror)
            CALL MPI_SEND(dosvv(1,i),nt,DAT_double,0,6, inverseheads_comm,MPIerror)
          endif
          if(leadsdos)then
            CALL MPI_SEND(dos2(1,i),nt,DAT_double,0,6,  inverseheads_comm,MPIerror)
          endif
          if(TransmissionChannels)then
            CALL MPI_SEND(tchannels(1,1,i),ntrcchannels * nt,DAT_double,0,6, inverseheads_comm,MPIerror)
            if(TransmOverk)then
              CALL MPI_SEND(tchannelsK(1,1,i,ik),ntrcchannels * nt, DAT_double,0,6, inverseheads_comm,MPIerror)
            endif
          endif

          if(GetT0S)then
            CALL MPI_SEND(Tr4G(1,1),2*nt*nTr4G,DAT_double,0,6, inverseheads_comm,MPIerror)
          endif

          if (TransmOverk) CALL MPI_SEND(TEk(1,i,ik),nt,DAT_double,0,14, inverseheads_comm,MPIerror)
        enddo
        do i=1,NspinBlocks
          if(emdos)then
            CALL MPI_SEND(emdostot(1,i),nt,DAT_double,0,6, inverseheads_comm,MPIerror)
            if(empdos)then
              CALL MPI_SEND(empdostot(1,1,i),em_nuo * nt, DAT_double,0,6, inverseheads_comm,MPIerror)
            endif
            if(empdosk)then
              CALL MPI_SEND(empdostotkpt(1,1,i),em_nuo * nt, DAT_double,0,6, inverseheads_comm,MPIerror)
            endif
          endif
        enddo
      else
        DO I=1,nheads-1
          do j=1,nspin
            CALL MPI_RECV(TE(I*nt+1,j), nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
            ! COOP
            if (coopinfo%ccoop) then
               do ibond=1,coopinfo%nbond
                 CALL MPI_RECV(coop1(I*nt+1,j,ibond), nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
                 CALL MPI_RECV(coop2(I*nt+1,j,ibond), nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
                 CALL MPI_RECV(cohp1(I*nt+1,j,ibond), nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
                 CALL MPI_RECV(cohp2(I*nt+1,j,ibond), nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
               enddo
            endif
            !
            if(TransmissionRL)then
              CALL MPI_RECV(TERL(I*nt+1,j), nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
              CALL MPI_RECV(TELL(I*nt+1,j), nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
              CALL MPI_RECV(TEL(I*nt+1,j), nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
            endif
            CALL MPI_RECV(NCHANL(I*nt+1,j),nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
            CALL MPI_RECV(NCHANR(I*nt+1,j),nt,DAT_double,I,6,inverseheads_comm,ISTATUS,  MPIerror)
            if(TransmOverk)then
              call MPI_RECV(nchannelsLK(I*nt+1,j,ik), nt,mpi_integer,I,6,inverseheads_comm,ISTATUS,MPIerror)
              call MPI_RECV(nchannelsRK(I*nt+1,j,ik), nt,mpi_integer,I,6,inverseheads_comm,ISTATUS,MPIerror)
            endif
            if(m_dosleads)then
              CALL MPI_RECV(dos(I*nt+1,j),nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
              CALL MPI_RECV(dosvv(I*nt+1,j), nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
            endif
            if(leadsdos)then
              CALL MPI_RECV(dos2(I*nt+1,j),nt,DAT_double,I,6,inverseheads_comm,ISTATUS,MPIerror)
            endif

            if(TransmissionChannels)then
              CALL MPI_RECV(tchannels(1,I*nt+1,j),ntrcchannels * nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
              if(TransmOverk)then
                CALL MPI_RECV(tchannelsK(1,I*nt+1,j,ik),ntrcchannels * nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
              endif
            endif

            if(GetT0S)then
              CALL MPI_RECV(Tr4G(1,I*nt+1),2*nt*nTr4G,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
            endif
            if (TransmOverk) CALL MPI_RECV(TEk(I*nt+1,j,ik), nt,DAT_double,I,14,inverseheads_comm,ISTATUS,MPIerror)
          enddo

          do j=1,NspinBlocks
            if(emdos)then
              CALL MPI_RECV(emdostot(I*nt+1,j),nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
              if(empdos)then
                CALL MPI_RECV(empdostot(1,I*nt+1,j),em_nuo *nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
              endif
              if(empdosk)then
                CALL MPI_RECV(empdostotkpt(1,I*nt+1,j),em_nuo *nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
              endif
            endif
          enddo
        ENDDO
      endif
#endif

      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
        write(12347,'(A,f12.6)')'sendrt',(sc_1-sc_0)*1.0d0/sc_r
        CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
      endif

      IF(empdosk) THEN
        IF(myhead .EQ. 0) THEN
         
          ALLOCATE(empdosbuf(ETransmGrid%nEnergiesGlobal,em_nuo))
          DO ISPIN=1,NspinBlocks
            DO I=1,ETransmGrid%nEnergiesGlobal
              indt=MOD(i-1,nheads) * nt+ (i-1)/nheads + 1
              empdosbuf(i,:)=empdostotkpt(:,indt,ISPIN)
            ENDDO
            call CopyReorderEMPOS(empdostotkpt(:,:,ispin),empdosbuf,em_nuo*ETransmGrid%nEnergiesGlobal)
          ENDDO
          deallocate(empdosbuf)
       
          call writePDOS(iv,ikpmod,'.TRC.EMPDOS',slabel,NspinBlocks,em_nuo,ETransmGrid%nEnergiesGlobal,DREAL(ETransmGrid%eGlobal),ef,empdostotkpt)
        endif
        deallocate(empdostotkpt)
      ENDIF


      if (TransmOverk) then
        IF (myhead .EQ. 0) THEN
          DO I=1,ETransmGrid%nEnergiesGlobal
            indt=MOD(i-1,nheads) * nt+ (i-1)/nheads + 1
            DO ISPIN = 1, NSPIN
              call io_assign(iukfiletr)
              IF (ISPIN .EQ. 1)filetrk = paste(slabel,'.TRC.k.up')
              IF (ISPIN .EQ. 2)filetrk = paste(slabel,'.TRC.k.down')
              if(nnodes_groupk > 1)then
                write( chivv, '(i7 )' ) mynode_groupk
                filetrk=TRIM(ADJUSTL(filetrk))// TRIM(ADJUSTL("_PartK_"))// TRIM(ADJUSTL(chivv))
              endif



              IF ((ik .EQ. 1).AND.(I.EQ.1)) THEN
                OPEN(UNIT=iukfiletr,FILE=filetrk,status='unknown')
              ELSE
                OPEN(UNIT=iukfiletr,FILE=filetrk,POSITION='append')
              ENDIF

              IF (TEk(indt,ISPIN,ik) .LT. 0.d0)TEk(indt,ISPIN,ik) = 0.d0
!              WRITE(iukfiletr,'(2F16.7,2d16.7)') kpoint,  TEk(indt,ISPIN,ik),DREAL(ETransmGrid%eGlobal(i))
!              WRITE(iukfiletr,'(3F16.7,2d16.7,2i4)') kpoint,wk, TEk(indt,ISPIN,ik),(DREAL(ETransmGrid%eGlobal(i))-ef)*RyToeV,nchannelsLK(indt,ISPIN,ik),nchannelsRK(indt,ISPIN,ik)
               
              if(.not.WriteIkTrcK)then
                WRITE(iukfiletr,'(2F12.6,3x, f20.14,2d16.7,2i4)') kpoint,wk, TEk(indt,ISPIN,ik),(DREAL(ETransmGrid%eGlobal(i))-ef)*RyToeV,nchannelsLK(indt,ISPIN,ik),nchannelsRK(indt,ISPIN,ik)
              else
                WRITE(iukfiletr,'(i9, 2F12.6,3x, f20.14,2d16.7,2i4)') ik, kpoint,wk, TEk(indt,ISPIN,ik),(DREAL(ETransmGrid%eGlobal(i))-ef)*RyToeV,nchannelsLK(indt,ISPIN,ik),nchannelsRK(indt,ISPIN,ik)
              endif
              call io_close(iukfiletr)
            ENDDO
          ENDDO
        ENDIF
      endif
!      CALL TIMER('TRANA3',2)
      IF (myhead .EQ. 0) THEN
        if(TransmissionChannels.and.TransmOverk)then
          call write_trc_channels_K(v,iv,slabel,nheads,nspin,ik,nk,ef,teK(:,:,ik),ntrcchannels,nchannelsLK(:,:,ik),nchannelsRK(:,:,ik),tchannelsK(:,:,:,ik),kpoint,wk)
        endif
      ENDIF

      IF (ik .EQ. nk) Then
!        CALL TIMER('TRANb3',1)
!output format
!E T_tot T_u T_d n_u n_d DOS_u DOS_d DOSvv_U DOSvv_d DOS2_u DOS2_d  (Right lead: n_u n_d DOS_u DOS_d DOSvv_U DOSvv_d DOS2_u DOS2_d  )
!        write(12347,*)"oikhead=",myhead,nheads,mynode_negf,nnodes_negf,mynode_negfo,nnodes_negfo,mynode_inverse,nnodes_inverse,mynode_groupk,nnodes_groupk

#ifdef MPI
        if(myhead == 0 .and. nnodes_groupk > 1)then
          if(.not.TransmissionRL) allocate(TERL(1,1),TELL(1,1),TEL(1,1))
          if(.not.m_dosleads)allocate(dos(1,1),dosvv(1,1))
          if(.not.leadsdos) allocate(dos2(1,1))
          if(.not.emdos) allocate(emdostot(1,1))
          call broadcast_k_TRC(TE,nchanl,terl,tell,tel,dos,dosvv,dos2,emdostot,nchanr,ETransmGrid%nEnergiesGlobal,NSPIN,NspinBlocks,nnodes_groupk,mynode_groupk,groupk_comm,TransmissionRL,m_dosleads,leadsdos,emdos)
          if(.not.TransmissionRL) deallocate(terl,tell,tel)
          if(.not.m_dosleads)deallocate(dos,dosvv)
          if(.not.leadsdos) deallocate(dos2)
          if(.not.emdos) deallocate(emdostot)
!          write( chivv, '(i7 )' ) mynode_groupk
!          filetr=TRIM(ADJUSTL(filetr))// TRIM(ADJUSTL("_PartK_"))// TRIM(ADJUSTL(chivv))
        endif
#endif

        IF (myhead == 0 .and. mynode_groupk == 0 ) THEN 
!          write(12347,*)"ikhead=",myhead,nheads,mynode_negf,nnodes_negf,mynode_negfo,nnodes_negfo,mynode_inverse,nnodes_inverse,mynode_groupk,nnodes_groupk

          call io_assign(iufiletr)
          write( chivv, '(i7 )' ) iv
          chivv = pasbias(chivv, '.')
          filetr = paste(slabel,'.TRC')
          filetr = paste( chivv, filetr )

          ! Meilin Bai, Dec 2012
          if (idyn .gt. 0 .and. idyn .lt. 6) then
            OPEN(UNIT=iufiletr,FILE=filetr,position='append',status='unknown',RECL=100000)
          else
            OPEN(UNIT=iufiletr,FILE=filetr,status='unknown',RECL=100000)
          endif
!          OPEN(UNIT=iufiletr,FILE=filetr,status='unknown',RECL=100000)
          ! End Meilin Bai

          if (frstme) then  !!
!              frstme = .false.  !! Meilin Bai
              if (idyn .gt. 0 .and. idyn .lt. 6) frstme = .false.
              WRITE(iufiletr,'(a6,f11.4,a14,i4)') '# V = ',V, '    k-points: ',nk
              WRITE(iufiletr,'(a18)',ADVANCE='NO') '# energy T_total, '
              DO ISPIN=1,NSPIN
                WRITE(iufiletr,'(a9,i1,a3)',ADVANCE='NO')'T_{ispin=',ispin,'} ,'
              ENDDO

              DO ISPIN=1,NSPIN
                WRITE(iufiletr,'(a18,i1,a3)',ADVANCE='NO') 'NCHAN_{left,ispin=',ispin,'} ,'
              ENDDO

              if(TransmissionRL)then
                DO ISPIN=1,NSPIN
                  WRITE(iufiletr,'(a9,i1,a3)',ADVANCE='NO')'TRL_{ispin=',ispin,'} ,'
                ENDDO
                DO ISPIN=1,NSPIN
                  WRITE(iufiletr,'(a9,i1,a3)',ADVANCE='NO')'TLL_{ispin=',ispin,'} ,'
                ENDDO
                DO ISPIN=1,NSPIN
                  WRITE(iufiletr,'(a9,i1,a3)',ADVANCE='NO')'TL_{ispin=',ispin,'} ,'
                ENDDO
                DO ISPIN=1,NSPIN
                  WRITE(iufiletr,'(a9,i1,a3)',ADVANCE='NO')'T+TLL-TL_{ispin=',ispin,'} ,'
                ENDDO
              endif

              if(m_dosleads)then
                DO ISPIN=1,NSPIN
                  WRITE(iufiletr,'(a20,i1,a3)',ADVANCE='NO')'DOS1_{leads,ispin=',ispin,'} ,'
                ENDDO
                DO ISPIN=1,NSPIN
                  WRITE(iufiletr,'(a18,i1,a3)',ADVANCE='NO')'DOSVV1_{leads,ispin=',ispin,'} ,'
                ENDDO
              endif

              if(leadsdos)then
                DO ISPIN=1,NSPIN
                  WRITE(iufiletr,'(a17,i1,a3)',ADVANCE='NO')  'DOS_{leads,ispin=',ispin,'} ,'
                ENDDO
              endif
              if(emdos)then
                DO ISPIN=1,NspinBlocks
                  WRITE(iufiletr,'(a14,i1,a3)',ADVANCE='NO')'DOS_{EM,ispin=',ispin,'} ,'
                ENDDO
              endif
!              if (GetT0S)then
!                WRITE(iufiletr,'(a24)',ADVANCE='NO') 'REAL(T_ud), IMAG(T_ud), '
!              endif
              if (different_eflr)then
                DO ISPIN=1,NSPIN
                  WRITE(iufiletr,'(a19,i1,a3)',ADVANCE='NO') 'NCHAN_{right,ispin=',ispin,'} ,'
                ENDDO
              endif
              WRITE(iufiletr,*)
          endif !! end frstme

          ! Meilin Bai, Dec 2012
          if (idyn .gt. 0 .and. idyn .lt. 6) write(iufiletr, '(a, i7)') '# istep = ', istep
          ! End Meilin Bai


          DO I=1,ETransmGrid%nEnergiesGlobal
            indt=MOD(i-1,nheads) * nt+ (i-1)/nheads + 1

              WRITE(iufiletr,'(2e16.7,A3)',ADVANCE='NO')(DREAL(ETransmGrid%eGlobal(i))-ef)* 13.6056981D0,SUM(TE(indt,:)), '   '
              DO ISPIN=1,NSPIN
                WRITE(iufiletr,'(e16.7,A3)',ADVANCE='NO')TE(indt,ISPIN),'   '
              ENDDO
              DO ISPIN=1,NSPIN
                WRITE(iufiletr,'(e16.7,A3)',ADVANCE='NO')NCHANL(indt,ISPIN),'   '
              ENDDO

              if(TransmissionRL)then
                DO ISPIN=1,NSPIN
                  WRITE(iufiletr,'(e16.7,A3)',ADVANCE='NO')TERL(indt,ISPIN),'   '
                ENDDO
                DO ISPIN=1,NSPIN
                  WRITE(iufiletr,'(e16.7,A3)',ADVANCE='NO')TELL(indt,ISPIN),'   '
                ENDDO
                DO ISPIN=1,NSPIN
                  WRITE(iufiletr,'(e16.7,A3)',ADVANCE='NO')TEL(indt,ISPIN),'   '
                ENDDO
                DO ISPIN=1,NSPIN
                  WRITE(iufiletr,'(e16.7,A3)',ADVANCE='NO')TE(indt,ISPIN)+TELL(indt,ISPIN)-TEL(indt,ISPIN),'   '
                ENDDO
              endif

              if(m_dosleads)then
                DO ISPIN=1,NSPIN
                  WRITE(iufiletr,'(d16.7,A3)',ADVANCE='NO')dos(indt,ISPIN),'   '
                ENDDO
                DO ISPIN=1,NSPIN
                  WRITE(iufiletr,'(d16.7,A3)',ADVANCE='NO')dosvv(indt,ISPIN),'   '
                ENDDO
              endif
              if(leadsdos)then
                DO ISPIN=1,NSPIN
                  WRITE(iufiletr,'(d16.7,A3)',ADVANCE='NO')dos2(indt,ISPIN),'   '
                ENDDO
              endif
              if(emdos)then
                DO ISPIN=1,NspinBlocks
                  WRITE(iufiletr,'(d16.7,A3)',ADVANCE='NO')emdostot(indt,ISPIN),'   '
                ENDDO
              endif
!              if(GetT0S)then
!                DO ISPIN=1,NSPIN
!                  WRITE(iufiletr,'(d16.7,A3)',ADVANCE='NO')TudTot(indt,ISPIN),'   '
!                ENDDO
!              endif
              if (different_eflr)then
                DO ISPIN=1,NSPIN
                  WRITE(iufiletr,'(e16.7,A3)',ADVANCE='NO')NCHANR(indt,ISPIN),'   '
                ENDDO
              endif
              WRITE(iufiletr,*)
          ENDDO
          WRITE(iufiletr,*)
          call io_close(iufiletr)
        ENDIF
!        CALL TIMER('TRANb3',2)
!        CALL TIMER('TRANA4',1)


!!!COOP/COHP
        if (coopinfo%ccoop .and. myhead .eq. 0) call output_coop(slabel, nheads, ETransmGrid%nEnergiesGlobal, nt, iv, nk, nspin, coopinfo%nbond, V, DREAL(ETransmGrid%eGlobal)-ef, coop1, coop2, cohp1, cohp2)
        if (coopinfo%ccoop)  deallocate(coop1, coop2, cohp1, cohp2)

!!!
        IF(empdos.and.myhead .EQ. 0) THEN

#ifdef MPI
          if(nnodes_groupk > 1)then
            call broadcast_k_PDOS(empdostot,em_nuo,ETransmGrid%nEnergiesGlobal,NspinBlocks,nnodes_groupk,mynode_groupk,groupk_comm)
          endif
#endif

          if(mynode_groupk==0)then
           
            ALLOCATE(empdosbuf(ETransmGrid%nEnergiesGlobal,em_nuo))
            DO ISPIN=1,NspinBlocks
              DO I=1,ETransmGrid%nEnergiesGlobal
                indt=MOD(i-1,nheads) * nt+ (i-1)/nheads + 1
              
                empdosbuf(i,:)=empdostot(:,indt,ISPIN)
              ENDDO
              call CopyReorderEMPOS(empdostot(:,:,ispin),empdosbuf,em_nuo*ETransmGrid%nEnergiesGlobal)
            ENDDO
            deallocate(empdosbuf)
         
            call writePDOS(iv,0,'.TRC.EMPDOS',slabel,NspinBlocks,em_nuo,ETransmGrid%nEnergiesGlobal,DREAL(ETransmGrid%eGlobal),ef,empdostot)
          endif
        ENDIF
!        CALL TIMER('TRANA4',2)
!        CALL TIMER('TRANA5',1)

        IF (myhead .EQ. 0.and.mynode_groupk==0) THEN
          if(.not.emdos) allocate(emdostot(1,1))
          if(.not.leadsdos) allocate(dos2(1,1))
          call write_trc_agr(v,iv,slabel,deltaenes, nheads,nspin,ef,te,emdostot, nchanl,nchanr,dos2)
          if(.not.emdos) deallocate(emdostot)
          if(.not.leadsdos) deallocate(dos2)
          if(TransmissionChannels)then
            call write_trc_channels(v,iv,slabel,nheads,nspin,ef,te,ntrcchannels,tchannels)
          endif
          if(GetT0S)then
            call write_trc_spindecomposition(v,iv,slabel,nheads,nspin,ef,te,nTr4G,Tr4G)
          endif
          if(emldos2)then
            call average_trc(te,nspin,ETransmGrid%nEnergiesGlobal)
          endif
        ENDIF
!        CALL TIMER('TRANSI',2)

        DEALLOCATE(TE,NCHANL,NCHANR)
        if(TransmissionRL)then
          deallocate(TERL)
          deallocate(TELL,TEL)
        endif

        if(m_dosleads)DEALLOCATE(dos,dosvv)
        if(allocated(dos2))DEALLOCATE(dos2)
        if(allocated(emdostot))DEALLOCATE(emdostot)
        if(allocated(Tr4G))DEALLOCATE(Tr4G)
        if(empdos)DEALLOCATE(empdostot)
        if (TransmOverk) DEALLOCATE(TEk)
        if (TransmOverk) DEALLOCATE(nchannelsLK,nchannelsRK)
        if(TransmissionChannels)then
          DEALLOCATE(tchannels)
          if(TransmOverk)then
            DEALLOCATE(tchannelsK)
          endif
        endif
        call deallocate_energygrid_selfenergies2(ETransmGrid)
!        CALL TIMER('TRANA5',2)
      ENDIF

      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
        write(12347,'(A,f12.6)')'writet',(sc_1-sc_0)*1.0d0/sc_r
        CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
      endif
      if(allocated(work))deallocate(work,GF_iter_dag)

      END SUBROUTINE transm



      SUBROUTINE write_trc_agr(v,iv,slabel,deltaenes, nnodes,nspin,ef,te,emdostot,nchanl,nchanr,dos2)

      USE precision
      use negfmod
      use mEnergyGrid

      IMPLICIT NONE

      integer iv,Nnodes,indt,nspin
      integer iufiletr,I,ispin,iplot,igraph,iset
      CHARACTER :: chivv*8
      CHARACTER :: slabel*20, paste*35, filetr*35
      DOUBLE PRECISION :: deltaenes,ef,v
      character pasbias*25
      double precision TE(ETransmGrid%nEnergiesGlobal,NSPIN), emdostot(ETransmGrid%nEnergiesGlobal,NSPIN),dos2(ETransmGrid%nEnergiesGlobal,NSPIN), nchanl(ETransmGrid%nEnergiesGlobal,NSPIN),nchanr(ETransmGrid%nEnergiesGlobal,NSPIN), dataplot(ETransmGrid%nEnergiesGlobal,NSPIN)

      call io_assign(iufiletr)
      write( chivv, '(i7 )' ) iv
      chivv = pasbias(chivv, '.')
      filetr = paste(slabel,'_TRC.agr')
      filetr = paste( chivv, filetr )
      OPEN(UNIT=iufiletr,FILE=filetr,status='unknown',RECL=100000)

      write(iufiletr,*)"@g0 on"
      write(iufiletr,*)"@g0 hidden false"
      write(iufiletr,*)"@with g0"
      write(iufiletr,*)"@   view 0.150000, 0.56, 0.55, 0.850000"
      write(iufiletr,*)'@   xaxis label "E-E\sF\N (eV)"'
      write(iufiletr,*)'@   yaxis label "Transmission"'
      write(iufiletr,*)'@   s0 line color 2'
      write(iufiletr,*)'@   s1 line color 3'
      write(iufiletr,*)"@g1 on"
      write(iufiletr,*)"@g1 hidden false"
      write(iufiletr,*)"@with g1"
      write(iufiletr,*)"@   view 0.692246, 0.56, 1.10, 0.850000"
      write(iufiletr,*)'@   xaxis label "E-E\sF\N (eV)"'
      write(iufiletr,*)'@   yaxis label "number of channels"'
      write(iufiletr,*)'@   s0 line color 2'
      write(iufiletr,*)'@   s1 line color 3'
      write(iufiletr,*)'@   s2 line color 4'
      write(iufiletr,*)'@   s3 line color 10'
      write(iufiletr,*)"@g2 on"
      write(iufiletr,*)"@g2 hidden false"
      write(iufiletr,*)"@with g2"
      write(iufiletr,*)"@   view 0.150000, 0.150000, 0.55, 0.44"
      write(iufiletr,*)'@   xaxis label "E-E\sF\N (eV)"'
      write(iufiletr,*)'@   yaxis label "EM DOS"'
      write(iufiletr,*)'@   s0 line color 2'
      write(iufiletr,*)'@   s1 line color 3'
      write(iufiletr,*)"@g3 on"
      write(iufiletr,*)"@g3 hidden false"
      write(iufiletr,*)"@with g3"
      write(iufiletr,*)"@   view 0.692246, 0.150000, 1.10, 0.44"
      write(iufiletr,*)'@   xaxis label "E-E\sF\N (eV)"'
      write(iufiletr,*)'@   yaxis label "Leads DOS"'
      write(iufiletr,*)'@   s0 line color 2'
      write(iufiletr,*)'@   s1 line color 3'

      igraph=0
      do iplot=1,5

        igraph=igraph+1
        iset=0
        if(iplot.eq.3)igraph=igraph-1
        if(iplot.eq.3)iset=2

        if(iplot.eq.1)dataplot=te
        if(iplot.eq.2)dataplot=nchanl
        if(iplot.eq.3)dataplot=nchanr
        if(iplot.eq.4)then
          if(emdos)then
            dataplot=emdostot
          else
            cycle
          endif
        endif
        if(iplot.eq.5)then
          if(leadsdos)then
            dataplot=dos2
          else
            cycle
          endif
        endif

        DO ISPIN=1,NSPIN
          write(iufiletr,'(a10,i1,a2,i1)')" @target G",igraph-1,".S",  ISPIN-1+iset
          write(iufiletr,*)"@type xy"

          DO I=1,ETransmGrid%nEnergiesGlobal
            indt=MOD(i-1,Nnodes) * ETransmGrid%nEnergies+(i-1)/Nnodes+1

            WRITE(iufiletr,'(e16.7,A3,e16.7)') (DREAL(ETransmGrid%eGlobal(i))-ef)* 13.6056981D0,'   ', dataplot(indt,ISPIN)

          ENDDO
          write(iufiletr,*)"&"
          write(iufiletr,*)"@autoscale"
        ENDDO
      enddo
      WRITE(iufiletr,*)
      call io_close(iufiletr)
      END SUBROUTINE write_trc_agr


      SUBROUTINE write_trc_spindecomposition(v,iv,slabel,nnodes,nspin,ef,te,ntrcchannels,tchannels)

      USE precision
      use negfmod
      use mEnergyGrid

      IMPLICIT NONE

      integer, intent(in) :: ntrcchannels,nspin
      double complex, intent(in) :: tchannels(ntrcchannels,ETransmGrid%nEnergiesGlobal)
      integer iv,Nnodes,indt
      integer iufiletr,I,ispin,iplot,igraph,iset,ic
      CHARACTER :: chivv*8
      CHARACTER :: slabel*20, paste*35, filetr*35
      DOUBLE PRECISION :: ef,v
      character pasbias*25
      double precision TE(ETransmGrid%nEnergiesGlobal,NSPIN), emdostot(ETransmGrid%nEnergiesGlobal,NSPIN),dos2(ETransmGrid%nEnergiesGlobal,NSPIN), nchanl(ETransmGrid%nEnergiesGlobal,NSPIN),nchanr(ETransmGrid%nEnergiesGlobal,NSPIN), dataplot(ETransmGrid%nEnergiesGlobal,NSPIN)

      call io_assign(iufiletr)
      write( chivv, '(i7 )' ) iv
      chivv = pasbias(chivv, '.')
      filetr = paste(slabel,'_TRC_Spin.dat')
      filetr = paste( chivv, filetr )
      OPEN(UNIT=iufiletr,FILE=filetr,status='unknown',RECL=100000)


      WRITE(iufiletr,'(a2)',ADVANCE='NO')'# '
      WRITE(iufiletr,'(a15)',ADVANCE='NO')'         E-E_F,'

      ispin=1
      do ic=1,ntrcchannels
        WRITE(iufiletr,'(a9,i6,a4)',ADVANCE='NO')'    Re[T(',ic,')], '
        WRITE(iufiletr,'(a9,i6,a4)',ADVANCE='NO')'    Im[T(',ic,')], '
      enddo

      DO ispin=1,nspin
          WRITE(iufiletr,'(a10,i6,a3)',ADVANCE='NO')' T_{ispin=',ispin,'}, '
      enddo

      WRITE(iufiletr,'(a1)')' ' 


      DO I=1,ETransmGrid%nEnergiesGlobal
        indt=MOD(i-1,Nnodes) * ETransmGrid%nEnergies+(i-1)/Nnodes+1
        WRITE(iufiletr,'(e16.7,A2)',ADVANCE='NO') (DREAL(ETransmGrid%eGlobal(i))-ef)* 13.6056981D0,'  '

        do ic=1,ntrcchannels
          WRITE(iufiletr,'(e16.7,A3)',ADVANCE='NO') DREAL(tchannels(ic,indt)),'   '
          WRITE(iufiletr,'(e16.7,A3)',ADVANCE='NO') DIMAG(tchannels(ic,indt)),'   '
        ENDDO

        DO ispin=1,nspin
          WRITE(iufiletr,'(e16.7,A3)',ADVANCE='NO') te(indt,ISPIN),'   '
        ENDDO
        WRITE(iufiletr,'(a1)')' ' 
      enddo
      WRITE(iufiletr,*)
      call io_close(iufiletr)

      END SUBROUTINE write_trc_spindecomposition



      SUBROUTINE write_trc_channels(v,iv,slabel,nnodes,nspin,ef,te,ntrcchannels,tchannels)

      USE precision
      use negfmod
      use mEnergyGrid

      IMPLICIT NONE

      integer, intent(in) :: ntrcchannels,nspin
      double precision, intent(in) :: tchannels(ntrcchannels,ETransmGrid%nEnergiesGlobal,nspin)
      integer iv,Nnodes,indt
      integer iufiletr,I,ispin,iplot,igraph,iset,ic
      CHARACTER :: chivv*8
      CHARACTER :: slabel*20, paste*35, filetr*35
      DOUBLE PRECISION :: ef,v
      character pasbias*25
      double precision TE(ETransmGrid%nEnergiesGlobal,NSPIN), emdostot(ETransmGrid%nEnergiesGlobal,NSPIN),dos2(ETransmGrid%nEnergiesGlobal,NSPIN), nchanl(ETransmGrid%nEnergiesGlobal,NSPIN),nchanr(ETransmGrid%nEnergiesGlobal,NSPIN), dataplot(ETransmGrid%nEnergiesGlobal,NSPIN)

      call io_assign(iufiletr)
      write( chivv, '(i7 )' ) iv
      chivv = pasbias(chivv, '.')
      filetr = paste(slabel,'_TRC_Channels.dat')
      filetr = paste( chivv, filetr )
      OPEN(UNIT=iufiletr,FILE=filetr,status='unknown',RECL=100000)


      WRITE(iufiletr,'(a2)',ADVANCE='NO')'# '
      WRITE(iufiletr,'(a7)',ADVANCE='NO')'E-E_F, '
      DO ispin=1,nspin
          WRITE(iufiletr,'(a9,i1,a3)',ADVANCE='NO')'T_{ispin=',ispin,'} ,'
      enddo

      do ic=1,ntrcchannels
        DO ispin=1,nspin
            WRITE(iufiletr,'(a11,i2,a7,i1,a3)',ADVANCE='NO')'T_{channel=',ic,',ispin=',ispin,'} ,'
        enddo
      enddo
      WRITE(iufiletr,'(a1)')' ' 


      DO I=1,ETransmGrid%nEnergiesGlobal
        indt=MOD(i-1,Nnodes) * ETransmGrid%nEnergies+(i-1)/Nnodes+1
        WRITE(iufiletr,'(e16.7,A2)',ADVANCE='NO') (DREAL(ETransmGrid%eGlobal(i))-ef)* 13.6056981D0,'  '

        DO ispin=1,nspin
          WRITE(iufiletr,'(e16.7,A3)',ADVANCE='NO') te(indt,ISPIN),'   '
        ENDDO
        do ic=1,ntrcchannels
          DO ispin=1,nspin
            WRITE(iufiletr,'(e16.7,A3)',ADVANCE='NO') tchannels(ic,indt,ISPIN),'   '
          ENDDO
        ENDDO
        WRITE(iufiletr,'(a1)')' ' 
      enddo
      WRITE(iufiletr,*)
      call io_close(iufiletr)

      END SUBROUTINE write_trc_channels



      SUBROUTINE write_trc_channels_K(v,iv,slabel,nnodes,nspin,ik,nk,ef,tK,ntrcchannels,nchannelsLK,nchannelsRK,tchannelsK,kpoint,wk)

      USE precision
      use negfmod
      use mEnergyGrid

      IMPLICIT NONE

      integer, intent(in) :: ntrcchannels,nk,nspin,ik
      double precision, intent(in) :: tchannelsK(ntrcchannels,ETransmGrid%nEnergiesGlobal,nspin)
      integer, intent(in) :: nchannelsLK(ETransmGrid%nEnergiesGlobal,nspin)
      integer, intent(in) :: nchannelsRK(ETransmGrid%nEnergiesGlobal,nspin)
      double precision, intent(in) :: tK(ETransmGrid%nEnergiesGlobal,nspin)
      double precision, intent(in) :: kpoint(2)
      double precision, intent(in) :: wk
      integer iv,Nnodes,indt
      integer iufiletr,I,ispin,iplot,igraph,iset,ic
      CHARACTER :: chivv*8
      CHARACTER :: slabel*20, paste*35, filetr*35
      DOUBLE PRECISION :: ef,v
      character pasbias*25

      call io_assign(iufiletr)
      write( chivv, '(i7 )' ) iv
      chivv = pasbias(chivv, '.')
      filetr = paste(slabel,'_TRC_Channels_K.dat')
      filetr = paste( chivv, filetr )

      if(ik.eq.1)then
        OPEN(UNIT=iufiletr,FILE=filetr,status='unknown',RECL=100000)
        WRITE(iufiletr,'(a2)',ADVANCE='NO')'# '
        WRITE(iufiletr,'(a21)',ADVANCE='NO')'E-E_F, k_x, k_y, wk, '
        DO ispin=1,nspin
            WRITE(iufiletr,'(a10,i1,a3)',ADVANCE='NO')'Tk_{ispin=',ispin,'} ,'
        enddo
  
        do ic=1,ntrcchannels
          DO ispin=1,nspin
              WRITE(iufiletr,'(a12,i2,a7,i1,a3)',ADVANCE='NO')'Tk_{channel=',ic,',ispin=',ispin,'} ,'
          enddo
        enddo
        DO ispin=1,nspin
            WRITE(iufiletr,'(a11,i2,a3)',ADVANCE='NO')'NLk_{ispin=',ispin,'} ,'
        enddo
        DO ispin=1,nspin
            WRITE(iufiletr,'(a11,i2,a3)',ADVANCE='NO')'NRk_{ispin=',ispin,'} ,'
        enddo
        WRITE(iufiletr,'(a1)')' ' 
      else
        OPEN(UNIT=iufiletr,FILE=filetr,position='append',RECL=100000)
      endif
      WRITE(iufiletr,'(a5,i6)')'# ik=',ik 


      DO I=1,ETransmGrid%nEnergiesGlobal
        indt=MOD(i-1,Nnodes) * ETransmGrid%nEnergies+(i-1)/Nnodes+1
        WRITE(iufiletr,'(e16.7,A2)',ADVANCE='NO') (DREAL(ETransmGrid%eGlobal(i))-ef)* 13.6056981D0,'  '
        WRITE(iufiletr,'(e16.7,A2,e16.7,A2)',ADVANCE='NO') kpoint(1),'  ',kpoint(2),'  '
        WRITE(iufiletr,'(e16.7,A2)',ADVANCE='NO') wk,'  '

        DO ispin=1,nspin
          WRITE(iufiletr,'(e16.7,A3)',ADVANCE='NO') tK(indt,ISPIN),'   '
        ENDDO
        do ic=1,ntrcchannels
          DO ispin=1,nspin
            WRITE(iufiletr,'(e16.7,A3)',ADVANCE='NO') tchannelsK(ic,indt,ISPIN),'   '
          ENDDO
        ENDDO
        DO ispin=1,nspin
          WRITE(iufiletr,'(i4,A3)',ADVANCE='NO') nchannelsLK(indt,ISPIN),'   '
        ENDDO
        DO ispin=1,nspin
          WRITE(iufiletr,'(i4,A3)',ADVANCE='NO') nchannelsRK(indt,ISPIN),'   '
        ENDDO
        WRITE(iufiletr,'(a1)')' ' 
      enddo
!      WRITE(iufiletr,*)
      call io_close(iufiletr)

      END SUBROUTINE write_trc_channels_K




  subroutine fermi_distribution_real(e,t,f)
 
    double precision,intent(in) :: e,t
    double precision,intent(out) :: f

    double precision eovert

    eovert=e/t

    IF ( eovert .GT. 50.D0) THEN
      f=0.D0
    ELSE
      f=1.D0/(1+DEXP(eovert))
    ENDIF
    
   end subroutine fermi_distribution_real




  subroutine writePDOS(iv,ik,suffix,slabel,nspin,nuo,nene,ei,ef,empdostot)

  implicit none

  integer, intent(in) :: iv,ik,nspin,nene,nuo
  CHARACTER(LEN=*), intent(in) :: suffix
  CHARACTER(LEN=*), intent(in) :: slabel
  double precision, intent(in) :: ei(nene),ef
  double precision, intent(in) :: empdostot(nene,nuo,nspin)

  CHARACTER :: chivv*8
  CHARACTER :: fileempdos*100
  INTEGER :: i,ii,iufileempdos,ispin

!  write(12347,*)"ninfo=",n,nspin,nene,nuo

  call io_assign(iufileempdos)
  if(ik>0)then
    write( chivv, '(i7 )' ) ik
    chivv = '_'//TRIM(ADJUSTL(chivv))
    fileempdos =   TRIM(ADJUSTL(chivv)) //  TRIM(ADJUSTL(suffix))
    fileempdos =   TRIM(ADJUSTL(slabel)) //  TRIM(ADJUSTL(fileempdos))
  elseif(ik<0)then
    write( chivv, '(i7 )' ) -ik-1
    fileempdos =  TRIM(ADJUSTL(slabel)) //  TRIM(ADJUSTL(suffix))// TRIM(ADJUSTL("_PartK_"))// TRIM(ADJUSTL(chivv))
  else
    fileempdos =  TRIM(ADJUSTL(slabel)) //  TRIM(ADJUSTL(suffix))
  endif
  if(iv >= 0 ) then
    write( chivv, '(i7 )' ) iv
    chivv = TRIM(ADJUSTL(chivv))// '.'
    fileempdos =  TRIM(ADJUSTL(chivv)) //  TRIM(ADJUSTL(fileempdos))
  endif
!  write(12347,*)"fileempdos=",fileempdos

  OPEN(UNIT=iufileempdos,FILE=fileempdos,status='unknown', RECL=100000)

  WRITE(iufileempdos,'(a6)') "<pdos>"
  WRITE(iufileempdos,'(a7,i1,a8)') "<nspin>",nspin,"</nspin>"
  WRITE(iufileempdos,'(a11,i6,a12)') "<norbitals>",nuo,"</norbitals>"
  WRITE(iufileempdos,'(a26)') "<energy_values units=""eV"">"
  DO I=1,nene
      WRITE(iufileempdos,'(f25.5)') (ei(i)-ef)* 13.6056981D0
  ENDDO
  WRITE(iufileempdos,'(a16)') "</energy_values>"

  do ii=1,nuo
    call xml_empdos(iufileempdos,ii)

    DO I=1,nene
      DO ISPIN=1,nspin
        WRITE(iufileempdos,'(f10.5,a3)',ADVANCE='NO')  empdostot(i,ii,ISPIN),'   '
      ENDDO
      WRITE(iufileempdos,*)
    ENDDO
    WRITE(iufileempdos,'(a7)') "</data>"
    WRITE(iufileempdos,'(a10)') "</orbital>"
  enddo
  WRITE(iufileempdos,'(a7)') "</pdos>"
  call io_close(iufileempdos)

  end subroutine writePDOS

  subroutine CopyReorderEMPOS(l1,l2,n)

  implicit none
  integer,intent(in) :: n
  double precision, intent(out):: l1(n)
  double precision, intent(in):: l2(n)

  integer i

  do i=1,n
    l1(i)=l2(i)
  enddo

  end subroutine CopyReorderEMPOS


  subroutine transm_wrap(ik,nk,N1,nl,nr,NspinBlocks, NspinComplexMatrix, &
  V,IV,slabel,kpoint,wk, Ef_Lead,kb,t, hgeneral,sgeneral,rhogeneral,ematgeneral,&
   istep, inicoor, idyn)
!
! Add inicoor, istep, and idyn for MD output
! Meilin Bai, Dec 2012
!

  use sigma
  use mTypes

  implicit none
  integer, intent(in) ::NspinComplexMatrix,NspinBlocks
  integer ik, nk,n1,nl,nr,iv,NspinLoops
  integer istep, inicoor, idyn ! Meilin Bai, Dec 2012
  double precision kbt,v,wk,Ef_Lead,kb,t
  DOUBLE PRECISION, DIMENSION (3) :: kpoint
  CHARACTER(LEN=*) :: slabel
  type(matrixTypeGeneral) :: hgeneral(NspinComplexMatrix),sgeneral,rhogeneral(NspinComplexMatrix),ematgeneral(NspinComplexMatrix)


  kbt=kb*T

  NspinLoops=NspinComplexMatrix
  if(NspinComplexMatrix.gt.2)NspinLoops=1

  CALL TRANSM(ik,nk,N1,nl,nr,NspinLoops,V,IV, NspinBlocks,NspinComplexMatrix, &
  slabel,kpoint(1:2),wk, H0_L,H0_R,S0_L,S0_R,S1_L,S1_R,H1_L,H1_R, Ef_Lead,kbt,&
  hgeneral,sgeneral,rhogeneral,ematgeneral, istep, inicoor, idyn)

  end subroutine transm_wrap




  subroutine average_trc(te,nspin,ntrc)

  use negfmod, only :TrcAverage
  use mMPI_NEGF

  implicit none
  integer, intent(in) :: nspin,ntrc
  double precision, intent(in) ::TE(ntrc,NSPIN)
 
  integer i,is,mpierror
  double precision, allocatable :: Tbuf(:)

  if(.not.allocated(TrcAverage))allocate(TrcAverage(nspin))
  
  TrcAverage=0.0D0

  if(ntrc==1)then
    TrcAverage(:)=te(1,:)
  else
    TrcAverage(:)=0.5D0 * te(1,:)
    do i=2,ntrc-1
      TrcAverage(:)=TrcAverage(:)+te(i,:)
    enddo
    TrcAverage(:)=TrcAverage(:)+0.5D0 * te(ntrc,:)

    TrcAverage=TrcAverage/(1.0D0*((ntrc-1)))
  endif

!  write(12347,*)"trcaverage_local=",trcaverage(1),te(1,1)

#ifdef MPI
  allocate(Tbuf(nspin))
  Tbuf=TrcAverage
  TrcAverage=0.0D0
  call MPI_REDUCE(Tbuf,TrcAverage,nspin,DAT_double,MPI_SUM,0,groupk_comm, MPIerror)
  deallocate(tbuf)
!  write(12347,*)"trcaverage=",trcaverage(1)
#endif


  end subroutine average_trc


  subroutine broadcast_k_TRC(transm,nchanl,terl,tell,tel,dos,dosvv,dos2,emdostot,nchanr,n,nspin,NspinBlocks,nnodes,node,comm,TransmissionRL,m_dosleads,leadsdos,emdos)

    use mConstants
    use mMPI_NEGF

    implicit none

    logical, intent(in)      :: TransmissionRL,m_dosleads,leadsdos,emdos
    integer, intent(in)      :: n,nspin,nnodes,node,comm,NspinBlocks
    real(kdp), intent(inout) :: transm(n,nspin)
    real(kdp), intent(inout) :: NCHANL(n,nspin)
    real(kdp), intent(inout) :: TERL(n,nspin)
    real(kdp), intent(inout) :: TELL(n,nspin)
    real(kdp), intent(inout) :: TEL(n,nspin)
    real(kdp), intent(inout) :: dos(n,nspin)
    real(kdp), intent(inout) :: dosvv(n,nspin)
    real(kdp), intent(inout) :: dos2(n,nspin)
    real(kdp), intent(inout) :: emdostot(n,NspinBlocks)
    real(kdp), intent(inout) :: NCHANR(n,nspin)

    real(kdp), allocatable   :: tbuf(:,:)
    integer mpierror

    allocate(tbuf(n,nspin))
    call MPI_REDUCE(transm,tbuf,n*nspin,DAT_double,MPI_SUM,0,comm, MPIerror)
    if(node==0) transm=tbuf
    call MPI_REDUCE(nchanl,tbuf,n*nspin,DAT_double,MPI_SUM,0,comm, MPIerror)
    if(node==0) nchanl=tbuf
    call MPI_REDUCE(nchanr,tbuf,n*nspin,DAT_double,MPI_SUM,0,comm, MPIerror)
    if(node==0) nchanr=tbuf

    if(TransmissionRL) then
      call MPI_REDUCE(TERL,tbuf,n*nspin,DAT_double,MPI_SUM,0,comm, MPIerror)
      if(node==0)TERL=tbuf
      call MPI_REDUCE(TELL,tbuf,n*nspin,DAT_double,MPI_SUM,0,comm, MPIerror)
      if(node==0)TELL=tbuf
      call MPI_REDUCE(TEL,tbuf,n*nspin,DAT_double,MPI_SUM,0,comm, MPIerror)
      if(node==0)TEL=tbuf
    endif

    if(m_dosleads) then
      call MPI_REDUCE(dos,tbuf,n*nspin,DAT_double,MPI_SUM,0,comm, MPIerror)
      if(node==0)dos=tbuf
      call MPI_REDUCE(dosvv,tbuf,n*nspin,DAT_double,MPI_SUM,0,comm, MPIerror)
      if(node==0)dosvv=tbuf
    endif

    if(leadsdos) then
      call MPI_REDUCE(dos2,tbuf,n*nspin,DAT_double,MPI_SUM,0,comm, MPIerror)
      if(node==0)dos2=tbuf
    endif
    deallocate(tbuf)
 
    if(emdos) then
      allocate(tbuf(n,NspinBlocks))
      call MPI_REDUCE(emdostot,tbuf,n*NspinBlocks,DAT_double,MPI_SUM,0,comm, MPIerror)
      if(node==0)emdostot=tbuf
      deallocate(tbuf)
    endif



  end subroutine broadcast_k_TRC

  subroutine broadcast_k_PDOS(empdostot,nuo,n,nspin,nodes,node,comm)

  use mConstants
  use mMPI_NEGF

  implicit none

  integer, intent(in)      :: n,nspin,nuo,nodes,node,comm
  real(kdp), intent(inout) :: empdostot(nuo,n,nspin)

  real(kdp), allocatable   :: tbuf(:,:,:)
  integer mpierror

  allocate(tbuf(nuo,n,nspin))
  call MPI_REDUCE(empdostot,tbuf,nuo*n*nspin,DAT_double,MPI_SUM,0,comm, MPIerror)
  if(node==0) empdostot=tbuf
  deallocate(tbuf)

  end subroutine broadcast_k_PDOS
