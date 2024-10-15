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
  SUBROUTINE loe(ik,nk,N1,NL,NR,NSPIN,V,IV, matSpin,slabel,kpoint,wk,H0L,H0R,S0L,S0R,S1Li,S1Ri,H1Li,H1Ri,ef,T,hgeneral,sgeneral,rhogeneral,ematgeneral)

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
! Altered by Tatsuhiko Ohto, March 2013
!   T function due to the e-ph coupling by LOE
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
! CHARACTER side_rankL(NSPIN)   : Side to which transformation must be
!                                 applied (left lead)
! CHARACTER side_rankR(NSPIN)   : Side to which transformaiont must be
!                                 applied (right lead)
! character*20  slabel          : System label
! real*8   wk                   : Weight for k points
! *********************** SAVED VARIABLES ***************************
! integer  NeneT              : Number of energy points for Transmission
! double precision TEnergI        : Initial Energy for Transmission
! double precision TEnergF        : Final Energy for Transmission
! ************************ OUTPUT/SAVED ******************************
! t0(Nenerg,NSPIN)              : Transmission Coefficients
! ********************************************************************

      USE precision
      USE parallel
      use sigma, only: sigma_method
      use negfmod
      use mTypes
      use mMatrixUtil
      use mSigmaMethod1
      use mMPI_NEGF
      use mEnergyGrid

      IMPLICIT NONE

#ifdef MPI
      INTEGER, DIMENSION (MPI_STATUS_SIZE) :: istatus
      INTEGER :: MPIerror
#endif

      INTEGER :: Nnodes
      INTEGER :: N1,NL,NR,NSPIN,IV,matSpin
      CHARACTER :: slabel*20, paste*35, filetr*35, filetrk*35,chivv*8

      character pasbias*25

      DOUBLE PRECISION, ALLOCATABLE, SAVE :: t0(:,:), t_in(:,:,:), t_inL(:,:,:), t_inR(:,:,:), t_eh(:,:,:), t_ec(:,:,:), t_ecL(:,:,:), t_ecR(:,:,:), t_AsymL(:,:,:), t_AsymR(:,:,:), t_inel(:,:,:), it_inel(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE   :: Afull(:,:,:)
      DOUBLE PRECISION, DIMENSION (2) :: kpoint
      DOUBLE COMPLEX, DIMENSION (NL,NL) :: S0L
      DOUBLE COMPLEX, DIMENSION (NL,NL) :: S1Li
      DOUBLE COMPLEX, DIMENSION (NL,NL,NSPIN) :: H0L,H1L,H1Li,S1L
      DOUBLE COMPLEX, DIMENSION (NR,NR) :: S0R
      DOUBLE COMPLEX, DIMENSION (NR,NR) :: S1Ri
      DOUBLE COMPLEX, DIMENSION (NR,NR,NSPIN) :: H0R,H1R,H1Ri,S1R
      DOUBLE COMPLEX, DIMENSION (N1,N1) :: Gamma_L_aux
      DOUBLE COMPLEX, DIMENSION (N1,N1) :: Gamma_R_aux
      DOUBLE COMPLEX, DIMENSION (N1,N1) :: GGammaLG
      DOUBLE COMPLEX, DIMENSION (N1,N1) :: GGammaRG
      DOUBLE COMPLEX, DIMENSION (N1,N1) :: Xmataux
      DOUBLE COMPLEX, DIMENSION (N1,N1) :: Ymataux
      DOUBLE COMPLEX, DIMENSION (N1,N1) :: G0, G1
      DOUBLE COMPLEX, DIMENSION (N1,N1) :: ephcoef
      DOUBLE COMPLEX, DIMENSION (N1**2) :: work
      double precision dHdQ(maxnh_eph)
      type(matrixTypeGeneral) :: gfgeneral,gfout,gfout2
      type(matrixTypeGeneral),intent(inout) :: ematgeneral(matSpin)
      type(matrixTypeGeneral) :: hgeneral(matSpin),rhogeneral(matSpin),sgeneral
      integer  nnzrow(n1),nnz
      type(ioType) :: io
      integer gfmattype

      INTEGER :: nk,ik,iufiletr,iueph, il, jl, iuo, juo, imode, jo
      INTEGER :: length_org
      double precision kxij, ckxij, skxij

      INTEGER, DIMENSION (N1) :: IPIV
      INTEGER     indt,NM
      DOUBLE COMPLEX,PARAMETER :: zone=(1.D0,0D0)
      DOUBLE COMPLEX,PARAMETER :: zzero=(0.D0,0D0)
      DOUBLE COMPLEX,PARAMETER :: zi=(0.D0,1.D0)
      DOUBLE PRECISION :: V,ef,T

      DOUBLE COMPLEX :: ei0,ener_sigma
      double precision dsigma
      double complex weightc,signweight,weightcurrent

      INTEGER :: I,II,JJ,ISPIN,L,J,k,MyNode,ind,opindex,nt,info,ispin2

      DOUBLE PRECISION :: wk,TEinterm, pi
      double PRECISION, save :: deltaene,deltaenes,deltaenebig
      double precision kbtoltrans,fl,fr,efl,efr

      external pasbias
      integer :: nleadslr=2
      integer, allocatable :: leadsdim(:)
      logical outputwfs
      double precision, allocatable :: LeadsVoltageShift(:)

#ifdef MPI
      CALL MPI_COMM_SIZE(negf_comm,Nnodes,MPIerror)
      CALL MPI_COMM_RANK(negf_comm,mynode,MPIerror)
#else
      MyNode=0
      Nnodes =1
#endif

      gfmattype=0 !for dense GF matrix

      pi=3.141592653589793D0
      IF (ik.eq.1) THEN

        if(mynode_inverse.eq.0)then
          allocate(leadsdim(nleadslr),LeadsVoltageShift(nleadslr))
          leadsdim(1)=nl
          leadsdim(2)=nr
          LeadsVoltageShift(1)=0.5D0 * V
          LeadsVoltageShift(2)=-0.5D0 * V

          call energygrid_transm(slabel,nspin,V,T,Ef,nleadslr,leadsdim,nk, LeadsVoltageShift,deltaene,deltaenes, deltaenebig)
          deallocate(leadsdim,LeadsVoltageShift)
        endif

      ENDIF

      IF (mynode_inverse.eq.0.and.ik .EQ. 1) Then

        if(MyNode.EQ.0) THEN
          allocate(t0(ETransmGrid%nEnergiesGlobal,nspin))
          allocate(t_in(ETransmGrid%nEnergiesGlobal,nmode,nspin))
          allocate(t_inL(ETransmGrid%nEnergiesGlobal,nmode,nspin))
          allocate(t_inR(ETransmGrid%nEnergiesGlobal,nmode,nspin))
          allocate(t_eh(ETransmGrid%nEnergiesGlobal,nmode,nspin))
          allocate(t_ec(ETransmGrid%nEnergiesGlobal,nmode,nspin))
          allocate(t_ecL(ETransmGrid%nEnergiesGlobal,nmode,nspin))
          allocate(t_ecR(ETransmGrid%nEnergiesGlobal,nmode,nspin))
          allocate(t_AsymL(ETransmGrid%nEnergiesGlobal,nmode,nspin))
          allocate(t_AsymR(ETransmGrid%nEnergiesGlobal,nmode,nspin))
          allocate(t_inel(ETransmGrid%nEnergiesGlobal,nmode,nspin))
          allocate(it_inel(ETransmGrid%nEnergiesGlobal,nmode,nspin))
        else
          allocate(t0(ETransmGrid%nEnergies,nspin))
          allocate(t_in(ETransmGrid%nEnergies,nmode,nspin))
          allocate(t_inL(ETransmGrid%nEnergies,nmode,nspin))
          allocate(t_inR(ETransmGrid%nEnergies,nmode,nspin))
          allocate(t_eh(ETransmGrid%nEnergies,nmode,nspin))
          allocate(t_ec(ETransmGrid%nEnergies,nmode,nspin))
          allocate(t_ecL(ETransmGrid%nEnergies,nmode,nspin))
          allocate(t_ecR(ETransmGrid%nEnergies,nmode,nspin))
          allocate(t_AsymR(ETransmGrid%nEnergies,nmode,nspin))
          allocate(t_AsymL(ETransmGrid%nEnergies,nmode,nspin))
          allocate(t_inel(ETransmGrid%nEnergies,nmode,nspin))
          allocate(it_inel(ETransmGrid%nEnergies,nmode,nspin))
        endif
        t0=0.d0
        t_in=0.d0
        t_inL=0.d0
        t_inR=0.d0
        t_eh=0.d0
        t_ec=0.d0
        t_ecL=0.d0
        t_ecR=0.d0
        t_AsymL=0.d0
        t_AsymR=0.d0
        t_inel=0.d0
        it_inel=0.d0
      ENDIF


      if((MyNode.eq.0).and.getsigma.and.(ik .EQ. 1)) call writesigmainfo(slabel,n1,nl,nr, Nnodes,ETransmGrid%nEnergies,NSPIN,nk,ETransmGrid%e(1), tenergi,deltaenebig,deltaenes, v,ef)

      weightc=1.D0/(1.0D0 * ETransmGrid%nEnergiesGlobal)

      eneloop: DO I=1,ETransmGrid%nEnergies
        spinloop: DO ISPIN=1,ETransmGrid%nSpin

          if(mynode_inverse.eq.0)then

#ifdef MPI
            CALL MPI_BARRIER(inverseheads_comm, MPIerror)
#endif

            allocate(LeadsVoltageShift(nleadslr))
            LeadsVoltageShift(1)=0.5D0 * V
            LeadsVoltageShift(2)=-0.5D0 * V

            call get_selfenergies(i,ispin,ik,ETransmGrid,v,0,LeadsVoltageShift,sigmatodisk)

            deallocate(LeadsVoltageShift)

            nnz=n1*n1
            call AllocateMatrixGeneral(n1,n1,nnz,gfmattype,gfgeneral, "keldyshreal", io)

            ei0=ETransmGrid%e(i)+zi*deltaimagtrc
            call setgfelementsgeneral_nc(ei0,matSpin,ispin, gfgeneral,nnz,n1,nl,nr,  ETransmGrid%sigma(1,1,1,1)%sigma, ETransmGrid%sigma(2,1,1,1)%sigma, hgeneral,sgeneral)
!end allocating sparse matrix


            if(mynode_inverse.eq.0)then
              CALL ZGETRF(N1,N1,gfgeneral%matdense%a,N1,IPIV,   INFO)
              CALL ZGETRI(N1,gfgeneral%matdense%a,  N1,IPIV,WORK,N1**2,INFO)
            else
              write(12347,*)"not calculating T"
            endif

            if(mynode_inverse.eq.0)then
              if(getsigma)then
                call writesigma(slabel,ispin,mynode,ik,i,n1,nl,nr, ETransmGrid%e(i),ef,   ETransmGrid%sigma(1,1,1,1)%sigma,ETransmGrid%sigma(2,1,1,1)%sigma, gfgeneral%matdense%a)
              endif

! calculate mode independent terms

              Gamma_L_aux(:,:)=(0.d0,0.d0)
              do il=1, NL
                iuo=il
                do jl=1, NL
                  juo=jl
                  Gamma_L_aux(juo,iuo)=zi*( ETransmGrid%sigma(1,1,1,1)%sigma(jl,il)-conjg(ETransmGrid%sigma(1,1,1,1)%sigma(il,jl)) )
                enddo
              enddo
              Gamma_R_aux(:,:)=(0.d0,0.d0)
              do il=1, NR
                iuo=N1-NR+il
                do jl=1, NR
                  juo=N1-NR+jl
                  Gamma_R_aux(juo,iuo)=zi*( ETransmGrid%sigma(2,1,1,1)%sigma(jl,il)-conjg(ETransmGrid%sigma(2,1,1,1)%sigma(il,jl)) )
                enddo
              enddo

! Xmataux=G*Gamma_L
              call zgemm('N','N',N1, N1, N1, zone, gfgeneral%matdense%a, N1, Gamma_L_aux, N1, zzero, Xmataux, N1)
! GGammaLG=G*Gamma_L*G^H
              call zgemm('N','C',N1, N1, N1, zone, Xmataux, N1, gfgeneral%matdense%a, N1, zzero, GGammaLG, N1)
! Xmataux=G*Gamma_R
              call zgemm('N','N',N1, N1, N1, zone, gfgeneral%matdense%a, N1, Gamma_R_aux, N1, zzero, Xmataux, N1)
! GGammaRG=G*Gamma_R*G^H
              call zgemm('N','C',N1, N1, N1, zone, Xmataux, N1, gfgeneral%matdense%a, N1, zzero, GGammaRG, N1)
! Xmataux=G*Gamma_R*G^H*Gamma_L
              call zgemm('N','N',N1, N1, N1, zone, GGammaRG, N1, Gamma_L_aux, N1, zzero, Xmataux, N1)
! G0=G*Gamma_R*G^H*Gamma_L*G
              call zgemm('N','N',N1, N1, N1, zone, Xmataux, N1, gfgeneral%matdense%a, N1, zzero, G0, N1)

              TEinterm = 0.d0
              do iuo = 1, N1
                TEinterm = TEinterm + Xmataux(iuo,iuo)
              enddo
              t0(I,ISPIN) = t0(I,ISPIN) + nk*wk*TEinterm

! loop over modes
              modeloop: do imode = 1, nmode

! read dHdQ from the file and construct e-ph coupling matrix (Approximated SSH)
                call io_assign(iueph)
                inquire(iolength=length_org) dHdQ
                open(iueph,file=fildhdq,access='direct',form='unformatted',recl=length_org)
                if ( ispin .eq. 1 ) then
                  read(iueph,rec=mode_index(imode)) dHdQ
                elseif ( ispin .eq. 2 ) then
                  read(iueph,rec=mode_index(imode)+3*na_fc_tot) dHdQ
                endif
                call io_close(iueph)

                allocate(Afull(2,N1,N1))
                Afull(:,:,:)=0.d0
                ephcoef(:,:)=zzero
                do iuo=1, N1
                  do j=1, numh_eph(iuo)
                    ind=listhptr_eph(iuo)+j
                    jo=listh_eph(ind)
                    juo=indxuo_eph(jo)
                    kxij=kpoint(1)*xij_eph(1,ind)+kpoint(2)*xij_eph(2,ind)
                    ckxij=cos(kxij)
                    skxij=sin(kxij)
                    Afull(1,juo,iuo)=Afull(1,juo,iuo)+dHdQ(ind)*ckxij
                    Afull(2,juo,iuo)=Afull(2,juo,iuo)-dHdQ(ind)*skxij
                  enddo
                enddo
                do iuo=ihead_eph, itail_eph
                  do juo=ihead_eph, itail_eph
                    ephcoef(juo,iuo)=Afull(1,juo,iuo)+zi*Afull(2,juo,iuo)
                  enddo
                enddo
                deallocate(Afull)
                ephcoef(:,:)=ephcoef(:,:)/sqrt(2.d0*mode_freq(imode))

! T_in and T_eh
! M*G*Gamma_L*G^H*M*G*Gamma_L*G^H
                call zgemm('N','N',N1, N1, N1, zone, ephcoef, N1, GGammaLG, N1, zzero, Ymataux, N1)
                call zgemm('N','N',N1, N1, N1, zone, Ymataux, N1, ephcoef, N1, zzero, Xmataux, N1)
                call zgemm('N','N',N1, N1, N1, zone, Xmataux, N1, GGammaLG, N1, zzero, Ymataux, N1)

! T_inL
                TEinterm = 0.d0
                do iuo=1, N1
                  TEinterm=TEinterm+Ymataux(iuo,iuo)
                enddo
                t_inL(i,imode,ispin) = t_inL(i,imode,ispin) + nk*wk*TEinterm

! M*G*Gamma_R*G^H*M*G*Gamma_L*G^H
                call zgemm('N','N',N1, N1, N1, zone, ephcoef, N1, GGammaRG, N1, zzero, Ymataux, N1)
                call zgemm('N','N',N1, N1, N1, zone, Ymataux, N1, ephcoef, N1, zzero, Xmataux, N1)
                call zgemm('N','N',N1, N1, N1, zone, Xmataux, N1, GGammaLG, N1, zzero, Ymataux, N1)

! T_in
                TEinterm = 0.d0
                do iuo=1, N1
                  TEinterm=TEinterm+Ymataux(iuo,iuo)
                enddo
                t_in(i,imode,ispin) = t_in(i,imode,ispin) + nk*wk*TEinterm

! M*G*Gamma_R*G^H*M*G*Gamma_R*G^H
                call zgemm('N','N',N1, N1, N1, zone, Xmataux, N1, GGammaRG, N1, zzero, Ymataux, N1)

! T_inR
                TEinterm = 0.d0
                do iuo=1, N1
                  TEinterm=TEinterm+Ymataux(iuo,iuo)
                enddo
                t_inR(i,imode,ispin) = t_inR(i,imode,ispin) + nk*wk*TEinterm

! M*ImG*M*ImG
                Ymataux(:,:)=imag(gfgeneral%matdense%a(:,:))
                call zgemm('N','N',N1, N1, N1, zone, ephcoef, N1, Ymataux, N1, zzero, Xmataux, N1)
                call zgemm('N','N',N1, N1, N1, zone, Xmataux, N1, Xmataux, N1, zzero, Ymataux, N1)

! T_eh
                TEinterm = 0.d0
                do iuo=1, N1
                  TEinterm=TEinterm+Ymataux(iuo,iuo)
                enddo
                t_eh(i,imode,ispin) = t_eh(i,imode,ispin) + nk*wk*TEinterm

! === T^ec = 2*Re Tr[M*G*M*G*Gamma_R*G^H*Gamma_L*G] ===
! G1=M*G*M
                call zgemm('N','N',N1, N1, N1, zone, ephcoef, N1, gfgeneral%matdense%a, N1, zzero, Xmataux, N1)
                call zgemm('N','N',N1, N1, N1, zone, Xmataux, N1, ephcoef, N1, zzero, G1, N1)
! T^ec
                call zgemm('N','N',N1, N1, N1, zone, G1, N1, G0, N1, zzero, Xmataux, N1)
                TEinterm = 0.d0
                do iuo=1, N1
                  TEinterm=TEinterm+2.d0*real(Xmataux(iuo,iuo))
                enddo
                t_ec(i,imode,ispin) = t_ec(i,imode,ispin) + nk*wk*TEinterm

! === T^ecL    = Im Tr[M*G*Gamma_L*G^H*M*G*Gamma_R*G^H*Gamma_L*G]===
! === T^AsymL  = Re Tr[M*G*Gamma_L*G^H*M*G*Gamma_R*G^H*Gamma_L*G]===
! G1=M*G*Gamma_L*G^H*M
                call zgemm('N','N',N1, N1, N1, zone, ephcoef, N1, GGammaLG, N1, zzero, Xmataux, N1)
                call zgemm('N','N',N1, N1, N1, zone, Xmataux, N1, ephcoef, N1, zzero, G1,N1)

! T^ecL
                call zgemm('N','N',N1, N1, N1, zone, G1, N1, G0, N1, zzero, Xmataux, N1)
                TEinterm = 0.d0
                do iuo=1, N1
                  TEinterm=TEinterm+imag(Xmataux(iuo,iuo))
                enddo
                t_ecL(i,imode,ispin) = t_ecL(i,imode,ispin) + nk*wk*TEinterm
! T^AsymL
                TEinterm = 0.d0
                do iuo=1, N1
                  TEinterm=TEinterm+real(Xmataux(iuo,iuo))
                enddo
                t_AsymL(i,imode,ispin) = t_AsymL(i,imode,ispin) + nk*wk*TEinterm

! === T^ecR    = Im Tr[M*G*Gamma_R*G^H*M*G*Gamma_R*G^H*Gamma_L*G]===
! === T^AsymR  = Re Tr[M*G*Gamma_R*G^H*M*G*Gamma_R*G^H*Gamma_L*G]===
! G1=M*G*Gamma_R*G^H*M
                call zgemm('N','N',N1, N1, N1, zone, ephcoef, N1, GGammaRG, N1, zzero, Xmataux, N1)
                call zgemm('N','N',N1, N1, N1, zone, Xmataux, N1, ephcoef, N1, zzero, G1, N1)
! T^ecR
                call zgemm('N','N',N1, N1, N1, zone, G1, N1, G0, N1, zzero, Xmataux, N1)
                TEinterm = 0.d0
                do iuo=1, N1
                  TEinterm=TEinterm+imag(Xmataux(iuo,iuo))
                enddo
                t_ecR(i,imode,ispin) = t_ecR(i,imode,ispin) + nk*wk*TEinterm

! T^AsymR
                TEinterm = 0.d0
                do iuo=1, N1
                  TEinterm=TEinterm+real(Xmataux(iuo,iuo))
                enddo
                t_AsymR(i,imode,ispin) = t_AsymR(i,imode,ispin) + nk*wk*TEinterm

! == Inelastic term:  Tr [ M*G*Gamma_R*G^H*M*G^H*Gamma_L*G ] ===

! G0=G^H*Gamma_L*G
                call zgemm('C','N',N1, N1, N1, zone, gfgeneral%matdense%a, N1, Gamma_L_aux, N1, zzero, Xmataux, N1)
                call zgemm('N','N',N1, N1, N1, zone, Xmataux, N1, gfgeneral%matdense%a, N1, zzero, Ymataux, N1)

! G1=M*G*GammaR*G^H*M
                call zgemm('N','N',N1, N1, N1, zone, ephcoef, N1, GGammaRG, N1, zzero, Xmataux,N1)
                call zgemm('N','N',N1, N1, N1, zone, Xmataux, N1, ephcoef, N1, zzero, G1, N1)
! T_inel
                call zgemm('N','N',N1, N1, N1, zone, G1, N1, Ymataux, N1, zzero, Xmataux, N1)
                TEinterm = 0.d0
                do iuo=1, N1
                  TEinterm=TEinterm+Xmataux(iuo,iuo)
                enddo
                t_inel(i,imode,ispin) = t_inel(i,imode,ispin) + nk*wk*TEinterm
                TEinterm = 0.d0
                do iuo=1, N1
                  TEinterm=TEinterm+imag(Xmataux(iuo,iuo))
                enddo
                it_inel(i,imode,ispin) = it_inel(i,imode,ispin) + nk*wk*TEinterm

              enddo modeloop

              call DestroyMatrixGeneral(gfgeneral,"keldyshreal",io)
            ENDIF

          endif

          call deallocate_selfenergies(i,ispin,ik,ETransmGrid)

        ENDDO spinloop
      ENDDO eneloop

      if(mynode_inverse.ne.0)return

      nt=ETransmGrid%nEnergies
#ifdef MPI
      if (myhead.NE.0) THEN
        do i=1,nspin
          CALL MPI_SEND(t0(1,i),nt,DAT_double,0,6, inverseheads_comm,MPIerror)
          do imode = 1, nmode
            CALL MPI_SEND(t_in(1,imode,i),nt,DAT_double,0,6,inverseheads_comm,MPIerror)
            CALL MPI_SEND(t_inL(1,imode,i),nt,DAT_double,0,6,inverseheads_comm,MPIerror)
            CALL MPI_SEND(t_inR(1,imode,i),nt,DAT_double,0,6,inverseheads_comm,MPIerror)
            CALL MPI_SEND(t_eh(1,imode,i),nt,DAT_double,0,6,inverseheads_comm,MPIerror)
            CALL MPI_SEND(t_ec(1,imode,i),nt,DAT_double,0,6,inverseheads_comm,MPIerror)
            CALL MPI_SEND(t_ecL(1,imode,i),nt,DAT_double,0,6,inverseheads_comm,MPIerror)
            CALL MPI_SEND(t_ecR(1,imode,i),nt,DAT_double,0,6,inverseheads_comm,MPIerror)
            CALL MPI_SEND(t_AsymL(1,imode,i),nt,DAT_double,0,6,inverseheads_comm,MPIerror)
            CALL MPI_SEND(t_AsymR(1,imode,i),nt,DAT_double,0,6,inverseheads_comm,MPIerror)
            CALL MPI_SEND(t_inel(1,imode,i),nt,DAT_double,0,6,inverseheads_comm,MPIerror)
            CALL MPI_SEND(it_inel(1,imode,i),nt,DAT_double,0,6,inverseheads_comm,MPIerror)
          enddo
        enddo
      else
        DO I=1,nheads-1
          do j=1,nspin
            CALL MPI_RECV(t0(I*nt+1,j), nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
            do imode = 1, nmode
              CALL MPI_RECV(t_in(I*nt+1,imode,j),nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
              CALL MPI_RECV(t_inL(I*nt+1,imode,j),nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
              CALL MPI_RECV(t_inR(I*nt+1,imode,j),nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
              CALL MPI_RECV(t_eh(I*nt+1,imode,j),nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
              CALL MPI_RECV(t_ec(I*nt+1,imode,j),nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
              CALL MPI_RECV(t_ecL(I*nt+1,imode,j),nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
              CALL MPI_RECV(t_ecR(I*nt+1,imode,j),nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
              CALL MPI_RECV(t_AsymL(I*nt+1,imode,j),nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
              CALL MPI_RECV(t_AsymR(I*nt+1,imode,j),nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
              CALL MPI_RECV(t_inel(I*nt+1,imode,j),nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
              CALL MPI_RECV(it_inel(I*nt+1,imode,j),nt,DAT_double,I,6,inverseheads_comm,ISTATUS, MPIerror)
            enddo
          enddo
        ENDDO
      endif
#endif


      IF (ik .EQ. nk) Then
        IF (myhead .EQ. 0) THEN

          call io_assign(iufiletr)
          write( chivv, '(i7 )' ) iv
          chivv = pasbias(chivv, '.')
          filetr = paste(slabel,'.INELTRC')
          filetr = paste( chivv, filetr )
          OPEN(UNIT=iufiletr,FILE=filetr,status='unknown',RECL=100000)

          WRITE(iufiletr,'(a6,f11.4,a14,i4)') '# V = ',V, '    k-points: ',nk
          WRITE(iufiletr,'(a18)',ADVANCE='NO') '# energy T_total, '

          DO ISPIN=1,NSPIN
            WRITE(iufiletr,'(a9,i1,a3)',ADVANCE='NO')'T_{ispin=',ispin,'} ,'
          ENDDO

          WRITE(iufiletr,*)

          DO I=1,1 ! ETransmGrid%nEnergiesGlobal
            indt=MOD(i-1,nheads) * nt+ (i-1)/nheads + 1

! T^bar = T(Ef)
            if ( i .eq. 1 ) then
              do ispin = 1, nspin
                do imode = 1, nmode
                  t_in_vib(imode,ispin,iv+1)  = t_in(indt,imode,ispin)/(1.d0*nk)
                  t_inL_vib(imode,ispin,iv+1) = t_inL(indt,imode,ispin)/(1.d0*nk)
                  t_inR_vib(imode,ispin,iv+1) = t_inR(indt,imode,ispin)/(1.d0*nk)
                  t_eh_vib(imode,ispin,iv+1)  = t_eh(indt,imode,ispin)/(1.d0*nk)
                  t_ec_Vb(imode,ispin,iv+1)   = t_ec(indt,imode,ispin)/(1.d0*nk)
                  t_ecL_Vb(imode,ispin,iv+1)  = t_ecL(indt,imode,ispin)/(1.d0*nk)
                  t_ecR_Vb(imode,ispin,iv+1)  = t_ecR(indt,imode,ispin)/(1.d0*nk)
                  t_AsymL_Vb(imode,ispin,iv+1)= t_AsymL(indt,imode,ispin)/(1.d0*nk)
                  t_AsymR_Vb(imode,ispin,iv+1)= t_AsymR(indt,imode,ispin)/(1.d0*nk)
                  t_inel_Vb(imode,ispin,iv+1) = t_inel(indt,imode,ispin)/(1.d0*nk)
                enddo
              enddo
            endif

            WRITE(iufiletr,'(2e16.7,A3)',ADVANCE='NO')(DREAL(ETransmGrid%eGlobal(i))-ef)* 13.6056981D0,SUM(t0(indt,:))/(1.D0 * nk), '   '
            DO ISPIN=1,NSPIN
              WRITE(iufiletr,'(e16.7,A3)',ADVANCE='NO')t0(indt,ISPIN)/(1.D0 * nk),'   '
              WRITE(iufiletr,*)
              do imode = 1, nmode
                WRITE(iufiletr,'(i16,A3)',ADVANCE='NO')imode,'   '
                WRITE(iufiletr,'(e16.7,A3)',ADVANCE='NO')t_in(indt,imode,ISPIN)/(1.D0 * nk),'   '
                WRITE(iufiletr,'(e16.7,A3)',ADVANCE='NO')t_inL(indt,imode,ISPIN)/(1.D0 * nk),'   '
                WRITE(iufiletr,'(e16.7,A3)',ADVANCE='NO')t_inR(indt,imode,ISPIN)/(1.D0 * nk),'   '
                WRITE(iufiletr,'(e16.7,A3)',ADVANCE='NO')t_eh(indt,imode,ISPIN)/(1.D0 * nk),'   '
                WRITE(iufiletr,'(e16.7,A3)',ADVANCE='NO')t_ec(indt,imode,ISPIN)/(1.D0 * nk),'   '
                WRITE(iufiletr,'(e16.7,A3)',ADVANCE='NO')t_ecL(indt,imode,ISPIN)/(1.D0 * nk),'   '
                WRITE(iufiletr,'(e16.7,A3)',ADVANCE='NO')t_ecR(indt,imode,ISPIN)/(1.D0 * nk),'   '
                WRITE(iufiletr,'(e16.7,A3)',ADVANCE='NO')t_AsymL(indt,imode,ISPIN)/(1.D0 * nk),'   '
                WRITE(iufiletr,'(e16.7,A3)',ADVANCE='NO')t_AsymR(indt,imode,ISPIN)/(1.D0 * nk),'   '
                WRITE(iufiletr,'(e16.7,A3)',ADVANCE='NO')t_inel(indt,imode,ISPIN)/(1.D0 * nk),'   '
                WRITE(iufiletr,*)
              enddo
            ENDDO
            WRITE(iufiletr,*)
          ENDDO
          WRITE(iufiletr,*)
          call io_close(iufiletr)
        ENDIF

        deallocate(t0)
        deallocate(t_in)
        deallocate(t_inL)
        deallocate(t_inR)
        deallocate(t_eh)
        deallocate(t_ec)
        deallocate(t_ecL)
        deallocate(t_ecR)
        deallocate(t_AsymL)
        deallocate(t_AsymR)
        deallocate(t_inel)
        deallocate(it_inel)

        call deallocate_energygrid_selfenergies2(ETransmGrid)
      ENDIF

      END SUBROUTINE

      subroutine loe_wrap(ik,nk,N1,nl,nr,NspinComplexMatrix,V,iv,&
         slabel,kpoint,wk, Ef_Lead,kb,t,&
         hgeneral,sgeneral,rhogeneral,ematgeneral)

      use sigma
      use mTypes

      implicit none
      integer ik, nk,n1,nl,nr,NspinComplexMatrix,iv,NspinLoops
      double precision kbt,v,wk,Ef_Lead,kb,t
      DOUBLE PRECISION, DIMENSION (3) :: kpoint
      CHARACTER(LEN=*) :: slabel
      type(matrixTypeGeneral) :: hgeneral(NspinComplexMatrix),sgeneral, rhogeneral(NspinComplexMatrix),ematgeneral(NspinComplexMatrix)

      kbt=kb*T

      NspinLoops=NspinComplexMatrix
      if(NspinComplexMatrix.gt.2)NspinLoops=1

      CALL loe(ik,nk,N1,nl,nr,NspinLoops,V,iv,&
          NspinComplexMatrix,slabel,kpoint(1:2),wk,&
          H0_L,H0_R,S0_L,S0_R,S1_L,S1_R,H1_L,H1_R,&
          Ef_Lead,kbt,hgeneral,sgeneral,rhogeneral,ematgeneral)

      end subroutine loe_wrap

