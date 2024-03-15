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
!                   KELDYSHREAL,
!                   KELDYSHIMAG  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!

      subroutine keldyshReal(N1,NSpinBlocks,NspinComplexMatrix,NL,NR,nk,
     &    ik,V,T,Delta,ef,rhobs_general,ematbs_general,
     &    hgeneralp,sgeneralp,rhogeneralp,ematgeneralp,
     &    iter,ldos,weightrho)

! *****************************************************************
! Calculates the Keldysh Green's function on the real axis
!
! Written by Alexandre Reily Rocha and Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: rochaa@tcd.ie, reilya@if.usp.br, runggeri@tcd.ie
! ********** HISTORY **********************************************
! Original version:      October 2007
! ********** INPUT ************************************************

      use precision
      use negfmod, only : bs_add, emtimings, emforces, rhobs_delta, 
     .    bs_writetrc, wkmod, bsskip, deltaimag, negfon, itermod,
     .    nbss, bssc, inversion_solver, nprocs_hs, geterrsigma,
     .    bseskip, bs_method
      use mTypes
      use mMatrixUtil
      use mONInterface
      use mBoundstates
      use mMPI_NEGF
      use sigma
      use mONBoundStates
      use mEnergyGrid
      use mSigmaMethod1, only: check_error_sigma
      use mMatrixUtilOMP
      use ScissorOperator, only : SCO_istart, SCO_nob, SCO_Hblock,
     .    SCOSetHamiltonianBlock
      
      implicit none
      
      include "const2.h"
      
      integer, intent(in) :: NSpinBlocks,NspinComplexMatrix
      double precision, intent(in) :: weightrho
      integer :: N1,NL,NR,nk,ik,II,JJ,ISPIN,is,
     &    INFO,I,inde,indemax,nw3,nw2,iter
        
      integer, dimension (:), allocatable :: IPIV
      
      double precision :: fermi_aux,Ei,V,T,ef,Delta
      
      double complex :: fR, fL,eicplx
        
      double complex, dimension (:,:), allocatable ::
     &    Tau1_aux
        
      DOUBLE COMPLEX,ALLOCATABLE ::  eiene(:,:,:)
      DOUBLE PRECISION,ALLOCATABLE ::  ene(:)
      DOUBLE COMPLEX, ALLOCATABLE, SAVE :: rhobstot(:,:,:,:)
      DOUBLE COMPLEX, ALLOCATABLE, DIMENSION (:,:) :: sdense_inv
      DOUBLE COMPLEX,ALLOCATABLE   ::  zv(:),
     &    WORK3(:),veigl(:,:),veigr(:,:)
      DOUBLE PRECISION,ALLOCATABLE ::  RWORK2(:)
      DOUBLE COMPLEX, ALLOCATABLE :: GF_iter1r(:,:),gammar(:,:),
     &    sigmar(:,:)
      DOUBLE COMPLEX, ALLOCATABLE :: GF_iter1l(:,:),gammal(:,:),
     &    sigmal(:,:)
      DOUBLE COMPLEX, ALLOCATABLE :: gf1(:,:),gf2(:,:)
      DOUBLE COMPLEX, ALLOCATABLE :: GF_iter2l(:,:),GF_iter2r(:,:),
     &    GF_iter3(:,:),gf_iter3b(:,:)
      double precision, allocatable :: ef_bssk(:,:)
      double precision, allocatable, save :: ef_bss(:,:)
      integer nleads,i1,i2,ind,j,il
      double complex, allocatable :: hdense(:,:,:),sdense(:,:)
      integer, allocatable :: nebss(:,:)
      double precision, allocatable :: deltabss(:)
      double precision, allocatable :: tij(:,:,:,:)
      type(matrixTypeGeneral) :: gfgeneral,gfP,
     &    rhobs_general(NspinComplexMatrix),
     &    ematbs_general(NspinComplexMatrix),gfout
      type(matrixTypeGeneral) :: hgeneralp(NspinComplexMatrix),
     &    sgeneralp,rhogeneralp(NspinComplexMatrix),
     &    ematgeneralp(NspinComplexMatrix)
      double complex, allocatable :: b(:)
      double complex, allocatable :: bnc(:,:)
      type(matrixTypeGeneral), allocatable :: sigmamp(:)
      type(matrixTypeGeneral), allocatable :: gammamp(:)
      type(ioType) :: io
      integer gfmattype,sizesigma
      integer  nnzrow(n1),nnz
      double complex gfadd
      integer*4:: sc_0,sc_1,sc_r,sc_m,scb_0,scb_1
      integer*4:: sca_0,sca_1
      INTEGER :: MPIerror
      double precision dsigma
      double complex ener_sigma
      logical ldos
      type(SelfEnergyType), allocatable :: sigmaleads(:)
      integer iend

      io%isDebug=.false.
      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
        CALL SYSTEM_CLOCK(sca_0,sc_r,sc_m)
      endif
     

      if(.not.negfon)then
        gfmattype=0 !for dense GF matrix
      else
        gfmattype=2 !for sparse GF matrix
      endif

      if(mynode_inverse.eq.0)then
        ALLOCATE(zv(N1))
        nw3=10*n1
        ALLOCATE(WORK3(nw3))
        nw2=2*n1
        ALLOCATE(RWORK2(nw2))
        ALLOCATE(GF_iter1r(N1,NR),gammar(NR,NR),sigmar(NR,NR))
        ALLOCATE(GF_iter1l(N1,NL),gammal(NL,NL),sigmal(NL,NL))
        ALLOCATE(GF_iter2l(N1,NL))
        ALLOCATE(GF_iter2r(N1,NR))
        ALLOCATE(GF_iter3(N1,NL+NR))
        ALLOCATE(GF_iter3b(N1,NL+NR))
        ALLOCATE(veigl(1,1),veigr(1,1))
        allocate(IPIV(N1))

        if(mynode_inverse.eq.0.and.bs_add.and.bs_method.eq.1)then
          ALLOCATE(Tau1_aux(N1,N1))
        endif

        if(mynode_inverse.eq.0.and.bs_add.and.bs_method.eq.0)then
          ALLOCATE(Tau1_aux(N1,N1))
          allocate(hdense(n1,n1,NspinComplexMatrix),sdense(n1,n1))
          hdense=0D0
          sdense=0D0

          do ii=1,n1
            do ind=sgeneralp%matSparse%q(ii),
     .          sgeneralp%matSparse%q(ii+1)-1
              do ispin=1,NspinComplexMatrix
                hdense(ii,hgeneralp(ispin)%matSparse%j(ind),ispin)=
     .              hgeneralp(ispin)%matSparse%b(ind)
              enddo
              sdense(ii,sgeneralp%matSparse%j(ind))=
     .            sgeneralp%matSparse%b(ind)
            enddo
          enddo

          if(SCOSetHamiltonianBlock)then
            write(12347,*)"Bound states: Applying scissor
     . operator to Hamiltonian",  SCO_istart,SCO_istart+SCO_nob-1,ei
            iend=SCO_istart+SCO_nob-1
            do ispin=1,NspinComplexMatrix !xxx : change for spiorbit or non-collinear spins
!xxx : SCO_Hblock is allocated with nspinrealinputmatrix, not NspinComplexMatrix
              hdense(SCO_istart:iend,SCO_istart:iend,ispin)=
     .            SCO_Hblock(:,:,ispin)
            enddo
          endif

        endif
                

        if(mynode_inverse.eq.0.and.(bs_add).and.(bs_method.eq.0)) then
          ALLOCATE(eiene(ERealGrid%nEnergies,N1,NspinComplexMatrix),
     &        ene(ERealGrid%nEnergiesGlobal))

          if(bsskip.gt.1)then
            IF ((MOD(itermod-1,bsskip) .EQ. 0) .AND. (ik .EQ. 1))
     .        then
              if(myhead.eq.0) write(*,*)"setting rhobstot=",itermod
              if(.not.ALLOCATED(rhobstot))
     &            ALLOCATE(rhobstot(N1,N1,NspinComplexMatrix,nk))
              rhobstot=0D0
            endif
          endif
 
          if(MOD(itermod-1,bsskip) .EQ. 0) then
            if(.not.ALLOCATED(sdense_inv))ALLOCATE(sdense_inv(N1,N1))
            sdense_inv=sdense
            CALL ZGETRF(N1,N1,sdense_inv,N1,IPIV,INFO)
            CALL ZGETRI(N1,sdense_inv,N1,IPIV,Tau1_aux,N1**2,INFO)
          endif
        endif



        if(mynode_inverse.eq.0.and.(bs_add).and.(bs_method.eq.1)) then

          nleads=nbss+2
          if(.not.allocated(ef_bss))then
            allocate(ef_bss(nleads,NspinComplexMatrix))
            ef_bss=ef
            ef_bss(nbss+1,:)=ef+e*V/2.D0
            ef_bss(nleads,:)=ef-e*V/2.D0
            ef_bss(1,:)=ef_bss(nbss+1,:)
            ef_bss(nbss,:)=ef_bss(nleads,:)
            call read_ef_bss(ef_bss,nleads,NspinComplexMatrix,ef)
            ef_bss(nbss+1,:)=ef+e*V/2.D0
            ef_bss(nleads,:)=ef-e*V/2.D0
!            do i1=1,nleads
!              if(nspin.eq.1)then
!                write(12347,*)"ef_ini=",i1,ef_bss(i1,1),myhead
!              else
!                write(12347,*)"ef_ini=",i1,ef_bss(i1,1),ef_bss(i1,2),
!     .              myhead
!              endif
!            enddo
          endif
          allocate(ef_bssk(nleads,NspinComplexMatrix))
          ef_bssk=ef_bss

          allocate(nebss(nleads,2),deltabss(nleads))
          allocate(tij(ERealGrid%nEnergies,nleads,nleads,
     &        NspinComplexMatrix))

          allocate(gammamp(nleads),sigmamp(nleads))

          call set_gammamp(n1,nl,nr,nleads,nbss,nebss,deltabss,
     .        gammamp,sigmamp,sgeneralp)

          if(bssc.eq.0)then
            call transmij_bs(ef_bssk,NspinComplexMatrix,
     .          ERealGrid%nEnergies, V,n1, delta,ik, 
     .          nl,nr,nleads,nbss,ef,kb * T,
     .          gammamp,sigmamp,nebss,deltabss,
     .          bs_writetrc,tij,hgeneralp,sgeneralp)
            if(.false.)then
               call find_ef_bss(tij,NspinComplexMatrix,nleads,nbss,
     .             ERealGrid%nEnergies,ef_bssk
     .             ,ef, kb * t,v)
            endif
          endif
        

        endif


        if(mynode_inverse.eq.0.and.bs_add.and.bs_method.eq.1)then
          nnz=n1*n1
          call AllocateMatrixGeneral(n1,n1,nnz,0,gfgeneral,
     .       "keldyshreal", io)
          
        endif

        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
          write(12347,'(A,f12.6)')
     $         'kr_uptoeneloop',(sc_1-sc_0)*1.0d0/sc_r
          CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
        endif
      endif

      if(emtimings)then
        CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
        write(12347,'(A,f12.6)')
     $       'kr_time1',(scb_1-sca_0)*1.0d0/sc_r
        CALL SYSTEM_CLOCK(sca_0,sc_r,sc_m)
      endif

  
      DO ISPIN=1,ERealGrid%nSpin

        inde=0
        DO I=1,ERealGrid%nEnergies
!          write(12347,*)"i_ene=",i,ERealGrid%nEnergies,ispin,
!     .         mynode_inverse
          if(emtimings)then
            CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
            write(12347,'(A,f12.6)')
     $           'kr_eneloop',(sc_1-sc_0)*1.0d0/sc_r
            CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
          endif

          Ei=ERealGrid%e(i)
      
          fermi_aux=(Ei-ef+e*V/2.D0)/(kB*T)
          IF ( fermi_aux .GT. 50.D0) THEN
            fR=0.D0
          ELSE
            fR=1.D0/(1+DEXP(fermi_aux))
          ENDIF
      
          fermi_aux=(Ei-ef-e*V/2.D0)/(kB*T)
          IF ( fermi_aux .GT. 50.D0) THEN
            fL=0.D0
          ELSE
            fL=1.D0/(1.D0+DEXP(fermi_aux))
          ENDIF
   
        
          if(mynode_inverse.eq.0.and.(bs_add).and.(bs_method.eq.0))
     &        then
            nnz=n1*n1
            call AllocateMatrixGeneral(n1,n1,nnz,0,gfgeneral,
     .         "keldyshreal", io)
            call bs_m0_eigenvalues
     &          (n1,NspinComplexMatrix,i,bseskip,bsskip,ispin,ik,inde,
     &          indemax,gfgeneral%matdense%a,hdense,sdense_inv,nl,nr,
     &          work3,  nw3,rwork2,
     &          nw2,ERealGrid%nEnergies,ERealGrid%nEnergiesGlobal,
     .          eiene,ene,ei)
            call DestroyMatrixGeneral(gfgeneral,"keldyshreal",io)
          endif

! *****************************************************************





!  Calculate the Keldysh Greens function
          

!          CALL TIMER('gfinv',1)
          if(emtimings)then
            CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
          endif

          allocate(sigmaleads(ERealGrid%nLeads))
          do il=1,ERealGrid%nLeads
            call allocate_sigma_single(sigmaleads(il),
     &          ERealGrid%leadsTotalDim(il),myhead,
     &          ERealGrid%InfoSigma,
     &          (0.0D0,0.0D0),0.0D0,0.0D0)
          enddo

          if(mynode_inverse.eq.0)then
            call extractsigma(ERealGrid,sigmaleads,i,ispin,ik)
          endif

          if(nnodes_inverse.gt.1)then
            do il=1,ERealGrid%nLeads
              sizeSigma=ERealGrid%leadsTotalDim(il)**2
#ifdef MPI
              call MPI_Bcast(sigmaleads(il)%sigma(1,1),sizeSigma,
     .            DAT_dcomplex,0, inverse_comm,MPIerror)
#endif
            enddo
          endif

          if(iter.eq.1.and.geterrsigma)then
            write(12347,*)"dsigmainfo=",ERealGrid%e(I),ispin,i,v
            ener_sigma=ERealGrid%e(I)-v*0.5D0+
     &          (0.0D0,1.0D0)*deltaimag
            call check_error_sigma(dsigma,.true.,ener_sigma,'L',nl,
     .          h0_L(:,:,ispin),H1_L(:,:,ispin),S0_L,S1_L,
     &          sigmaleads(1)%sigma)
            ener_sigma=ERealGrid%e(I)+v*0.5D0+
     &          (0.0D0,1.0D0)*deltaimag
            call check_error_sigma(dsigma,.true.,ener_sigma,'R',nr,
     .          h0_r(:,:,ispin),H1_r(:,:,ispin),S0_r,S1_r,
     &          sigmaleads(2)%sigma)
          endif


          if(bs_add.and.bs_method.eq.1)then

            call set_gammamp_lr(n1,nl,nr,nleads,nbss,nebss,
     &          gammamp,sigmamp,
     &          sigmaleads(1)%sigma,sigmaleads(2)%sigma)         
            
            eicplx=ei
!            eicplx=ei+zi*delta !remove, just for test
            
!            call setgfelementsdense_bs(eicplx,ispin,gfgeneral,
!     &          nleads,sigmamp,nebss,hgeneralp(ispin),sgeneralp)
            call setgfelementsgeneral_nc_bs(eicplx,NspinComplexMatrix,
     &            ispin,
     &            gfgeneral,nnz,n1,nl,nr,nbss,nleads,sigmamp,nebss,
     &            hgeneralp,sgeneralp)

            do i1=nleads-1,nleads
              call DestroyMatrixGeneral(gammamp(i1),"keldyshreal",io)
              call DestroyMatrixGeneral(sigmamp(i1),"keldyshreal",io)
            enddo

          else

!start allocating sparse matrix
         
            if(mynode_inverse==0)then

              sgeneralp%mattype=2
              hgeneralp(ispin)%mattype=2
              if(gfmattype.eq.0)then
                nnz=n1*n1
              elseif(gfmattype.eq.2)then
                call findnnzgf2(nnz,nnzrow,n1,n1,nl,nr,
     &              sigmaleads(1)%sigma,
     &              sigmaleads(2)%sigma,
     .              hgeneralp(ispin))
              endif

              call AllocateMatrixGeneral(n1,n1,nnz,gfmattype,gfgeneral,
     .            "keldyshreal", io)

              eicplx=Ei+zi*delta
!              write(12347,*)"ereal=",dreal(eicplx),dimag(eicplx),1
              call setgfelementsgeneral_nc(eicplx,NspinComplexMatrix,
     &            ispin, gfgeneral,nnz,n1,nl,nr,
     &            sigmaleads(1)%sigma,
     &            sigmaleads(2)%sigma,
     &            hgeneralp,sgeneralp)
            endif
 

            if(nprocs_hs.ne.1)then
              sgeneralp%mattype=3
              hgeneralp(ispin)%mattype=3
              call findnnzgf2P(nnz,nnzrow,
     &            hgeneralp(ispin)%matSparseP%matSparse%iRows,n1,nl,nr,
     &            sigmaleads(1)%sigma, sigmaleads(2)%sigma,
     &            hgeneralp(ispin))

              call AllocateMatrixGeneral(-1,inverse_comm,
     .            nnodes_inverse,mynode_inverse,n1,n1,
     .            hgeneralp(ispin)%matSparseP%matSparse%iRows,n1,1,
     .            hgeneralp(ispin)%matSparseP%matSparse%iVert,nnz,
     .            3,gfP,"gf_keldyshreal",io)

              eicplx=Ei+zi*delta
              call setgfelementsgeneral_nc(eicplx,NspinComplexMatrix,
     &            ispin, gfP,nnz,n1,nl,nr,
     &            sigmaleads(1)%sigma, sigmaleads(2)%sigma,
     &            hgeneralp,sgeneralp)
              sgeneralp%mattype=2
              hgeneralp(ispin)%mattype=2

              call MaxDifferenceMatrixCRS_CRSP(gfgeneral,gfP,
     .            mynode_inverse,"(gfer) ",io)
            endif



!end allocating sparse matrix


          endif

          if(emtimings)then
            CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
            write(12347,'(A,f12.6)')
     $           'kr_setgf',(scb_1-scb_0)*1.0d0/sc_r
          endif
! *****************************************************************
!  Calculate the Keldysh Greens function

          if(emtimings)then
            CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
          endif

          if(mynode_inverse.eq.0)then
            if(bs_add.and.bs_method.eq.1)then
              CALL ZGETRF(N1,N1,gfgeneral%matdense%a,N1,IPIV,INFO)
              CALL ZGETRI(N1,gfgeneral%matdense%a,N1,IPIV,Tau1_aux,
     .            N1**2,INFO)
              if(NspinComplexMatrix<=2)then
                GF_iter3(:,1:NL)=gfgeneral%matdense%a(:,1:NL)
                GF_iter3(:,NL+1:NL+NR)=
     .               gfgeneral%matdense%a(:,N1-NR+1:N1)
              else
                GF_iter3(:,1:NL/2)=gfgeneral%matdense%a(:,1:NL/2)
                GF_iter3(:,nl/2+1:NL)=
     .              gfgeneral%matdense%a(:,n1/2+1:n1/2+NL/2)
                GF_iter3(:,NL+1:NL+NR/2)=
     .                 gfgeneral%matdense%a(:,N1/2-NR/2+1:N1/2)
                GF_iter3(:,NL+nr/2+1:NL+NR)=
     .                 gfgeneral%matdense%a(:,N1-NR/2+1:N1)
              endif

              GF_iter1l=GF_iter3(:,1:NL)
              GF_iter1r=GF_iter3(:,NL+1:NL+NR)
            else

              if(.not.negfon)then
                GF_iter1l=0D0
                GF_iter1r=0D0

                if(NspinComplexMatrix<=2)then
                  do ii=1,NL
                    GF_iter1l(ii,ii)=1D0
                  enddo
                  do ii=1,NR
                    GF_iter1r(N1-NR+ii,ii)=1D0
                  enddo
                else
                  do ii=1,NL/2
                    GF_iter1l(ii,ii)=1D0
                  enddo
                  do ii=1,NL/2
                    GF_iter1l(n1/2+ii,nl/2+ii)=1D0
                  enddo

                  do ii=1,NR/2
                    GF_iter1r(N1/2-NR/2+ii,ii)=1D0
                  enddo
                  do ii=1,NR/2
                    GF_iter1r(N1-NR/2+ii,nr/2+ii)=1D0
                  enddo
                endif

                GF_iter3(:,1:NL)=GF_iter1l
                GF_iter3(:,NL+1:NL+NR)=GF_iter1r

                call ZGESV(N1,NL+NR,gfgeneral%matdense%a,
     .              N1,IPIV,GF_iter3,N1,INFO)

                GF_iter1l=GF_iter3(:,1:NL)
                GF_iter1r=GF_iter3(:,NL+1:NL+NR)
              else

                call AllocateMatrixGeneral(n1,nl+nr,n1*(nl+nr),0,gfout,
     .              "keldyshreal", io)

                call InvertONGeneral(N1,gfgeneral,nl,nr,gfout,2,
     .              inversion_solver)

                gf_iter1l=gfout%matdense%a(:,1:nl)
                gf_iter1r=gfout%matdense%a(:,nl+1:nl+nr)

                call DestroyMatrixGeneral(gfout,"keldyshreal",io)
     
              endif
              call DestroyMatrixGeneral(gfgeneral,"keldyshreal",io)

            endif
            if(emtimings)then
              CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')
     $             'kr_gfinv',(scb_1-scb_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
            endif

          else
            write(12347,*)"not inverting real"
          endif
          if(nprocs_hs.ne.1)
     &      call DestroyMatrixGeneral(gfP,"keldyshreal",io)
!          CALL TIMER('gfinv',2)
          
!          CALL TIMER('bs_m1',1)
          if(mynode_inverse.eq.0)then
            if((bs_add).and.(bs_method.eq.1))then


!              do i1=1,nleads
!                if(NspinComplexMatrix.eq.1)then
!                  write(12347,*)"ef_bssuse=",i1,13.6057D0 * 
!     .                (ef_bss(i1,1)-ef)
!                else
!                  write(12347,*)"ef_bssuse=",i1,13.6057D0 * 
!     .                (ef_bss(i1,1)-ef),13.6057D0 * 
!     .                (ef_bss(i1,2)-ef)
!                endif
!              enddo


              if(NspinComplexMatrix.le.2)then
              call BS_M1(N1,gfgeneral%matdense%a,
     &            rhobs_general(ISPIN),
     &            ematbs_general(ISPIN),
     &            Delta,ERealGrid%w(I),fl,fr,
     &            bs_add,bs_method, NspinComplexMatrix,ispin,ei,
     &            gammamp,sigmamp,
     &            ef_bss,nbss,nleads,nebss,
     &            v,kb*t,ef,emforces)
              else

              call BS_M1_nc2(N1,gfgeneral%matdense%a,
     &            rhobs_general,
     &            ematbs_general,
     &            Delta,ERealGrid%w(I),fl,fr,
     &            bs_add,bs_method, NspinComplexMatrix,ispin,ei,
     &            gammamp,sigmamp,
     &            ef_bss,nbss,nleads,nebss,
     &            v,kb*t,ef,emforces)

              endif

            endif
 
!#ifdef MP  I
!        c  all MPI_Barrier(MPI_Comm_World,MPIerror)
!        c  all MPI_Finalize( MPIerror )
!#endif
!        s  top
!



! ********  *********************************************************
! Calculate the left Gamma matrix for the left lead in order to calculate
! the non-  equilibrium Greens function

            if(ldos)then
              fl=1.0D0
              fr=0.0D0
            endif
            if(emtimings)then
              CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')
     $             'kr_copyrhoystart',(scb_1-scb_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
            endif

            call ConjugateTransposeMatrixCOMP(sigmal,
     &          sigmaleads(1)%sigma,nl,nl)
            if(emtimings)then
              CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')
     $             'kr_copyrhoy0a',(scb_1-scb_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
            endif


            call SubstractMultiplyCMatrixCOMP(gammal,
     &          sigmaleads(1)%sigma,sigmal,nl,nl,-(fR-fL) * zi)

            if(emtimings)then
              CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')
     $             'kr_copyrhoy0b',(scb_1-scb_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
            endif



            CALL ZGEMM('N','N',N1,NL,NL,(1.D0,0.D0),GF_iter1l,N1,
     &          gammal,NL,(0.D0,0.D0),GF_iter2l,N1)
            


            if(emtimings)then
              CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')
     $             'kr_copyrhoy2',(scb_1-scb_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
            endif



            ALLOCATE(gf1(NL,N1),gf2(NL,N1))

            call ConjugateTransposeMatrixCOMP(gf1,
     &          GF_iter1l,nl,n1)

            if(emtimings)then
              CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')
     $             'kr_copyrhoy3a',(scb_1-scb_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
            endif



            call TransposeMatrixCOMP(gf2,
     &          GF_iter2l,nl,n1)
            if(emtimings)then
              CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')
     $             'kr_copyrhoy3b',(scb_1-scb_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
            endif




            if(NspinComplexMatrix<=2)then
              allocate(b(rhogeneralp(1)%matSparse%nnz))
              b=0.0D0
              if(emtimings)then
                CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
                write(12347,'(A,f12.6)')
     $               'kr_copyrhoy4',(scb_1-scb_0)*1.0d0/sc_r
                CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
              endif

              call UpdateRhoNEQ(rhogeneralp(1)%matSparse%nnz,n1,nl,nr,
     &            nl,
     &            rhogeneralp(1)%matSparse%q,rhogeneralp(1)%matSparse%j,
     &            b(:),gf1,gf2,
     &            weightrho * (1.D0,0D0)/(2.D0*PI)*ERealGrid%w(I),
     &            .false.)

              if(emtimings)then
                CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
                write(12347,'(A,f12.6)')
     $               'kr_copyrhoy5',(scb_1-scb_0)*1.0d0/sc_r
                CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
              endif


              rhogeneralp(ispin)%matSparse%b(:)=
     &            rhogeneralp(ispin)%matSparse%b(:)+b(:)
          
              if(emtimings)then
                CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
                write(12347,'(A,f12.6)')
     $               'kr_copyrhoy6',(scb_1-scb_0)*1.0d0/sc_r
                CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
              endif
              if(emforces)
     &            ematgeneralp(ispin)%matSparse%b(:)=
     &            ematgeneralp(ispin)%matSparse%b(:)+Ei * b(:)

              if(emtimings)then
                CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
                write(12347,'(A,f12.6)')
     $               'kr_copyrhoy7',(scb_1-scb_0)*1.0d0/sc_r
                CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
              endif
              deallocate(b)

            else
              allocate(bnc(rhogeneralp(1)%matSparse%nnz,
     .            NspinComplexMatrix))
              bnc=0.0D0
              call UpdateRhoNEQ_nc(rhogeneralp(1)%matSparse%nnz,n1,
     &          nl/2,nr/2,nl,
     &          rhogeneralp(1)%matSparse%q,rhogeneralp(1)%matSparse%j,
     &          bnc(:,1),bnc(:,2), bnc(:,3),bnc(:,4), gf1,gf2,
     &          weightrho * (1.D0,0D0)/(2.D0*PI)*ERealGrid%w(I),
     &          .false.)

              do is=1,NspinComplexMatrix
                rhogeneralp(is)%matSparse%b=
     &              rhogeneralp(is)%matSparse%b+bnc(:,is)
              enddo

              if(emforces)then
                do is=1,NspinComplexMatrix
                  ematgeneralp(is)%matSparse%b=
     &                ematgeneralp(is)%matSparse%b+ei* bnc(:,is)
                enddo
              endif

              deallocate(bnc)
            endif
 
!            call MPI_Barrier(inverse_comm,MPIerror)
            if(emtimings)then
              CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')
     $             'kr_copyrho',(scb_1-scb_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
            endif
            deallocate(gf1,gf2)

! ********  *********************************************************
! Calculate the right Gamma matrix for the right lead in order to calculate
! the non-  equilibrium Greens function

            if(ldos)then
              fl=0.0D0
              fr=1.0D0
            endif
            if(emtimings)then
              CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')
     $             'kr_copyrhoxstart',(scb_1-scb_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
            endif


            call ConjugateTransposeMatrixCOMP(sigmar,
     &          sigmaleads(2)%sigma,nr,nr)
            if(emtimings)then
              CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')
     $             'kr_copyrhox0a',(scb_1-scb_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
            endif



            call SubstractMultiplyCMatrixCOMP(gammar,
     &          sigmaleads(2)%sigma,sigmar,nr,nr,(fR-fL) * zi)

            if(emtimings)then
              CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')
     $             'kr_copyrhox0b',(scb_1-scb_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
            endif




            CALL ZGEMM('N','N',N1,NR,NR,(1.D0,0.D0),GF_iter1r,N1,
     &          gammar,NR,(0.D0,0.D0),GF_iter2r,N1)
 

            if(emtimings)then
              CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')
     $             'kr_copyrhox2',(scb_1-scb_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
            endif


! ********  **********************************************************
! Calculate the "out of equilibrium" contribution to the charge for the
! left lea  d

            ALLOCATE(gf1(NR,N1),gf2(NR,N1))

            call ConjugateTransposeMatrixCOMP(gf1,
     &          GF_iter1r,nr,n1)

            if(emtimings)then
              CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')
     $             'kr_copyrhox3a',(scb_1-scb_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
            endif



            call TransposeMatrixCOMP(gf2,
     &          GF_iter2r,nr,n1)

            if(emtimings)then
              CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')
     $             'kr_copyrhox3b',(scb_1-scb_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
            endif


            if(NspinComplexMatrix<=2)then
              allocate(b(rhogeneralp(1)%matSparse%nnz))
              b=0.0D0

              if(emtimings)then
                CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
                write(12347,'(A,f12.6)')
     $               'kr_copyrhox4',(scb_1-scb_0)*1.0d0/sc_r
                CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
              endif
              call UpdateRhoNEQ(rhogeneralp(1)%matSparse%nnz,n1,
     &            nl,nr,nr,
     &            rhogeneralp(1)%matSparse%q,rhogeneralp(1)%matSparse%j,
     &            b(:),gf1,gf2,
     &            (1D0 - weightrho) *(1.D0,0D0)/(2.D0*PI)*ERealGrid%w(I)
     &            ,.false.)
              if(emtimings)then
                CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
                write(12347,'(A,f12.6)')
     $               'kr_copyrhox5',(scb_1-scb_0)*1.0d0/sc_r
                CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
              endif


              rhogeneralp(ispin)%matSparse%b(:)=
     &            rhogeneralp(ispin)%matSparse%b(:)+b(:)
 
              if(emtimings)then
                CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
                write(12347,'(A,f12.6)')
     $               'kr_copyrhox6',(scb_1-scb_0)*1.0d0/sc_r
                CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
              endif
              if(emforces)
     &            ematgeneralp(ispin)%matSparse%b(:)=
     &            ematgeneralp(ispin)%matSparse%b(:)+Ei * b(:)

              if(emtimings)then
                CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
                write(12347,'(A,f12.6)')
     $               'kr_copyrhox7',(scb_1-scb_0)*1.0d0/sc_r
                CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
              endif

              deallocate(b)



            else
              allocate(bnc(rhogeneralp(1)%matSparse%nnz,
     .            NspinComplexMatrix))
              bnc=0.0D0
              call UpdateRhoNEQ_nc(rhogeneralp(1)%matSparse%nnz,n1,
     &            nl/2,nr/2,nr,
     &            rhogeneralp(1)%matSparse%q,rhogeneralp(1)%matSparse%j,
     &            bnc(:,1),bnc(:,2), bnc(:,3),bnc(:,4), gf1,gf2,
     &            (1D0 - weightrho) *(1.D0,0D0)/(2.D0*PI)*ERealGrid%w(I)
     &            ,.false.)

              do is=1,NspinComplexMatrix
                rhogeneralp(is)%matSparse%b=
     &              rhogeneralp(is)%matSparse%b+bnc(:,is)
              enddo

              if(emforces) then
                do is=1,NspinComplexMatrix
                   ematgeneralp(is)%matSparse%b=
     &                ematgeneralp(is)%matSparse%b+ei* bnc(:,is)
                enddo
              endif

              deallocate(bnc)
            endif

            deallocate(gf1,gf2)

            if(emtimings)then
              CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')
     $             'kr_copyrho2',(scb_1-scb_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
            endif

          endif
!          CALL TIMER('bs_m1',2)
          do il=1,ERealGrid%nLeads
            deallocate(sigmaleads(il)%sigma)
          enddo
          deallocate(sigmaleads)

        enddo

      enddo

      if(emtimings)then
        CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
        write(12347,'(A,f12.6)')
     $       'kr_time2',(scb_1-sca_0)*1.0d0/sc_r
        CALL SYSTEM_CLOCK(sca_0,sc_r,sc_m)
      endif

 

      if(mynode_inverse.eq.0.and.bs_add.and.bs_method.eq.1)then
        call DestroyMatrixGeneral(gfgeneral,"keldyshreal",io)
      endif

      if(emtimings)then
        CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
        write(12347,'(A,f12.6)')
     $       'kr_outeneloop',(sc_1-sc_0)*1.0d0/sc_r
        CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
      endif


      if(mynode_inverse.eq.0)then
        if((bs_add).and.(bs_method.eq.1))then

          do i1=1,nbss
            call DestroyMatrixGeneral(gammamp(i1),"keldyshreal",io)
            call DestroyMatrixGeneral(sigmamp(i1),"keldyshreal",io)
          enddo
          deallocate(gammamp,sigmamp)

          if(bssc.eq.0)then

            if(ik.eq.nk)then
              do i1=1,nleads
                if(NspinComplexMatrix.eq.1)then
                  write(12347,*)"efib_bssuse=",i1,13.6057D0 * 
     .                (ef_bss(i1,1)-ef),ik
                else
                  write(12347,*)"efib_bssuse=",i1,13.6057D0 * 
     .                (ef_bss(i1,1)-ef),13.6057D0 * 
     .                (ef_bss(i1,2)-ef),ik
                endif
              enddo
            endif

            call transm_probes_tot(ef_bss,tij,wkmod,nleads,nbss,
     .      NspinComplexMatrix,
     .      ERealGrid%nEnergies,ik,nk,v,kb * t,ef)

            if(ik.eq.nk)then

              if(myhead.eq.0)then
                open(12349,file='DM.EF.dat')
                do ispin=1,NspinComplexMatrix
                  do i1=1,nleads
                    write(12349,'(f15.10,A3)',ADVANCE='NO')
     .                 13.6057D0 *(ef_bss(i1,ispin)-ef), "   "
                  enddo
                  write(12349,*)
                enddo
                close(12349)
              endif
              
              do i1=1,nleads
                if(NspinComplexMatrix.eq.1)then
                  write(12347,*)"efi_bssuse=",i1,13.6057D0 * 
     .                (ef_bss(i1,1)-ef),ik
                else
                  write(12347,*)"efi_bssuse=",i1,13.6057D0 * 
     .                (ef_bss(i1,1)-ef),13.6057D0 * 
     .                (ef_bss(i1,2)-ef),ik
                endif
              enddo
            endif



            if(bs_writetrc)
     .          call bs_transm(ef_bssk,ef_bss,NspinComplexMatrix,
     .          ERealGrid%nEnergies,V,ik,nk,nleads,nbss,ef,kb * T,
     .          nebss,tij,wkmod)
          endif
          deallocate(ef_bssk,nebss,deltabss,tij)
        endif

        if(bs_add) deallocate(Tau1_aux)
        deallocate(IPIV)
        deallocate(zv,work3,rwork2)
        deallocate(GF_iter1r,gammar,sigmar)
        deallocate(GF_iter1l,gammal,sigmal)
        deallocate(GF_iter2l)
        deallocate(GF_iter2r)
        deallocate(GF_iter3)

! ####  ### calculate bound states from eigenvalues #########

        if((bs_add).and.(bs_method.eq.0))then

          if(.not.ALLOCATED(rhobs_delta))
     &        ALLOCATE(rhobs_delta(N1,N1,NspinComplexMatrix))
          if(MOD(itermod-1,bsskip) .EQ. 0) then
            CALL TIMER('bssel',1)

            rhobs_delta=0D0
            call bs_m0(ERealGrid%nEnergiesGlobal,N1,NspinComplexMatrix,
     &          ERealGrid%nEnergies,ik,nk,NL,NR,
     &          V,ef,T,H0_L,
     &          H1_L,S0_L,S1_L,H0_R, H1_R,S0_R,S1_R, deltaimag,
     &          hdense,sdense,sdense_inv,
     &          delta,rhobs_delta,eiene,ene,indemax,wkmod)
            if(bsskip.gt.1) rhobstot(:,:,:,ik)=rhobs_delta

            CALL TIMER('bssel',2)
          else
            rhobs_delta=rhobstot(:,:,:,ik)
          endif

        endif
        if(bs_add.and.bs_method.eq.0)then
          deallocate(hdense,sdense)
        endif
        if(allocated(ef_bss))deallocate(ef_bss)
        if(emtimings)then
          CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
          write(12347,'(A,f12.6)')
     $         'kr_finalizing',(sc_1-sc_0)*1.0d0/sc_r
        endif
      endif
      if(emtimings)then
        CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
        write(12347,'(A,f12.6)')
     $       'kr_time3',(scb_1-sca_0)*1.0d0/sc_r
        CALL SYSTEM_CLOCK(sca_0,sc_r,sc_m)
      endif



! ####### end calculate bound states from eigenvalues #########

      end subroutine keldyshReal


      subroutine keldyshImag(N1,NSpinBlocks,NspinComplexMatrix,
     &    NL,NR,ik,iter,
     &    NEnerg1,NEnerg2,NPOLES,R0,
     &    V,T,ef, hgeneralp,sgeneralp,rhogeneralp,ematgeneralp,
     &    weightrho)

! *****************************************************************
! Calculates the Keldysh Green's function on the complex plane
! along a predefined contour.
!
! Written by Alexandre Reily Rocha and Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: rochaa@tcd.ie, reilya@if.usp.br, runggeri@tcd.ie
! ********** HISTORY **********************************************
! Original version:      October 2007
! *****************************************************************

      use negfmod, only: gamma_negf,negfon,inversion_solver,
     & geterrsigma,emtimings,nprocs_hs,emforces,
     & ComputeImpurityGfMatsubara,PrintImpurityGfMatsubara,
     & CallImpuritySolver
      use mTypes
      use mMatrixUtil
      use mONInterface
      use mMPI_NEGF
      use mEnergyGrid
      use sigma, only: h0_L,H1_L,S0_L,S1_L,h0_r,H1_r,S0_r,S1_r
      use mSigmaMethod1, only: check_error_sigma
      use mImpuritySolver,only: 
     .    PrintGFHSMatsubaraGeneral,AddRhoTilde,CalculateSinv,
     .    RegularizeGFPolesGeneral,CTQMCHyb_ImpuritySolverInterface
      
      IMPLICIT NONE

      include "const2.h" 
                  
      double precision, intent(in) :: weightrho
      integer, intent(in) :: NSpinBlocks,NspinComplexMatrix
      INTEGER :: N1,NL,NR,Iprime,Nenerg1,NEnerg2,iter,
     &    NPOLES,ISPIN,NenergImNode,I,JJ,II,ik,
     &    INFO,i1,i2,nnz,ind,j,indold,ind2,indnew,indadd,il,n1h

      INTEGER, DIMENSION (:), ALLOCATABLE :: IPIV 
      
      DOUBLE PRECISION :: ef,V,T,R0
      
      DOUBLE COMPLEX :: CONSTR,CONSTL,drhoij
      
      
      DOUBLE COMPLEX, DIMENSION (:), ALLOCATABLE :: WORK
      
      type(matrixTypeGeneral) :: gf,gfout,gfP,gfserial
      type(matrixTypeGeneral) :: hgeneralp(NspinComplexMatrix),
     .    sgeneralp,rhogeneralp(NspinComplexMatrix),
     .    ematgeneralp(NspinComplexMatrix)
      type(ioType) :: io
      integer  nnzrow(n1),sizeSigma
      double precision deltag,dsigma
      double complex gfij,gfji,gfij2,gfji2,ener_sigma
      double complex, allocatable :: rhobuf(:)

      type(SelfEnergyType), allocatable :: sigmaleads(:)

      integer gfmattype
      integer*4:: sc_0,sc_1,sc_r,sc_m
      INTEGER :: MPIerror
      logical writeGFHSheader
     
      

      if(.not.negfon)then
        gfmattype=0 !for dense GF matrix
      else
        gfmattype=2 !for sparse GF matrix
      endif

      io%isDebug=.false.


      if(ComputeImpurityGfMatsubara.and.iter.eq.1)then
!if we use k-points this needs to be re-run for each kpoints, or else
!stored in memory for all k-points
        call CalculateSinv(sgeneralp,gfmattype)
      endif

      if(.not.negfon)allocate(IPIV(N1),WORK(N1*N1)) 

      NenergImNode=EImagGrid%nEnergies
      DO ISPIN=1,EImagGrid%nSpin


        DO I=1,EImagGrid%nEnergies
!          write(12347,*)"ie=",i,ispin,n1
            
          if(emtimings)CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)

          if(EImagGrid%InfoSigma==0.or.EImagGrid%InfoSigma==1)then
            Iprime=I+myhead*NenergImNode
          else
            Iprime=(I-1)*nheads+myhead+1
          endif
              
          IF ((Iprime).LE.Nenerg1+NEnerg2) THEN
            IF (DREAL((EImagGrid%e(I)-ef-e*V/2.D0)/(kB*T)) .GT.
     &          40.D0) THEN
              CONSTL=0.D0
            ELSE
              IF ((Iprime) .LE. Nenerg1) THEN
                CONSTL=1.D0/(CDEXP(DCMPLX(DREAL(EImagGrid%e(I))-
     &              ef-e*V/2.D0)/(kB*T))+1.D0)
              ELSE
                CONSTL=zi*(EImagGrid%e(I)-DCMPLX(R0))/(CDEXP((
     &              EImagGrid%e(I)-ef-e*V/2.D0)/(kB*T))+1.D0)
              ENDIF
            ENDIF
          ELSE IF ((Iprime)
     &       .LE.NEnerg1+NEnerg2+NPOLES) THEN
            CONSTL=1.D0
          ELSE
            CONSTL=0.D0
          ENDIF

          IF ((Iprime) .LE. Nenerg1+NEnerg2) 
     &        THEN
            IF (DREAL((EImagGrid%e(I)-ef+e*V/2.D0)/(kB*T)) .GT.
     &         40.D0) THEN
              CONSTR=0.D0
            ELSE
              IF ((Iprime) .LE. Nenerg1) THEN
                CONSTR=1.D0/(CDEXP(DCMPLX(DREAL(EImagGrid%e(I))-
     &              ef+e*V/2.D0)/(kB*T))+1.D0)
              ELSE
                CONSTR=zi*(EImagGrid%e(I)-DCMPLX(R0))/(CDEXP((
     &              EImagGrid%e(I)-ef+e*V/2.D0)/(kB*T))+1.D0)
              ENDIF
            ENDIF
          ELSE IF ((Iprime).LE.
     &        NEnerg1+NEnerg2+NPOLES)THEN
            CONSTR=0.D0
          ELSE
            CONSTR=1.D0
          ENDIF

!start allocating sparse matrix

             
          if(emtimings)then
            CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
            write(12347,'(A,f12.6)')
     $           'setconst',(sc_1-sc_0)*1.0d0/sc_r
            CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
          endif

          allocate(sigmaleads(EImagGrid%nLeads))

          do il=1,EImagGrid%nLeads
            call allocate_sigma_single(sigmaleads(il),
     &          EImagGrid%leadsTotalDim(il),myhead,
     &          EImagGrid%InfoSigma,
     &          (0.0D0,0.0D0),0.0D0,0.0D0)
          enddo

          if(mynode_inverse.eq.0)then
            call extractsigma(EImagGrid,sigmaleads,i,ispin,ik)
          endif

          if(nnodes_inverse.gt.1)then
            do il=1,EImagGrid%nLeads
              sizeSigma=EImagGrid%leadsTotalDim(il)**2
#ifdef MPI
              call MPI_Bcast(sigmaleads(il)%sigma(1,1),sizeSigma,
     .            DAT_dcomplex,0, inverse_comm,MPIerror)
#endif
            enddo
          endif


          if(emtimings)then
            CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
            write(12347,'(A,f12.6)')
     $           't_extractsigma',(sc_1-sc_0)*1.0d0/sc_r
            CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
          endif

          if(iter.eq.1.and.geterrsigma)then
            write(12347,*)"dsigmainfo=",EImagGrid%e(I),ispin,i
            ener_sigma=EImagGrid%e(I)-v*0.5D0
            call check_error_sigma(dsigma,.true.,ener_sigma,'L',nl,
     .          h0_L(:,:,ispin),H1_L(:,:,ispin),S0_L,S1_L,
     &          sigmaleads(1)%sigma)
            ener_sigma=EImagGrid%e(I)+v*0.5D0
            call check_error_sigma(dsigma,.true.,ener_sigma,'R',nr,
     .          h0_r(:,:,ispin),H1_r(:,:,ispin),S0_r,S1_r,
     &          sigmaleads(2)%sigma)


            if(emtimings)then
              CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')
     $             't_checkerrorsigma',(sc_1-sc_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
            endif
          endif


          if(mynode_inverse.eq.0)then

            sgeneralp%mattype=2
            hgeneralp(ispin)%mattype=2
            if(gfmattype.eq.0)then
              nnz=n1*n1
            else
              call findnnzgf2(nnz,nnzrow,n1,n1,nl,nr,
     &           sigmaleads(1)%sigma, sigmaleads(2)%sigma,
     .           hgeneralp(ispin))
            endif

            if(emtimings)then
              CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')
     $             'findnnzgf2',(sc_1-sc_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
            endif

!            write(12347,*)"nnz=",n1,nnz,ispin
!            write(*,*)"nnzrow=",nnzrow

            call AllocateMatrixGeneral(n1,n1,nnz,gfmattype,gf,
     .          "keldyshimag", io)

            call setgfelementsgeneral_nc(EImagGrid%e(I),
     &          NspinComplexMatrix,ispin,
     &          gf,nnz,n1,nl,nr,
     &          sigmaleads(1)%sigma, sigmaleads(2)%sigma,
     &          hgeneralp,sgeneralp)

          endif

          if(nprocs_hs.ne.1)then
            sgeneralp%mattype=3
            hgeneralp(ispin)%mattype=3
            call findnnzgf2P(nnz,nnzrow,
     *          hgeneralp(ispin)%matSparseP%matSparse%iRows,n1,nl,nr,
     &          sigmaleads(1)%sigma, sigmaleads(2)%sigma,
     .          hgeneralp(ispin))


            call AllocateMatrixGeneral(-1,inverse_comm,
     .          nnodes_inverse,mynode_inverse,n1,n1,
     .          hgeneralp(ispin)%matSparseP%matSparse%iRows,n1,1,
     .          hgeneralp(ispin)%matSparseP%matSparse%iVert,nnz,
     .          3,gfP,"gf_keldyshimag",io)

!firs run:hprime=h
!second run
!       h-> hprime=h+sigmamb
!xxx double counting: for half filling -U/2
            call setgfelementsgeneral_nc(EImagGrid%e(I),
     .          NspinComplexMatrix,ispin,
     &          gfP,nnz,n1,nl,nr,
     &          sigmaleads(1)%sigma, sigmaleads(2)%sigma,
     &          hgeneralp,sgeneralp)

            sgeneralp%mattype=2
            hgeneralp(ispin)%mattype=2

            call MaxDifferenceMatrixCRS_CRSP(gf,gfP,
     .          mynode_inverse,"(gf) ",io)
          endif


          if(emtimings)then
            CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
            write(12347,'(A,f12.6)')
     $           'setgfelementsgeneral',(sc_1-sc_0)*1.0d0/sc_r
            CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
          endif

!          CALL TIMER('gfinvi',1)

      if(.false.)then
        call PrintSparse2DenseReorderedNC(hgeneralp,n1/2,
     .      NspinComplexMatrix, NSpinBlocks, NspinComplexMatrix,1,
     .      "hk_3")
        call PrintSparse2DenseReorderedNC(sgeneralp,n1/2, 1,
     .      NSpinBlocks,NspinComplexMatrix,2,"sk_3")
        call PrintSparse2DenseReorderedNC(gf,n1/2, 1, NSpinBlocks,
     .      NspinComplexMatrix,3,"gk_3")
      endif


          if(.not.negfon) then        
            if(mynode_inverse.eq.0)then
              CALL ZGETRF(N1,N1,gf%matdense%a,N1,IPIV,INFO)
              CALL ZGETRI(N1,gf%matdense%a,N1,IPIV,WORK,N1**2,INFO)
            endif
          else

            if(mynode_inverse.eq.0)then
!Uncomment the commented out lines below to print out the GF before and
!after inversion
!              if(myhead==0) call PrintMatrixCRS(gf%matSparse,"gfinv",io)
              call InvertONGeneral(N1,gf,nl,nr,gfout,1,inversion_solver)
!              if(myhead==0) call PrintMatrixCRS(gf%matSparse,"gf",io)
!              call stopnegf
            endif

            if(nprocs_hs.ne.1)then

              if(emtimings)then
                CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
                write(12347,'(A,f12.6)')
     $               't_serial_inversion',(sc_1-sc_0)*1.0d0/sc_r
                CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
              endif



              call InvertONGeneral2(N1,gfp,gfserial,nl,nr,gfout,1,
     .            inversion_solver)


              if(emtimings)then
                CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
                write(12347,'(A,f12.6)')
     $               't_parallel_inversion',(sc_1-sc_0)*1.0d0/sc_r
                CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
              endif



              if(mynode_inverse.eq.0)then
                call MaxDifferenceMatrixCRS_CRS(gf,gfserial,
     .              mynode_inverse,"(gfeiinv) ",io)
              endif

              if(emtimings)then
                CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
                write(12347,'(A,f12.6)')
     $               't_difference_matrices',(sc_1-sc_0)*1.0d0/sc_r
                CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
              endif

            endif


          endif

          if(.false.)then
            call PrintSparse2DenseReorderedNC(gf,n1/2, 1, NSpinBlocks,
     .          NspinComplexMatrix,3,"gk_4")
          endif

          if(emtimings)then
            CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
            write(12347,'(A,f12.6)')
     $           'invertki',(sc_1-sc_0)*1.0d0/sc_r
            CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
          endif

          if(ComputeImpurityGfMatsubara)then
            if(i==1.and.mynode_inverse==0)then
              writeGFHSheader=.true.
!              write(12347,*)"setting writeheader T"
            else
              writeGFHSheader=.false.
!              write(12347,*)"setting writeheader F"
            endif

            call PrintGFHSMatsubaraGeneral(kb * T,EImagGrid%e(I),ef,gf,
     .          hgeneralp(ispin),sgeneralp,sigmaleads(1)%sigma,nl,
     .          sigmaleads(2)%sigma,nr,writeGFHSheader,
     .          myhead,PrintImpurityGfMatsubara,CallImpuritySolver,
     .          EImagGrid%ig(i),
     .          EImagGrid%nEnergiesGlobal,ispin,EImagGrid%nSpin)


            if(emtimings)then
              CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')
     $             'printgf',(sc_1-sc_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
            endif



            call RegularizeGFPolesGeneral(EImagGrid%e(I),ef,gf)

            if(emtimings)then
              CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')
     $             'RegularizeDelta',(sc_1-sc_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
            endif



          endif

          do il=1,EImagGrid%nLeads
            deallocate(sigmaleads(il)%sigma)
          enddo
          deallocate(sigmaleads)



!          CALL TIMER('gfinvi',2)

          if(mynode_inverse.eq.0)then
          
            call updaterho_nc(rhogeneralp,ematgeneralp,emforces,ispin,
     .          NspinComplexMatrix,gf, nl,nr,gfmattype, 
     .          EImagGrid%w(i),constl,constr,weightrho,EImagGrid%e(i),
     .          .false.)

      if(.false.)then
        call PrintSparse2DenseReorderedNC(rhogeneralp,n1/2,
     .      NspinComplexMatrix, NSpinBlocks, NspinComplexMatrix,1,
     .      "rhok_5")
      endif






            if(emtimings)then
              CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
              write(12347,'(A,f12.6)')
     $             'setrhoki2',(sc_1-sc_0)*1.0d0/sc_r
              CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
            endif

            call DestroyMatrixGeneral(gf,"keldyshimag2",io)
            call DestroyMatrixGeneral(gfserial,"keldyshimag2",io)
! end allocating sparse matrix
          endif
          if(nprocs_hs.ne.1)
     &        call DestroyMatrixGeneral(gfP,"keldyshimag",io)

!          write(12347,*)"iee=",i,ispin,n1
        ENDDO

        if(ComputeImpurityGfMatsubara)then

!          if(.true.)then
!            allocate(rhobuf(rhogeneralp(ispin)%matSparse%nnz))
!            CALL MPI_REDUCE(rhogeneralp(ispin)%matSparse%b(1),
!     &          rhobuf,
!     &          rhogeneralp(ispin)%matSparse%nnz, DAT_dcomplex,MPI_SUM,
!     &          0,inverseheads_comm,MPIerror)
!            rhogeneralp(ispin)%matSparse%b=rhobuf
!            deallocate(rhobuf)
!          endif
          if(myhead==0)
     &     call AddRhoTilde(rhogeneralp(ispin),gfmattype, nl,nr,.false.)

!          call MPI_Barrier(inverseheads_comm,MPIerror)
!          call stopnegf
        endif

      ENDDO

      if(CallImpuritySolver)then
        call CTQMCHyb_ImpuritySolverInterface(kb*T,myhead,nheads,
     .      inverseheads_comm)
!xxx: for charge self-consistency        call calculateRhogeneralManyBody
!xxx: propagate to other full GF
      endif

      if(.not.negfon)deallocate(IPIV,WORK) 
      
      end subroutine keldyshImag


