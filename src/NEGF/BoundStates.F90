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
!                   BS_M0_EIGENVALUES,
!                   BS_M1_NC2,
!                   BS_M1_NC,
!                   BS_M1,
!                   BS_COLLECTRHO_M1,
!                   BS_COLLECTRHO,
!                   EIOFE,
!                   REDUCEENE,
!                   BS_M0,
!                   EIOFE,
!                   REDUCEENE,
!                   GETEBS,
!                   GETRHOBS,
!                   TRANSM_PROBES_TOT,
!                   BS_TRANSM,
!                   TRANSMIJ_BS,
!                   SET_GAMMA_BS,
!                   FIND_EF_BSS,
!                   FIND_EFOUT  
! AND
! THE MODULE
!                   MBOUNDSTATES  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!
module mBoundStates
  use mConstants
  use mTypes
  use mMatrixUtil
    
  implicit none
  private
  public :: transmij_bs
  public :: find_ef_bss
  public :: bs_m0_eigenvalues
  public :: bs_m1
  public :: bs_m1_nc
  public :: bs_m1_nc2
  public :: transm_probes_tot
  public :: bs_transm
  public :: bs_m0
  public :: bs_collectrho_m1
  public :: bs_collectrho
  public :: getenebs
    
  contains

  subroutine bs_m0_eigenvalues(n1,nspin,i,bseskip,bsskip,ispin,ik,inde,indemax,gf_iter,H_chain,S0_chain_inv,nl,nr,work3,nw3,rwork2,nw2, Nenerg_div,nenerg_div_nodes,eiene,ene,ei)

! *****************************************************************
! Calculate the Eigenvalues of the effective Hamiltonian of the scattering
! region
!
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

    use negfmod, only: itermod, outinfo
    use mEnergyGrid

    implicit none
    integer n1,nspin,i,bseskip,bsskip,ispin,inde,indemax,ik,nl,nr, izv,nw3,nw2,info,Nenerg_div,nenerg_div_nodes
    double complex gf_iter(n1,n1),H_chain(n1,n1,NSPIN), S0_chain_inv(n1,n1)
    DOUBLE COMPLEX  zv(n1),veigl(1,1),veigr(1,1),work3(nw3)
    double precision rwork2(nw2),ei
    double complex eiene(Nenerg_div,N1,NSPIN)
    double precision ene(nenerg_div_nodes)

    if((MOD(i-1,bseskip) .EQ. 0))then
      CALL TIMER('bsev',1)
      IF (MOD(itermod-1,bsskip) .EQ. 0)then
    
        inde=inde+1
        indemax=inde
        GF_iter=H_Chain(:,:,ISPIN)

        GF_iter(1:NL,1:NL)=GF_iter(1:NL,1:NL)+ERealGrid%sigma(1,i,ispin,ik)%sigma
        GF_iter(N1-NR+1:N1,N1-NR+1:N1)= GF_iter(N1-NR+1:N1,N1-NR+1:N1)+ ERealGrid%sigma(2,i,ispin,ik)%sigma
        GF_iter=matmul(S0_Chain_inv,GF_iter)

        call ZGEEV( 'N', 'N', N1, GF_iter, N1, zv, veigl, 1, veigr, 1, WORK3, nw3 , RWORK2, INFO )

        do izv=1,N1
!     write(*,*)"inde=",zv(izv),inde
          eiene(inde,izv,ISPIN)=zv(izv)
          if (outinfo)  write(12347,*)"eneslr=",ei,DREAL(zv(izv)),DIMAG(zv(izv)),ispin,ik
        enddo
        ene(inde)=ei
      endif
      CALL TIMER('bsev',2)
    endif

  end subroutine bs_m0_eigenvalues

  SUBROUTINE BS_M1_NC2(N1,GF_iter,rhobs_general,ematbs_general,Delta,const,fl,fr, bs_add,bs_method, nspin,ispin,ei, gammamp,sigmamp,ef_bss,nbss,nleads,nebss, v,t,ef,emforces)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

    use negfmod, only : bs_nmid,bs_nmid1,bs_nmid2,emtimings

    use mTypes

    IMPLICIT NONE
    
    INCLUDE "const2.h"

    INTEGER N1,bs_method,bs_nmide, II,JJ,nleads,nbss,ispin,nspin,jindex,kk
    LOGICAL bs_add
    logical, intent(in) :: emforces
    DOUBLE COMPLEX GF_iter(N1,N1),const,fl,fr, const2
    double complex, allocatable :: gpart(:,:),gpartaux(:,:), gammadense(:,:),gammadense2(:,:)
    double complex, allocatable :: gf1(:,:),gf2(:,:)
    DOUBLE PRECISION aln(N1),al(N1),alr(N1),Delta, ei,ef,v,t,fi,fe(n1),fmid
    double precision ef_bss(nleads,nspin),faux
    integer nebss(nleads,2),i1,i2
    type(matrixTypeGeneral) :: gammamp(nleads)
    type(matrixTypeGeneral) :: sigmamp(nleads)
    type(matrixTypeGeneral) :: rhobs_general(nspin)
    type(matrixTypeGeneral) :: ematbs_general(nspin)
    integer irows,ni,ni2
    integer*4:: sc_0,sc_1,sc_r,sc_m,scb_0,scb_1

    if(emtimings)then
      CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
    endif
 

    fmid=(fl+fr)* 0.5D0
    const2=(1.D0/(2.D0*PI))*CONST
    irows=n1/2

    do i1=1,nbss

      ni=gammamp(i1)%iRows
      ni2=2 * ni

      fi=(Ei-ef_bss(i1,ispin))/t
      if(fi.gt.0)then
        fi=EXP(-fi)/(1D0+EXP(-fi))
      else
        fi=1.D0/(1D0+EXP(fi))
      endif

      allocate(gpart(n1,ni2))
      allocate(gpartaux(n1,ni2))
      allocate(gammadense(ni,ni))
      allocate(gammadense2(ni2,ni2))

      if(emtimings)then
        CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
      endif

      call sparsetodense(gammamp(i1)%MatSparse,gammadense,ni)
      gammadense=(fi-fmid) * gammadense
      gammadense2=0.0D0
      gammadense2(1:ni,1:ni)=gammadense
      gammadense2(ni+1:ni2,ni+1:ni2)=gammadense

!      call writemat7c(ei,GF_iter,n1,n1,"gftotal")
      
      gpart(:,1:ni)=GF_iter(:,nebss(i1,1):nebss(i1,2))
      gpart(:,ni+1:ni2)=GF_iter(:,irows+nebss(i1,1):irows+nebss(i1,2))

      if(emtimings)then
        CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
        write(12347,'(A,f12.6)')'bs_m1_nc2_copymatrix',(scb_1-scb_0)*1.0d0/sc_r
        CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
      endif




      CALL ZGEMM('N','N',n1,ni2,ni2,(1.D0,0.D0),gpart,n1,gammadense2,ni2,(0.D0,0.D0),gpartaux,n1)
!      write(12347,*)"maxg1=",maxval(abs(gammadense)),maxval(abs(gpart)),maxval(abs(gpartaux))
      if(emtimings)then
        CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
        write(12347,'(A,f12.6)')'bs_m1_nc2_matmul',(scb_1-scb_0)*1.0d0/sc_r
        CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
      endif



       ALLOCATE(gf1(ni2,n1),gf2(ni2,n1))
       gf1=DCONJG(transpose(gpart))
       gf2=transpose(gpartaux)

       if(emtimings)then
         CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
         write(12347,'(A,f12.6)')'bs_m1_nc2_copymatrix2',(scb_1-scb_0)*1.0d0/sc_r
         CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
       endif


       call UpdateRhoNEQ_nc(rhobs_general(1)%matSparse%nnz,n1, ni,ni,ni2, rhobs_general(1)%matSparse%q,rhobs_general(1)%matSparse%j,rhobs_general(1)%matSparse%b,rhobs_general(2)%matSparse%b, rhobs_general(3)%matSparse%b,rhobs_general(4)%matSparse%b, gf1,gf2, const2,.true.)

        if(emtimings)then
         CALL SYSTEM_CLOCK(scb_1,sc_r,sc_m)
         write(12347,'(A,f12.6)')'bs_m1_nc2_updaterho',(scb_1-scb_0)*1.0d0/sc_r
         CALL SYSTEM_CLOCK(scb_0,sc_r,sc_m)
       endif



       deallocate(gf1,gf2)


      deallocate(gpart)
      deallocate(gpartaux)
      deallocate(gammadense)
      deallocate(gammadense2)


    enddo
    if(emtimings)then
      CALL SYSTEM_CLOCK(sc_1,sc_r,sc_m)
      write(12347,'(A,f12.6)')'bs_m1_nc2',(sc_1-sc_0)*1.0d0/sc_r
      CALL SYSTEM_CLOCK(sc_0,sc_r,sc_m)
    endif

  END SUBROUTINE BS_M1_NC2



  SUBROUTINE BS_M1_NC(N1,GF_iter,rhobs_general1,rhobs_general2,ematbs_general1,ematbs_general2,Delta,const,fl,fr, bs_add,bs_method, nspin,ispin,ei, gammamp,sigmamp,ef_bss,nbss,nleads,nebss, v,t,ef,emforces)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

    use negfmod, only : bs_nmid,bs_nmid1,bs_nmid2, outinfo

    use mTypes

    IMPLICIT NONE
    
    INCLUDE "const2.h"

    INTEGER N1,bs_method,bs_nmide, II,JJ,nleads,nbss,ispin,nspin,jindex,kk
    LOGICAL bs_add
    logical, intent(in) :: emforces
    DOUBLE COMPLEX GF_iter(N1,N1),const,fl,fr, const2
    double complex, allocatable :: gpart(:,:),gpartaux(:,:), gammadense(:,:)
    DOUBLE PRECISION aln(N1),al(N1),alr(N1),Delta, ei,ef,v,t,fi,fe(n1),fmid
    double precision ef_bss(nleads,nspin),faux
    integer nebss(nleads,2),i1,i2
    type(matrixTypeGeneral) :: gammamp(nleads)
    type(matrixTypeGeneral) :: sigmamp(nleads)
    type(matrixTypeGeneral) :: rhobs_general1
    type(matrixTypeGeneral) :: rhobs_general2
    type(matrixTypeGeneral) :: ematbs_general1
    type(matrixTypeGeneral) :: ematbs_general2
    integer irows



    fmid=(fl+fr)* 0.5D0
    const2=(1.D0/(2.D0*PI))*CONST

    do i1=1,nbss

      fi=(Ei-ef_bss(i1,ispin))/t
      if(fi.gt.0)then
        fi=EXP(-fi)/(1D0+EXP(-fi))
      else
        fi=1.D0/(1D0+EXP(fi))
      endif


      allocate(gpart(n1,gammamp(i1)%iRows))
      allocate(gpartaux(n1,gammamp(i1)%iRows))
      allocate(gammadense(gammamp(i1)%iRows,gammamp(i1)%iRows))

      call sparsetodense(gammamp(i1)%MatSparse,gammadense,gammamp(i1)%iRows)
      gammadense=(fi-fmid) * gammadense

!      call writemat7c(ei,GF_iter,n1,n1,"gftotal")

      gpart(:,:)=GF_iter(:,nebss(i1,1):nebss(i1,2))
!      call writemat7c(ei,gpart,n1,gammamp(i1)%iRows,"gpart1")

      CALL ZGEMM('N','N',N1,gammamp(i1)%iRows,gammamp(i1)%iRows, (1.D0,0.D0),gpart,N1, gammadense,gammamp(i1)%iRows,(0.D0,0.D0),gpartaux,N1)
      if (outinfo) write(12347,*)"maxg1=",maxval(abs(gammadense)),maxval(abs(gpart)),maxval(abs(gpartaux))

      do ii=1,rhobs_general1%iRows
        do jj=rhobs_general1%matSparse%q(ii), rhobs_general1%matSparse%q(ii+1)-1

          jindex=rhobs_general1%matSparse%j(jj)

          do kk=1,gammamp(i1)%iRows
            rhobs_general1%matSparse%b(jj)= rhobs_general1%matSparse%b(jj)+const2 * gpartaux(ii,kk)*DCONJG(gpart(jindex,kk))
!            if(kk<4) write(12347,*)"rhe1=",ii,jj,kk,rhobs_general1%matSparse%b(jj),gpartaux(ii,kk),DCONJG(gpart(jindex,kk))
          enddo

        enddo
      enddo

      if(emforces)then
        do ii=1,rhobs_general1%iRows
          do jj=rhobs_general1%matSparse%q(ii), rhobs_general1%matSparse%q(ii+1)-1

            jindex=rhobs_general1%matSparse%j(jj)
            do kk=1,gammamp(i1)%iRows
              ematbs_general1%matSparse%b(jj)= ematbs_general1%matSparse%b(jj)+ ei * const2 * gpartaux(ii,kk)*DCONJG(gpart(jindex,kk))
            enddo

          enddo
        enddo

      endif

 


      irows=n1/2
      if (outinfo) write(12347,*)"maxbsrho2=",maxval(abs(rhobs_general1%matSparse%b)),n1,irows

      gpart(:,:)=GF_iter(:,irows+nebss(i1,1):irows+nebss(i1,2))

!      call writemat7c(ei,gpart,n1,gammamp(i1)%iRows,"gpart2")
      CALL ZGEMM('N','N',N1,gammamp(i1)%iRows,gammamp(i1)%iRows, (1.D0,0.D0),gpart,N1, gammadense,gammamp(i1)%iRows,(0.D0,0.D0),gpartaux,N1)
      if (outinfo) write(12347,*)"maxg2=",maxval(abs(gammadense)),maxval(abs(gpart)),maxval(abs(gpartaux))

      do ii=1,rhobs_general1%iRows
        do jj=rhobs_general1%matSparse%q(ii), rhobs_general1%matSparse%q(ii+1)-1

          jindex=rhobs_general1%matSparse%j(jj)

          do kk=1,gammamp(i1)%iRows
            rhobs_general2%matSparse%b(jj)= rhobs_general2%matSparse%b(jj)+const2 * gpartaux(irows+ii,kk)*DCONJG(gpart(irows+jindex,kk))
!            if(kk<4) write(12347,*)"rhe2=",ii,jj,kk,rhobs_general2%matSparse%b(jj),gpartaux(ii,kk),DCONJG(gpart(jindex,kk))
          enddo

        enddo
      enddo
      if (outinfo) write(12347,*)"maxbsrho2=",maxval(abs(rhobs_general2%matSparse%b))

      if(emforces)then
        do ii=1,rhobs_general1%iRows
          do jj=rhobs_general1%matSparse%q(ii), rhobs_general1%matSparse%q(ii+1)-1

            jindex=rhobs_general1%matSparse%j(jj)
            do kk=1,gammamp(i1)%iRows
              ematbs_general2%matSparse%b(jj)= ematbs_general2%matSparse%b(jj)+ ei * const2 * gpartaux(ii,kk)*DCONJG(gpart(jindex,kk))
            enddo

          enddo
        enddo

      endif

      deallocate(gpart)
      deallocate(gpartaux)
      deallocate(gammadense)


    enddo


  END SUBROUTINE BS_M1_NC



  SUBROUTINE BS_M1 (N1,GF_iter,rhobs_generals,ematbs_generals,Delta,const,fl,fr, bs_add,bs_method, nspin,ispin,ei, gammamp,sigmamp,ef_bss,nbss,nleads,nebss, v,t,ef,emforces)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

    use negfmod, only : bs_nmid,bs_nmid1,bs_nmid2

    use mTypes

    IMPLICIT NONE
    
    INCLUDE "const2.h"

    INTEGER N1,bs_method,bs_nmide, II,JJ,nleads,nbss,ispin,nspin,jindex,kk
    LOGICAL bs_add
    logical, intent(in) :: emforces
    DOUBLE COMPLEX GF_iter(N1,N1),const,fl,fr, const2
    double complex, allocatable :: gpart(:,:),gpartaux(:,:), gammadense(:,:)
    DOUBLE PRECISION aln(N1),al(N1),alr(N1),Delta, ei,ef,v,t,fi,fe(n1),fmid
    double precision ef_bss(nleads,nspin),faux
    integer nebss(nleads,2),i1,i2
    type(matrixTypeGeneral) :: gammamp(nleads)
    type(matrixTypeGeneral) :: sigmamp(nleads)
    type(matrixTypeGeneral) :: rhobs_generals
    type(matrixTypeGeneral) :: ematbs_generals


    fmid=(fl+fr)* 0.5D0

    if((bs_add).and.(bs_method.eq.0)) then
      do II=1,bs_nmid
        al(II)=1.0D0
      enddo
      do II=bs_nmid+1,N1
        al(II)=0.0D0
      enddo
    elseif(.not.bs_add)then
      al=0.5D0
    endif


    if((bs_add).and.(bs_method.eq.1)) then

      const2=(1.D0/(2.D0*PI))*CONST

      do i1=1,nbss

        fi=(Ei-ef_bss(i1,ispin))/t
        if(fi.gt.0)then
          fi=EXP(-fi)/(1D0+EXP(-fi))
        else
          fi=1.D0/(1D0+EXP(fi))
        endif


        allocate(gpart(n1,gammamp(i1)%iRows))
        allocate(gpartaux(n1,gammamp(i1)%iRows))
        allocate(gammadense(gammamp(i1)%iRows,gammamp(i1)%iRows))

        gammadense=0D0
        do ii=1,gammamp(i1)%iRows
          do jj=gammamp(i1)%matSparse%q(ii), gammamp(i1)%matSparse%q(ii+1)-1
            gammadense(ii, gammamp(i1)%matSparse%j(jj))= gammadense(ii, gammamp(i1)%matSparse%j(jj))+ (fi-fmid) * gammamp(i1)%matSparse%b(jj)
          enddo
        enddo

        gpart(:,:)=GF_iter(:,nebss(i1,1):nebss(i1,2))

        CALL ZGEMM('N','N',N1,gammamp(i1)%iRows,gammamp(i1)%iRows, (1.D0,0.D0),gpart,N1, gammadense,gammamp(i1)%iRows,(0.D0,0.D0),gpartaux,N1)


        do ii=1,rhobs_generals%iRows
          do jj=rhobs_generals%matSparse%q(ii), rhobs_generals%matSparse%q(ii+1)-1

            jindex=rhobs_generals%matSparse%j(jj)

            do kk=1,gammamp(i1)%iRows
              rhobs_generals%matSparse%b(jj)= rhobs_generals%matSparse%b(jj)+const2 * gpartaux(ii,kk)*DCONJG(gpart(jindex,kk))
            enddo

          enddo
        enddo

        if(emforces)then
          do ii=1,rhobs_generals%iRows
            do jj=rhobs_generals%matSparse%q(ii), rhobs_generals%matSparse%q(ii+1)-1

              jindex=rhobs_generals%matSparse%j(jj)
              do kk=1,gammamp(i1)%iRows
                ematbs_generals%matSparse%b(jj)= ematbs_generals%matSparse%b(jj)+ ei * const2 * gpartaux(ii,kk)*DCONJG(gpart(jindex,kk))
              enddo

            enddo
          enddo

        endif

        deallocate(gpart)
        deallocate(gpartaux)
        deallocate(gammadense)


      enddo

    endif

  END SUBROUTINE BS_M1


  SUBROUTINE bs_collectrho_m1(N1,NSPIN, mynode,rhobs_general,sgeneral)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

    use negfmod
    use mTypes
    use mMPI_NEGF
    IMPLICIT NONE
    
    INTEGER N1,NSPIN,II,JJ,mynode,ISPIN,jindex
    double complex, allocatable :: aux_sparse(:)
    type(matrixTypeGeneral) :: rhobs_general(nspin)
    type(matrixTypeGeneral) :: sgeneral
    DOUBLE COMPLEX  trrhos
    INTEGER :: MPIerror


#ifdef MPI
#ifdef NODAT

    allocate(aux_sparse(rhobs_general(1)%matSparse%nnz))
    do ispin=1,nspin
      CALL MPI_REDUCE(rhobs_general(ispin)%matSparse%b, aux_sparse,rhobs_general(1)%matSparse%nnz, MPI_DOUBLE_COMPLEX,MPI_SUM,0,inverseheads_comm,MPIerror)
      rhobs_general(ispin)%matSparse%b(:)=aux_sparse(:)
    enddo
    deallocate(aux_sparse)
#else
    allocate(aux_sparse(rhobs_general(1)%matSparse%nnz))
    do ispin=1,nspin
      CALL MPI_REDUCE(rhobs_general(ispin)%matSparse%b(1), aux_sparse(1),rhobs_general(1)%matSparse%nnz, DAT_dcomplex,MPI_SUM,0,inverseheads_comm,MPIerror)
      rhobs_general(ispin)%matSparse%b(:)=aux_sparse(:)
    enddo
    deallocate(aux_sparse)

#endif
#endif

    if(MyNode.ne.0)then
      do ispin=1,nspin
        rhobs_general(ispin)%matSparse%b(:)=0D0
      enddo
    endif


    if(MyNode.eq.0) then

      do ispin=1,nspin ! if bs_add true this should only go to Min(2,nspin)
        trrhos=0D0
        do ii=1,rhobs_general(ispin)%iRows
          do jj=rhobs_general(ispin)%matSparse%q(ii), rhobs_general(ispin)%matSparse%q(ii+1)-1

            jindex=rhobs_general(ispin)%matSparse%j(jj)
            trrhos=trrhos+ rhobs_general(ispin)%matSparse%b(jj) * DCONJG(sgeneral%matSparse%b(jj))
!     .            S0_Chain(jindex,II)
          enddo
        enddo

        if (outinfo) write(12347,*)"tr_psi[rho s]=",trrhos, ispin,ikpmod,itermod
      enddo

    endif

  END SUBROUTINE bs_collectrho_m1


  SUBROUTINE bs_collectrho(N1,NSPIN,S0_Chain, mynode)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

    use negfmod
    use mMPI_NEGF
    IMPLICIT NONE
    
    INTEGER N1,NSPIN,II,JJ,mynode,ISPIN
    DOUBLE COMPLEX  S0_Chain(N1,N1),trrhos
    INTEGER :: MPIerror

#ifdef MPI
#ifdef NODAT
      CALL MPI_REDUCE(rhobs_delta,aux_par,N1*N1*NSPIN, MPI_DOUBLE_COMPLEX,MPI_SUM,0, inverseheads_comm,MPIerror)
      rhobs_delta=aux_par
#else
      CALL MPI_REDUCE(rhobs_delta(1,1,1),aux_par(1,1,1), N1*N1*NSPIN, DAT_dcomplex,MPI_SUM,0,inverseheads_comm,MPIerror)
      rhobs_delta=aux_par
#endif
#endif

    if(MyNode.ne.0)then
      rhobs_delta=0D0
    endif

    if(MyNode.eq.0) then
        do ispin=1,NSPIN
          trrhos=0D0
          DO II=1,N1
            DO JJ=1,N1
              trrhos=trrhos+rhobs_delta(II,JJ,ispin) * S0_Chain(JJ,II)
            ENDDO
          ENDDO
          if (outinfo) write(12347,*)"tr_psi[rho s]=",trrhos, ispin,ikpmod,itermod
        enddo

    endif


  END SUBROUTINE bs_collectrho

  recursive SUBROUTINE getenebs(ebs,V,Ef_Lead,T,NL,H0_L, H1_L,S0_L,S1_L, NR,H0_R, H1_R,S0_R,S1_R, deltaene,ik, H_chain,S0_chain,S0_chain_inv,N1,bs_tol,bs_min,bs_nmid, nbs,nebs,delta, eizerobs,ezerobs,nbsin,nbsout,nbs2,psibs,psibst, glpsipst,grpsipst)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

#ifdef MPI
    use mMPI_NEGF
#endif

    IMPLICIT NONE

    INTEGER :: NL,NR,N1,nbs,I,J,ik,JJ, nzero,indz,indz2,nebs,nbs2,nbsin,nbsin2,nbsout,i1

    INTEGER bs_nmid, nener

    DOUBLE COMPLEX, DIMENSION (NL,NL) :: H0_L,H1_L,S0_L,S1_L, gamma_l,gamma_lb,gamma_r,gamma_rb
    DOUBLE COMPLEX, DIMENSION (NR,NR) :: H0_R,H1_R,S0_R,S1_R
    DOUBLE COMPLEX, DIMENSION (N1,N1) :: H_chain, S0_chain,S0_chain_inv,veigr, veigr2,veigrb,veigr2b
    DOUBLE COMPLEX zv(N1),zv2(N1), eiene1d(2 * N1)
    DOUBLE PRECISION ebs(2),ebs2(2)
    DOUBLE PRECISION V,deltaene, bs_tol,bs_min,delta
    DOUBLE COMPLEX :: psib1(N1),psib2(N1),psib3(N1),psib4(N1)
    DOUBLE PRECISION :: Ef_Lead,T, ene1d(2 * N1),eimin,releimag
    double complex, pointer :: eizero(:,:)
    double precision, pointer :: ezero(:,:)
    double precision, pointer :: ene1dr(:)
    double complex, pointer :: eiene1dr(:)
    double complex, pointer :: eizerobs(:,:)
    double precision, pointer :: ezerobs(:,:)
    double complex, pointer :: eizerobs2(:,:)
    double precision, pointer :: ezerobs2(:,:)
    DOUBLE COMPLEX,pointer :: psibs(:,:),psibst(:,:),glpsipst(:,:), grpsipst(:,:)
    DOUBLE COMPLEX,pointer :: psibs2(:,:),psibst2(:,:), glpsipst2(:,:),grpsipst2(:,:)

    INTEGER :: MPIerror,Mynode

    interface eiofe
      subroutine eiofe(e,ei,neitot,eizero,ezero,nzero)
        integer neitot,nzero
        double complex ei(neitot)
        double precision e(neitot)
        double complex, pointer :: eizero(:,:)
        double precision, pointer :: ezero(:,:)
      end subroutine eiofe
      subroutine reduceene(eimin, ene1d,eiene1d,nene, ene1dr,eiene1dr,nener)

        INTEGER nene,nener
        double precision ene1d(nene),eimin
        double complex eiene1d(nene)
        double precision, pointer :: ene1dr(:)
        double complex, pointer :: eiene1dr(:)

      end subroutine reduceene
    end interface eiofe

#ifdef MPI
    CALL MPI_COMM_RANK(negf_comm,mynode,MPIerror)
#else
    MyNode=0
#endif
    write(*,*)"in getenebs",ik
    write(*,*)"grebs1=",ebs(1),ebs(1)
    write(*,*)"grebs2=",ebs(2),ebs(2)


    nbsin2=nbsin
    call getebs(ebs(1),v, NL,H0_L, H1_L,S0_L,S1_L, NR,H0_R, H1_R,S0_R,S1_R, N1,H_chain,S0_chain_inv, veigr,veigr2,zv,gamma_l,gamma_r,deltaene )

    do jj=1,N1
      eiene1d(jj)=zv(JJ)
      ene1d(jj)=ebs(1)
    enddo

    call getebs(ebs(2),v, NL,H0_L, H1_L,S0_L,S1_L, NR,H0_R, H1_R,S0_R,S1_R, N1,H_chain,S0_chain_inv, veigrb,veigr2b,zv2,gamma_lb,gamma_rb,deltaene )

    do jj=1,N1
      eiene1d(N1 + jj)=zv2(JJ)
      ene1d(N1 + jj)=ebs(2)
    enddo

!      do jj=1,2 * N1
!        write(*,*)"ene2=",ene1d(jj),ene1d(jj),DREAL(eiene1d(jj))
!      enddo

    if(bs_tol.lt.1d-3)then
      eimin=1d-3
    else
      eimin=bs_tol
    endif
    call reduceene(eimin, ene1d,eiene1d,2 * N1, ene1dr,eiene1dr,nener)
    write(*,*)"nenerbs=",nener
    call f77flush
    call eiofe(ene1dr,eiene1dr,nener,eizero, ezero, nzero)
    write(*,*)"nzerobs=",nzero

    do i=1,nzero

      if(((ABS(DIMAG(eizero(i,1))).le.bs_tol).and. (ABS(DIMAG(eizero(i,1))).ge.bs_min)).NEQV. ((ABS(DIMAG(eizero(i,2))).le.bs_tol).and. (ABS(DIMAG(eizero(i,2))).ge.bs_min)))then

        write(*,*)"calling subroutine recursively"
        write(*,*)"ebs_recursive:",ebs(1),ebs(2), (ebs(1)+ebs(2))* 0.5D0,DREAL(eizero(i,3))
        write(*,*)"ebs_recursive_imag:",DIMAG(eizero(i,1)), DIMAG(eizero(i,2))

        ebs2(1)=ebs(1)
!          ebs2(2)=DREAL(eizero(i,3))
        ebs2(2)=(ebs(1)+ebs(2))* 0.5D0
        call getenebs(ebs2,V,Ef_Lead,T,NL,H0_L, H1_L,S0_L,S1_L, NR,H0_R, H1_R,S0_R,S1_R, deltaene,ik, H_chain,S0_chain,S0_chain_inv,N1,bs_tol,bs_min,bs_nmid, nbs,nebs,delta, eizerobs,ezerobs,nbsin2,nbsout,nbs2,psibs,psibst,glpsipst, grpsipst)
        nbsin2=nbsout

        ebs2(1)=(ebs(1)+ebs(2))* 0.5D0
!          ebs2(1)=DREAL(eizero(i,3))
        ebs2(2)=ebs(2)
        call getenebs(ebs2,V,Ef_Lead,T,NL,H0_L, H1_L,S0_L,S1_L, NR,H0_R, H1_R,S0_R,S1_R, deltaene,ik, H_chain,S0_chain,S0_chain_inv,N1,bs_tol,bs_min,bs_nmid, nbs,nebs,delta, eizerobs,ezerobs,nbsin2,nbsout,nbs2,psibs,psibst,glpsipst, grpsipst)

        return
      endif

    enddo

    do i=1,nzero

      if(((ABS(DIMAG(eizero(i,1))).le.bs_tol).and. (ABS(DIMAG(eizero(i,1))).ge.bs_min)).and. ((ABS(DIMAG(eizero(i,2))).le.bs_tol).and. (ABS(DIMAG(eizero(i,2))).ge.bs_min)))then

        if(ABS(DIMAG(eizero(i,1))).eq.0D0.or. ABS(DIMAG(eizero(i,2))).eq.0D0) then
          write(*,*)"warning, one of the imaginary parts of the eigenvalues is 0."
          if(ABS(DIMAG(eizero(i,1))).eq.0D0.and. ABS(DIMAG(eizero(i,2))).eq.0D0) then
            releimag=1D0
          else
            releimag=0D0
          endif
        else
          releimag=ABS(DIMAG(eizero(i,1)))/ABS(DIMAG(eizero(i,2)))
        endif

        if(releimag.gt.1D1.or.releimag.lt.1D-1)then

          write(*,*)"calling subroutine recursively due to too big difference in eimag"
          write(*,*)"ebs_recursive:",ebs(1),ebs(2), (ebs(1)+ebs(2))* 0.5D0,DREAL(eizero(i,3))
          write(*,*)"ebs_recursive_imag:",DIMAG(eizero(i,1)), DIMAG(eizero(i,2))

          ebs2(1)=ebs(1)
!            ebs2(2)=DREAL(eizero(i,3))
          ebs2(2)=(ebs(1)+ebs(2))* 0.5D0
          call getenebs(ebs2,V,Ef_Lead,T,NL,H0_L, H1_L,S0_L,S1_L, NR,H0_R, H1_R,S0_R,S1_R, deltaene,ik, H_chain,S0_chain,S0_chain_inv,N1,bs_tol,bs_min, bs_nmid, nbs,nebs,delta, eizerobs,ezerobs,nbsin2,nbsout,nbs2,psibs,psibst,glpsipst, grpsipst)
          nbsin2=nbsout

          ebs2(1)=(ebs(1)+ebs(2))* 0.5D0
!            ebs2(1)=DREAL(eizero(i,3))
          ebs2(2)=ebs(2)
          call getenebs(ebs2,V,Ef_Lead,T,NL,H0_L, H1_L,S0_L,S1_L, NR,H0_R, H1_R,S0_R,S1_R, deltaene,ik, H_chain,S0_chain,S0_chain_inv,N1,bs_tol,bs_min, bs_nmid, nbs,nebs,delta, eizerobs,ezerobs,nbsin2,nbsout,nbs2,psibs,psibst,glpsipst, grpsipst)

          return
        endif
      endif

    enddo

    nbs2=0
    do i=1,nzero

      if(((ABS(DIMAG(eizero(i,1))).le.bs_tol).and. (ABS(DIMAG(eizero(i,1))).ge.bs_min)).and. ((ABS(DIMAG(eizero(i,2))).le.bs_tol).and. (ABS(DIMAG(eizero(i,2))).ge.bs_min)))then
        nbs2=nbs2+1
      endif

    enddo
    write(*,*)"numberi of bs for kp=",ik,nbs2

    nbsout=nbsin+nbs2
    if(nbs2.eq.0)then
      return
    elseif(nbsin.eq.0)then
      allocate(eizerobs(nbs2,3),ezerobs(nbs2,2))
      allocate(psibs(nbs2,N1),psibst(nbs2,N1), glpsipst(nbs2,N1),grpsipst(nbs2,N1))
    else
      write(*,*)"ninfo",ik,nbsin,nbs2,nbsout
      call f77flush
      allocate(eizerobs2(nbsin,3),ezerobs2(nbsin,2))
      allocate(psibs2(nbsin,N1),psibst2(nbsin,N1), glpsipst2(nbsin,N1),grpsipst2(nbsin,N1))

      eizerobs2=eizerobs
      ezerobs2=ezerobs
      psibs2=psibs
      psibst2=psibst
      glpsipst2=glpsipst
      grpsipst2=grpsipst
      deallocate(eizerobs,ezerobs,psibs,psibst,glpsipst,grpsipst)

      write(*,*)"ninfo2",ik,nbsin,nbs2,nbsout
      call f77flush
      allocate(eizerobs(nbsin+nbs2,3),ezerobs(nbsin+nbs2,2))
      allocate(psibs(nbsin+nbs2,N1),psibst(nbsin+nbs2,N1) ,glpsipst(nbsin+nbs2,N1),grpsipst(nbsin+nbs2,N1))
      write(*,*)"ninfo3",ik,nbsin,nbs2,nbsout
      call f77flush
      eizerobs(1:nbsin,:)=eizerobs2(1:nbsin,:)
      ezerobs(1:nbsin,:)=ezerobs2(1:nbsin,:)
      psibs(1:nbsin,:)=psibs2(1:nbsin,:)
      psibst(1:nbsin,:)=psibst2(1:nbsin,:)
      glpsipst(1:nbsin,:)=glpsipst2(1:nbsin,:)
      grpsipst(1:nbsin,:)=grpsipst2(1:nbsin,:)
      deallocate(eizerobs2,ezerobs2,psibs2,psibst2)
      deallocate(glpsipst2,grpsipst2)
      write(*,*)"ninfo4",ik,nbsin,nbs2,nbsout
      call f77flush
    endif

    nbs2=0
    do i=1,nzero

      if(((ABS(DIMAG(eizero(i,1))).le.bs_tol).and. (ABS(DIMAG(eizero(i,1))).ge.bs_min)).and. ((ABS(DIMAG(eizero(i,2))).le.bs_tol).and. (ABS(DIMAG(eizero(i,2))).ge.bs_min)))then
        nbs2=nbs2+1
        eizerobs(nbsin+nbs2,:)=eizero(i,:)
        ezerobs(nbsin+nbs2,:)=ezero(i,:)
      endif

    enddo
      write(*,*)"ninfo4enddo",ik,nbsin,nbs2,nbsout
      call f77flush

    do i=nbsin+1,nbsout

        write(*,*)"ebs_1=",ezerobs(i,1),DREAL(eizerobs(i,1)), DIMAG(eizerobs(i,1))
        write(*,*)"ebs_2=",ezerobs(i,2),DREAL(eizerobs(i,2)), DIMAG(eizerobs(i,2))
        write(*,*)"ebs_a=",DREAL(eizerobs(i,3)),DREAL(eizerobs(i,3)), DIMAG(eizerobs(i,3))

    enddo

    do i1=1,nbs2

        i=i1+nbsin

        do j=1,n1
          if(eizerobs(i,1).eq.zv(j))then
            indz=j 
!              exit
            write(*,*)"indzinloop=",indz
          endif
        enddo
        write(*,*)"indz=",indz
        do j=1,n1
          if(eizerobs(i,2).eq.zv2(j))then
            indz2=j 
!              exit
            write(*,*)"indzinloop=",indz2
          endif
        enddo
        write(*,*)"indz2=",indz2

     
        if(ABS(DIMAG(eizerobs(i,1))).lt.ABS(DIMAG(eizerobs(i,2))))then
          write(*,*)"using energy 1"
          psib1(:)=veigr(:,indz)
          psib2(:)=veigr2(:,indz)
        else
          write(*,*)"using energy 2"
          psib1(:)=veigrb(:,indz2)
          psib2(:)=veigr2b(:,indz2)
          gamma_l=gamma_lb
        endif
        psibs(i,:)=psib1
        psibst(i,:)=psib2

        psib3=MATMUL(S0_chain_inv,psib2)
        psib4=0D0
        psib4(1:NL)=MATMUL(gamma_l,psib3(1:NL))
        psib3=MATMUL(S0_chain_inv,psib4)
        glpsipst(i,:)=psib3

        psib3=MATMUL(S0_chain_inv,psib2)
        psib4=0D0
        psib4(N1-NR+1:N1)=MATMUL(gamma_r,psib3(N1-NR+1:N1))
        psib3=MATMUL(S0_chain_inv,psib4)
        grpsipst(i,:)=psib3

    enddo

  end SUBROUTINE getenebs



  SUBROUTINE bs_m0(nenerg_div_nodes,N1,NSPIN,Nenerg_div,ik, nk, NL,NR, V,Ef_Lead,T,H0_L, H1_L,S0_L,S1_L,H0_R, H1_R,S0_R,S1_R, deltaene, H_chain,S0_chain,S0_chain_inv, delta,rhobstot,eiene,ene,indemax,wk )

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

    use negfmod
#ifdef MPI
    use mMPI_NEGF
#endif

    IMPLICIT NONE

!startdefs
    INTEGER nenerg_div_nodes,N1,NSPIN,Nenerg_div,I,J,ik,nebssave, nk,ISPIN,ii,jj,nener,ie,iz,nbstot,nbsnode,nbsrest,ino, nbsin,nbsout,NL,NR,nbs,nbs2,nzero,indemax,indetotal
    DOUBLE COMPLEX  eiene(Nenerg_div,N1,NSPIN), eienea(nenerg_div_nodes,N1,NSPIN)
    DOUBLE PRECISION ene(nenerg_div_nodes)
    double precision,ALLOCATABLE, SAVE:: ebssave(:,:,:,:), ebssaven(:,:,:,:),ebssaveno(:,:,:,:,:)
    DOUBLE COMPLEX,ALLOCATABLE      ::  eiene1d(:)
    DOUBLE PRECISION,ALLOCATABLE      ::  ene1d(:)
    double precision, pointer :: ene1dr(:)
    double complex, pointer :: eiene1dr(:)
    double complex, pointer :: eizero(:,:)
    double precision, pointer :: ezero(:,:)
    DOUBLE PRECISION ebs
    logical  printdebug,printdebug2
    DOUBLE PRECISION kpoint(3)
    integer iebs(NSPIN,nk),iebsn(NSPIN,nk)
    INTEGER, ALLOCATABLE    ::   nbspnode(:,:,:)
    DOUBLE COMPLEX, DIMENSION (N1,N1) :: S0_chain,S0_chain_inv
    DOUBLE COMPLEX  H_chain(N1,N1,NSPIN), rhobstot(N1,N1,NSPIN)
    DOUBLE PRECISION :: Ef_Lead,T,V,eimin, deltaene,delta,wk
    DOUBLE COMPLEX, DIMENSION (NL,NL,NSPIN) :: H0_L,H1_L
    DOUBLE COMPLEX, DIMENSION (NL,NL) :: S0_L,S1_L
    DOUBLE COMPLEX, DIMENSION (NR,NR,NSPIN) :: H0_R,H1_R
    DOUBLE COMPLEX, DIMENSION (NR,NR) :: S0_R,S1_R
    DOUBLE COMPLEX,pointer :: psibs(:,:),psibst(:,:), glpsibst(:,:),grpsibst(:,:)
    double complex, pointer :: eizerobs(:,:)
    double precision, pointer :: ezerobs(:,:)

!enddefs

    INTEGER :: MPIerror,Mynode,Nnodes
#ifdef MPI
    INTEGER, DIMENSION (MPI_STATUS_SIZE) :: istatus
#endif

    interface eiofe
      subroutine eiofe(e,ei,neitot,eizero,ezero,nzero)
        integer neitot,nzero
        double complex ei(neitot)
        double precision e(neitot)
        double complex, pointer :: eizero(:,:)
        double precision, pointer :: ezero(:,:)
      end subroutine eiofe
      subroutine reduceene(eimin, ene1d,eiene1d,nene, ene1dr,eiene1dr,nener)

        INTEGER nene,nener
        double precision ene1d(nene),eimin
        double complex eiene1d(nene)
        double precision, pointer :: ene1dr(:)
        double complex, pointer :: eiene1dr(:)

      end subroutine reduceene
    end interface eiofe





#ifdef MPI
    CALL MPI_COMM_RANK(negf_comm,mynode,MPIerror)
    CALL MPI_COMM_SIZE(negf_comm,Nnodes,MPIerror)
#else
    Nnodes=1
    MyNode=0
#endif
      printdebug=.false.
      printdebug2=.false.
      nebssave=100
      allocate(nbspnode(Nnodes,NSPIN,nk))
!        rhobstot=0D0

      write(*,*)"in getbsmo", MyNode,ik
        IF (MyNode.EQ.0) THEN
          eienea(1:indemax,1:N1,1:NSPIN)= eiene(1:indemax,1:N1,1:NSPIN)
#ifdef MPI
          DO I=1,Nnodes-1
            CALL MPI_RECV(ene(I*indemax+1:(I+1)*indemax), Nenerg_div,DAT_double,I,1, negf_comm,istatus,MPIerror)
            CALL MPI_RECV(eiene(1,1,1), Nenerg_div * N1 * NSPIN,DAT_dcomplex,I,1, negf_comm,istatus,MPIerror)
            eienea(I*indemax+1:(I+1)*indemax,1:N1,1:NSPIN)= eiene(1:indemax,1:N1,1:NSPIN)
          ENDDO
        ELSE
          CALL MPI_SEND(ene(1:indemax),Nenerg_div,DAT_double,0,1, negf_comm,MPIerror)
          CALL MPI_SEND(eiene(1,1,1),Nenerg_div * N1 * NSPIN, DAT_dcomplex,0,1, negf_comm,MPIerror)
#endif
        ENDIF
        indetotal=indemax * Nnodes


!          IF (MyNode.EQ.0) THEN
!          do i=1,indetotal
!            DO ISPIN=1,NSPIN
!              write(12346,'(e," ")',ADVANCE='NO')ene(i)
!              do izv=1,N1
!                write(12346,'(" ",e," ")',ADVANCE='NO')
!     .             DREAL(eienea(i,izv,ISPIN))
!              enddo
!              write(12346,'(" ")')
!            enddo
!          enddo
!          endif

        if(ik.eq.1)then
          if(allocated(ebssaven)) DEALLOCATE(ebssaven,ebssaveno)
          allocate(ebssaven(nebssave,2,NSPIN,nk))
          allocate(ebssaveno(Nnodes,nebssave,2,NSPIN,nk))
        endif

        IF (MyNode.EQ.0)then

          if(ik.eq.1)then
            if(allocated(ebssave)) DEALLOCATE(ebssave)
            allocate(ebssave(nebssave,2,NSPIN,nk))
          endif


          allocate(eiene1d(indetotal * N1))
          allocate(ene1d(indetotal * N1))
          DO ISPIN=1,NSPIN
            do ii=1,indetotal
              do jj=1,N1
                eiene1d((ii-1) * N1+jj)=eienea(ii,jj,ISPIN)
                ene1d((ii-1) * N1+jj)=ene(ii)
              enddo
            enddo

            if(bs_tol.lt.1d-3)then
              eimin=5d-3
            else
              eimin=5d0 * bs_tol
            endif
            call reduceene(eimin, ene1d,eiene1d,indetotal * N1, ene1dr,eiene1dr,nener)
            write(*,*)"nenerng,total=",nener,indetotal * N1


            call eiofe(ene1dr,eiene1dr,nener,eizero, ezero, nzero)
            deallocate(ene1dr,eiene1dr)
            write(*,*)"nzero=",nzero

            if(nzero.lt.100 * nebssave)then
            ie=0
            do iz=1,nzero

            if (outinfo) write(12347,'("bsinfo: ", 2(e12.5," ")," ik=",i7," ", 3(e12.5," "), " ispin=",i7," e1: ", 3(e12.5," ")," e2: ", 3(e12.5," ")," ei: ",2(e12.5," "))') DREAL(eizero(iz,3)),DREAL(eizero(iz,3)),ik, kpoint(1),kpoint(2),wk,ispin,ezero(iz,1), DREAL(eizero(iz,1)),-DIMAG(eizero(iz,1)), ezero(iz,2), DREAL(eizero(iz,2)),-DIMAG(eizero(iz,2)), DREAL(eizero(iz,3)),-DIMAG(eizero(iz,3))

              if(((ABS(DIMAG(eizero(iz,1))).le.bs_tol).and. (ABS(DIMAG(eizero(iz,1))).ge.bs_min)).or. ((ABS(DIMAG(eizero(iz,2))).le.bs_tol).and. (ABS(DIMAG(eizero(iz,2))).ge.bs_min)))then

                ebs=DREAL(eizero(iz,3))
                if (outinfo) write(12347,'("ebsinfo: ", 2(e12.5," ")," ik=",i7," ", 3(e12.5," "), " ispin=",i7," e1: ", 3(e12.5," ")," e2: ", 3(e12.5," ")," ei: ",2(e12.5," "))') DREAL(eizero(iz,3)),DREAL(eizero(iz,3)),ik, kpoint(1),kpoint(2),wk,ispin,ezero(iz,1), DREAL(eizero(iz,1)),-DIMAG(eizero(iz,1)), ezero(iz,2), DREAL(eizero(iz,2)),-DIMAG(eizero(iz,2)), DREAL(eizero(iz,3)),-DIMAG(eizero(iz,3))

!! --> deallocate eizero and ezero!
                if(ie.eq.0)then
                  ie=ie+1
                  ebssave(ie,:,ISPIN,ik)=ezero(iz,:)
                else
                  if(ebssave(ie,1,ISPIN,ik).ne.ezero(iz,1))then
                    ie=ie+1
                    ebssave(ie,:,ISPIN,ik)=ezero(iz,:)
                  else
                      if (outinfo) write(12347,*)"multiple ei in energy interval"
                  endif
                endif


              endif
            enddo
            iebs(ISPIN,ik)=ie
            
            deallocate(eizero,ezero)
            endif
          
          ENDDO
          deallocate(eiene1d,ene1d)

          if (outinfo) then
          if(nspin.eq.2)then
            write(12347,*)"nbsall:",ik,iebs(1,ik),iebs(2,ik)
          else
            write(12347,*)"nbsall:",ik,iebs(1,ik)
          endif
          endif

          if (outinfo) then
          do ISPIN=1,NSPIN
            do j=1,iebs(ISPIN,ik)
              write(12347,*) "ebs1:=",ebssave(j,:,ISPIN,ik),ISPIN,ik
            enddo
          enddo
          endif


        endif

#ifdef MPI
        CALL MPI_BARRIER( negf_comm, MPIerror )
#endif

        if(bsrun.eq.2)then
          write(*,*)"bsrun=",bsrun
#ifdef MPI
          call MPI_Finalize( MPIerror )
#endif
          stop
        endif


!startnew


!          iebs(1,ik)=4
!          iebs(2,ik)=2
        IF (MyNode.EQ.0) THEN
          nbstot=0D0
          do ispin=1,nspin
            nbstot=nbstot+iebs(ispin,ik)
          enddo
          nbsnode=nbstot/Nnodes
          nbsrest=nbstot-nbsnode * Nnodes 
          write(*,*)"nbstot=",nbstot,nbsnode,nbsrest
          nbspnode(:,:,ik)=0
          ino=1
          do ISPIN=1,NSPIN
            do j=1,iebs(ISPIN,ik)

              write(*,*) "ebssapar:=",ebssave(j,:,ISPIN,ik),ISPIN,ik
              write(*,*)"adding to node",ino

              nbspnode(ino,ispin,ik)=nbspnode(ino,ispin,ik)+1
              ebssaveno(ino,nbspnode(ino,ispin,ik),:,ISPIN,ik)= ebssave(j,:,ISPIN,ik)
              ino=MODULO(ino,Nnodes)+1
            enddo
          enddo
          do ino=1,Nnodes
            if(nspin.eq.2)then
              write(*,*)"nnodebs=",ino,nbspnode(ino,1,ik), nbspnode(ino,2,ik)
            else
              write(*,*)"nnodebs=",ino,nbspnode(ino,1,ik)
            endif
            do ispin=1,NSPIN
              do j=1,nbspnode(ino,ISPIN,ik)
                write(*,*) "nebssapar:=",ino, ebssaveno(ino,j,:,ISPIN,ik),ISPIN,ik
                
              enddo
            enddo
          enddo
        endif


#ifdef MPI
        call MPI_Bcast(nbspnode(1,1,1),Nnodes * NSPIN * nk, MPI_integer,0,negf_comm,MPIerror)
        call MPI_Bcast(ebssaveno(1,1,1,1,1), Nnodes * nebssave * 2 * NSPIN * nk, DAT_double,0,negf_comm,MPIerror)
#endif

        do ino=1,Nnodes
          if(nspin.eq.2)then
            write(*,*)"nno=",MyNode,"sep:",ino,nbspnode(ino,1,ik), nbspnode(ino,2,ik)
          else
            write(*,*)"nno=",MyNode,"sep:",ino,nbspnode(ino,1,ik)
          endif
          do ispin=1,NSPIN
            do j=1,nbspnode(ino,ISPIN,ik)
              write(*,*) "nebssapaar:=",mynode,ino, ebssaveno(ino,j,:,ISPIN,ik),ISPIN,ik
              
            enddo
          enddo
        enddo
        iebsn(:,ik)=nbspnode(MyNode+1,:,ik)
        ebssaven(:,:,:,ik)=ebssaveno(MyNode+1,:,:,:,ik)

        write(*,*)"iebsn=",MyNode,iebsn(:,ik)

        do ispin=1,NSPIN
          do j=1,iebsn(ISPIN,ik)
            write(*,*) "nbsinfop:=",MyNode, ebssaven(j,:,ISPIN,ik),ISPIN,ik
          enddo
        enddo

!      if(ik.eq.2)then
!#ifdef MPI
!      call MPI_Finalize( MPIerror )
!#endif
!      stop
!      endif



        




          if (outinfo) then
          if(nspin.eq.2)then
            write(12347,*)"nbsnode:",ik,iebsn(1,ik),iebsn(2,ik),mynode
          else
            write(12347,*)"nbsnode:",ik,iebsn(1,ik),mynode
          endif
          write(12347,*)"enei_bss="
      endif
          DO ISPIN=1,NSPIN
            nbsin=0
            nbsout=0
            if(iebsn(ISPIN,ik).eq.0 .and. outinfo)write(12347,*)"enei_bss=nostate"
            do ie=1,iebsn(ISPIN,ik)
              write(12347,*) "lebs_saved:=",ebssaven(ie,:,ISPIN,ik),ISPIN,ik

              nbsin=nbsout
              write(*,*)"ie_nbs=",ie,nbsin
              CALL getenebs(ebssaven(ie,:,ISPIN,ik), V,Ef_Lead,T,NL,H0_L(:,:,ISPIN), H1_L(:,:,ISPIN),S0_L,S1_L, NR,H0_R(:,:,ISPIN), H1_R(:,:,ISPIN),S0_R,S1_R, deltaene,ik, H_chain(:,:,ISPIN),S0_chain,S0_chain_inv,N1, bs_tol, bs_min,bs_nmid,nbs,iebsn(ISPIN,ik),Delta, eizerobs,ezerobs,nbsin,nbsout,nbs2,psibs,psibst, glpsibst,grpsibst)
              write(*,*)"ie_nbsout=",ik,ISPIN,iebsn(ISPIN,ik),nbsout
                

            enddo

            if (outinfo) write(12347,*)"nbs_r=",nbsout
!              do i=1,nbsout
!                write(12347,*)"ebsn_1=",ezerobs(i,1),
! DREAL(eizerobs(i,1)),DIMAG(eizerobs(i,1))
!                write(12347,*)"ebsn_2=",ezerobs(i,2),
! DREAL(eizerobs(i,2)), DIMAG(eizerobs(i,2))
!                write(12347,*)"ebsn_a=",DREAL(eizerobs(i,3)),
! DREAL(eizerobs(i,3)), DIMAG(eizerobs(i,3))
!
!                psi1=MATMUL(S0_chain_inv,psibst(i,:))
!                normvt=DOT_PRODUCT(psibst(i,:),psi1)
!
!                zevL=DOT_PRODUCT(psibst(i,:),glpsibst(i,:))
!                zevR=DOT_PRODUCT(psibst(i,:),grpsibst(i,:))
!
!                zevL=-0.5D0 * zevL/normvt
!                zevR=-0.5D0 * zevR/normvt
!                zev=zevL+zevR
!                write(*,*)"zev=",zev,zevL,zevR
!                write(12347,*)"zev=",zev,zevL,zevR
!                write(12347,*)"zbs=",DREAL(eizerobs(i,3)),
! DREAL(eizerobs(i,3)), DIMAG(eizerobs(i,3)),
! DREAL(zev), DIMAG(zev),
! DREAL(zevL), DIMAG(zevL),
! DREAL(zevR), DIMAG(zevR)
!
!              enddo

            if(nbsout.gt.0)then
                call getrhobs(V,Ef_Lead,T,ik,N1, S0_chain,S0_chain_inv,rhobstot(:,:,ISPIN), Delta,eizerobs,ezerobs,nbsout,psibs,psibst, glpsibst,grpsibst,ispin)
            endif


            if(nbsout.gt.0) deallocate(eizerobs,ezerobs,psibs,psibst, glpsibst,grpsibst)


          ENDDO
          

  end SUBROUTINE bs_m0



  SUBROUTINE getebs(ebs,v, NL,H0_L, H1_L,S0_L,S1_L, NR,H0_R, H1_R,S0_R,S1_R, N1,H_chain,S0_chain_inv, veig,veigt,zv,gamma_l,gamma_r,deltaimag )

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************
    use mSelfenergies, only : SelfEnergyGeneral

    IMPLICIT NONE
    INTEGER :: NL,NR,N1,nrchanL,nrchanR,II,ilo,ihi,INFO

    INTEGER IPIV(N1)

    DOUBLE COMPLEX, DIMENSION (NL,NL) :: H0_L,H1_L,S0_L,S1_L,Sigma_L,gamma_l
    DOUBLE COMPLEX, DIMENSION (NR,NR) :: H0_R,H1_R,S0_R,S1_R,Sigma_R, gamma_r
    DOUBLE COMPLEX, DIMENSION (N1,N1) :: Sigma_aux,H_chain, S0_chain_inv,veig,veigt
    DOUBLE COMPLEX Ei,zv(N1),work(2 * N1 * N1), work2(N1**2),zvimag(N1),normvt
    DOUBLE PRECISION ebs
    DOUBLE PRECISION V,scaleev(N1),abnrm,rconde(N1), rcondv(N1),rwork(2 * N1),deltaimag
    DOUBLE COMPLEX :: psib1(N1)
    DOUBLE COMPLEX :: vbufL(NL),vbufR(NR)
    DOUBLE COMPLEX :: zi
    PARAMETER(zi=(0.0d0,1.0d0))


    Ei=ebs-v * 0.5D0 
    call SelfEnergyGeneral('L',NL,Ei,H0_L, H1_L,S0_L,S1_L, Sigma_L,nrchanL, deltaimag,.true.)
    Ei=ebs+v * 0.5D0 
    call SelfEnergyGeneral('R',NR,Ei,H0_R, H1_R,S0_R,S1_R, Sigma_R,nrchanR, deltaimag,.true.)

    Sigma_aux=(0.D0,0.D0)
    Sigma_aux(1:NL,1:NL)=Sigma_L(:,:)
    Sigma_aux(N1-NR+1:N1,N1-NR+1:N1)= Sigma_R(:,:)
    gamma_l=zi * (Sigma_L-DCONJG(TRANSPOSE(Sigma_L)))
    gamma_r=zi * (Sigma_R-DCONJG(TRANSPOSE(Sigma_R)))

    Sigma_aux=H_Chain(:,:)+Sigma_aux
!      S0_chain_inv=0D0
!      DO II=1,N1
!      S0_chain_inv(II,II)=1D0
!      ENDDO
    Sigma_aux=matmul(S0_Chain_inv,Sigma_aux)

    call ZGEEVX( 'B', 'N', 'V','N', N1, Sigma_aux, N1, zv, veigt, N1, veig, N1,ilo,ihi,scaleev,abnrm, rconde,rcondv, WORK, 2 * N1 * N1 , RWORK, INFO )

    veigt=veig
    CALL ZGETRF(N1,N1,veigt,N1,IPIV,INFO)
    CALL ZGETRI(N1,veigt,N1,IPIV,WORK2,N1**2,INFO)

    veigt=DCONJG(TRANSPOSE(veigt))
 
    if(.false.)then

      DO II=1,N1

        psib1=MATMUL(S0_chain_inv, veigt(:,II))

        vbufL=MATMUL(gamma_l,psib1(1:NL))
        vbufR=MATMUL(gamma_r,psib1(N1-NR+1:N1))
        zvimag(II)=DOT_PRODUCT(psib1(1:NL),vbufL)+ DOT_PRODUCT(psib1(N1-NR+1:N1),vbufR)
        normvt=DOT_PRODUCT(veigt(:,II),psib1)
        zvimag(II)=-0.5D0 * zvimag(II)/normvt
        write(*,*)"zvgamma",zv(II),zvimag(II)
        write(*,*)"relat",DIMAG(zv(II))/DREAL(zvimag(II))

      ENDDO
      call f77flush

    endif

  end SUBROUTINE getebs


  SUBROUTINE getrhobs(V,Ef_Lead,T,ik,N1, S0_chain,S0_chain_inv,rhobs2, delta,eizerobs,ezerobs,nbs,psibs,psibst,glpsipst,grpsipst, ispin)

! *****************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! *****************************************************************

    use negfmod
#ifdef MPI
    use mMPI_NEGF
#endif

    IMPLICIT NONE

    INTEGER :: N1,nbs,I,J,ik,II,JJ,indmaxn,indmaxm, infoeq,ispin
    DOUBLE PRECISION alpha
    DOUBLE COMPLEX   alphasc,pn,pm

    DOUBLE COMPLEX, DIMENSION (N1,N1) :: rhobs, S0_chain,S0_chain_inv,rhobs2
    DOUBLE PRECISION dei_den,dei_dem,de,dei_de
    DOUBLE PRECISION V,delta
    DOUBLE COMPLEX :: ppt,trrhos,zi,wmesh,wdelta
    DOUBLE PRECISION :: fermi_aux,Ef_Lead,T,fL,fR,cbuf

    double complex  :: eizerobs(nbs,3)
    double precision :: ezerobs(nbs,2)
    DOUBLE COMPLEX :: psibs(nbs,N1),psibst(nbs,N1), glpsipst(nbs,N1),grpsipst(nbs,N1)
    DOUBLE COMPLEX :: psin(N1),psint(N1),psib3(N1)
    DOUBLE COMPLEX :: psim(N1),psimt(N1),normvt,zevL,zevR,zev
    DOUBLE PRECISION :: gammaL,gammaR
    character fname*33
    character*20  charnode

    INTEGER :: MPIerror,Mynode
    PARAMETER(zi=(0.0d0,1.0d0))

#ifdef MPI
    CALL MPI_COMM_RANK(negf_comm,mynode,MPIerror)
#else
    MyNode=0
#endif
    write(*,*)"in getrhobs",ik

    if (outinfo) write(12347,*)" Number of k-points =        1 Number of Spins =        1 Number of basis orbs =",     n1

    If (itermod .EQ. 1) THEN
      write(charnode,'(i10)')MyNode
      charnode="wfs_" // TRIM(ADJUSTL(charnode))

      if(ispin.eq.1)then
        fname=TRIM(ADJUSTL(charnode))// ".up"
      else
        fname=TRIM(ADJUSTL(charnode))// ".down"
      endif
    endif

    If ((itermod .EQ. 1).and.ik.eq.1) THEN
      open(12348, file=fname, form='unformatted', status='unknown' )
      write(12348) 1
      write(12348) 2
      write(12348) n1

      endfile (12348)
      backspace (12348)
      close (12348)
    endif

    if (outinfo) then
    do i=1,nbs
        write(12347,*)"ebsg_1=",ezerobs(i,1),DREAL(eizerobs(i,1)), DIMAG(eizerobs(i,1))
        write(12347,*)"ebsg_2=",ezerobs(i,2),DREAL(eizerobs(i,2)), DIMAG(eizerobs(i,2))
        write(12347,*)"ebsg_a=",DREAL(eizerobs(i,3)), DREAL(eizerobs(i,3)),DIMAG(eizerobs(i,3))
        write(*,*)"ebsg_a=",DREAL(eizerobs(i,3)), DREAL(eizerobs(i,3)),DIMAG(eizerobs(i,3))
        write(12347,*)"eibsg_a=",DREAL(eizerobs(i,3)), DREAL(eizerobs(i,3)),DIMAG(eizerobs(i,3)), ezerobs(i,1),DREAL(eizerobs(i,1)), DIMAG(eizerobs(i,1)), ezerobs(i,2),DREAL(eizerobs(i,2)), DIMAG(eizerobs(i,2)), ispin,ik,itermod
    enddo

    if(ik.eq.1.and.itermod.eq.1)then
    write(12347,*)"start of wf output",ispin,ik,itermod
    do i=1,nbs
        write(12347,*)"ei_wf=",DREAL(eizerobs(i,3)), DREAL(eizerobs(i,3)),DIMAG(eizerobs(i,3)), ezerobs(i,1),DREAL(eizerobs(i,1)), DIMAG(eizerobs(i,1)), ezerobs(i,2),DREAL(eizerobs(i,2)), DIMAG(eizerobs(i,2)), ispin,ik,itermod

!         Printout of wave-function



        write(12347,*)
    write(12347,'(a72)')     ' *************************************** ********************************'
        write(12347,'(a22,2x,i6,2x,3f10.6)') 'k-point = ',ik, kpointmod(1),kpointmod(2),kpointmod(3)
        write(12347,'(a22,2x,i6)') 'Spin component = ',ispin
        write(12347,'(a22,2x,i6)') 'Num. wavefunctions = ',nbs

        psin=psibs(i,:)
        psint=psibst(i,:)

        call writephiascii(N1,psin,psint,eizerobs(i,3),i)

    enddo
    write(12347,*)"end of wf output",ispin,ik,itermod
endif

 
      if(itermod.eq.1.and.ik.eq.1)then ! xxx: in principle this can be output for all k-points
        open(12348, file=fname, form='unformatted', position='append', status='old' )
        write(12348) ik,kpointmod(1),kpointmod(2),kpointmod(3)
        write(12348) 1
        write(12348) nbs

          do i=1,nbs

          psin=psibs(i,:)
          psint=psibst(i,:)

          call writephibin(N1,psin,eizerobs(i,3),i)

          enddo

        write(12348) ik,kpointmod(1),kpointmod(2),kpointmod(3)
        write(12348) 2
        write(12348) nbs
      endif

      do i=1,nbs

        psin=psibs(i,:)
        psint=psibst(i,:)

        call writephibin(N1,psin,eizerobs(i,3),i)

      enddo

    close(12348)

    endif
    
    do i=1,nbs
    if (outinfo) write(12347,*)"ei=",DREAL(eizerobs(i,3))
      do j=1,nbs
        if(j.ne.i)then
            if (outinfo) write(12347,*)"cycling loop",i,j
          cycle
        else
            if (outinfo) write(12347,*)"j=i",i,j
        endif
        de=DREAL(eizerobs(i,3))-DREAL(eizerobs(j,3))
        if(ABS(de).lt.10D0 * delta)then
          if(i.ne.j .and. outinfo) write(12347,*)"warning: closely spaced levels, de=",DREAL(eizerobs(i,3)),abs(de),delta,i,j
        endif

        dei_den=(DREAL(eizerobs(i,2))-DREAL(eizerobs(i,1)))/ (ezerobs(i,2)-ezerobs(i,1))
        dei_dem=(DREAL(eizerobs(j,2))-DREAL(eizerobs(j,1)))/ (ezerobs(j,2)-ezerobs(j,1))
        dei_den=ABS(1D0-dei_den)
        dei_dem=ABS(1D0-dei_dem)
        dei_de=SQRT(dei_den * dei_dem)
!          dei_de=1D0
if (outinfo) write(12347,*)"|1-dei_de|=",dei_de,dei_den,dei_dem

        if(.true.)then
          wmesh=(1D0/ABS(dei_den)+1D0/ABS(dei_dem)) * delta/ (zi * (eizerobs(i,3)-DCONJG(eizerobs(j,3)))+2D0 * delta)
        elseif(.false.)then
          wmesh= 2D0 * zi * delta/ (DREAL(eizerobs(j,3))-DREAL(eizerobs(i,3))+ zi * ((-DIMAG(eizerobs(j,3))+delta)/ dei_dem+ (-DIMAG(eizerobs(i,3))+delta)/ dei_den))
          wmesh= wmesh/(dei_den * dei_dem)

!         old wrong wmesh:
!            wmesh= 2D0 * zi * delta/
! (dei_den * DCONJG(eizerobs(j,3))-
! dei_dem * eizerobs(i,3)+
! zi * delta * (dei_den+dei_dem))
        else
          wmesh=1D0
        endif
        if (outinfo) write(12347,*)"wmesh=",ABS(wmesh),wmesh

        psin=psibs(i,:)
        psint=psibst(i,:)
        psim=psibs(j,:)
        psimt=psibst(j,:)

        psib3=0D0
        DO II=1,N1
          DO JJ=1,N1
            psib3(II)=psib3(ii)+ S0_Chain_inv(ii,jj) * psimt(jj)
          ENDDO
        ENDDO
        ppt=0D0
        do II=1,N1
          ppt=ppt+DCONJG(psint(II)) * psib3(II)
        enddo

        if (outinfo) then
        if(i.eq.j)write(12347,*)"ppte=",i,j,DREAL(ppt),DIMAG(ppt)
        if(i.ne.j)write(12347,*)"pptn=",i,j,DREAL(ppt),DIMAG(ppt)
    endif

        rhobs=0D0
        DO II=1,N1
          DO JJ=1,N1
            rhobs(II,JJ)=ppt * psin(ii) * DCONJG(psim(jj))
          ENDDO
        ENDDO

        fermi_aux=(DREAL(eizerobs(i,3))-Ef_Lead-V/2.D0)/T
        fL=1.D0/(1.D0+DEXP(fermi_aux))
        fermi_aux=(DREAL(eizerobs(i,3))-Ef_Lead+V/2.D0)/T
        fR=1.D0/(1.D0+DEXP(fermi_aux))

        
        cbuf=Maxval(ABS(psin))
        do II=1,N1
          if(cbuf.eq.ABS(psin(II)))indmaxn=II
        enddo
 
        cbuf=Maxval(ABS(psim))
        do II=1,N1
          if(cbuf.eq.ABS(psim(II)))indmaxm=II
        enddo
        if (outinfo) write(12347,*)"maxval=",indmaxn,indmaxm

        psib3=MATMUL(S0_chain_inv,psimt)
        normvt=DOT_PRODUCT(psint,psib3)

        zevL=DOT_PRODUCT(psint,glpsipst(i,:))
        zevR=DOT_PRODUCT(psint,grpsipst(i,:))

        zevL=-0.5D0 * zevL/normvt
        zevR=-0.5D0 * zevR/normvt
        zev=zevL+zevR
        if (outinfo) write(12347,*)"zevinbs=",zev,zevL,zevR

        gammaL=2D0 * ABS(-1D0 * DREAL(zevL))
        gammaR=2D0 * ABS(-1D0 * DREAL(zevR))


        alphasc=0.5D0 * (gammaL-gammaR)/(gammaL+gammaR)

        alpha=-0.5D0
        if(indmaxn.le.bs_nmid) then
          alpha=alpha+0.5D0
        endif
        if(indmaxm.le.bs_nmid) then
          alpha=alpha+0.5D0
        endif
        if (outinfo) write(12347,*)"alphac=",alpha-DREAL(alphasc),alpha, DREAL(alphasc),DREAL(zevL),DREAL(zevR)

        write(*,*)"bssc=",bssc
        if(bssc.eq.0)then
          alpha=DREAL(alphasc)
        elseif(bssc.eq.2)then
          alpha=0D0
        endif

        trrhos=0D0
        DO II=1,N1
          DO JJ=1,N1
            trrhos=trrhos+rhobs(II,JJ) * S0_Chain(JJ,II)
          ENDDO
        ENDDO

        infoeq=0
        if(i.ne.j)infoeq=1
        
        if (outinfo) write(12347,*)"allinfobs=",DREAL(eizerobs(i,3)),infoeq, DREAL(eizerobs(i,3)), DIMAG(eizerobs(i,3)), DREAL(eizerobs(j,3)), DIMAG(eizerobs(j,3)), DREAL(zev), DIMAG(zev), DREAL(zevL), DIMAG(zevL), DREAL(zevR), DIMAG(zevR),alpha,DREAL(alphasc), DIMAG(alphasc), trrhos * alpha * (fL-fR) * wmesh,trrhos,ppt,wmesh, dei_de,dei_den,dei_dem,(fL-fR), vmod,itermod, ikpmod,kpointmod(1),kpointmod(2)



!          gammaL=-2D0 * DREAL(zevL)
!          gammaR=-2D0 * DREAL(zevR)
        pn=0.5D0 * (fl+fr) + alphasc * (fl - fr)
        pm=0.5D0 * (fl+fr) + alpha * (fl - fr)

        wdelta=delta/ (zi * (eizerobs(i,3)-DCONJG(eizerobs(j,3)))+2D0 * delta)

        if (outinfo ) write(12347,*)"bsgamma=",DREAL(eizerobs(i,3)), DREAL(eizerobs(i,3)), DIMAG(eizerobs(i,3)), DREAL(eizerobs(j,3)), DIMAG(eizerobs(j,3)), gammaL,gammaR,gammaL/(gammaL+gammaR), DREAL(alphasc),DIMAG(alphasc),(fL-fR), DREAL(wmesh),DIMAG(wmesh), DREAL((fL-fR) * wmesh * gammaL/(gammaL+gammaR)), DIMAG((fL-fR) * wmesh * gammaL/(gammaL+gammaR)), DREAL((fL-fR) * wmesh * gammaR/(gammaL+gammaR)), DIMAG((fL-fR) * wmesh * gammaR/(gammaL+gammaR)), DREAL((fL-fR) * wmesh), DIMAG((fL-fR) * wmesh), DREAL(pn * wmesh),DIMAG(pn * wmesh), DREAL(wdelta),DIMAG(wdelta), DREAL(pn),DIMAG(pn), DREAL(pm),DIMAG(pm), vmod,itermod, ikpmod,kpointmod(1),kpointmod(2), dei_de,dei_den,dei_dem,ispin

!0          write(12347,*)"ing= 1  ",DREAL(eizerobs(i,3))
!          write(12347,*)"ing= 2  ",DREAL(eizerobs(i,3))
!1          write(12347,*)"ing= 3  ",DIMAG(eizerobs(i,3))
!          write(12347,*)"ing= 4  ",DREAL(eizerobs(j,3))
!          write(12347,*)"ing= 5  ",DIMAG(eizerobs(j,3))
!2          write(12347,*)"ing= 6  ",gammaL
!3          write(12347,*)"ing= 7  ",gammaR
!4          write(12347,*)"ing= 8  ",gammaL/(gammaL+gammaR)
!5          write(12347,*)"ing= 9  ",DREAL(alphasc)
!          write(12347,*)"ing= 10 ",DIMAG(alphasc)
!6          write(12347,*)"ing= 11 ",(fL-fR)
!7          write(12347,*)"ing= 12 ",DREAL(wmesh)
!          write(12347,*)"ing= 13 ",DIMAG(wmesh)
!          write(12347,*)"ing= 14 ",DREAL((fL-fR) * wmesh * gammaL/(gammaL+gammaR))
!          write(12347,*)"ing= 15 ",DIMAG((fL-fR) * wmesh * gammaL/(gammaL+gammaR))
!          write(12347,*)"ing= 16 ",DREAL((fL-fR) * wmesh * gammaR/(gammaL+gammaR))
!          write(12347,*)"ing= 17 ",DIMAG((fL-fR) * wmesh * gammaR/(gammaL+gammaR))
!          write(12347,*)"ing= 18 ",DREAL((fL-fR) * wmesh)
!          write(12347,*)"ing= 19 ",DIMAG((fL-fR) * wmesh)
!          write(12347,*)"ing= 20 ",DREAL(pn * wmesh)
!          write(12347,*)"ing= 21 ",DIMAG(pn * wmesh)
!8          write(12347,*)"ing= 22 ",DREAL(wdelta)
!          write(12347,*)"ing= 23 ",DIMAG(wdelta)
!9          write(12347,*)"ing= 24 ",DREAL(pn)
!          write(12347,*)"ing= 25 ",DIMAG(pn)
!10          write(12347,*)"ing= 26 ",DREAL(pm)
!          write(12347,*)"ing= 27 ",DIMAG(pm)
!x          write(12347,*)"ing= 28 ",vmod
!          write(12347,*)"ing= 29 ",itermod 
!          write(12347,*)"ing= 30 ",ikpmod
!          write(12347,*)"ing= 31 ",kpointmod(1)
!          write(12347,*)"ing= 32 ",kpointmod(2)
!to add
!          write(12347,*)"ing= 33 ",dei_den
!          write(12347,*)"ing= 34 ",dei_dem
!          write(12347,*)"ing= 33 ",1D0/abs(dei_den)
!          write(12347,*)"ing= 34 ",1D0/abs(dei_dem)



        if(DREAL(trrhos).lt.10.0D0) then
          rhobs2=rhobs2+alpha * (fL-fR) * rhobs * wmesh
        else
            if (outinfo) write(12347,*)"warning: tr[rho] is too big"
        endif
      enddo

    enddo

  end SUBROUTINE getrhobs

  subroutine transm_probes_tot(ef_bss,tij,wk,nleads,nbss,nspin, Nenerg_div,ik,nk,v,t,ef_lead)
  use negfmod, only: outinfo

    IMPLICIT NONE

    integer nleads,nbss,nspin,Nenerg_div,ik,nk,i1
    double precision tij(Nenerg_div,nleads,nleads,nspin)
    double precision V,ef_bss(nleads,nspin),ef_lead,t,wk
    double precision, allocatable, save :: tijallk(:,:,:,:)
    double precision, allocatable, save :: tijallktot(:,:,:,:)

    if(ik.eq.1)then
      allocate(tijallk(Nenerg_div,nleads,nleads,nspin))
      allocate(tijallktot(Nenerg_div,nleads,nleads,nspin))
      tijallk=0D0
    endif

    tijallk=tijallk+tij * wk

    if(ik.eq.nk)then

        if (outinfo) write(12347,*)"ef_total"
      
      if(nspin.eq.2)then
        tijallktot(:,:,:,1)=tijallk(:,:,:,1)+tijallk(:,:,:,2)
        tijallktot(:,:,:,2)=tijallk(:,:,:,1)+tijallk(:,:,:,2)
      else
        tijallktot(:,:,:,1)=tijallk(:,:,:,1)
      endif

      call find_ef_bss(tijallktot,nspin,nleads,nbss,Nenerg_div, ef_bss,ef_lead,t,v)

      deallocate(tijallk,tijallktot)
    endif
    
  end subroutine transm_probes_tot



  subroutine bs_transm(ef_bssk,ef_bss,nspin,Nenerg_div,V, ik,nk,nleads,nbss,ef_lead,t, nebss,tij,wk)

    use mEnergyGrid
    use negfmod,only:outinfo
#ifdef MPI
    use mMPI_NEGF
#endif
    
    IMPLICIT NONE

    integer nspin,Nenerg_div,ik,nl,nr,nk
    double precision V,delta,ef_bssk(nleads,nspin), ef_bss(nleads,nspin),ef_lead,t,deltabss(nleads),wk
    double precision iclocmat(nleads,nleads,nspin), ictotmat(nleads,nleads,nspin),ictot(nleads,nspin)
    double precision, allocatable, save :: ictotmatallk(:,:,:),ictotallk(:,:)
    integer nebss(nleads,2)

    integer ispin,i,j,i1,i2,i3,i4,MPIerror,mynode,Nnodes
    integer nbss,ibs,nleads,nbsmax
    double complex ei,fi,fj
    DOUBLE COMPLEX, PARAMETER :: zi=(0.D0,1.D0)
    double precision, allocatable :: ti(:,:,:)
    double precision tij(Nenerg_div,nleads,nleads,nspin)
    double precision, allocatable, save :: tijallk(:,:,:,:), tijallk_cur(:,:,:,:),tiallk_cur(:,:,:)
    double precision, allocatable :: tij_cur(:,:,:,:),ti_cur(:,:,:)
    logical writetrc
    DOUBLE PRECISION, PARAMETER :: eh_const=1.6022D-4*13.6058D0/4.1357D0

#ifdef MPI
    CALL MPI_COMM_SIZE(negf_comm,Nnodes,MPIerror)
    CALL MPI_COMM_RANK(negf_comm,mynode,MPIerror)
#else
    MyNode=0
    Nnodes=1
#endif
 
    allocate(ti(Nenerg_div,nleads,nspin))
    allocate(tij_cur(Nenerg_div,nleads,nleads,nspin), ti_cur(Nenerg_div,nleads,nspin))

    if(ik.eq.1)then
      allocate(tijallk(Nenerg_div,nleads,nleads,nspin))
      allocate(tijallk_cur(Nenerg_div,nleads,nleads,nspin))
      allocate(tiallk_cur(Nenerg_div,nleads,nspin))
      tijallk=0D0
      tijallk_cur=0D0
      tiallk_cur=0D0
      if(MyNode.eq.0)then
        allocate( ictotmatallk(nleads,nleads,nspin),ictotallk(nleads,nspin))
        ictotmatallk=0D0
        ictotallk=0D0
      endif
    endif


    ti_cur=0D0
    do ispin=1,nspin
      do i1=1,nleads
        do i=1,Nenerg_div
          Ei=ERealGrid%e(i)

          fi=(Ei-ef_bss(i1,ispin))/t
          if(DREAL(fi).gt.0)then
            fi=EXP(-fi)/(1D0+EXP(-fi))
          else
            fi=1.D0/(1D0+EXP(fi))
          endif

          do i2=1,nleads

            fj=(Ei-ef_bss(i2,ispin))/t
            if(DREAL(fj).gt.0)then
              fj=EXP(-fj)/(1D0+EXP(-fj))
            else
              fj=1.D0/(1D0+EXP(fj))
            endif
            tij_cur(i,i1,i2,ispin)=tij(i,i1,i2,ispin) * (fi-fj)

            if(i1.ne.i2) ti_cur(i,i1,ispin)=ti_cur(i,i1,ispin)+ tij_cur(i,i1,i2,ispin)
          enddo
        enddo

      enddo
    enddo
 
    iclocmat=0D0
    do i=1,Nenerg_div
      iclocmat(:,:,:)=iclocmat(:,:,:)+ DREAL(ERealGrid%w(i))* tij_cur(i,:,:,:)
    enddo

#ifdef MPI
#ifdef NODAT
    CALL MPI_REDUCE(iclocmat(1,1,1),ictotmat(1,1,1), nleads * nleads*NSPIN, MPI_DOUBLE_PRECISION,MPI_SUM,0,negf_comm,MPIerror)
#else
    CALL MPI_REDUCE(iclocmat(1,1,1),ictotmat(1,1,1), nleads*nleads *NSPIN, DAT_double,MPI_SUM,0,negf_comm,MPIerror)
#endif
#else
    ictotmat=iclocmat
#endif
    ictotmat=eh_const * ictotmat
    ictot=0D0
    if(MyNode.eq.0)then
      if (outinfo) then
      do ispin=1,nspin
        do i1=1,nleads
         write(12347,'("icmatk= ",1e16.7," ")',ADVANCE='NO') 13.6058D0 * v
          do i2=1,nleads
           write(12347,'(1e16.7," ")',ADVANCE='NO') ictotmat(i1,i2,ispin)
            ictot(i1,ispin)=ictot(i1,ispin)+ictotmat(i1,i2,ispin)
          enddo
           write(12347,*)ispin,ik
        enddo
      enddo
      do i1=1,nleads
        if(nspin.eq.2)then
          write(12347,'("ick= ",3e16.7," ",i7)') 13.6058D0 * v, ictot(i1,1), ictot(i1,2),ik
        else
          write(12347,'("ick= ",2e16.7," ",i7)') 13.6058D0 * v, ictot(i1,1),ik
        endif
      enddo
      endif
     
      if(ik.eq.nk )then

      ictotallk=ictotallk+ictot * wk
      ictotmatallk=ictotmatallk+ictotmat * wk

      if (outinfo) then
        do ispin=1,nspin
          do i1=1,nleads
            write(12347,'("icmat= ",1e16.7," ")',ADVANCE='NO') 13.6058D0 * v
            do i2=1,nleads
              write(12347,'(1e16.7," ")',ADVANCE='NO') ictotmatallk(i1,i2,ispin)
            enddo
            write(12347,*)ispin
          enddo
        enddo
        do i1=1,nleads
          if(nspin.eq.2)then
            write(12347,'("ic= ",4e16.7)') 13.6058D0 * v, ictotallk(i1,1)+ ictotallk(i1,2), ictotallk(i1,1), ictotallk(i1,2)
          else
            write(12347,'("ic= ",3e16.7)') 13.6058D0 * v, 2D0 * ictotallk(i1,1), ictotallk(i1,1)
          endif
        enddo
       endif
 
        deallocate(ictotallk,ictotmatallk)

      endif
    endif


    ti=0D0
    do i1=1,nleads
      do i2=1,nleads
        if(i1.ne.i2) ti(:,i1,:)=ti(:,i1,:)+tij(:,i1,i2,:)
      enddo
    enddo

    if (outinfo) then
    do i1=1,nleads
      do ispin=1,nspin
        write(12347,*)"tijk=  ",ispin,i1,ik
        write(12347,*)"tijk= ",ispin,i1," # ",i1,ik
        write(12347,*)"tij_curk=  ",ispin,i1,ik
        write(12347,*)"tij_curk= ",ispin,i1," # ",i1,ik
      enddo
      write(12347,*)"ti_curk=  "
      write(12347,*)"tik=  "

      do i=1,Nenerg_div
        Ei=ERealGrid%e(i)
        write(12347,*)"tik=  ",i1,13.6057D0 * (DREAL(ei)-ef_lead), ti(i,i1,1)
        write(12347,*)"ti_curk=  ",i1, 13.6057D0 * (DREAL(ei)-ef_lead), ti_cur(i,i1,1)


        do ispin=1,nspin
          write(12347,'("tijk= ",2i5," ",1e16.7," ")',ADVANCE='NO') ispin,i1,13.6057D0 * (DREAL(ei)-ef_lead)
          do i2=nleads,1,-1
            write(12347,'(1e16.7,A3)',ADVANCE='NO') tij(i,i1,i2,ispin),'   '
          enddo
          write(12347,*)ik
        enddo

        do ispin=1,nspin
          write(12347,'("tij_curk= ",2i5," ",1e16.7," ")', ADVANCE='NO') ispin,i1,13.6057D0 * (DREAL(ei)-ef_lead)
          do i2=nleads,1,-1
            write(12347,'(1e16.7,A3)',ADVANCE='NO') tij_cur(i,i1,i2,ispin),'   '
          enddo
          write(12347,*)ik
        enddo

      enddo
    enddo
endif

    tijallk=tijallk+tij * wk
    tijallk_cur=tijallk_cur+tij_cur * wk
    tiallk_cur=tiallk_cur+ti_cur * wk

    if(ik.eq.nk)then

      ti=0D0
      do i1=1,nleads
        do i2=1,nleads
          if(i1.ne.i2) ti(:,i1,:)=ti(:,i1,:)+tijallk(:,i1,i2,:)
        enddo
      enddo

      if (outinfo) then
      do i1=1,nleads

        do ispin=1,nspin
          write(12347,*)"tij=  ",ispin,i1
          write(12347,*)"tij= ",ispin,i1," # ",i1
          write(12347,*)"tij_cur=  ",ispin,i1
          write(12347,*)"tij_cur= ",ispin,i1," # ",i1
        enddo
        write(12347,*)"ti_cur=  "
        write(12347,*)"ti=  "

        do i=1,Nenerg_div
          Ei=ERealGrid%e(i)
          write(12347,*)"ti=  ",i1,13.6057D0 * (DREAL(ei)-ef_lead), ti(i,i1,1)
          write(12347,*)"ti_cur=  ",i1, 13.6057D0 * (DREAL(ei)-ef_lead), tiallk_cur(i,i1,1)

          do ispin=1,nspin
            write(12347,'("tij= ",2i5," ",1e16.7," ")',ADVANCE='NO') ispin,i1,13.6057D0 * (DREAL(ei)-ef_lead)
!              do i2=1,nleads
            do i2=nleads,1,-1
              write(12347,'(1e16.7,A3)',ADVANCE='NO') tijallk(i,i1,i2,ispin),'   '
            enddo
            write(12347,*)
          enddo

          do ispin=1,nspin
            write(12347,'("tij_cur= ",2i5," ",1e16.7," ")', ADVANCE='NO') ispin,i1,13.6057D0 * (DREAL(ei)-ef_lead)
            do i2=nleads,1,-1
              write(12347,'(1e16.7,A3)',ADVANCE='NO') tijallk_cur(i,i1,i2,ispin),'   '
            enddo
            write(12347,*)
          enddo

        enddo
      enddo
  endif


      deallocate(tijallk,tijallk_cur,tiallk_cur)
    endif
    
!#ifdef MPI
!        call MPI_Barrier(negf_comm,MPIerror)
!        call MPI_Finalize( MPIerror )
!#endif
!        stop

    deallocate(ti,tij_cur,ti_cur)

  end subroutine bs_transm


  subroutine transmij_bs(ef_bss,nspin,Nenerg_div,V, n1,delta,ik, nl,nr,nleads,nbss,ef_lead,t,gammamp,sigmamp, nebss,deltabss, writetrc,tij,hgeneral,sgeneral)

    use mTypes
    use mMatrixUtil
    use mEnergyGrid
    use mONBoundstates
#ifdef MPI
    use mMPI_NEGF
#endif
    
    IMPLICIT NONE

    integer nspin,Nenerg_div,n1,ik,nl,nr
    double precision V,delta,ef_bss(nleads,nspin),ef_lead,t, deltabss(nleads)
    double precision iclocmat(nleads,nleads,nspin), ictotmat(nleads,nleads,nspin),ictot(nleads,nspin)
    integer nebss(nleads,2)

    integer ispin,i,j,i1,i2,i3,i4,MPIerror,mynode,Nnodes,ii,jj
    integer nbss,ibs,nleads,nbsmax
    INTEGER IPIV(N1),info
    double complex ei,fi,fj
    DOUBLE COMPLEX, PARAMETER :: zi=(0.D0,1.D0)
    double complex aux(n1,n1)
    double precision tij(Nenerg_div,nleads,nleads,nspin)
    double precision tij2(Nenerg_div,nleads,nleads,nspin)
    logical writetrc
    DOUBLE PRECISION, PARAMETER :: eh_const=1.6022D-4*13.6058D0/4.1357D0
    type(matrixTypeGeneral) :: gammamp(nleads)
    type(matrixTypeGeneral) :: sigmamp(nleads)
    type(matrixTypeGeneral) :: gfgeneral
    type(ioType) :: io
    double complex, allocatable :: gpart(:,:),gpart2(:,:), gammadense1(:,:),gammadense2(:,:),gpart3(:,:),gtransm(:,:)
    type(matrixTypeGeneral) :: hgeneral(nspin),sgeneral

#ifdef MPI
    CALL MPI_COMM_SIZE(negf_comm,Nnodes,MPIerror)
    CALL MPI_COMM_RANK(negf_comm,mynode,MPIerror)
#else
    MyNode=0
    Nnodes=1
#endif
    io%isDebug=.false.
    call AllocateMatrixGeneral(n1,n1,n1*n1,0,gfgeneral, "set_bs_ef", io)

    DO ISPIN=1,NSPIN
         
      DO I=1,Nenerg_div
        Ei=ERealGrid%e(i)
        
        call set_gammamp_lr(n1,nl,nr,nleads,nbss,nebss, gammamp,sigmamp,ERealGrid%sigma(1,i,ispin,ik)%sigma,ERealGrid%sigma(2,i,ispin,ik)%sigma)         

        call setgfelementsdense_bs(ei,ispin,gfgeneral, nleads,sigmamp,nebss,hgeneral(ispin),sgeneral)

!       write(12347,*)"deltagfg=",maxval(abs(gfgeneral%matdense%a-
! GF_iter))

        CALL ZGETRF(N1,N1,gfgeneral%matdense%a,N1,IPIV,INFO)
        CALL ZGETRI(N1,gfgeneral%matdense%a,N1,IPIV,aux,N1**2,INFO)

        if(writetrc)then
          nbsmax=nleads
        else
          nbsmax=nbss
        endif


        do i1=1,nbsmax

          allocate(gammadense1(gammamp(i1)%iRows,gammamp(i1)%iRows))

          gammadense1=0D0
          do ii=1,gammamp(i1)%iRows
            do jj=gammamp(i1)%matSparse%q(ii), gammamp(i1)%matSparse%q(ii+1)-1
              gammadense1(ii, gammamp(i1)%matSparse%j(jj))= gammamp(i1)%matSparse%b(jj)
            enddo
          enddo


          do i2=1,nleads

            allocate(gpart(gammamp(i2)%iRows,gammamp(i1)%iRows))
            allocate(gpart2(gammamp(i2)%iRows,gammamp(i1)%iRows))
            allocate(gpart3(gammamp(i1)%iRows,gammamp(i2)%iRows))
            allocate(gammadense2(gammamp(i2)%iRows,gammamp(i2)%iRows))
            allocate(gtransm(gammamp(i2)%iRows,gammamp(i2)%iRows))


            gammadense2=0D0
            do ii=1,gammamp(i2)%iRows
              do jj=gammamp(i2)%matSparse%q(ii), gammamp(i2)%matSparse%q(ii+1)-1
                gammadense2(ii, gammamp(i2)%matSparse%j(jj))= gammamp(i2)%matSparse%b(jj)
              enddo
            enddo

            gpart(:,:)=gfgeneral%matdense%a(nebss(i2,1):nebss(i2,2), nebss(i1,1):nebss(i1,2))


            call ZGEMM('N','N',gammamp(i2)%iRows,gammamp(i1)%iRows, gammamp(i1)%iRows,(1D0,0D0),gpart, gammamp(i2)%iRows,gammadense1,gammamp(i1)%iRows, (0D0,0D0),gpart2,gammamp(i2)%iRows)

            call ZGEMM('C','N',gammamp(i1)%iRows,gammamp(i2)%iRows, gammamp(i2)%iRows,(1D0,0D0),gpart, gammamp(i2)%iRows,gammadense2,gammamp(i2)%iRows, (0D0,0D0),gpart3,gammamp(i1)%iRows)

            tij(i,i1,i2,ispin)=0D0
            do i3=1,gammamp(i2)%iRows
              do i4=1,gammamp(i1)%iRows
                tij(i,i1,i2,ispin)= tij(i,i1,i2,ispin)+ gpart2(i3,i4)*gpart3(i4,i3)
              enddo
            enddo


            deallocate(gpart)
            deallocate(gpart2)
            deallocate(gpart3)
            deallocate(gammadense2)
            deallocate(gtransm)

          enddo

          deallocate(gammadense1)

        enddo


        do i1=nleads-1,nleads
          call DestroyMatrixGeneral(gammamp(i1),"set_bs_ef",io)
          call DestroyMatrixGeneral(sigmamp(i1),"set_bs_ef",io)
        enddo

      enddo

    enddo
    call DestroyMatrixGeneral(gfgeneral,"set_bs_ef",io)

  end subroutine transmij_bs


  subroutine set_gamma_bs(nspin,n1,S, nl,nr,nleads,nbss,nebss,deltabss,sigmabstot)

    use negfmod, only: nebss_bs,deltabss_bs
    
    IMPLICIT NONE

    integer nspin,n1,nl,nr,nleads,nbss
    double complex s(n1,n1)!,GF_iter(n1,n1)
    integer nebss(nleads,2)
    double precision deltabss(nleads)
    double complex gammacsbs(n1,n1,nleads),sigmabstot(n1,n1)
    DOUBLE COMPLEX, PARAMETER :: zi=(0.D0,1.D0)

    integer ispin,i,j,i1,i2,i3,i4

    deltabss=0D0
    do i=1,nbss
      nebss(i,:)=nebss_bs(i,:)
      deltabss(i)=deltabss_bs(i)
!        write(*,*)"deltabss=",i,deltabss(i)
!        write(*,*)"nebss=",i,nebss(i,1),nebss(i,2)
    enddo

    nebss(nbss+1,1)=1
    nebss(nbss+1,2)=nl
    nebss(nleads,1)=n1-nr+1
    nebss(nleads,2)=n1
!
!      do i=1,nleads
!        write(*,*)"nebss=",nebss(i,1),nebss(i,2)
!      enddo

    gammacsbs=0D0

    do i1=1,nbss
      do i2=nebss(i1,1),nebss(i1,2)
        do i3=1,n1
          gammacsbs(i2,i3,i1)=gammacsbs(i2,i3,i1)+ deltabss(i1) * s(i2,i3)
        enddo
      enddo

      do i2=1,n1
        do i3=nebss(i1,1),nebss(i1,2)
          gammacsbs(i2,i3,i1)=gammacsbs(i2,i3,i1)+ deltabss(i1) * s(i2,i3)
        enddo
      enddo
!        write(11111,*)"mati=",i1
!        call writemat7(DREAL(ei),gammacsbs(:,:,i1),n1,n1,"gammai")
    enddo

    sigmabstot=0D0
    do i1=1,nbss
      sigmabstot(:,:)=sigmabstot(:,:)-0.5D0 * zi * gammacsbs(:,:,i1)
    enddo


!      gf_iter=0D0
!      do i1=1,nbss
!        gf_iter=gf_iter+gammacsbs(:,:,i1)
!      enddo
!      write(*,*)"igammatotal-s",maxval(abs(gf_iter-2D0 * delta * s))

  end subroutine set_gamma_bs




  subroutine find_ef_bss(tij,nspin,nleads,nbss,Nenerg_div,ef_bss, ef_lead,t,v)

  use negfmod,only: outinfo
    IMPLICIT NONE

    integer nspin,nleads,nbss,Nenerg_div
    double precision tij(Nenerg_div,nleads,nleads,nspin), ef_bss(nleads,nspin),ef_bssout(nleads,nspin),ef_lead,t,v

    integer i1,i2,i3
    double precision tolef,alpha

!      ef_bss(1:2,:)=ef_bss(nbss+1,1)
!      ef_bss(3,:)=ef_lead
!      ef_bss(4:5,:)=ef_bss(nleads,1)

    do i1=1,2000
      write(*,*)"i1=",i1

      call find_efout(tij,nspin,nleads,nbss,Nenerg_div,ef_bss, ef_lead,ef_bssout,t,v)

!        do i2=1,nleads
!          if(nspin.eq.2)then
!          write(12347,*)"ef_bss iter=",i1,i2,13.6057D0 * 
!     .        (ef_bss(i2,1)-ef_lead),13.6057D0 * 
!     .        (ef_bssout(i2,1)-ef_lead),
!     .        abs(ef_bss(i2,1)-ef_bssout(i2,1)),13.6057D0 * 
!     .        (ef_bss(i2,2)-ef_lead),13.6057D0 * 
!     .        (ef_bssout(i2,2)-ef_lead),
!     .        abs(ef_bss(i2,2)-ef_bssout(i2,2))
!          else
!          write(12347,*)"ef_bss iter=",i1,i2,13.6057D0 * 
!     .        (ef_bss(i2,1)-ef_lead),13.6057D0 * 
!     .        (ef_bssout(i2,1)-ef_lead),
!     .        abs(ef_bss(i2,1)-ef_bssout(i2,1))
!          endif
!        enddo

      tolef=0D-7
      if (outinfo) then
      write(12347,*)"deltaef=",maxval(abs(ef_bss-ef_bssout)),tolef
      write(12347,'("ef_converge= ",1i5," ")',ADVANCE='NO') i1
!        do i2=1,nleads
      do i2=1,nleads
        write(12347,'(2e16.7,A3)',ADVANCE='NO') 13.6057D0 * (ef_bssout(i2,1)-ef_lead), 13.6057D0 * (ef_bssout(i2,2)-ef_lead),'   '
      enddo
      write(12347,*)
  endif


      if(maxval(abs(ef_bss-ef_bssout)).le.tolef)then
!          write(*,*)"found ef_bss"

         if (outinfo) then
        do i2=1,nleads
          if(nspin.eq.1)then
            write(12347,*)"ef_bssfinal=",i1,i2,13.6057D0 * (ef_bss(i2,1)-ef_lead)
          else
            write(12347,*)"ef_bssfinal=",i1,i2,13.6057D0 * (ef_bss(i2,1)-ef_lead),13.6057D0 * (ef_bss(i2,2)-ef_lead)
          endif
        enddo
    endif

        return
      endif
      alpha=0.9D0
      ef_bss(1:nbss,:)=(1D0-alpha)*ef_bss(1:nbss,:)+alpha * ef_bssout(1:nbss,:) 
    enddo

  end subroutine find_ef_bss


  subroutine find_efout(tij,nspin,nleads,nbss,Nenerg_div,ef_bss, ef_lead,ef_bssout,t,V)

    use mEnergyGrid
#ifdef MPI
    use mMPI_NEGF
#endif


    IMPLICIT NONE

    integer nspin,nleads,nbss,Nenerg_div
    double precision tij(Nenerg_div,nleads,nleads,nspin), ef_bss(nleads,nspin),ef_bssout(nleads,nspin),ef_lead,t,v
    double complex ei,fi,fj

    integer i1,i2,i3,ie,ispin,ief,mpierror,mynode,nnodes
    double precision tidf(Nenerg_div,nspin),ictot(nspin),icloc(nspin), ic(nleads,nspin),tolic,ef_high(nspin),ef_low(nspin)

    
#ifdef MPI
    CALL MPI_COMM_SIZE(negf_comm,Nnodes,MPIerror)
    CALL MPI_COMM_RANK(negf_comm,mynode,MPIerror)
#else
    MyNode=0
    Nnodes=1
#endif

    ef_bssout=ef_bss
    do i1=1,nbss
      do ispin=1,nspin
        if(ef_bss(nbss+1,ispin).gt.ef_bss(nleads,ispin))then
          ef_high(ispin)=ef_bss(nbss+1,ispin)
          ef_low(ispin)=ef_bss(nleads,ispin)
        else
          ef_high(ispin)=ef_bss(nleads,ispin)
          ef_low(ispin)=ef_bss(nbss+1,ispin)
        endif
      enddo
      do ief=1,1000
        tidf=0D0
!          write(12347,*)"tidf"
        write(*,*)"ief=",ief
        do ie=1,Nenerg_div
          Ei=ERealGrid%e(ie)
          do ispin =1,nspin
            fi=(Ei-ef_bssout(i1,ispin))/t
            if(DREAL(fi).gt.0)then
              fi=EXP(-fi)/(1D0+EXP(-fi))
            else
              fi=1.D0/(1D0+EXP(fi))
            endif
            do i2=1,nleads
              if(i2.eq.i1)cycle
              fj=(Ei-ef_bss(i2,ispin))/t
              if(DREAL(fj).gt.0)then
                fj=EXP(-fj)/(1D0+EXP(-fj))
              else
                fj=1.D0/(1D0+EXP(fj))
              endif
              tidf(ie,ispin)=tidf(ie,ispin)+ tij(ie,i1,i2,ispin) *(fi-fj)
            enddo
          enddo
!            write(12347,*)"tidf",13.6057D0 *(dreal(ei)-ef_lead),
!     .          tidf(ie,1),DREAL(fi)
        enddo

        
        icloc=0D0
        do ie=1,Nenerg_div
          do ispin=1,nspin
            icloc(ISPIN)=icloc(ISPIN)+ DREAL(ERealGrid%w(ie))*tidf(ie,ispin)
          enddo
        enddo
!          write(12347,*)"icloc",icloc(1),icloc(2)
!          write(*,*)"icloc",icloc(1)


#ifdef MPI
#ifdef NODAT
        CALL MPI_REDUCE(icloc(1:NSPIN),ictot,NSPIN, MPI_DOUBLE_PRECISION,MPI_SUM,0,negf_comm,MPIerror)
        CALL MPI_BCAST(ictot,NSPIN,MPI_DOUBLE_PRECISION,0, negf_comm,MPIerror)
#else
        CALL MPI_REDUCE(icloc(1:NSPIN),ictot,NSPIN,DAT_double, MPI_SUM,0,negf_comm,MPIerror)
        CALL MPI_BCAST(ictot(1),NSPIN,DAT_double,0, negf_comm,MPIerror)
#endif
#else
        ictot=icloc
#endif

        tolic=1d-13
        if(maxval(abs(ictot)).lt.tolic)exit
        do ispin=1,nspin
          if(ictot(ispin).gt.0)then
            ef_high(ispin)=ef_bssout(i1,ispin)
            ef_bssout(i1,ispin)=0.5D0 * (ef_bssout(i1,ispin)+ ef_low(ispin))
          else
            ef_low(ispin)=ef_bssout(i1,ispin)
            ef_bssout(i1,ispin)=0.5D0 * (ef_bssout(i1,ispin)+ ef_high(ispin))
          endif
        enddo
!         write(12347,*)"ef_bssout1",ef_low(1),ef_high(1),ef_bssout(i1,1)
!         write(12347,*)"ef_bssout2",ef_low(2),ef_high(2),ef_bssout(i1,2)



      enddo
    enddo




  end subroutine find_efout




end module mBoundStates
