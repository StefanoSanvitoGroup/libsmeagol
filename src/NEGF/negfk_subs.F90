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
!                   LEADS_INFO,
!                   OPTIONS_NEGFK,
!                   BROADCAST_LEADS,
!                   ENERGY_POINTS,
!                   CONVERTCOPY_HSLEADS_SIGMA  
! IN THIS FILE ARE LICENSED FOR DISTRIBUTION TO THE SMEAGOL
! COPYRIGHT HOLDERS AND AUTHORS UNDER AND ONLY THE "SMEAGOL 
! ACADEMIC LICENSE" (www.smeagol.tcd.ie). 
!
! FOR INFORMATION OR QUERIES PLEASE CONTACT THE E-MAIL: smeagol@tcd.ie
!

  SUBROUTINE leads_info(slabelL, NL, nspinL, maxnhL, Ef_LeadL, slabelR, NR, nspinR, maxnhR, Ef_LeadR,Ef_Lead,nspin)

! ********************************************************************
! Subroutine to read information about the electrodes
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! ********** HISTORY **********************************************
! Original version:	October 2007
! ********** INPUT ************************************************

  use sigma
  use negfmod, only: ef_em,ef_l,ef_r
#ifdef MPI
  use mMPI_NEGF
#endif
  IMPLICIT NONE

  INTEGER   iu,Mynode,NL,NR,maxnhL,maxnhR,nspinL,nspinR,Nnodes, nspin,nspinloop,nlLoop,nrLoop
  CHARACTER slabelL*20,slabelR*20
  DOUBLE PRECISION Ef_LeadR,Ef_LeadL,Ef_Lead
#ifdef MPI
  INTEGER :: MPIerror
#endif

#ifdef MPI
  CALL MPI_COMM_RANK(negf_comm,mynode,MPIerror)
  CALL MPI_COMM_SIZE(negf_comm,Nnodes,MPIerror)
#else
  Nnodes=1
  MyNode=0
#endif


  if (Mynode .EQ. 0) then
    call io_assign(iu)
    open(iu,file='bulklft.DAT',status='old')
    read(iu,*) slabelL, NL, nspinL, maxnhL, Ef_LeadL
    call io_close(iu)

    call io_assign(iu)
    open(iu,file='bulkrgt.DAT',status='old')
    read(iu,*) slabelR, NR, nspinR, maxnhR, Ef_LeadR
    call io_close(iu)
  endif

#ifdef MPI
#ifdef NODAT
  CALL MPI_BCAST(Ef_LeadL,1,MPI_DOUBLE_PRECISION,0, negf_comm,MPIerror)
  CALL MPI_BCAST(Ef_LeadR,1,MPI_DOUBLE_PRECISION,0, negf_comm,MPIerror)
#else
  CALL MPI_BCAST(Ef_LeadL,1,DAT_double,0,negf_comm,MPIerror)
  CALL MPI_BCAST(Ef_LeadR,1,DAT_double,0,negf_comm,MPIerror)
#endif
  CALL MPI_BCAST(NL,1,MPI_INTEGER,0,negf_comm,MPIerror)
  CALL MPI_BCAST(NR,1,MPI_INTEGER,0,negf_comm,MPIerror)
#endif

  Ef_Lead = Ef_LeadL

  ef_em=Ef_Lead
  ef_l=Ef_LeadL
  ef_r=Ef_LeadR

  nspinloop=nspin
  if(nspin>2)nspinloop=1
  nlLoop=NL
  nrLoop=NR
  if(nspin>2)nlLoop=2 * NL
  if(nspin>2)nrLoop=2 * NR

  SigmaNl=nlLoop
  SigmaNr=nrLoop
  allocate(H0_L(nlLoop,nlLoop,nspinloop), H1_L(nlLoop,nlLoop,nspinloop))
  allocate(S0_L(nlLoop,nlLoop), S1_L(nlLoop,nlLoop))
  allocate(H0_R(nrLoop,nrLoop,nspinloop), H1_R(nrLoop,nrLoop,nspinloop))
  allocate(S0_R(nrLoop,nrLoop), S1_R(nrLoop,nrLoop))

  SigmaNl2=nl
  SigmaNr2=nr
  allocate(H0_L2(NL,NL,NSPIN), H1_L2(NL,NL,NSPIN))
  allocate(S0_L2(NL,NL), S1_L2(NL,NL))
  allocate(H0_R2(NR,NR,NSPIN), H1_R2(NR,NR,NSPIN))
  allocate(S0_R2(NR,NR), S1_R2(NR,NR))
   

  end SUBROUTINE leads_info


  SUBROUTINE OPTIONS_NEGFK(ICASE, PERIODIC)

  IMPLICIT NONE

  INTEGER  ICASE
  LOGICAL  PERIODIC 
  
  If (ICASE .EQ. 1) Then
    PERIODIC=.FALSE.
  Else
    PERIODIC=.TRUE.
  EndIf


  end SUBROUTINE OPTIONS_NEGFK
  
  SUBROUTINE broadcast_leads(NSPIN,NL,H0_L,S0_L,H1_L,S1_L, NR,H0_R,S0_R,H1_R,S1_R)

! ********************************************************************
! Subroutine to broadcast the Hamiltonian and Overlap matrix elements
! to all the compute nodes in a parallel environment
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! ********** HISTORY **********************************************
! Original version:	October 2007
! ********** INPUT ************************************************


#ifdef MPI
  use mMPI_NEGF
#endif

  IMPLICIT NONE
  INTEGER  NSPIN,NL,NR
  DOUBLE COMPLEX  H0_L(NL,NL,NSPIN),H1_L(NL,NL,NSPIN)
  DOUBLE COMPLEX  S0_L(NL,NL),S1_L(NL,NL)
  DOUBLE COMPLEX  H0_R(NR,NR,NSPIN),H1_R(NR,NR,NSPIN)
  DOUBLE COMPLEX  S0_R(NR,NR),S1_R(NR,NR)
#ifdef MPI
  INTEGER :: MPIerror
#endif

#ifdef MPI
#ifdef NODAT
  CALL MPI_BCAST(H0_L,NL*NL*NSPIN,MPI_DOUBLE_COMPLEX,0, negf_comm,MPIerror)
  CALL MPI_BCAST(S0_L,NL*NL,MPI_DOUBLE_COMPLEX,0, negf_comm,MPIerror)
  CALL MPI_BCAST(H1_L,NL*NL*NSPIN,MPI_DOUBLE_COMPLEX,0, negf_comm,MPIerror)
  CALL MPI_BCAST(S1_L,NL*NL,MPI_DOUBLE_COMPLEX,0, negf_comm,MPIerror)
  CALL MPI_BCAST(H0_R,NR*NR*NSPIN,MPI_DOUBLE_COMPLEX,0, negf_comm,MPIerror)
  CALL MPI_BCAST(S0_R,NR*NR,MPI_DOUBLE_COMPLEX,0, negf_comm,MPIerror)
  CALL MPI_BCAST(H1_R,NR*NR*NSPIN,MPI_DOUBLE_COMPLEX,0, negf_comm,MPIerror)
  CALL MPI_BCAST(S1_R,NR*NR,MPI_DOUBLE_COMPLEX,0, negf_comm,MPIerror)
#else
  CALL MPI_BCAST(H0_L(1,1,1),NL*NL*NSPIN,DAT_dcomplex,0, negf_comm,MPIerror)
  CALL MPI_BCAST(S0_L(1,1),NL*NL,DAT_dcomplex,0, negf_comm,MPIerror)
  CALL MPI_BCAST(H1_L(1,1,1),NL*NL*NSPIN,DAT_dcomplex,0, negf_comm,MPIerror)
  CALL MPI_BCAST(S1_L(1,1),NL*NL,DAT_dcomplex,0, negf_comm,MPIerror)
  CALL MPI_BCAST(H0_R(1,1,1),NR*NR*NSPIN,DAT_dcomplex,0, negf_comm,MPIerror)
  CALL MPI_BCAST(S0_R(1,1),NR*NR,DAT_dcomplex,0, negf_comm,MPIerror)
  CALL MPI_BCAST(H1_R(1,1,1),NR*NR*NSPIN,DAT_dcomplex,0, negf_comm,MPIerror)
  CALL MPI_BCAST(S1_R(1,1),NR*NR,DAT_dcomplex,0, negf_comm,MPIerror)
#endif
#else
  write(*,*)"bcast_negf called in serial version, exiting..."
  stop
#endif

  end SUBROUTINE broadcast_leads
  



  SUBROUTINE energy_points(integraltype,EnergI,EnergF, V,Nenerg_div,nenerg_div_nodes, EiCompL,EiCompR,CONST,mynode,nnodes)

! ********************************************************************
! Written by Ivan Rungger, October 2008
! Computational Spintronics Group
! Trinity College Dublin
! e-mail: runggeri@tcd.ie
! ********** HISTORY **********************************************
! Original version:	October 2007
! ********************************************************************

#ifdef MPI
  use mMPI_NEGF
#endif

  IMPLICIT NONE
  
  include "const2.h"
  
  INTEGER Mynode,Nenerg_div,Nnodes, nenerg_div_nodes,I
  DOUBLE PRECISION EnergI,EnergF,V
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION (:) ::  EX,EW
  DOUBLE COMPLEX EiCompL(Nenerg_div),EiCompR(Nenerg_div), CONST(Nenerg_div)
  CHARACTER(LEN=15) :: integraltype
#ifdef MPI
  INTEGER :: MPIerror
  INTEGER, DIMENSION (MPI_STATUS_SIZE) :: istatus
#endif

      IF (MyNode.EQ.0) THEN
	   ALLOCATE(EX(nenerg_div_nodes),EW(nenerg_div_nodes))
       if (integraltype.eq.'gauss-legendre') then
        CALL GAULEG(EnergI,EnergF,EX(1: nenerg_div_nodes),EW(1: nenerg_div_nodes), nenerg_div_nodes)
!           write(*,*)"ex=",ex
       elseif (integraltype.eq.'gauss-chebyshev') then
        CALL GAUCHEB(EnergI,EnergF,EX,EW,Nenerg_div*Nnodes)
       else
        write(6,'(a)') "negf : Wrong type of Gaussian integration"
        write(6,'(a)') "negf : ",integraltype
#ifdef MPI
        CALL MPI_Abort(negf_comm, 1, MPIerror)
#else
        stop
#endif
       endif
#ifdef MPI
       DO I=1,Nnodes-1
#ifdef NODAT
     CALL MPI_SEND(EX(I*Nenerg_div+1:(I+1)*Nenerg_div),Nenerg_div, MPI_DOUBLE_PRECISION,I,1,inverseheads_comm,MPIerror)
     CALL MPI_SEND(EW(I*Nenerg_div+1:(I+1)*Nenerg_div),Nenerg_div, MPI_DOUBLE_PRECISION,I,2,inverseheads_comm,MPIerror)
#else
     CALL MPI_SEND(EX(I*Nenerg_div+1:(I+1)*Nenerg_div),Nenerg_div, DAT_double,I,1,inverseheads_comm,MPIerror)
     CALL MPI_SEND(EW(I*Nenerg_div+1:(I+1)*Nenerg_div),Nenerg_div, DAT_double,I,2,inverseheads_comm,MPIerror)
#endif
       ENDDO
      ELSE
	   ALLOCATE(EX(nenerg_div),EW(nenerg_div))
#ifdef NODAT
       CALL MPI_RECV(EX,Nenerg_div,MPI_DOUBLE_PRECISION,0,1,inverseheads_comm,istatus,MPIerror)
       CALL MPI_RECV(EW,Nenerg_div,MPI_DOUBLE_PRECISION,0,2,inverseheads_comm,istatus,MPIerror)
#else
       CALL MPI_RECV(EX,Nenerg_div,DAT_double,0,1, inverseheads_comm,istatus,MPIerror)
       CALL MPI_RECV(EW,Nenerg_div,DAT_double,0,2, inverseheads_comm,istatus,MPIerror)
#endif
#endif
      ENDIF

      DO I=1,Nenerg_div
       EiCompL(I)=DCMPLX(EX(I)-V/2.D0)
       EiCompR(I)=DCMPLX(EX(I)+V/2.D0)
       CONST(I)=DCMPLX(EW(I))
      ENDDO
	  
	  DEALLOCATE(EX,EW)

  END SUBROUTINE energy_points


subroutine convertcopy_hsleads_sigma(nspin,NspinComplexMatrix,nl,nr)

 use sigma

 implicit none
 integer, intent(in)::nspin,nl,nr,NspinComplexMatrix
 integer i1,i2

 s0_l=0.0D0
 s1_l=0.0D0
 s0_r=0.0D0
 s1_r=0.0D0

 if(nspin<=2)then
   S0_L=S0_L2
   S0_R=S0_R2
   H0_L=H0_L2
   H0_R=H0_R2
   S1_L=S1_L2
   S1_R=S1_R2
   H1_L=H1_L2
   H1_R=H1_R2
 else
   if(NspinComplexMatrix==4)then
     S0_L(1:nl,1:nl)=S0_L2(:,:)
     S0_L(nl+1:2*nl,nl+1:2*nl)=S0_L2(:,:)

     H0_L(1:nl,1:nl,1)=H0_L2(:,:,1)
     H0_L(nl+1:2*nl,nl+1:2*nl,1)=H0_L2(:,:,2)
     H0_L(1:nl,nl+1:2*nl,1)=H0_L2(:,:,3)
     H0_L(nl+1:2* nl,1:nl,1)=H0_L2(:,:,4)

     S1_L(1:nl,1:nl)=S1_L2(:,:)
     S1_L(nl+1:2*nl,nl+1:2*nl)=S1_L2(:,:)

     H1_L(1:nl,1:nl,1)=H1_L2(:,:,1)
     H1_L(nl+1:2*nl,nl+1:2*nl,1)=H1_L2(:,:,2)
     H1_L(1:nl,nl+1:2*nl,1)=H1_L2(:,:,3)
     H1_L(nl+1:2* nl,1:nl,1)=H1_L2(:,:,4)

     S0_r(1:nr,1:nr)=S0_r2(:,:)
     S0_r(nr+1:2*nr,nr+1:2*nr)=S0_r2(:,:)

     H0_r(1:nr,1:nr,1)=H0_r2(:,:,1)
     H0_r(nr+1:2*nr,nr+1:2*nr,1)=H0_r2(:,:,2)
     H0_r(1:nr,nr+1:2*nr,1)=H0_r2(:,:,3)
     H0_r(nr+1:2* nr,1:nr,1)=H0_r2(:,:,4)

     S1_r(1:nr,1:nr)=S1_r2(:,:)
     S1_r(nr+1:2*nr,nr+1:2*nr)=S1_r2(:,:)

     H1_r(1:nr,1:nr,1)=H1_r2(:,:,1)
     H1_r(nr+1:2*nr,nr+1:2*nr,1)=H1_r2(:,:,2)
     H1_r(1:nr,nr+1:2*nr,1)=H1_r2(:,:,3)
     H1_r(nr+1:2* nr,1:nr,1)=H1_r2(:,:,4)
   else

     do i1=1,nl
   do i2=1,nl
     S0_L(2*i1-1,2*i2-1)=S0_L2(i1,i2)
     S0_L(2*i1,2*i2)=S0_L2(i1,i2)
   enddo
     enddo

     do i1=1,nl
   do i2=1,nl
     H0_L(2*i1-1,2*i2-1,1)=H0_L2(i1,i2,1)
     H0_L(2*i1-1,2*i2,1)=H0_L2(i1,i2,3)
     H0_L(2*i1,2*i2,1)=H0_L2(i1,i2,2)
     H0_L(2*i1,2*i2-1,1)=H0_L2(i1,i2,4)
   enddo
     enddo

     do i1=1,nl
   do i2=1,nl
     S1_L(2*i1-1,2*i2-1)=S1_L2(i1,i2)
     S1_L(2*i1,2*i2)=S1_L2(i1,i2)
   enddo
     enddo

     do i1=1,nl
   do i2=1,nl
     H1_L(2*i1-1,2*i2-1,1)=H1_L2(i1,i2,1)
     H1_L(2*i1-1,2*i2,1)=H1_L2(i1,i2,3)
     H1_L(2*i1,2*i2,1)=H1_L2(i1,i2,2)
     H1_L(2*i1,2*i2-1,1)=H1_L2(i1,i2,4)
   enddo
     enddo

     do i1=1,nr
   do i2=1,nr
     S0_r(2*i1-1,2*i2-1)=S0_r2(i1,i2)
     S0_r(2*i1,2*i2)=S0_r2(i1,i2)
   enddo
     enddo

     do i1=1,nr
   do i2=1,nr
     H0_r(2*i1-1,2*i2-1,1)=H0_r2(i1,i2,1)
     H0_r(2*i1-1,2*i2,1)=H0_r2(i1,i2,3)
     H0_r(2*i1,2*i2,1)=H0_r2(i1,i2,2)
     H0_r(2*i1,2*i2-1,1)=H0_r2(i1,i2,4)
   enddo
     enddo

     do i1=1,nr
   do i2=1,nr
     S1_r(2*i1-1,2*i2-1)=S1_r2(i1,i2)
     S1_r(2*i1,2*i2)=S1_r2(i1,i2)
   enddo
     enddo

     do i1=1,nr
   do i2=1,nr
     H1_r(2*i1-1,2*i2-1,1)=H1_r2(i1,i2,1)
     H1_r(2*i1-1,2*i2,1)=H1_r2(i1,i2,3)
     H1_r(2*i1,2*i2,1)=H1_r2(i1,i2,2)
     H1_r(2*i1,2*i2-1,1)=H1_r2(i1,i2,4)
   enddo
     enddo
   endif

 endif

end subroutine convertcopy_hsleads_sigma


